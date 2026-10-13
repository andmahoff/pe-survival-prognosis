"""18j_arch_neural.py
Neural networks with one encoder per modality and a fusion layer, the
intermediate fusion arms of A10 (PREREGISTRATION_ADDENDUM_2.md, A26).

Fitted for death in every frame, cohort and modality set of A10, and every
seed and fold of 04_frames.py. In the network (p2_nn.py, numpy) each
modality's features pass through their own encoder (dense + ReLU: EHR 16
units, ECG 32, CTPA 16). The encodings are concatenated and passed through
a fusion layer (32 units, ReLU) and a head. Multimodal dropout of 0.25
removes each modality's encoding for a patient with probability 0.25
during training, never all of them; ordinary dropout is tuned.

Intermediate PH has one output, the Cox partial likelihood (DeepSurv head),
with risk at every horizon from the training fold's Breslow baseline.
Intermediate TV has one output per interval between consecutive horizons
(0 to 10.5, 10.5 to 30.5, ...), the discrete-time logistic hazard (Kvamme
and Borgan 2021), so each modality's effect may differ between intervals;
risk is 1 minus the product of (1 - hazard).

Each head has two candidates (dropout 0.1 or 0.4), chosen on a validation
split of 20% of the training fold's patients (grouped, stratified by
death) with early stopping (patience 15, at most 300 epochs), then refitted
on the whole training fold for the best number of epochs. Adam (1e-3),
batches of 256, weight decay 1e-4; preprocessing as in 18a_early_linear.R,
per modality, fitted on the training fold. Risks are averaged over seeds
and written in the layout of 18a (hadm_id, model, set, horizon, risk) for
18k_arch_eval2.R.

Checks:  analytic against numerical gradients (FAIL beyond 1e-5); every
         patient predicted, risks in [0, 1] and never falling with the
         horizon, D2 per-fold rule (FAIL); share of fits choosing each
         dropout (INFO).
Notes:   each finished task is saved in fusion_arch/neural_ckpt/. A timing
         probe fits one fold of the largest set on synthetic data, projects
         the run time and pauses 5 s.

Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
         feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
         ctpa46.csv
Outputs: data/processed/fusion_arch/arch_preds_neural_*.csv, with
         params_arch_neural.csv and neural_ckpt/ (not in the repo);
         results/arch_neural_checks.csv, results/arch_neural_log.txt
Usage:   python scripts/18j_arch_neural.py [--probe-only] [--workers=N]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ["PYTHONPATH"] = HERE + os.pathsep + os.environ.get(
    "PYTHONPATH", "")
import p2_nn as nn

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
IN = os.path.join(P2, "data", "inputs")
PR = os.path.join(P2, "data", "processed")
FA = os.path.join(PR, "fusion_arch")
RS = os.path.join(P2, "results")
CKPT = os.path.join(FA, "neural_ckpt")
os.makedirs(CKPT, exist_ok=True)
HZ = {"A": [10, 30], "B": [10, 30, 90, 180, 365]}
MATS = {"ehr": "ehr28.csv", "ecg": "ecg_opt.csv", "ctpa": "ctpa46.csv"}
SETS = {"EHR+ECG": {"EHR+ECG": ["ehr", "ecg"]},
        "3-mod": {"EHR+ECG": ["ehr", "ecg"], "EHR+CTPA": ["ehr", "ctpa"],
                  "EHR+ECG+CTPA": ["ehr", "ecg", "ctpa"]}}
ENC = {"ehr": 16, "ecg": 32, "ctpa": 16}
TRUNK = 32
MD = 0.25
CANDS = [0.1, 0.4]
LR, WD, BS, EPOCHS, PATIENCE = 1e-3, 1e-4, 256, 300, 15
VALFRAC = 0.2
args = sys.argv[1:]
PROBE_ONLY = "--probe-only" in args
wk = [a for a in args if a.startswith("--workers=")]
WORKERS = (int(wk[0].split("=")[1]) if wk
           else max(1, (os.cpu_count() or 2) // 2 - 1))
T0 = time.time()
CHK = []


def chk(item, value, status):
    CHK.append({"item": item, "value": str(value), "status": status})
    print("  %-5s %-62s %s" % (status, item, value), flush=True)


def prep(X, tr):
    X = X.copy()
    for j in range(X.shape[1]):
        col = X[tr, j]
        m = np.nanmedian(col) if np.isfinite(col).any() else 0.0
        X[np.isnan(X[:, j]), j] = m
        if len(np.unique(X[tr, j])) > 2:
            lo, hi = np.quantile(X[tr, j], [0.005, 0.995])
            X[:, j] = np.clip(X[:, j], lo, hi)
    mu = X[tr].mean(0)
    sd = X[tr].std(0, ddof=1)
    keep = np.isfinite(sd) & (sd > 0)
    return (X[:, keep] - mu[keep]) / sd[keep]


def inner_split(subj, ev, rows, seed):
    rng = np.random.RandomState(seed)
    s = subj[rows]
    us, inv = np.unique(s, return_inverse=True)
    pos = np.zeros(len(us))
    np.maximum.at(pos, inv, ev[rows].astype(float))
    val = np.zeros(len(us), bool)
    for g in (0, 1):
        ii = np.where(pos == g)[0]
        rng.shuffle(ii)
        val[ii[:int(round(VALFRAC * len(ii)))]] = True
    return rows[~val[inv]], rows[val[inv]]


def run_task(tk, probe_folds=None):
    t_task = time.time()
    H = HZ[tk["fr"]]
    cuts = np.array(H, float) + 0.5
    t, e = tk["t"], tk["e"]
    n = len(t)
    out = {"Intermediate PH": np.full((n, len(H)), np.nan),
           "Intermediate TV": np.full((n, len(H)), np.nan)}
    prm = []
    folds = sorted(np.unique(tk["fold"])) if probe_folds is None \
        else probe_folds
    for k in folds:
        te = np.where(tk["fold"] == k)[0]
        trall = np.where(tk["fold"] != k)[0]
        tr, va = inner_split(tk["subj"], e, trall,
                             tk["seed"] * 100 + int(k))
        parts = [prep(X, trall) for X in tk["Xs"]]
        Xs = np.hstack(parts)
        blocks, s0 = [], 0
        for P_ in parts:
            blocks.append(np.arange(s0, s0 + P_.shape[1]))
            s0 += P_.shape[1]
        enc = [ENC[m] for m in tk["mods"]]
        heads = {"Intermediate PH": ("cox", 1, {"t": t, "e": e}),
                 "Intermediate TV": ("hazard", len(H),
                                     nn.hazard_targets(t, e, cuts))}
        for hi, (am, (loss, k_out, tgt)) in enumerate(heads.items()):
            sd = tk["seed"] * 10000 + int(k) * 100 + hi * 10
            best = None
            for ci, dr in enumerate(CANDS):
                net = nn.Net(blocks, enc, TRUNK, k_out, dropout=dr, md=MD,
                             rng=np.random.RandomState(sd + ci))
                ep, vl = nn.train(net, Xs[tr], nn.subset_target(tgt, tr),
                                  loss, lr=LR, wd=WD, bs=BS, epochs=EPOCHS,
                                  Xv=Xs[va],
                                  tv=nn.subset_target(tgt, va),
                                  patience=PATIENCE, seed=sd + ci)
                if best is None or vl < best[2]:
                    best = (dr, max(ep, 5), vl)
            dr, ep, vl = best
            net = nn.Net(blocks, enc, TRUNK, k_out, dropout=dr, md=MD,
                         rng=np.random.RandomState(sd + 7))
            nn.train(net, Xs[trall], nn.subset_target(tgt, trall), loss,
                     lr=LR, wd=WD, bs=BS, epochs=ep, seed=sd + 7)
            o_te = net.predict(Xs[te])
            if loss == "cox":
                o_tr = net.predict(Xs[trall])
                out[am][te] = nn.breslow_risk(o_tr[:, 0], t[trall],
                                              e[trall], o_te[:, 0], cuts)
            else:
                out[am][te] = nn.hazard_risk(o_te)
            prm.append({"fold": int(k), "arm": am, "dropout": dr,
                        "epochs": ep, "val_loss": vl})
    res = {"fr": tk["fr"], "co": tk["co"], "set": tk["set"],
           "seed": tk["seed"], "R": out, "prm": prm,
           "min": (time.time() - t_task) / 60}
    if probe_folds is None:
        tmp = tk["ckpt"] + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(res, f)
        os.replace(tmp, tk["ckpt"])
    return res


def cindex(t, e, s):
    ev = np.where(e == 1)[0]
    num = den = 0.0
    for i in ev:
        later = t > t[i]
        den += later.sum()
        num += (s[later] < s[i]).sum() + 0.5 * (s[later] == s[i]).sum()
    return num / den if den > 0 else np.nan


def main():
    print("=== GRADIENT CHECK ===")
    gc = nn.grad_check(seed=1)
    worst = max(gc.values())
    chk("analytic vs numerical gradients (worst, all losses)",
        "%.1e" % worst, "PASS" if worst < 1e-5 else "FAIL")
    if worst >= 1e-5:
        pd.DataFrame(CHK).to_csv(os.path.join(
            RS, "arch_neural_checks.csv"), index=False)
        sys.exit(1)
    print("\n=== BUILD ===")
    fr_d = {f: pd.read_csv(os.path.join(PR, "frame_%s.csv" % f))
            for f in ("A", "B")}
    folds = pd.read_csv(os.path.join(PR, "folds.csv"))
    wl = pd.read_csv(os.path.join(PR, "feature_whitelist.csv"))
    mats = {m: pd.read_csv(os.path.join(IN, f)).drop_duplicates("hadm_id")
            for m, f in MATS.items()}
    tasks, META = [], {}
    for fr in ("A", "B"):
        for co in SETS:
            sub = fr_d[fr] if co == "EHR+ECG" else \
                fr_d[fr][fr_d[fr]["in_three_mod"] == 1]
            sub = sub.reset_index(drop=True)
            code = "ehrecg" if co == "EHR+ECG" else "3mod"
            fo = folds[(folds["frame"] == fr) & (folds["cohort"] == co)]
            seeds = sorted(fo["seed"].unique())
            fmap = {s: fo[fo["seed"] == s].set_index("hadm_id")["fold"]
                    .reindex(sub["hadm_id"]).to_numpy() for s in seeds}
            XM = {}
            for m in MATS:
                cc = [c for c in wl.loc[wl["modality"] == m, "column"]
                      if c in mats[m].columns]
                XM[m] = (mats[m].set_index("hadm_id")[cc]
                         .reindex(sub["hadm_id"]).to_numpy(float))
            META[(fr, co)] = {"sub": sub, "code": code, "seeds": seeds,
                              "fmap": fmap, "sets": list(SETS[co])}
            for sn, ms in SETS[co].items():
                for s in seeds:
                    if np.isnan(fmap[s].astype(float)).any():
                        sys.exit("fold missing: %s %s seed %s"
                                 % (fr, co, s))
                    tasks.append({
                        "fr": fr, "co": co, "set": sn, "mods": ms,
                        "seed": int(s), "Xs": [XM[m] for m in ms],
                        "fold": fmap[s].astype(int),
                        "subj": sub["subject_id"].to_numpy(),
                        "t": sub["time_death"].to_numpy(float),
                        "e": sub["event_death"].to_numpy(float),
                        "ckpt": os.path.join(CKPT, "%s_%s_%s_%d.pkl" % (
                            fr, code, sn.replace("+", ""), s))})
            print("  %s %-7s n=%d  sets %s  seeds %d"
                  % (fr, co, len(sub), ", ".join(SETS[co]), len(seeds)))

    print("\n=== TIMING PROBE (synthetic data, largest task's shape) ===")
    size = lambda t: len(t["t"]) * sum(X.shape[1] for X in t["Xs"])
    big = max(tasks, key=size)
    rng = np.random.RandomState(1)
    syn = dict(big)
    syn["Xs"] = [rng.randn(*X.shape) for X in big["Xs"]]
    t0 = time.time()
    run_task(syn, probe_folds=[0])
    per_fold = time.time() - t0
    todo = [t for t in tasks if not os.path.exists(t["ckpt"])]
    cost = sum(per_fold * 5 * size(t) / size(big) for t in todo)
    print("  largest task: %.1f s per fold (both heads, 2 candidates each)"
          % per_fold)
    print("  %d of %d tasks to fit: about %.1f hours serial, %.1f hours on"
          " %d workers" % (len(todo), len(tasks), cost / 3600,
                           cost / 3600 / WORKERS * 1.15, WORKERS))
    if PROBE_ONLY:
        print("--probe-only: stopping here.")
        return
    print("  starting in 5 s (Ctrl+C to stop)...", flush=True)
    time.sleep(5)
    if todo:
        todo.sort(key=lambda t: -size(t))
        Parallel(n_jobs=WORKERS, backend="loky", verbose=5)(
            delayed(run_task)(t) for t in todo)

    print("\n=== ASSEMBLE AND CHECK ===")
    res = []
    for t in tasks:
        with open(t["ckpt"], "rb") as f:
            res.append(pickle.load(f))
    PRM = []
    for (fr, co), M in META.items():
        sub = M["sub"]
        H = HZ[fr]
        n = len(sub)
        long, nbad, mono = [], 0, 0
        for sn in M["sets"]:
            rr = [r for r in res if r["fr"] == fr and r["co"] == co
                  and r["set"] == sn]
            for r in rr:
                for q in r["prm"]:
                    PRM.append(dict(q, frame=fr, cohort=co, set=sn,
                                    seed=r["seed"]))
            for am in ("Intermediate PH", "Intermediate TV"):
                P = np.mean([r["R"][am] for r in rr], axis=0)
                mono += int((np.diff(P, axis=1) < -1e-9).any(1).sum())
                for r in rr:
                    f = M["fmap"][r["seed"]].astype(int)
                    v = [cindex(sub["time_death"].to_numpy()[f == k],
                                sub["event_death"].to_numpy()[f == k],
                                r["R"][am][f == k, -1])
                         for k in sorted(np.unique(f))]
                    if np.nanmean(v) < 0.40:
                        nbad += 1
                long.append(pd.DataFrame({
                    "hadm_id": np.tile(sub["hadm_id"].to_numpy(), len(H)),
                    "model": am, "set": sn,
                    "horizon": np.repeat(H, n), "risk": P.T.reshape(-1)}))
        L = pd.concat(long, ignore_index=True)
        L.to_csv(os.path.join(FA, "arch_preds_neural_%s_%s.csv"
                              % (fr, M["code"])), index=False)
        bad = int((~np.isfinite(L["risk"])).sum())
        ok = bad == 0 and L["risk"].min() >= 0 and L["risk"].max() <= 1
        chk("%s %s: missing / risk range" % (fr, co),
            "%d / %.4f to %.4f" % (bad, L["risk"].min(), L["risk"].max()),
            "PASS" if ok else "FAIL")
        chk("%s %s: patients whose risk falls with the horizon"
            % (fr, co), mono, "PASS" if mono == 0 else "FAIL")
        chk("%s %s: seed-level scores failing the per-fold rule (D2)"
            % (fr, co), nbad, "PASS" if nbad == 0 else "FAIL")
    P = pd.DataFrame(PRM)
    P.to_csv(os.path.join(FA, "params_arch_neural.csv"), index=False)
    for (am, dr), g in P.groupby(["arm", "dropout"]):
        chk("%s: share of fits choosing dropout %.1f" % (am, dr),
            "%.0f%%" % (100 * len(g) / (P["arm"] == am).sum()), "INFO")
    ck = pd.DataFrame(CHK)
    ck.to_csv(os.path.join(RS, "arch_neural_checks.csv"), index=False)
    print("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours"
          % ((ck.status == "PASS").sum(), (ck.status == "WARN").sum(),
             (ck.status == "FAIL").sum(), (time.time() - T0) / 3600))
    print("\nOUTPUT FILES")
    for (fr, co), M in META.items():
        print("  " + os.path.join(FA, "arch_preds_neural_%s_%s.csv"
                                  % (fr, M["code"])), " (local only)")
    print("  " + os.path.join(FA, "params_arch_neural.csv"),
          " (local only)")
    print("  " + os.path.join(RS, "arch_neural_checks.csv"))
    print("  " + os.path.join(RS, "arch_neural_log.txt"),
          " (terminal log)")
    print("\nOutputs: arch_neural_checks.csv")


if __name__ == "__main__":
    main()