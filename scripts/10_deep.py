"""10_deep.py
Registered neural level-1 learners, out of fold over the stage 4 folds,
per frame, cohort, modality and seed (PREREGISTRATION.md, sections 4 and
6, S4; PREREGISTRATION_ADDENDUM_2.md, A26).

Learners (numpy, in p2_nn.py; one network per learner and endpoint):
  DeepSurv  Cox partial-likelihood network (Katzman et al. 2018) for
            death, composite and cause-specific CV readmission (deaths
            first censored): a linear predictor and, for death and
            composite, the risk at every horizon from the Breslow baseline.
  DeepHit   competing risks (Lee et al. 2018), CV readmission first
            (cause 1) and death first (cause 2): cumulative incidence of
            CV readmission at every horizon. 20 intervals whose edges
            include every horizon; alpha 0.5 (pycox's default of 0.2
            inflated the incidence twofold in a simulation), sigma 0.1.
  MLP twin  the binary twin: one output per horizon with binary
            cross-entropy, for death, composite and unconditional CV
            readmission by h (the labels of stage 6a's logistic twins).

Two hidden layers (64 and 32 units, ReLU), dropout, L2 weight decay 1e-4,
Adam (learning rate 1e-3), batches of 256. The tuning budget is matched
across learners (section 4): two candidates each (dropout 0.1 or 0.5),
chosen on a validation split of 20% of the training fold's patients
(grouped by patient, stratified by any event) with early stopping
(patience 15, at most 300 epochs). The chosen network is refitted on the
whole training fold for its best number of epochs. Preprocessing is as in
stage 18a (median imputation, winsorising at the 0.5th and 99.5th
percentiles, standardisation), fitted on the training fold.

Checks (FAIL): analytic gradients match finite differences to 1e-5, every
patient is scored, risks lie in [0, 1] and never fall with the horizon,
and the death scores pass the D2 per-fold rule.

Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
         feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
         ctpa46.csv; scripts/p2_nn.py
Outputs: results/level1_deep_checks.csv, results/level1_deep_log.txt,
         data/processed/level1/oof_deep_*.csv (patient level, not in the
         repo; the layout of stages 6a and 6b, read by stage 11e)
Usage:   python scripts/10_deep.py [--probe-only] [--workers=N]
Notes:   A timing probe on synthetic data of the largest task's shape
         pauses 5 s; --probe-only stops after it. Tasks are checkpointed.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import re
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
L1 = os.path.join(PR, "level1")
RS = os.path.join(P2, "results")
CKPT = os.path.join(L1, "deep_ckpt")
os.makedirs(CKPT, exist_ok=True)
HZ = {"A": [10, 30], "B": [10, 30, 90, 180, 365]}
MATS = {"ehr": "ehr28.csv", "ecg": "ecg_opt.csv", "ctpa": "ctpa46.csv"}
CANDS = [0.1, 0.5]
UNITS = 64
LR, WD, BS, EPOCHS, PATIENCE = 1e-3, 1e-4, 256, 300, 15
VALFRAC = 0.2
NGRID = 20
ALPHA, SIGMA = 0.5, 0.1
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
    """Validation patients (VALFRAC), grouped and stratified by event."""
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


def specs_for(fr, t_d, e_d, t_f, e_c, e_r):
    H = HZ[fr]
    cut = np.array(H) + 0.5
    lab = lambda e, t: np.column_stack([(e == 1) & (t <= c)
                                        for c in cut]).astype(float)
    S = []
    S.append(("dsurv", "death", "cox", 1, {"t": t_d, "e": e_d}))
    S.append(("dsurv", "composite", "cox", 1, {"t": t_f, "e": e_c}))
    S.append(("dsurv", "cs_cv", "cox", 1,
              {"t": t_f, "e": (e_r == 1).astype(float)}))
    S.append(("dhit", "cv", "deephit", None, None))
    for nm, e, t in (("death", e_d, t_d), ("composite", e_c, t_f),
                     ("cvu", (e_r == 1).astype(float), t_f)):
        Y = lab(e, t)
        S.append(("mlp", nm, "bce", len(H),
                  {"Y": Y, "M": np.ones_like(Y)}))
    return S


def fit_one(Xs, tr, va, te, kind, loss, k_out, tgt, fr, seed):
    """Tune dropout on (tr, va), refit on tr + va (the whole training
    fold), predict te. Returns (test output, training output, chosen
    dropout, epochs, validation loss, training rows)."""
    best = None
    full = np.concatenate([tr, va])
    for ci, dr in enumerate(CANDS):
        net = nn.Net([np.arange(Xs.shape[1])], [UNITS], UNITS // 2, k_out,
                     dropout=dr, rng=np.random.RandomState(seed + ci))
        ep, vl = nn.train(net, Xs[tr], nn.subset_target(tgt, tr), loss,
                          lr=LR, wd=WD, bs=BS, epochs=EPOCHS,
                          Xv=Xs[va], tv=nn.subset_target(tgt, va),
                          patience=PATIENCE, seed=seed + ci)
        if best is None or vl < best[2]:
            best = (dr, max(ep, 5), vl)
    dr, ep, vl = best
    net = nn.Net([np.arange(Xs.shape[1])], [UNITS], UNITS // 2, k_out,
                 dropout=dr, rng=np.random.RandomState(seed + 7))
    nn.train(net, Xs[full], nn.subset_target(tgt, full), loss, lr=LR,
             wd=WD, bs=BS, epochs=ep, seed=seed + 7)
    return net.predict(Xs[te]), net.predict(Xs[full]), dr, ep, vl, full


def run_task(tk, probe_folds=None):
    t_task = time.time()
    fr = tk["fr"]
    H = HZ[fr]
    cuts = np.array(H) + 0.5
    n = len(tk["t_d"])
    m = tk["mod"]
    S = specs_for(fr, tk["t_d"], tk["e_d"], tk["t_f"], tk["e_c"],
                  tk["e_r"])
    cols = {}
    prm = []
    folds = sorted(np.unique(tk["fold"])) if probe_folds is None \
        else probe_folds
    for k in folds:
        te = np.where(tk["fold"] == k)[0]
        trall = np.where(tk["fold"] != k)[0]
        ev_any = (tk["e_r"] > 0).astype(float)
        tr, va = inner_split(tk["subj"], ev_any, trall,
                             tk["seed"] * 100 + int(k))
        Xs = prep(tk["X"], trall)
        for si, (kind, en, loss, k_out, tgt) in enumerate(S):
            sd = tk["seed"] * 10000 + int(k) * 100 + si * 10
            if kind == "dhit":
                ev_t = tk["t_f"][trall][tk["e_r"][trall] > 0]
                q = np.quantile(ev_t, np.linspace(0, 1, NGRID)) \
                    if len(ev_t) else np.array([])
                edges = np.unique(np.concatenate([
                    cuts, np.ceil(q) + 0.5]))
                edges = edges[(edges > 0) & (edges <= cuts[-1])]
                K = len(edges)
                tgt = nn.deephit_targets(tk["t_f"], tk["e_r"], edges, 2)
                tgt["alpha"], tgt["sigma"] = ALPHA, SIGMA
                loss = "deephit"
                k_out = 2 * K + 1
            o_te, o_tr, dr, ep, vl, full = fit_one(
                Xs, tr, va, te, kind, loss, k_out, tgt, fr, sd)
            prm.append({"fold": int(k), "learner": kind, "endpoint": en,
                        "dropout": dr, "epochs": ep, "val_loss": vl})
            if kind == "dsurv":
                c = "%s.dsurv.%s.lp" % (m, en)
                cols.setdefault(c, np.full(n, np.nan))[te] = o_te[:, 0]
                if en != "cs_cv":
                    rk = nn.breslow_risk(o_tr[:, 0], tgt["t"][full],
                                         tgt["e"][full], o_te[:, 0], cuts)
                    for j, h in enumerate(H):
                        c = "%s.dsurv.%s.r%d" % (m, en, h)
                        cols.setdefault(c, np.full(n, np.nan))[te] = \
                            rk[:, j]
            elif kind == "dhit":
                cif = nn.deephit_cif(o_te, 2, K)
                for j, h in enumerate(H):
                    jj = int(np.where(np.isclose(edges, h + 0.5))[0][0])
                    c = "%s.dhit.cv.r%d" % (m, h)
                    cols.setdefault(c, np.full(n, np.nan))[te] = \
                        cif[:, 0, jj]
            else:
                p = nn.sigmoid(o_te)
                for j, h in enumerate(H):
                    c = "%s.mlp.%s_h%d" % (m, en, h)
                    cols.setdefault(c, np.full(n, np.nan))[te] = p[:, j]
    out = {"fr": fr, "co": tk["co"], "mod": m, "seed": tk["seed"],
           "cols": cols, "prm": prm, "min": (time.time() - t_task) / 60}
    if probe_folds is None:
        tmp = tk["ckpt"] + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(out, f)
        os.replace(tmp, tk["ckpt"])
    return out


def main():
    print("=== GRADIENT CHECK ===")
    gc = nn.grad_check(seed=1)
    worst = max(gc.values())
    chk("analytic vs numerical gradients (worst, all losses)",
        "%.1e (%s)" % (worst, ", ".join("%s %.0e" % kv
                                        for kv in gc.items())),
        "PASS" if worst < 1e-5 else "FAIL")
    if worst >= 1e-5:
        pd.DataFrame(CHK).to_csv(os.path.join(
            RS, "level1_deep_checks.csv"), index=False)
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
        for co in ("EHR+ECG", "3-mod"):
            sub = fr_d[fr] if co == "EHR+ECG" else \
                fr_d[fr][fr_d[fr]["in_three_mod"] == 1]
            sub = sub.reset_index(drop=True)
            code = "ehrecg" if co == "EHR+ECG" else "3mod"
            mods = ["ehr", "ecg"] + (["ctpa"] if co == "3-mod" else [])
            fo = folds[(folds["frame"] == fr) & (folds["cohort"] == co)]
            seeds = sorted(fo["seed"].unique())
            fmap = {}
            for s in seeds:
                f_ = fo[fo["seed"] == s].set_index("hadm_id")["fold"]
                fmap[s] = f_.reindex(sub["hadm_id"]).to_numpy()
                if np.isnan(fmap[s].astype(float)).any():
                    sys.exit("fold missing: %s %s seed %s" % (fr, co, s))
            META[(fr, co)] = {"sub": sub, "code": code, "mods": mods,
                              "seeds": seeds, "fmap": fmap}
            for m in mods:
                cc = [c for c in wl.loc[wl["modality"] == m, "column"]
                      if c in mats[m].columns]
                X = (mats[m].set_index("hadm_id")[cc]
                     .reindex(sub["hadm_id"]).to_numpy(float))
                for s in seeds:
                    tasks.append({
                        "fr": fr, "co": co, "mod": m, "seed": int(s),
                        "X": X, "fold": fmap[s].astype(int),
                        "subj": sub["subject_id"].to_numpy(),
                        "t_d": sub["time_death"].to_numpy(float),
                        "e_d": sub["event_death"].to_numpy(float),
                        "t_f": sub["time_first"].to_numpy(float),
                        "e_c": sub["event_comp"].to_numpy(float),
                        "e_r": sub["event_cr"].to_numpy(int),
                        "ckpt": os.path.join(CKPT, "%s_%s_%s_%d.pkl"
                                             % (fr, code, m, s))})
            print("  %s %-7s n=%d  modalities %s  seeds %d"
                  % (fr, co, len(sub), "/".join(mods), len(seeds)))

    print("\n=== TIMING PROBE (synthetic data, largest task's shape) ===")
    big = max(tasks, key=lambda t: t["X"].shape[0] * t["X"].shape[1])
    rng = np.random.RandomState(1)
    n, p = big["X"].shape
    syn = dict(big)
    syn["X"] = rng.randn(n, p)
    t0 = time.time()
    run_task(syn, probe_folds=[0])
    per_fold = time.time() - t0
    cost = sum(per_fold * 5 * (t["X"].shape[0] * t["X"].shape[1])
               / (n * p) for t in tasks)
    print("  largest task: %.1f s per fold (all 7 learners, 2 candidates"
          " each)" % per_fold)
    todo = [t for t in tasks if not os.path.exists(t["ckpt"])]
    frac = len(todo) / max(1, len(tasks))
    print("  %d of %d tasks to fit: about %.1f hours serial, %.1f hours on"
          " %d workers" % (len(todo), len(tasks), cost * frac / 3600,
                           cost * frac / 3600 / WORKERS * 1.15, WORKERS))
    if PROBE_ONLY:
        print("--probe-only: stopping here.")
        return
    print("  starting in 5 s (Ctrl+C to stop)...", flush=True)
    time.sleep(5)
    if todo:
        todo.sort(key=lambda t: -t["X"].shape[0] * t["X"].shape[1])
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
        parts = []
        for s in M["seeds"]:
            rr = [r for r in res if r["fr"] == fr and r["co"] == co
                  and r["seed"] == s]
            d = pd.DataFrame({"seed": s, "subject_id": sub["subject_id"],
                              "hadm_id": sub["hadm_id"],
                              "fold": M["fmap"][s].astype(int)})
            for r in rr:
                for c, v in r["cols"].items():
                    d[c] = v
                for q in r["prm"]:
                    PRM.append(dict(q, frame=fr, cohort=co,
                                    modality=r["mod"], seed=s))
            parts.append(d)
        out = pd.concat(parts, ignore_index=True)
        pc = [c for c in out.columns
              if c not in ("seed", "subject_id", "hadm_id", "fold")]
        out.to_csv(os.path.join(L1, "oof_deep_%s_%s.csv"
                                % (fr, M["code"])), index=False)
        nmiss = int(out[pc].isna().sum().sum())
        chk("%s %s: rows / columns / missing cells" % (fr, co),
            "%d / %d / %d" % (len(out), len(pc), nmiss),
            "PASS" if nmiss == 0 and len(out) == len(sub) *
            len(M["seeds"]) else "FAIL")
        rc = [c for c in pc if re.search(r"\.r\d+$", c)
              or ".mlp." in c]
        rk = out[rc].to_numpy(float)
        ok = bool(np.nanmin(rk) >= 0 and np.nanmax(rk) <= 1)
        mono = 0
        for m in M["mods"]:
            for stem in ("dsurv.death", "dsurv.composite", "dhit.cv"):
                cs = ["%s.%s.r%d" % (m, stem, h) for h in HZ[fr]]
                if all(c in out for c in cs):
                    mono += int((np.diff(out[cs].to_numpy(), axis=1)
                                 < -1e-9).any(axis=1).sum())
        chk("%s %s: risks in [0, 1]; rows with a risk falling over time"
            % (fr, co), "%s; %d" % (ok, mono),
            "PASS" if ok and mono == 0 else "FAIL")
        nbad = 0
        tdead = sub["time_death"].to_numpy()
        edead = sub["event_death"].to_numpy()
        for c in [c for c in pc if c.endswith(".lp")]:
            v = []
            for s in M["seeds"]:
                os_ = out[out["seed"] == s]
                for k in sorted(os_["fold"].unique()):
                    ii = np.where(os_["fold"].to_numpy() == k)[0]
                    sc = os_[c].to_numpy()[ii]
                    v.append(cindex(tdead[ii], edead[ii], sc)
                             if ".death." in c else np.nan)
            v = [x for x in v if np.isfinite(x)]
            if v and np.mean(v) < 0.40:
                nbad += 1
        chk("%s %s: death scores failing the per-fold rule (D2)"
            % (fr, co), nbad, "PASS" if nbad == 0 else "FAIL")
    P = pd.DataFrame(PRM)
    P.to_csv(os.path.join(L1, "params_deep.csv"), index=False)
    for (kind, dr), g in P.groupby(["learner", "dropout"]):
        chk("%s: share of fits choosing dropout %.1f" % (kind, dr),
            "%.0f%%" % (100 * len(g) / (P["learner"] == kind).sum()),
            "INFO")
    ck = pd.DataFrame(CHK)
    ck.to_csv(os.path.join(RS, "level1_deep_checks.csv"), index=False)
    print("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours"
          % ((ck.status == "PASS").sum(), (ck.status == "WARN").sum(),
             (ck.status == "FAIL").sum(), (time.time() - T0) / 3600))
    print("\nOUTPUT FILES")
    for (fr, co), M in META.items():
        print("  " + os.path.join(L1, "oof_deep_%s_%s.csv"
                                  % (fr, M["code"])), " (local only)")
    print("  " + os.path.join(L1, "params_deep.csv"), " (local only)")
    print("  " + os.path.join(RS, "level1_deep_checks.csv"))
    print("  " + os.path.join(RS, "level1_deep_log.txt"),
          " (terminal log)")
    print("\nOutputs: level1_deep_checks.csv")


def cindex(t, e, s):
    """Harrell's C (higher score = earlier event), ties in score count
    one half."""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    s = np.asarray(s, float)
    ev = np.where(e == 1)[0]
    num = den = 0.0
    for i in ev:
        later = t > t[i]
        den += later.sum()
        num += (s[later] < s[i]).sum() + 0.5 * (s[later] == s[i]).sum()
    return num / den if den > 0 else np.nan


if __name__ == "__main__":
    main()