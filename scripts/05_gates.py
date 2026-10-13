"""05_gates.py
Pre-registered gates G1 to G5, Python half (stage 5a).

No analysis result is produced here; the gates decide whether the analysis
starts. G1 (leakage): no label or timing column appears in any feature
matrix. G2 (reproduction): the learner of the side project
(pe-fusion-transfer), script 88, is copied verbatim (oof(), BEST_C, CT_CS,
seeds, folds, row order) and run on this study's frozen matrices. It must
reproduce that project's fusion_final.csv within 0.005 AUROC for every
MIMIC-trained single-modality reference: ecg_opt in both cohorts and ctpa
in the three-modality cohort. The EHR has no MIMIC-trained reference
there, because its score is zero-shot INSPECT, which stage 1 reproduced
exactly.

On the ECG death score in Frame A, G3 (equivalence) requires the
time-dependent AUC at 30.5 days to equal the ordinary death_30d AUROC to
within 5e-4. G4 (sign) requires Harrell's C above 0.5, its negation to
give 1 - C, and 200 permuted scores to average 0.50 +- 0.01. G5 requires
Harrell's C to equal Uno's C truncated at day 31, to 1e-4. 05b_gates.R
repeats G3 to G5 with the R evaluator.

Inputs:  data/inputs/mimic_labels_harmonised.csv, cohorts.csv,
         ecg_opt.csv, ctpa46.csv, side_config.csv;
         data/processed/feature_whitelist.csv, frames_config.csv,
         frame_A.csv; fusion_final.csv of the side project, copied to
         data/inputs/reference/ with its md5 recorded in gates.csv
Outputs: results/gates.csv, results/gates_log.txt;
         data/processed/gate_scores.csv (patient level, not in the repo)
Usage:   python scripts/05_gates.py [--probe-only]
Notes:   --probe-only times fits on synthetic data of the real shapes,
         projects the run time and stops. Otherwise the fits run in
         parallel (joblib) after a 5 s pause.
"""
import os
import sys
for v_ in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
import time
import shutil
import hashlib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import (
    StratifiedGroupKFold)
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sksurv.util import Surv
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc)

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
SIDE = os.path.join(DS, "SideProject")
INP = os.path.join(P2, "data", "inputs")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")
REFD = os.path.join(INP, "reference")
KEY = ["subject_id", "hadm_id"]
TOL = 0.005
NFOLD = 5
NJOBS = max(1, min(10, (os.cpu_count() or 4) - 2))

# ======= copied from SIDE PROJECT script 88 =======
CT_CS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0,
         1e1, 1e2]
BEST_C = {"cv_first": 1e-4,
          "composite_30d": 1e-3,
          "death_30d": 1e-3,
          "death_30d_inhosp": 3e-4}
OUTS = ["death_30d", "composite_30d",
        "cv_first", "death_30d_inhosp"]


def oof(X, y, grp, seed, C, cw="balanced",
        tune=None):
    p = np.zeros(len(y))
    cv = StratifiedGroupKFold(
        n_splits=NFOLD, shuffle=True,
        random_state=seed)
    for tr, te in cv.split(X, y, grp):
        im = SimpleImputer(strategy="median")
        a = im.fit_transform(X[tr])
        b = im.transform(X[te])
        sc = StandardScaler()
        a, b = sc.fit_transform(a), sc.transform(b)
        c = C
        if tune:
            best = -1.0
            k = min(3, max(2,
                           int(y[tr].sum())
                           // 10))
            try:
                icv = StratifiedGroupKFold(
                    n_splits=k, shuffle=True,
                    random_state=42)
                for cand in tune:
                    qq = np.zeros(len(tr))
                    for t2, v2 in icv.split(
                            a, y[tr], grp[tr]):
                        mm = LogisticRegression(
                            C=cand,
                            max_iter=8000,
                            class_weight=cw)
                        mm.fit(a[t2],
                               y[tr][t2])
                        qq[v2] = \
                            mm.predict_proba(
                                a[v2])[:, 1]
                    s = roc_auc_score(
                        y[tr], qq)
                    if s > best:
                        best, c = s, cand
            except Exception:
                c = C
        m = LogisticRegression(
            C=c, max_iter=8000,
            class_weight=cw)
        m.fit(a, y[tr])
        p[te] = m.predict_proba(b)[:, 1]
    return p


def labels(df, outcome):
    # fusion_lib.labels, copied
    d = df
    if outcome == "cv_first":
        d = df[df["death_first"] == 0]
    y = pd.to_numeric(d[outcome],
                      errors="coerce")
    y = y.fillna(0).astype(int).values
    return d, y
# ============== end of copied code ==============


GATES = []


def rec(gate, cohort, outcome, model, item, ref,
        val, status, note=""):
    GATES.append({"gate": gate, "cohort": cohort,
                  "outcome": outcome,
                  "model": model, "item": item,
                  "reference": ref, "value": val,
                  "status": status, "note": note})
    print("  %-4s %-3s %-7s %-16s %-7s %-34s"
          " ref %-9s got %s"
          % (status, gate, cohort, outcome, model,
             item, ref, val), flush=True)


def banner(t):
    print("", flush=True)
    print("=" * 78)
    print(t)
    print("=" * 78, flush=True)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20),
                      b""):
            h.update(b)
    return h.hexdigest()


def synth(n, p, prev, rng):
    X = rng.normal(size=(n, p))
    X[rng.random((n, p)) < 0.04] = np.nan
    y = (rng.random(n) < prev).astype(int)
    grp = rng.integers(0, int(n * 0.9), n)
    return X, y, grp


def probe(n_full, n_3e, n_ct):
    banner("TIMING PROBE (synthetic data, real"
           " shapes)")
    rng = np.random.default_rng(0)
    out = []
    for lab_, (n, p, cw, tune, C) in (
            ("ECG, EHR+ECG cohort",
             (3501, 264, "balanced", None, 1e-3)),
            ("ECG, three-modality",
             (1636, 264, "balanced", None, 1e-3)),
            ("CTPA with inner C tuning",
             (1636, 46, None, CT_CS, 1.0))):
        X, y, g = synth(n, p, 0.11, rng)
        t0 = time.time()
        oof(X, y, g, 42, C, cw=cw, tune=tune)
        dt = time.time() - t0
        out.append(dt)
        print("  %-26s %6.1f s per 5-fold fit"
              % (lab_, dt), flush=True)
    serial = (n_full * out[0] + n_3e * out[1]
              + n_ct * out[2])
    par = serial / NJOBS * 1.2
    print("  fits: %d + %d + %d" % (n_full, n_3e,
                                    n_ct))
    print("  projected: serial %.1f min, parallel"
          " (%d workers) about %.1f min"
          % (serial / 60, NJOBS, par / 60))
    print("  (rough: real data can take longer"
          " to converge than synthetic)")
    return par


def main():
    t_start = time.time()
    probe_only = "--probe-only" in sys.argv
    os.makedirs(REFD, exist_ok=True)

    banner("LOAD (script-88 row order)")
    lab = pd.read_csv(os.path.join(
        INP, "mimic_labels_harmonised.csv"))
    coh = pd.read_csv(os.path.join(
        INP, "cohorts.csv"))
    ecg = pd.read_csv(os.path.join(
        INP, "ecg_opt.csv"))
    ctp = pd.read_csv(os.path.join(
        INP, "ctpa46.csv"))
    wl = pd.read_csv(os.path.join(
        PROC, "feature_whitelist.csv"))
    cfg = pd.read_csv(os.path.join(
        INP, "side_config.csv")).set_index("key")
    seeds = [int(s) for s in str(
        cfg.loc["seeds", "value"]).split(",")]
    LG = wl.loc[wl["block"] == "ecg_logit",
                "column"].tolist()
    MC = wl.loc[wl["block"] == "ecg_meas",
                "column"].tolist()
    CT = wl.loc[wl["modality"] == "ctpa",
                "column"].tolist()

    k1 = set(coh.loc[coh["in_ehr_ecg"] == 1,
                     "hadm_id"])
    full = ecg[ecg["hadm_id"].isin(k1)].merge(
        lab, on=KEY, how="inner")
    k3 = set(coh.loc[coh["in_three_mod"] == 1,
                     "hadm_id"])
    ctc = ctp[ctp["hadm_id"].isin(k3)].merge(
        lab, on=KEY, how="inner").merge(
        ecg.drop(columns=["subject_id"]),
        on="hadm_id", how="inner")
    print("  EHR+ECG rows %d, three-modality rows"
          " %d, seeds %s" % (len(full), len(ctc),
                             seeds))

    banner("G1 LEAKAGE")
    banned = (set(lab.columns) - set(KEY)) | {
        "time_death", "event_death", "time_first",
        "event_cr", "event_comp", "first_adm",
        "hospice", "final_year_possible",
        "in_three_mod", "in_ehr_ecg"}
    for nm, cols, want in (("ecg_logit", LG, 71),
                           ("ecg_meas", MC, 193),
                           ("ctpa", CT, 46)):
        bad = [c for c in cols if c in banned]
        rec("G1", "-", "-", nm,
            "banned columns / feature count", "0/%d"
            % want, "%d/%d" % (len(bad), len(cols)),
            "PASS" if not bad and len(cols) == want
            else "FAIL")
    for nm, n_, want in (("EHR+ECG", len(full), 3501),
                         ("3-mod", len(ctc), 1636)):
        rec("G1", nm, "-", "-", "cohort rows", want,
            n_, "PASS" if n_ == want else "FAIL")

    tasks = []
    for cohort, D in (("ehr_ecg", full),
                      ("3mod", ctc)):
        for oc in OUTS:
            d, y = labels(D, oc)
            if y.sum() < 30:
                continue
            grp = d["subject_id"].values
            mods = [("ecg_opt", LG + MC, BEST_C[oc],
                     "balanced", None)]
            if cohort == "3mod":
                mods.append(("ctpa", CT, 1.0, None,
                             CT_CS))
            for m, cols, C, cw, tune in mods:
                X = d[cols].values.astype(float)
                for s in seeds:
                    tasks.append((cohort, oc, m, s, X,
                                  y, grp, C, cw, tune,
                                  d["hadm_id"].values))
    n_full = sum(1 for t in tasks
                 if t[0] == "ehr_ecg")
    n_3e = sum(1 for t in tasks if t[0] == "3mod"
               and t[2] == "ecg_opt")
    n_ct = sum(1 for t in tasks if t[2] == "ctpa")
    probe(n_full, n_3e, n_ct)
    if probe_only:
        print("\n--probe-only: stopping here.")
        return
    print("\n  starting in 5 s (Ctrl+C to stop)...",
          flush=True)
    time.sleep(5)

    banner("G2 REPRODUCTION: %d fits on %d workers"
           % (len(tasks), NJOBS))
    res = Parallel(n_jobs=NJOBS, verbose=5)(
        delayed(oof)(t[4], t[5], t[6], t[3], t[7],
                     t[8], t[9]) for t in tasks)

    src = os.path.join(SIDE, "data", "processed",
                       "fusion_final.csv")
    if not os.path.exists(src):
        print("  SIDE PROJECT fusion_final.csv not"
              " found: " + src)
        sys.exit(1)
    dst = os.path.join(REFD, "fusion_final.csv")
    shutil.copy2(src, dst)
    ff = pd.read_csv(dst)
    note = "md5 " + md5(dst)

    rows = []
    for t, p in zip(tasks, res):
        rows.append({"cohort": t[0], "outcome": t[1],
                     "model": t[2], "seed": t[3],
                     "auc": roc_auc_score(t[5], p)})
    ra = pd.DataFrame(rows)
    for (c_, o_, m_), g_ in ra.groupby(
            ["cohort", "outcome", "model"],
            sort=False):
        r = ff[(ff["model_type"] == c_)
               & (ff["outcome"] == o_)
               & (ff["model"] == m_)]
        mean_ = float(g_["auc"].mean())
        s42 = float(g_.loc[g_["seed"] == seeds[0],
                           "auc"].iloc[0])
        if not len(r):
            rec("G2", c_, o_, m_,
                "mean AUROC over seeds", "none",
                round(mean_, 4), "FAIL",
                "no reference row")
            continue
        ref = float(r["mean"].iloc[0])
        rec("G2", c_, o_, m_, "mean AUROC over seeds",
            round(ref, 4), round(mean_, 4),
            "PASS" if abs(mean_ - ref) <= TOL
            else "FAIL", note)
        if "seed42" in r.columns:
            r42 = float(r["seed42"].iloc[0])
            rec("G2", c_, o_, m_,
                "seed %d AUROC (information)"
                % seeds[0], round(r42, 4),
                round(s42, 4), "INFO")

    banner("G3-G5 ON FRAME A (ECG death score)")
    conf = pd.read_csv(os.path.join(
        PROC, "frames_config.csv")).set_index("frame")
    cut = 30 + float(conf.loc["A", "auc_cut_offset"])
    tau_c = 30 + float(conf.loc["A", "c_tau_offset"])
    sel = [(t, p) for t, p in zip(tasks, res)
           if t[0] == "ehr_ecg"
           and t[1] == "death_30d"
           and t[2] == "ecg_opt"]
    score = np.mean([p for _, p in sel], axis=0)
    hadm = sel[0][0][10]
    fa = pd.read_csv(os.path.join(
        PROC, "frame_A.csv"))
    sc = pd.DataFrame({"hadm_id": hadm,
                       "score_ecg_death": score}).merge(
        fa[["subject_id", "hadm_id", "time_death",
            "event_death"]], on="hadm_id",
        how="inner")
    rec("G3", "EHR+ECG", "death", "ecg_opt",
        "rows matched to Frame A", len(full),
        len(sc), "PASS" if len(sc) == len(full)
        else "FAIL")
    ev = sc["event_death"].astype(bool).values
    tt = sc["time_death"].values.astype(float)
    s_ = sc["score_ecg_death"].values
    yb = ev & (tt <= cut)
    auc_b = roc_auc_score(yb, s_)
    ys = Surv.from_arrays(event=ev, time=tt)
    tauc = cumulative_dynamic_auc(
        ys, ys, s_, [cut])[0][0]
    rec("G3", "EHR+ECG", "death", "ecg_opt",
        "binary AUROC death by 30.5 d", "-",
        round(auc_b, 6), "INFO")
    rec("G3", "EHR+ECG", "death", "ecg_opt",
        "time-dependent AUC at 30.5 d",
        round(auc_b, 6), round(tauc, 6),
        "PASS" if abs(tauc - auc_b) < 5e-4
        else "FAIL")
    h = concordance_index_censored(ev, tt, s_)[0]
    hn = concordance_index_censored(ev, tt, -s_)[0]
    rng = np.random.default_rng(42)
    cp = [concordance_index_censored(
        ev, tt, rng.permutation(s_))[0]
        for _ in range(200)]
    rec("G4", "EHR+ECG", "death", "ecg_opt",
        "Harrell C", "> 0.5", round(h, 6),
        "PASS" if h > 0.5 else "FAIL")
    rec("G4", "EHR+ECG", "death", "ecg_opt",
        "C of negated score", round(1 - h, 6),
        round(hn, 6),
        "PASS" if abs(h + hn - 1) < 1e-6
        else "FAIL")
    rec("G4", "EHR+ECG", "death", "ecg_opt",
        "mean C of 200 permutations", "0.50+-0.01",
        round(float(np.mean(cp)), 6),
        "PASS" if abs(np.mean(cp) - 0.5) <= 0.01
        else "FAIL")
    u = concordance_index_ipcw(ys, ys, s_,
                               tau=tau_c)[0]
    rec("G5", "EHR+ECG", "death", "ecg_opt",
        "Uno C tau 31", round(h, 6), round(u, 6),
        "PASS" if abs(u - h) < 1e-4 else "FAIL")

    banner("WRITE")
    gt = pd.DataFrame(GATES)
    gt.to_csv(os.path.join(RES, "gates.csv"),
              index=False)
    sc.to_csv(os.path.join(PROC, "gate_scores.csv"),
              index=False)
    nf = int((gt["status"] == "FAIL").sum())
    npass = int((gt["status"] == "PASS").sum())
    print("  PASS %d   FAIL %d" % (npass, nf))
    if nf:
        print(gt[gt["status"] == "FAIL"][
            ["gate", "cohort", "outcome", "model",
             "item", "reference", "value"]]
            .to_string(index=False))
    print("  elapsed %.1f min"
          % ((time.time() - t_start) / 60))
    print("")
    print("OUTPUT FILES")
    print("  " + os.path.join(RES, "gates.csv"))
    print("  " + os.path.join(RES, "gates_log.txt"),
          "(terminal log)")
    print("  " + os.path.join(PROC, "gate_scores.csv"),
          " (patient-level, local only)")
    print("  " + dst)
    print("")
    print("NEXT: run 05b_gates.R, then check"
          " gates.csv and gates_r.csv")


if __name__ == "__main__":
    main()