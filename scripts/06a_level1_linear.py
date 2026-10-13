"""06a_level1_linear.py
Level-1 linear learners, out of fold, for every frame, cohort, modality
and seed (stage 6a).

Each modality gets a survival learner and a binary learner. The survival
learner is a ridge-penalised Cox model (scikit-survival
CoxPHSurvivalAnalysis, Efron ties) for four targets: death (time_death,
event_death), composite (time_first, event_comp), cs_cv (CV readmission
first; a death first is censored at its time) and cs_death (death first; a
CV readmission first is censored at its time). It stores the linear
predictor (lp) and, for death and composite, the risk at every horizon,
1 - exp(-H0(h + 0.5) * exp(lp)). The binary learner is L2 logistic
regression at every horizon h (Frame A: 10, 30; Frame B: 10, 30, 90, 180,
365) for death_h, composite_h, cvu_h (CV readmission first by h; a death
first counts as no event; all rows used) and cvc_h (the same label, but
patients dying first by h are left out of training, as in the dissertation
(pe-multimodal-prognosis-msc); predicted for everyone).

Budgets are matched: 5 candidates each on mirror-image grids
(alpha = 1 / C, alpha 1 to 10,000), the stage 4 outer folds, the same
inner 3-fold grouped CV (seeded by the outer seed), median imputation and
scaling fitted inside each training fold, and no class weights (the Cox
likelihood has none). Tuning uses Harrell's C (Cox) and AUROC (logistic).
No performance numbers are written: the checks file holds completeness and
pass/fail sanity checks only. Forests are fitted in stage 6b; the
registered primary (stage 7) uses only these scores.

Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
         feature_whitelist.csv; data/inputs/side_config.csv (seeds),
         ehr28.csv, ecg_opt.csv, ctpa46.csv
Outputs: results/level1_linear_checks.csv, results/level1_linear_log.txt;
         data/processed/level1/oof_linear_{A,B}_{ehrecg,3mod}.csv, a row per
         seed and admission and a column per modality.learner.target
         (patient level, not in the repo);
         data/processed/level1/params_linear.csv (tuning choices, not in the
         repo)
Usage:   python scripts/06a_level1_linear.py [--probe-only]
Notes:   --probe-only times synthetic fits of the real shapes, projects the
         run time and stops; otherwise the run starts after a 5 s pause.
"""
import os
import sys
for v_ in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
import time
import warnings
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
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")
OUTD = os.path.join(PROC, "level1")
NFOLD = 5
NJOBS = max(1, min(10, (os.cpu_count() or 4) - 2))
ALPHAS = [1.0, 10.0, 100.0, 1000.0, 10000.0]
CS = [1.0 / a for a in ALPHAS]
HZ = {"A": [10, 30], "B": [10, 30, 90, 180, 365]}
SURV = ["death", "composite", "cs_cv", "cs_death"]
CODE = {"EHR+ECG": "ehrecg", "3-mod": "3mod"}
MATS = {"ehr": "ehr28.csv", "ecg": "ecg_opt.csv",
        "ctpa": "ctpa46.csv"}


def surv_targets(fr):
    tf = fr["time_first"].values.astype(float)
    er = fr["event_cr"].values
    return {
        "death": (fr["time_death"].values
                  .astype(float),
                  fr["event_death"].values == 1),
        "composite": (tf, fr["event_comp"].values
                      == 1),
        "cs_cv": (tf, er == 1),
        "cs_death": (tf, er == 2)}


def bin_targets(fr, hs):
    td = fr["time_death"].values.astype(float)
    tf = fr["time_first"].values.astype(float)
    er = fr["event_cr"].values
    ed = fr["event_death"].values == 1
    ec = fr["event_comp"].values == 1
    allr = np.ones(len(fr), bool)
    out = {}
    for h in hs:
        cut = h + 0.5
        cvb = (er == 1) & (tf <= cut)
        dfb = (er == 2) & (tf <= cut)
        out["death_h%d" % h] = (ed & (td <= cut),
                                allr)
        out["composite_h%d" % h] = (ec & (tf <= cut),
                                    allr)
        out["cvu_h%d" % h] = (cvb, allr)
        out["cvc_h%d" % h] = (cvb, ~dfb)
    return out


def inner(y, grp, seed):
    cv = StratifiedGroupKFold(
        n_splits=3, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(y)),
                         y.astype(int), grp))


def cox_fit(a, t, e, alpha):
    m = CoxPHSurvivalAnalysis(alpha=alpha,
                              ties="efron")
    m.fit(a, Surv.from_arrays(event=e, time=t))
    return m


def tune_cox(a, t, e, grp, seed, W):
    sc = []
    for al in ALPHAS:
        vals = []
        for t2, v2 in inner(e, grp, seed):
            if e[t2].sum() < 5 or e[v2].sum() < 2:
                continue
            try:
                with warnings.catch_warnings(
                        record=True) as w:
                    warnings.simplefilter("always")
                    m = cox_fit(a[t2], t[t2], e[t2], al)
                W[0] += len(w)
                vals.append(concordance_index_censored(
                    e[v2], t[v2], m.predict(a[v2]))[0])
            except Exception:
                W[1] += 1
        sc.append(float(np.mean(vals)) if vals
                  else -1.0)
    return ALPHAS[int(np.argmax(sc))], sc


def tune_lr(a, y, grp, seed, W):
    sc = []
    for c in CS:
        vals = []
        for t2, v2 in inner(y, grp, seed):
            if (y[t2].sum() < 5 or y[v2].sum() < 2
                    or y[v2].all()):
                continue
            try:
                with warnings.catch_warnings(
                        record=True) as w:
                    warnings.simplefilter("always")
                    m = LogisticRegression(
                        C=c, max_iter=8000).fit(
                        a[t2], y[t2])
                W[0] += len(w)
                vals.append(roc_auc_score(
                    y[v2], m.predict_proba(a[v2])[:, 1]))
            except Exception:
                W[1] += 1
        sc.append(float(np.mean(vals)) if vals
                  else -1.0)
    return CS[int(np.argmax(sc))], sc


def cox_risk(m, lp, hs):
    H = m.cum_baseline_hazard_
    out = {}
    for h in hs:
        i = np.searchsorted(H.x, h + 0.5,
                            side="right") - 1
        h0 = float(H.y[i]) if i >= 0 else 0.0
        out[h] = 1.0 - np.exp(-h0 * np.exp(lp))
    return out


def run_task(frame, cohort, mod, seed, X, ST, BT,
             fold, grp, hs):
    t_start = time.time()
    n = X.shape[0]
    cols, prm, W = {}, [], [0, 0]

    def col(nm):
        if nm not in cols:
            cols[nm] = np.full(n, np.nan)
        return cols[nm]

    for k in range(NFOLD):
        tr = np.where(fold != k)[0]
        te = np.where(fold == k)[0]
        im = SimpleImputer(strategy="median")
        a = im.fit_transform(X[tr])
        b = im.transform(X[te])
        ss = StandardScaler()
        a, b = ss.fit_transform(a), ss.transform(b)
        g_tr = grp[tr]
        for tg in SURV:
            t_, e_ = ST[tg]
            al, sc = tune_cox(a, t_[tr], e_[tr], g_tr,
                              seed, W)
            try:
                with warnings.catch_warnings(
                        record=True) as w:
                    warnings.simplefilter("always")
                    m = cox_fit(a, t_[tr], e_[tr], al)
                W[0] += len(w)
                lp = np.clip(m.predict(b), -50, 50)
                col("%s.cox.%s.lp" % (mod, tg))[te] = lp
                if tg in ("death", "composite"):
                    for h, r in cox_risk(m, lp,
                                         hs).items():
                        col("%s.cox.%s.r%d"
                            % (mod, tg, h))[te] = r
            except Exception:
                W[1] += 1
            prm.append({"fold": k, "learner": "cox",
                        "target": tg, "param": al,
                        "edge": int(al in (ALPHAS[0],
                                           ALPHAS[-1])),
                        "inner": ";".join(
                            "%.4f" % s_ for s_ in sc)})
        for tg, (y, keep) in BT.items():
            sel = keep[tr]
            a2, y2 = a[sel], y[tr][sel].astype(int)
            c, sc = tune_lr(a2, y2, g_tr[sel], seed, W)
            try:
                with warnings.catch_warnings(
                        record=True) as w:
                    warnings.simplefilter("always")
                    m = LogisticRegression(
                        C=c, max_iter=8000).fit(a2, y2)
                W[0] += len(w)
                col("%s.lr.%s" % (mod, tg))[te] = \
                    m.predict_proba(b)[:, 1]
            except Exception:
                W[1] += 1
            prm.append({"fold": k, "learner": "lr",
                        "target": tg, "param": c,
                        "edge": int(c in (CS[0], CS[-1])),
                        "inner": ";".join(
                            "%.4f" % s_ for s_ in sc)})
    return (frame, cohort, mod, seed,
            pd.DataFrame(cols), prm, W,
            time.time() - t_start)


def synth_surv(n, p, rng):
    X = rng.normal(size=(n, p))
    lp = X[:, :3].sum(axis=1) * 0.3
    day = np.floor(rng.exponential(
        300 * np.exp(-lp)))
    e = day <= 365
    return X, np.where(e, day + 0.5, 366.0), e


def probe(jobs):
    print("", flush=True)
    print("=" * 70)
    print("TIMING PROBE (synthetic data, real shapes)")
    print("=" * 70, flush=True)
    rng = np.random.default_rng(0)
    unit = {}
    for p in (264, 46, 28):
        X, t, e = synth_surv(2800, p, rng)
        t0 = time.time()
        cox_fit(X, t, e, 100.0)
        tc = time.time() - t0
        t0 = time.time()
        LogisticRegression(C=0.01, max_iter=8000).fit(
            X, e.astype(int))
        tl = time.time() - t0
        unit[p] = (tc, tl)
        print("  p=%-3d  Cox fit %.2f s   logistic fit"
              " %.2f s (n=2800)" % (p, tc, tl),
              flush=True)
    for j in jobs:
        n, p = j["X"].shape
        pk = 264 if p > 100 else (46 if p > 30 else 28)
        tc, tl = unit[pk]
        s = (0.8 * n) / 2800.0
        n_cox = len(SURV) * (len(ALPHAS) * 3 + 1)
        n_lr = len(j["BT"]) * (len(CS) * 3 + 1)
        j["cost"] = NFOLD * s * (
            n_cox * tc * 0.85 + n_lr * tl * 0.85)
    serial = sum(j["cost"] for j in jobs)
    big = max(j["cost"] for j in jobs)
    par = max(serial / NJOBS * 1.15, big)
    print("  tasks: %d   slowest single task about"
          " %.0f min" % (len(jobs), big / 60))
    print("  projected: serial %.0f min, parallel"
          " (%d workers) about %.0f min"
          % (serial / 60, NJOBS, par / 60))
    print("  (rough: real data may converge more"
          " slowly than synthetic)", flush=True)


def main():
    t_start = time.time()
    probe_only = "--probe-only" in sys.argv
    os.makedirs(OUTD, exist_ok=True)

    fa = pd.read_csv(os.path.join(PROC, "frame_A.csv"))
    fb = pd.read_csv(os.path.join(PROC, "frame_B.csv"))
    folds = pd.read_csv(os.path.join(PROC,
                                     "folds.csv"))
    wl = pd.read_csv(os.path.join(
        PROC, "feature_whitelist.csv"))
    cfg = pd.read_csv(os.path.join(
        INP, "side_config.csv")).set_index("key")
    seeds = [int(s) for s in str(
        cfg.loc["seeds", "value"]).split(",")]
    mats = {m: pd.read_csv(os.path.join(INP, f_))
            .drop_duplicates("hadm_id")
            for m, f_ in MATS.items()}

    jobs, meta = [], {}
    for frame, F in (("A", fa), ("B", fb)):
        for cohort in ("EHR+ECG", "3-mod"):
            sub = (F if cohort == "EHR+ECG"
                   else F[F["in_three_mod"] == 1])
            sub = sub.reset_index(drop=True)
            ST = surv_targets(sub)
            BT = bin_targets(sub, HZ[frame])
            grp = sub["subject_id"].values
            fo = folds[(folds["frame"] == frame)
                       & (folds["cohort"] == cohort)]
            fmap = {}
            for s in seeds:
                fs = sub[["hadm_id"]].merge(
                    fo[fo["seed"] == s][
                        ["hadm_id", "fold"]],
                    on="hadm_id", how="left")["fold"]
                if fs.isna().any():
                    sys.exit("fold missing for %s %s"
                             " seed %d" % (frame, cohort,
                                           s))
                fmap[s] = fs.values.astype(int)
            meta[(frame, cohort)] = (sub, ST, BT, fmap)
            mods = ["ehr", "ecg"] + (
                ["ctpa"] if cohort == "3-mod" else [])
            for mod in mods:
                cc = wl.loc[wl["modality"] == mod,
                            "column"].tolist()
                X = sub[["hadm_id"]].merge(
                    mats[mod][["hadm_id"] + cc],
                    on="hadm_id", how="left")[cc]
                X = X.values.astype(float)
                for s in seeds:
                    jobs.append(dict(
                        frame=frame, cohort=cohort,
                        mod=mod, seed=s, X=X, ST=ST,
                        BT=BT, fold=fmap[s], grp=grp,
                        hs=HZ[frame]))
            print("  %s %-7s n=%d  modalities %s"
                  % (frame, cohort, len(sub), mods),
                  flush=True)

    probe(jobs)
    if probe_only:
        print("\n--probe-only: stopping here.")
        return
    print("\n  starting in 5 s (Ctrl+C to stop)...",
          flush=True)
    time.sleep(5)
    jobs.sort(key=lambda j: -j["cost"])
    res = Parallel(n_jobs=NJOBS, verbose=5)(
        delayed(run_task)(
            j["frame"], j["cohort"], j["mod"],
            j["seed"], j["X"], j["ST"], j["BT"],
            j["fold"], j["grp"], j["hs"])
        for j in jobs)

    print("", flush=True)
    print("=" * 70)
    print("ASSEMBLE AND CHECK")
    print("=" * 70, flush=True)
    CHK, PRM = [], []
    n_warn = n_fail = 0
    blocks = {}
    for fr, co, mod, s, df, prm, W, dt in res:
        blocks.setdefault((fr, co, s), []).append(df)
        for p_ in prm:
            PRM.append(dict(frame=fr, cohort=co,
                            modality=mod, seed=s, **p_))
        n_warn += W[0]
        n_fail += W[1]

    def chk(fr, co, item, value, status):
        CHK.append({"frame": fr, "cohort": co,
                    "item": item, "value": value,
                    "status": status})
        print("  %-4s %s %-7s %-44s %s"
              % (status, fr, co, item, value),
              flush=True)

    for (fr, co), (sub, ST, BT, fmap) in meta.items():
        parts = []
        for s in seeds:
            dfs = blocks[(fr, co, s)]
            w_ = pd.concat(
                [sub[["subject_id", "hadm_id"]]
                 .reset_index(drop=True)]
                + [d.reset_index(drop=True)
                   for d in dfs], axis=1)
            w_.insert(0, "seed", s)
            w_.insert(3, "fold", fmap[s])
            parts.append(w_)
        out = pd.concat(parts, ignore_index=True)
        pc = sorted(c for c in out.columns
                    if "." in c)
        out = out[["seed", "subject_id", "hadm_id",
                   "fold"] + pc]
        fn = "oof_linear_%s_%s.csv" % (fr, CODE[co])
        out.to_csv(os.path.join(OUTD, fn), index=False)
        nan = int(out[pc].isna().sum().sum())
        chk(fr, co, "rows (admissions x seeds)",
            "%d = %d x %d" % (len(out), len(sub),
                              len(seeds)),
            "PASS" if len(out) == len(sub) * len(seeds)
            else "FAIL")
        chk(fr, co, "prediction columns / missing"
            " cells", "%d / %d" % (len(pc), nan),
            "PASS" if nan == 0 else "FAIL")
        avg = out.groupby("hadm_id", sort=False)[
            pc].mean().reindex(sub["hadm_id"].values)
        n_bad = 0
        for c in pc:
            mod, lrn, tg = c.split(".")[:3]
            if lrn == "cox" and not c.endswith(".lp"):
                continue
            v = avg[c].values
            if np.isnan(v).any():
                n_bad += 1
                continue
            if lrn == "cox":
                t_, e_ = ST[tg]
                ok = concordance_index_censored(
                    e_, t_, v)[0] > 0.5
            else:
                y_ = BT[tg][0]
                ok = roc_auc_score(y_, v) > 0.5
            n_bad += int(not ok)
        chk(fr, co, "scores ranking no better than"
            " chance", n_bad,
            "PASS" if n_bad == 0 else "FAIL")

    pr = pd.DataFrame(PRM)
    pr.to_csv(os.path.join(OUTD, "params_linear.csv"),
              index=False)
    chk("-", "-", "final fits that failed", n_fail,
        "PASS" if n_fail == 0 else "FAIL")
    chk("-", "-", "warnings during fitting", n_warn,
        "INFO")
    for lrn, grid in (("cox", ALPHAS), ("lr", CS)):
        s_ = pr[pr["learner"] == lrn]
        for mod in ("ehr", "ecg", "ctpa"):
            q = s_[s_["modality"] == mod]
            if not len(q):
                continue
            lo = float((q["param"] == grid[0]).mean())
            hi = float((q["param"] == grid[-1]).mean())
            chk("-", "-", "%s %s: share at grid ends"
                " (first/last)" % (lrn, mod),
                "%.2f / %.2f" % (lo, hi),
                "WARN" if max(lo, hi) > 0.5
                else "INFO")

    ck = pd.DataFrame(CHK)
    ck.to_csv(os.path.join(
        RES, "level1_linear_checks.csv"), index=False)
    nf = int((ck["status"] == "FAIL").sum())
    print("")
    print("  FAIL %d   WARN %d   elapsed %.0f min"
          % (nf, int((ck["status"] == "WARN").sum()),
             (time.time() - t_start) / 60))
    print("")
    print("OUTPUT FILES")
    for (fr, co) in meta:
        print("  " + os.path.join(
            OUTD, "oof_linear_%s_%s.csv"
            % (fr, CODE[co])) + "  (local only)")
    print("  " + os.path.join(OUTD,
                              "params_linear.csv"))
    print("  " + os.path.join(
        RES, "level1_linear_checks.csv"))
    print("  " + os.path.join(
        RES, "level1_linear_log.txt"),
        "(terminal log)")
    print("")
    print("Outputs: level1_linear_checks.csv")


if __name__ == "__main__":
    main()