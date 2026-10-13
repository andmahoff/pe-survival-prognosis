"""08_late_fusion_linear.py
Late fusion of the linear level-1 scores by weighted rank averaging (S1,
S3 and their three-modality versions, S5).

Each modality's out-of-fold score becomes a percentile rank within its
seed. Ranks are combined with weights on a 0.05 grid that sum to 1. The
weights are chosen inside each outer training fold and applied to the
held-out fold, with the same folds and seeds as stage 6a.

S1, formulation factorial (Frame A, 30 days; death, composite and CV
readmission as competing risks). Level 1 is either the binary twin at 30
days (lr: death_h30, composite_h30, cvu_h30) or the ridge-Cox linear
predictor (cox: death, composite, cs_cv). The weights are chosen by either
AUROC for an event by day 30 (auc) or Uno's C truncated at day 31 (uno).
The four cells give four fused scores, each scored on both measures beside
every single modality.
S3, horizon curve (Frame B, death): weights chosen by AUROC at 10, 30, 90,
180 and 365 days, from the Cox scores (one model for every horizon) and
from the binary twins trained at each horizon, plus the Cox scores fused
by Uno's C to day 366.
S5: the same in the three-modality cohorts, with CTPA as a third modality
(three-way weight grid).

The numbers are descriptive: pooled over seed-averaged scores, plus the
per-fold average (deviation D2). Formal contrasts with intervals come from
the R evaluator (stage 11).

Inputs:  data/processed/frame_A.csv, frame_B.csv, frames_config.csv,
         data/processed/level1/oof_linear_<frame>_<cohort>.csv,
         data/inputs/side_config.csv
Outputs: results/fusion_weights.csv, results/fusion_s1.csv,
         results/fusion_s3.csv, results/fusion_log.txt,
         data/processed/fusion/fused_*.csv (patient level, not in the repo)
Usage:   python scripts/08_late_fusion_linear.py [--probe-only]
Notes:   Parallel with joblib over (frame, cohort, seed) on physical
         cores minus 1, at most 10. A timing probe on real tuning units
         projects the run time and pauses 5 s; --probe-only stops after
         it. On Windows the machine is kept awake during the run.
"""
import os
import sys
for v_ in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
import time
import ctypes
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score
from sksurv.util import Surv
from sksurv.metrics import concordance_index_ipcw

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")
L1D = os.path.join(PROC, "level1")
FUD = os.path.join(PROC, "fusion")
HB = [10, 30, 90, 180, 365]
STEP = 0.05
CODE = {"EHR+ECG": "ehrecg", "3-mod": "3mod"}


def physical_cores():
    try:
        from joblib import cpu_count
        n = cpu_count(only_physical_cores=True)
    except Exception:
        n = max(1, (os.cpu_count() or 2) // 2)
    return max(1, n)


NJOBS = max(1, min(10, physical_cores() - 1))


def keep_awake(on):
    if os.name != "nt":
        return
    ctypes.windll.kernel32.SetThreadExecutionState(
        0x80000000 | (0x00000001 if on else 0))


def grid(k):
    n = int(round(1 / STEP))
    if k == 2:
        return [np.array([1 - i / n, i / n])
                for i in range(n + 1)]
    out = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            out.append(np.array([i / n, j / n,
                                 (n - i - j) / n]))
    return out


def pct_rank(v):
    return pd.Series(v).rank(pct=True).values


def uno(t, e, s, tau):
    y = Surv.from_arrays(event=e, time=t)
    return concordance_index_ipcw(y, y, s, tau=tau)[0]


def make_specs(fr, sub, conf):
    td = sub["time_death"].values.astype(float)
    tf = sub["time_first"].values.astype(float)
    ed = sub["event_death"].values == 1
    ec = sub["event_comp"].values == 1
    er = sub["event_cr"].values
    off = float(conf.loc[fr, "auc_cut_offset"])
    ctau = float(conf.loc[fr, "c_tau_offset"])
    S = []
    if fr == "A":
        cut, tau = 30 + off, 30 + ctau
        ends = {"death": ("cox.death.lp",
                          "lr.death_h30",
                          ed & (td <= cut), td, ed),
                "composite": ("cox.composite.lp",
                              "lr.composite_h30",
                              ec & (tf <= cut), tf, ec),
                "cv": ("cox.cs_cv.lp", "lr.cvu_h30",
                       (er == 1) & (tf <= cut), tf,
                       er == 1)}
        for en, (cx, lr, lab, t_, e_) in ends.items():
            for l1, col in (("cox", cx), ("lr", lr)):
                for crit in ("auc", "uno"):
                    S.append(dict(
                        name="fuse.%s.%s.%s"
                        % (l1, crit, en),
                        l1col=col, crit=crit, label=lab,
                        t=t_, e=e_, tau=tau, endpoint=en,
                        level1=l1, horizon=30))
    else:
        tau = 365 + ctau
        for h in HB:
            lab = ed & (td <= h + off)
            for l1, col in (("cox", "cox.death.lp"),
                            ("lr", "lr.death_h%d" % h)):
                S.append(dict(
                    name="fuse.%s.auc%d.death" % (l1, h),
                    l1col=col, crit="auc", label=lab,
                    t=td, e=ed, tau=tau, endpoint="death",
                    level1=l1, horizon=h))
        S.append(dict(
            name="fuse.cox.uno.death",
            l1col="cox.death.lp", crit="uno",
            label=ed & (td <= 365 + off), t=td, e=ed,
            tau=tau, endpoint="death", level1="cox",
            horizon=365))
    return S


def run_task(fr, co, seed, mods, SC, fold, specs):
    t0 = time.time()
    gw = grid(len(mods))
    n = len(fold)
    fused, W = {}, []
    for sp in specs:
        M = SC[sp["l1col"]]
        R = np.column_stack([pct_rank(M[:, j])
                             for j in range(M.shape[1])])
        F = np.stack([R @ w for w in gw], axis=1)
        out = np.full(n, np.nan)
        for k in range(5):
            tr = fold != k
            te = ~tr
            if sp["crit"] == "auc":
                y = sp["label"][tr]
                sc = [roc_auc_score(y, F[tr, g])
                      for g in range(F.shape[1])]
            else:
                t_, e_ = sp["t"][tr], sp["e"][tr]
                sc = [uno(t_, e_, F[tr, g], sp["tau"])
                      for g in range(F.shape[1])]
            g = int(np.argmax(sc))
            out[te] = F[te, g]
            rec = dict(frame=fr, cohort=co, seed=seed,
                       fold=k, fused=sp["name"],
                       endpoint=sp["endpoint"],
                       level1=sp["level1"],
                       criterion=sp["crit"],
                       horizon=sp["horizon"])
            for j, m in enumerate(mods):
                rec["w_" + m] = float(gw[g][j])
            W.append(rec)
        fused[sp["name"]] = out
    return fr, co, seed, pd.DataFrame(fused), W, \
        time.time() - t0


def metric_pair(scores, folds, kind, y=None, t=None,
                e=None, tau=None):
    def met(idx, s):
        try:
            if kind == "auc":
                return roc_auc_score(y[idx], s[idx])
            return uno(t[idx], e[idx], s[idx], tau)
        except Exception:
            return np.nan
    n = len(scores[0])
    pooled = met(np.arange(n), np.mean(scores, axis=0))
    pf = []
    for s_, f_ in zip(scores, folds):
        for k in np.unique(f_):
            v = met(np.where(f_ == k)[0], s_)
            if not np.isnan(v):
                pf.append(v)
    return pooled, (float(np.mean(pf)) if pf
                    else np.nan)


def main():
    t_start = time.time()
    probe_only = "--probe-only" in sys.argv
    os.makedirs(FUD, exist_ok=True)
    keep_awake(True)
    conf = pd.read_csv(os.path.join(
        PROC, "frames_config.csv")).set_index("frame")
    cfg = pd.read_csv(os.path.join(
        INP, "side_config.csv")).set_index("key")
    seeds = [int(s) for s in str(
        cfg.loc["seeds", "value"]).split(",")]
    frames = {f_: pd.read_csv(os.path.join(
        PROC, "frame_%s.csv" % f_)) for f_ in "AB"}

    tasks, META = [], {}
    for fr in ("A", "B"):
        F = frames[fr]
        for co in ("EHR+ECG", "3-mod"):
            sub = (F if co == "EHR+ECG"
                   else F[F["in_three_mod"] == 1])
            sub = sub.reset_index(drop=True)
            mods = ["ehr", "ecg"] + (
                ["ctpa"] if co == "3-mod" else [])
            specs = make_specs(fr, sub, conf)
            oof = pd.read_csv(os.path.join(
                L1D, "oof_linear_%s_%s.csv"
                % (fr, CODE[co])))
            need = sorted({sp["l1col"] for sp in specs})
            per = {}
            for s in seeds:
                g = oof[oof["seed"] == s].set_index(
                    "hadm_id").reindex(
                    sub["hadm_id"].values)
                SC = {c: np.column_stack(
                    [g["%s.%s" % (m, c)].values
                     for m in mods]) for c in need}
                if any(np.isnan(v).any()
                       for v in SC.values()):
                    sys.exit("missing level-1 score in"
                             " %s %s seed %d"
                             % (fr, co, s))
                fold = g["fold"].values.astype(int)
                per[s] = (g, fold)
                tasks.append((fr, co, s, mods, SC, fold,
                              specs))
            META[(fr, co)] = (sub, mods, specs, per)
            print("  %s %-7s n=%d modalities %s fused"
                  " scores %d" % (fr, co, len(sub), mods,
                                  len(specs)), flush=True)

    # ---- timing probe: real units, run in parallel
    print("\n" + "=" * 70)
    print("TIMING PROBE (real units, %d in parallel)"
          % NJOBS)
    print("=" * 70, flush=True)
    unit_t = {}
    for key in (("3-mod", "uno"), ("3-mod", "auc"),
                ("EHR+ECG", "uno"), ("EHR+ECG", "auc")):
        tk = next(t for t in tasks
                  if t[0] == "A" and t[1] == key[0])
        sp = next(s for s in tk[6]
                  if s["crit"] == key[1])
        t0 = time.time()
        Parallel(n_jobs=NJOBS, max_nbytes=None)(
            delayed(run_task)(*tk[:6], [sp])
            for _ in range(NJOBS))
        unit_t[key] = time.time() - t0
        print("  %-8s %-4s unit: %.1f s under full load"
              % (key[0], key[1], unit_t[key]), flush=True)
    total = 0.0
    for tk in tasks:
        for sp in tk[6]:
            total += unit_t[(tk[1], sp["crit"])]
    print("  projected: about %.1f min on %d workers"
          % (total / NJOBS / 60, NJOBS), flush=True)
    if probe_only:
        keep_awake(False)
        print("\n--probe-only: stopping here.")
        return
    print("  starting in 5 s (Ctrl+C to stop)...",
          flush=True)
    time.sleep(5)

    tasks.sort(key=lambda t: -len(t[3]))
    res = Parallel(n_jobs=NJOBS, verbose=5,
                   max_nbytes=None)(
        delayed(run_task)(*t) for t in tasks)

    print("\n" + "=" * 70)
    print("ASSEMBLE AND SUMMARISE")
    print("=" * 70, flush=True)
    Wall, S1, S3 = [], [], []
    for (fr, co), (sub, mods, specs, per) in META.items():
        rr = {r[2]: r for r in res
              if r[0] == fr and r[1] == co}
        parts = []
        for s in seeds:
            df = rr[s][3].copy()
            g, fold = per[s]
            df.insert(0, "fold", fold)
            df.insert(0, "hadm_id", sub["hadm_id"].values)
            df.insert(0, "subject_id",
                      sub["subject_id"].values)
            df.insert(0, "seed", s)
            parts.append(df)
            Wall.extend(rr[s][4])
        fu = pd.concat(parts, ignore_index=True)
        nan = int(fu.filter(like="fuse.").isna()
                  .sum().sum())
        fu.to_csv(os.path.join(
            FUD, "fused_%s_%s.csv" % (fr, CODE[co])),
            index=False)
        print("  %s %-7s fused file written, missing"
              " cells %d" % (fr, co, nan), flush=True)
        folds = [per[s][1] for s in seeds]

        def fused_scores(name):
            return [rr[s][3][name].values for s in seeds]

        def uni_scores(m, col):
            return [per[s][0]["%s.%s" % (m, col)].values
                    for s in seeds]

        wdf = pd.DataFrame([w for w in Wall
                            if w["frame"] == fr
                            and w["cohort"] == co])
        if fr == "A":
            for en in ("death", "composite", "cv"):
                sp0 = next(sp for sp in specs
                           if sp["endpoint"] == en)
                y, t_, e_, tau = (sp0["label"], sp0["t"],
                                  sp0["e"], sp0["tau"])
                cands = []
                for l1 in ("cox", "lr"):
                    for crit in ("auc", "uno"):
                        nm = "fuse.%s.%s.%s" % (l1, crit,
                                                en)
                        cands.append(("fused", l1, crit,
                                      nm, fused_scores(nm)))
                    col = next(sp["l1col"] for sp in specs
                               if sp["endpoint"] == en
                               and sp["level1"] == l1)
                    for m in mods:
                        cands.append((m.upper(), l1, "-",
                                      "%s.%s" % (m, col),
                                      uni_scores(m, col)))
                for kind_, l1, crit, nm, sc in cands:
                    a_p, a_f = metric_pair(sc, folds, "auc",
                                           y=y)
                    u_p, u_f = metric_pair(sc, folds, "uno",
                                           t=t_, e=e_,
                                           tau=tau)
                    row = dict(frame=fr, cohort=co,
                               endpoint=en, score=kind_,
                               level1=l1, criterion=crit,
                               column=nm,
                               events=int(y.sum()),
                               auc30_pooled=a_p,
                               auc30_perfold=a_f,
                               unoC_pooled=u_p,
                               unoC_perfold=u_f)
                    if kind_ == "fused":
                        q = wdf[wdf["fused"] == nm]
                        for m in mods:
                            row["w_" + m] = q["w_" + m].mean()
                    S1.append(row)
        else:
            off = float(conf.loc["B", "auc_cut_offset"])
            td = sub["time_death"].values.astype(float)
            ed = sub["event_death"].values == 1
            for h in HB:
                y = ed & (td <= h + off)
                for l1 in ("cox", "lr"):
                    nm = "fuse.%s.auc%d.death" % (l1, h)
                    col = ("cox.death.lp" if l1 == "cox"
                           else "lr.death_h%d" % h)
                    a_p, a_f = metric_pair(
                        fused_scores(nm), folds, "auc", y=y)
                    q = wdf[wdf["fused"] == nm]
                    row = dict(frame=fr, cohort=co,
                               horizon=h, level1=l1,
                               events=int(y.sum()),
                               fused_auc_pooled=a_p,
                               fused_auc_perfold=a_f)
                    for m in mods:
                        row["w_%s_mean" % m] = \
                            q["w_" + m].mean()
                        row["w_%s_sd" % m] = q["w_" + m].std()
                        u_p, u_f = metric_pair(
                            uni_scores(m, col), folds, "auc",
                            y=y)
                        row["%s_auc_pooled" % m] = u_p
                        row["%s_auc_perfold" % m] = u_f
                    S3.append(row)

    wt = pd.DataFrame(Wall)
    wcols = [c for c in wt.columns if c.startswith("w_")]
    ws = wt.groupby(["frame", "cohort", "fused",
                     "level1", "criterion", "endpoint",
                     "horizon"], sort=False)[wcols].agg(
        ["mean", "std", "min", "max"])
    ws.columns = ["%s_%s" % c for c in ws.columns]
    ws = ws.reset_index()
    ws.to_csv(os.path.join(RES, "fusion_weights.csv"),
              index=False)
    s1 = pd.DataFrame(S1)
    s1.to_csv(os.path.join(RES, "fusion_s1.csv"),
              index=False)
    s3 = pd.DataFrame(S3)
    s3.to_csv(os.path.join(RES, "fusion_s3.csv"),
              index=False)
    keep_awake(False)

    print("")
    print("  S3 ECG weight by horizon (mean over 30"
          " folds):")
    print(s3.pivot_table(index=["cohort", "level1"],
                         columns="horizon",
                         values="w_ecg_mean")
          .round(3).to_string())
    print("")
    print("  elapsed %.1f min"
          % ((time.time() - t_start) / 60))
    print("")
    print("OUTPUT FILES")
    for (fr, co) in META:
        print("  " + os.path.join(
            FUD, "fused_%s_%s.csv" % (fr, CODE[co]))
            + "  (local only)")
    for f_ in ("fusion_weights.csv", "fusion_s1.csv",
               "fusion_s3.csv"):
        print("  " + os.path.join(RES, f_))
    print("  " + os.path.join(RES, "fusion_log.txt"),
          "(terminal log)")
    print("")
    print("Outputs: fusion_weights.csv, fusion_s1.csv"
          " and fusion_s3.csv")


if __name__ == "__main__":
    main()