"""06a_diag.py
Diagnostics for the stage 6a level-1 linear scores, following up
level1_linear_checks.csv.

Part 1 asks which scores rank no better than chance, and whether that is a
sign error or a modality with no signal for a small endpoint. For every
seed-averaged score column at or below 0.5 (Harrell's C for Cox linear
predictors, AUROC for logistic columns), it reports target, events, the
value with a 95% patient-level bootstrap interval (1,000 resamples, seed
42) and the range across seeds. A sign error gives values clearly below
0.5 (interval excluding 0.5); no signal gives values near 0.5 (interval
spanning it). Only these columns are written, so the primary analysis
stays blind as long as its own cohort (Frame B EHR+ECG) has none.

Part 2 asks where the stage 6a fitting warnings came from and whether the
final fits converged. It refits a representative sample on the real data
(seed 42, outer fold 0, Frames A and B, EHR+ECG cohort, EHR and ECG
modalities, every target and penalty value), recording warning types,
convergence and seconds per fit. The fit chosen by tuning in that fold is
marked final = 1.

Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
         feature_whitelist.csv, level1/params_linear.csv,
         level1/oof_linear_*.csv; data/inputs/ehr28.csv, ecg_opt.csv,
         ctpa46.csv; target and grid definitions imported from
         scripts/06a_level1_linear.py
Outputs: results/level1_diag_scores.csv, results/level1_diag_fits.csv,
         results/level1_diag_log.txt
Usage:   python scripts/06a_diag.py [--probe-only]
Notes:   --probe-only times one real ECG Cox fit and one real EHR Cox fit,
         projects the parallel run (workers = physical cores - 1, at most
         10) and stops; otherwise the fits start after a 5 s pause. On
         Windows the machine is kept awake while the script runs.
"""
import os
import sys
for v_ in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
import time
import ctypes
import warnings
import importlib.util
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
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
L1D = os.path.join(PROC, "level1")
NBOOT = 1000


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
    es_cont, es_sys = 0x80000000, 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(
        es_cont | (es_sys if on else 0))


def load_l1():
    p = os.path.join(P2, "scripts",
                     "06a_level1_linear.py")
    spec = importlib.util.spec_from_file_location(
        "l1_defs", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fit_one(kind, a, t, e, y, param):
    t0 = time.time()
    err = ""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            if kind == "cox":
                CoxPHSurvivalAnalysis(
                    alpha=param, ties="efron").fit(
                    a, Surv.from_arrays(event=e,
                                        time=t))
            else:
                LogisticRegression(
                    C=param, max_iter=8000).fit(a, y)
        except Exception as ex:
            err = str(ex)[:120]
    msgs = sorted(set("%s: %s" % (
        wi.category.__name__,
        str(wi.message).replace("\n", " ")[:100])
        for wi in w))
    conv = not any("Convergence" in m_ for m_ in msgs)
    return (time.time() - t0, len(w), " | ".join(msgs),
            conv, err)


def boot_ci(metric, idx_groups, rng):
    vals = []
    for _ in range(NBOOT):
        pick = rng.integers(0, len(idx_groups),
                            len(idx_groups))
        ii = np.concatenate([idx_groups[j]
                             for j in pick])
        try:
            vals.append(metric(ii))
        except Exception:
            pass
    return (np.percentile(vals, 2.5),
            np.percentile(vals, 97.5))


def main():
    t_start = time.time()
    probe_only = "--probe-only" in sys.argv
    keep_awake(True)
    l1 = load_l1()
    fa = pd.read_csv(os.path.join(PROC, "frame_A.csv"))
    fb = pd.read_csv(os.path.join(PROC, "frame_B.csv"))
    folds = pd.read_csv(os.path.join(PROC,
                                     "folds.csv"))
    wl = pd.read_csv(os.path.join(
        PROC, "feature_whitelist.csv"))
    prm = pd.read_csv(os.path.join(
        L1D, "params_linear.csv"))
    mats = {m: pd.read_csv(os.path.join(INP, f_))
            .drop_duplicates("hadm_id")
            for m, f_ in l1.MATS.items()}
    frames = {"A": fa, "B": fb}

    # ---------- PART 1: scores at or below 0.5 ----
    print("=" * 70)
    print("PART 1: scores ranking no better than chance")
    print("=" * 70, flush=True)
    rng = np.random.default_rng(42)
    rows1 = []
    for fr in ("A", "B"):
        for co, code in (("EHR+ECG", "ehrecg"),
                         ("3-mod", "3mod")):
            F = frames[fr]
            sub = (F if co == "EHR+ECG"
                   else F[F["in_three_mod"] == 1])
            sub = sub.reset_index(drop=True)
            ST = l1.surv_targets(sub)
            BT = l1.bin_targets(sub, l1.HZ[fr])
            out = pd.read_csv(os.path.join(
                L1D, "oof_linear_%s_%s.csv"
                % (fr, code)))
            pc = [c for c in out.columns if "." in c]
            avg = out.groupby("hadm_id", sort=False)[
                pc].mean().reindex(sub["hadm_id"].values)
            per = {s: g_.set_index("hadm_id")[pc]
                   .reindex(sub["hadm_id"].values)
                   for s, g_ in out.groupby("seed")}
            grp_idx = [np.where(sub["subject_id"].values
                                == s_)[0] for s_ in
                       pd.unique(sub["subject_id"])]
            for c in pc:
                mod, lrn, tg = c.split(".")[:3]
                if lrn == "cox" and not c.endswith(".lp"):
                    continue
                if lrn == "cox":
                    t_, e_ = ST[tg]

                    def met(ii, v=None, t_=t_, e_=e_):
                        return concordance_index_censored(
                            e_[ii], t_[ii], v[ii])[0]
                    ev = int(e_.sum())
                else:
                    y_ = BT[tg][0]

                    def met(ii, v=None, y_=y_):
                        return roc_auc_score(y_[ii], v[ii])
                    ev = int(y_.sum())
                v = avg[c].values
                allidx = np.arange(len(v))
                val = met(allidx, v)
                if val > 0.5:
                    continue
                lo, hi = boot_ci(
                    lambda ii: met(ii, v), grp_idx, rng)
                seeds_v = [met(allidx, per[s][c].values)
                           for s in sorted(per)]
                verdict = ("BELOW chance: possible sign"
                           " error" if hi < 0.5 else
                           "consistent with no signal")
                rows1.append({
                    "frame": fr, "cohort": co,
                    "column": c, "target": tg,
                    "n": len(v), "events": ev,
                    "value": round(val, 4),
                    "ci_lo": round(lo, 4),
                    "ci_hi": round(hi, 4),
                    "seed_min": round(min(seeds_v), 4),
                    "seed_max": round(max(seeds_v), 4),
                    "verdict": verdict})
                print("  %s %-7s %-26s ev=%-4d %.4f"
                      " [%.4f, %.4f]  %s"
                      % (fr, co, c, ev, val, lo, hi,
                         verdict), flush=True)
    s1 = pd.DataFrame(rows1)
    s1.to_csv(os.path.join(RES,
                           "level1_diag_scores.csv"),
              index=False)
    print("  columns at or below 0.5:", len(s1))

    # ---------- PART 2: warnings and convergence --
    print("")
    print("=" * 70)
    print("PART 2: warnings and convergence (seed 42,"
          " fold 0)")
    print("=" * 70, flush=True)
    jobs = []
    for fr in ("A", "B"):
        sub = frames[fr].reset_index(drop=True)
        ST = l1.surv_targets(sub)
        BT = l1.bin_targets(sub, l1.HZ[fr])
        fo = folds[(folds["frame"] == fr)
                   & (folds["cohort"] == "EHR+ECG")
                   & (folds["seed"] == 42)]
        fold = sub[["hadm_id"]].merge(
            fo[["hadm_id", "fold"]], on="hadm_id",
            how="left")["fold"].values
        tr = np.where(fold != 0)[0]
        for mod in ("ehr", "ecg"):
            cc = wl.loc[wl["modality"] == mod,
                        "column"].tolist()
            X = sub[["hadm_id"]].merge(
                mats[mod][["hadm_id"] + cc],
                on="hadm_id", how="left")[cc]
            X = X.values.astype(float)[tr]
            a = SimpleImputer(
                strategy="median").fit_transform(X)
            a = StandardScaler().fit_transform(a)
            q = prm[(prm["frame"] == fr)
                    & (prm["cohort"] == "EHR+ECG")
                    & (prm["modality"] == mod)
                    & (prm["seed"] == 42)
                    & (prm["fold"] == 0)]
            chosen = dict(zip(zip(q["learner"],
                                  q["target"]),
                              q["param"]))
            for tg in l1.SURV:
                t_, e_ = ST[tg]
                for al in l1.ALPHAS:
                    jobs.append(dict(
                        frame=fr, modality=mod,
                        learner="cox", target=tg,
                        param=al, final=int(
                            chosen.get(("cox", tg))
                            == al),
                        args=("cox", a, t_[tr], e_[tr],
                              None, al)))
            for tg, (y_, keep) in BT.items():
                sel = keep[tr]
                for c in l1.CS:
                    jobs.append(dict(
                        frame=fr, modality=mod,
                        learner="lr", target=tg,
                        param=c, final=int(
                            chosen.get(("lr", tg)) == c),
                        args=("lr", a[sel], None, None,
                              y_[tr][sel].astype(int),
                              c)))

    # timing probe on real data
    ecg_job = next(j for j in jobs
                   if j["modality"] == "ecg"
                   and j["learner"] == "cox")
    ehr_job = next(j for j in jobs
                   if j["modality"] == "ehr"
                   and j["learner"] == "cox")
    t_ecg = fit_one(*ecg_job["args"])[0]
    t_ehr = fit_one(*ehr_job["args"])[0]
    n_ecg = sum(1 for j in jobs if j["modality"]
                == "ecg" and j["learner"] == "cox")
    n_ehr = sum(1 for j in jobs if j["modality"]
                == "ehr" and j["learner"] == "cox")
    serial = n_ecg * t_ecg + n_ehr * t_ehr + 0.05 * len(
        [j for j in jobs if j["learner"] == "lr"])
    print("  REAL-data probe: ECG Cox %.1f s, EHR Cox"
          " %.1f s per fit" % (t_ecg, t_ehr))
    print("  %d fits; serial about %.1f min, parallel"
          " (%d workers = physical cores - 1) about"
          " %.1f min" % (len(jobs), serial / 60, NJOBS,
                         serial / 60 / NJOBS * 1.3))
    if probe_only:
        keep_awake(False)
        print("\n--probe-only: stopping here.")
        return
    print("  starting in 5 s (Ctrl+C to stop)...",
          flush=True)
    time.sleep(5)
    jobs.sort(key=lambda j: -(j["modality"] == "ecg"
                              and j["learner"] == "cox"))
    res = Parallel(n_jobs=NJOBS, verbose=5)(
        delayed(fit_one)(*j["args"]) for j in jobs)

    rows2 = []
    for j, (sec, nw, msg, conv, err) in zip(jobs, res):
        rows2.append({k: j[k] for k in (
            "frame", "modality", "learner", "target",
            "param", "final")} | {
            "seconds": round(sec, 2),
            "n_warnings": nw, "converged": int(conv),
            "warnings": msg, "error": err})
    s2 = pd.DataFrame(rows2)
    s2.to_csv(os.path.join(RES, "level1_diag_fits.csv"),
              index=False)
    print("")
    print(s2.groupby(["learner", "modality"]).agg(
        fits=("seconds", "size"),
        mean_s=("seconds", "mean"),
        with_warnings=("n_warnings",
                       lambda x: (x > 0).mean()),
        converged=("converged", "mean")).round(3)
        .to_string())
    fin = s2[s2["final"] == 1]
    print("")
    print("  FINAL fits: %d, converged %d, with errors %d"
          % (len(fin), int(fin["converged"].sum()),
             int((fin["error"] != "").sum())))
    print("  warning types seen:")
    for m_ in sorted(set(" | ".join(
            s2["warnings"]).split(" | ")) - {""}):
        print("    " + m_)
    keep_awake(False)
    print("")
    print("  elapsed %.1f min"
          % ((time.time() - t_start) / 60))
    print("")
    print("OUTPUT FILES")
    print("  " + os.path.join(RES,
                              "level1_diag_scores.csv"))
    print("  " + os.path.join(RES,
                              "level1_diag_fits.csv"))
    print("  " + os.path.join(RES,
                              "level1_diag_log.txt"),
          "(terminal log)")
    print("")
    print("Outputs: level1_diag_scores.csv and"
          " level1_diag_fits.csv")


if __name__ == "__main__":
    main()