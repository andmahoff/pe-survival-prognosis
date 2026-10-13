"""06a_foldbias.py
Check whether below-chance CTPA scores from stage 6a are a pooling
artefact.

When out-of-fold predictions from different folds are pooled into one
AUROC and the signal is weak, the AUROC is biased downward, because
training and test class proportions are negatively correlated (Parker,
Gunter and Bedo, 2007, BMC Bioinformatics 8:326). A fold with more events
in its test set was trained on fewer, so its logistic intercept (baseline
risk) is lower and its events are scored low. Their remedy is to average
the per-fold AUROCs instead of pooling.

The columns checked, in Frames A and B (three-modality cohort), are
ctpa.lr.cvc_h10 and ctpa.lr.cvu_h10, with the matching CTPA Cox score
ctpa.cox.cs_cv.lp (no intercept, so expected to be largely immune). All
are judged on CV readmission as first event by day 10. The script reports
the pooled AUROC of the seed-averaged predictions (as in stage 6a), the
mean of the per-fold AUROCs over 6 seeds x 5 folds (Parker's unbiased
version), and the correlation across the 30 test sets between a test
fold's event rate and its mean prediction (strongly negative points to the
intercept mechanism). The verdict is "below chance even per fold: SIGN" if
the mean per-fold AUROC is under 0.40, "pooling artefact" if the pooled
AUROC is more than 0.05 below it, and "no pooling gap" otherwise.

Only CTPA columns in the three-modality cohorts are read, not the
primary's cohort.

Inputs:  data/processed/frame_A.csv, frame_B.csv,
         level1/oof_linear_A_3mod.csv, level1/oof_linear_B_3mod.csv
Outputs: results/level1_foldbias.csv, results/level1_foldbias_log.txt
Usage:   python scripts/06a_foldbias.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")
L1D = os.path.join(PROC, "level1")
CUT = 10.5

CHECK = [("A", "ctpa.lr.cvc_h10"),
         ("A", "ctpa.lr.cvu_h10"),
         ("A", "ctpa.cox.cs_cv.lp"),
         ("B", "ctpa.lr.cvc_h10"),
         ("B", "ctpa.lr.cvu_h10"),
         ("B", "ctpa.cox.cs_cv.lp")]

rows = []
for fr, col in CHECK:
    F = pd.read_csv(os.path.join(
        PROC, "frame_%s.csv" % fr))
    sub = F[F["in_three_mod"] == 1].reset_index(
        drop=True)
    y = ((sub["event_cr"].values == 1)
         & (sub["time_first"].values <= CUT)
         ).astype(int)
    out = pd.read_csv(os.path.join(
        L1D, "oof_linear_%s_3mod.csv" % fr),
        usecols=["seed", "hadm_id", "fold", col])
    avg = out.groupby("hadm_id")[col].mean(
        ).reindex(sub["hadm_id"].values).values
    pooled_avg = roc_auc_score(y, avg)
    pooled_seed, fold_auc, rate, meanp = [], [], [], []
    for s, g in out.groupby("seed"):
        g = g.set_index("hadm_id").reindex(
            sub["hadm_id"].values)
        p = g[col].values
        f = g["fold"].values
        pooled_seed.append(roc_auc_score(y, p))
        for k in np.unique(f):
            m = f == k
            if 0 < y[m].sum() < m.sum():
                fold_auc.append(roc_auc_score(y[m],
                                              p[m]))
            rate.append(y[m].mean())
            meanp.append(p[m].mean())
    r = float(np.corrcoef(rate, meanp)[0, 1])
    pf = float(np.mean(fold_auc))
    if pf < 0.40:
        verdict = "below chance even per fold: SIGN"
    elif pooled_avg < pf - 0.05:
        verdict = "pooling artefact"
    else:
        verdict = "no pooling gap"
    rows.append({
        "frame": fr, "cohort": "3-mod",
        "column": col, "events_by_day10": int(y.sum()),
        "pooled_auc_seedavg_pred": round(pooled_avg, 4),
        "pooled_auc_mean_over_seeds":
        round(float(np.mean(pooled_seed)), 4),
        "perfold_auc_mean": round(pf, 4),
        "perfold_auc_sd": round(float(np.std(
            fold_auc)), 4),
        "test_sets_used": len(fold_auc),
        "corr_eventrate_vs_meanpred": round(r, 3),
        "verdict": verdict})
    print("  %s %-20s pooled %.3f  per-fold %.3f"
          "  r(rate, mean pred) %+.2f  %s"
          % (fr, col, pooled_avg, pf, r, verdict),
          flush=True)

res = pd.DataFrame(rows)
fp = os.path.join(RES, "level1_foldbias.csv")
res.to_csv(fp, index=False)
print("")
print("OUTPUT FILES")
print("  " + fp)
print("  " + os.path.join(RES,
                          "level1_foldbias_log.txt"),
      "(terminal log)")
print("")
print("Outputs: level1_foldbias.csv")