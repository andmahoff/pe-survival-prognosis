"""13e_verified_results.py
Reference file of verified results (PREREGISTRATION_ADDENDUM_2.md A38).

Every headline number of the study is read by code straight from the saved
aggregate CSVs and tagged with its source file and the rows it came from,
so that the write-up can quote this file. Rerun it whenever a results file
changes. Nothing is refitted and no patient-level file is read.

The reference file has 13 sections:
  1  Cohorts and events (stages 4r, 4p, 13a)
  2  The registered primary analysis, H1 (stage 7) and its sensitivity
     analyses (stages 7b, 7c, 9a, 19a); Frame B
  3  Final models: discrimination and calibration (Tables 2 and S2)
  4  Added value over clinical scores and the EHR alone (Tables 3, S3),
     the scores' own AUROC, and Frame A v4 against the at-diagnosis model
  5  Decision curves (stage 13b)
  6  Survival against matched binary models (stage 18s)
  7  Temporal validation (stage 21b)
  8  Model choices and hospice discharges (stages 18o, 18p, 18r)
  9  At-diagnosis sensitivity analyses: strict history, vital signs
  10 Interpretability at diagnosis (stage 20p)
  11 CTPA blocks at diagnosis (stage 16p)
  12 Subgroups at discharge (stage 22a)
  13 Checks across the final stages

Counts from 1 to 9 are written as "<10" everywhere. Checks: core inputs
found (FAIL) and optional ones (WARN, section skipped); death events by 30
and 365 days agree between stage 4r and Table 1 (FAIL); no FAIL in any
final stage's checks file (FAIL).

Inputs:  results/*.csv (aggregates only; the files read are listed in the
         last table of the output)
Outputs: results/verified_results.md (the reference file),
         results/verified_results.csv (every number, long format),
         results/verified_checks.csv, results/verified_log.txt (terminal
         log)
Usage:   python scripts/13e_verified_results.py
"""
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

P2 = os.path.join(os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd())), "SurvivalProject")
RS = os.path.join(P2, "results")
HZ = [10, 30, 90, 180, 365]
CHK = []
ROWS = []
MD = []
USED = {}

CORE = ["optionB_outcome_events.csv", "table1.csv", "primary.csv",
        "tables23_long.csv", "slP_summary.csv", "slP_metrics.csv",
        "fin_dca_key.csv", "fin_metrics.csv"]
OPTIONAL = ["presentation_summary.csv", "primary_robust.csv",
            "primary_hospice.csv", "s5_primary.csv", "ab_ecg_ratio.csv",
            "binsurv_summary.csv", "tv2_P_summary.csv", "tv2_B_summary.csv",
            "tv2_P_contrasts.csv", "tv2_B_contrasts.csv",
            "tv2_P_drift.csv", "tv2_B_drift.csv", "sl7_decision.csv",
            "optB2_decision.csv", "optB2_core_check.csv",
            "sl7nh_summary.csv", "optB2_sl_summary.csv",
            "interpP_spline_terms.csv", "interpP_shap_importance.csv",
            "interpP_shap_domains.csv", "interpP_agreement.csv",
            "interpP_fig_twins.csv", "ctpa_blocks_P_time.csv",
            "ctpa_blocks_P_contrasts.csv", "subgroup_metrics.csv"]
CHECKS = ["presentation_checks.csv", "optionB_outcome_checks.csv",
          "table1_checks.csv", "slP_checks.csv", "optB2_sl_checks.csv",
          "sl7_checks.csv", "binsurv_checks.csv", "tv2_P_checks.csv",
          "tv2_B_checks.csv", "interpP_checks.csv",
          "ctpa_blocks_P_checks.csv", "subgroup_checks.csv",
          "fin_checks.csv", "tables23_checks.csv"]


def chk(item, value, status):
    CHK.append({"item": item, "value": str(value), "status": status})
    print("  %-5s %-62s %s" % (status, item, value), flush=True)


def stop_if_fail():
    ck = pd.DataFrame(CHK)
    if (ck.status == "FAIL").any():
        ck.to_csv(os.path.join(RS, "verified_checks.csv"), index=False)
        print("\nSTOP: a check failed. see: verified_checks.csv")
        sys.exit(1)


def rd(f):
    p = os.path.join(RS, f)
    if not os.path.exists(p):
        return None
    if f not in USED:
        USED[f] = datetime.fromtimestamp(os.path.getmtime(p))
    return pd.read_csv(p)


# ---------------- formatting ----------------
def num(v, d=3):
    if v is None or not np.isfinite(v):
        return "n/a"
    s = ("{:.%df}" % d).format(v)
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    return s.replace("-", "\u2212")


def ci(e, lo, hi, d=3):
    if lo is None or not np.isfinite(lo):
        return num(e, d)
    return "%s (%s to %s)" % (num(e, d), num(lo, d), num(hi, d))


def cnt(n):
    if n is None or not np.isfinite(n):
        return "n/a"
    n = int(round(n))
    return "<10" if 1 <= n <= 9 else "{:,}".format(n)


def pct(v, d=1):
    return ("{:.%df}%%" % d).format(100 * v)


def rec(sec, item, est, lo=np.nan, hi=np.nan, d=3, src="", filt="",
        horizon=None, text=None):
    """record one number for the CSV and return its display text."""
    t = text if text is not None else ci(est, lo, hi, d)
    ROWS.append({"section": sec, "item": item, "horizon": horizon,
                 "estimate": est, "lo": lo, "hi": hi, "text": t,
                 "source_file": src, "source_rows": filt})
    return t


def h(s, level=2):
    MD.append("")
    MD.append("#" * level + " " + s)
    MD.append("")


def p(s):
    MD.append(s)


def table(head, rows):
    MD.append("| " + " | ".join(head) + " |")
    MD.append("|" + "|".join(["---"] * len(head)) + "|")
    for r in rows:
        MD.append("| " + " | ".join(str(x) for x in r) + " |")
    MD.append("")


def src(f, filt):
    p("*Source: `results/%s` (%s).*" % (f, filt))
    p("")


# clearer names for MIMIC items (as stage 20q)
RELAB = {
    "lab_50820_ph_first": "Blood pH", "lab_51491_ph_first": "Urine pH",
    "lab_50861_alanine_aminotransferase_alt_first": "ALT",
    "lab_50878_asparate_aminotransferase_as_first": "AST",
    "lab_50911_creatine_kinase_mb_isoenzyme_first": "CK-MB",
    "lab_50910_creatine_kinase_ck_first": "Creatine kinase",
    "lab_50954_lactate_dehydrogenase_ld_first": "LDH",
    "lab_51516_wbc_first": "Urine white cells",
    "lab_51493_rbc_first": "Urine red cells",
    "lab_51498_specific_gravity_first": "Urine specific gravity",
    "lab_51476_epithelial_cells_first": "Urine epithelial cells",
    "lab_51082_creatinine_urine_first": "Urine creatinine",
    "lab_51492_protein_first": "Urine protein",
    "lab_52172_rdw_sd_first": "RDW-SD",
    "lab_51214_fibrinogen_functional_first": "Fibrinogen",
    "lab_50804_calculated_total_co2_first": "Total CO2 (blood gas)",
    "lab_50808_free_calcium_first": "Ionised calcium",
    "lab_50953_iron_binding_capacity_total_first":
        "Total iron binding capacity"}


def lab_of(feature, label):
    return RELAB.get(feature, label)


def tidy(s):
    s = s.replace("auc: ", "AUROC: ").replace("ipa: ", "IPA: ")
    s = s.replace("Super learner, recalibrated", "SL")
    s = s.replace("Super learner-EHR, recalibrated", "SL-EHR")
    s = s.replace(", recalibrated", "").replace(" - ", " minus ")
    return s


HH = ["%d d" % x for x in HZ]
FRL = {"P": "At diagnosis", "B": "At discharge"}
OCL = {"death": "Death", "comp": "Death or CV readmission",
       "cv": "CV readmission (death competing)"}

# ---------------- inputs ----------------
print("=== INPUTS ===")
miss = [f for f in CORE if not os.path.exists(os.path.join(RS, f))]
chk("core inputs found", miss[0] if miss else "all %d" % len(CORE),
    "FAIL" if miss else "PASS")
stop_if_fail()
mo = [f for f in OPTIONAL if not os.path.exists(os.path.join(RS, f))]
chk("optional inputs found (sections skipped if missing)",
    "all %d" % len(OPTIONAL) if not mo else "missing: " + ", ".join(mo),
    "PASS" if not mo else "WARN")

# ================= 1. cohorts and events =================
print("\n=== BUILD ===")
S = "1 Cohorts and events"
h("1. Cohorts and events")
EV = rd("optionB_outcome_events.csv")
T1 = rd("table1.csv")
first = T1.columns[1]
pat = T1[T1.iloc[:, 0] == "Patients, n"][first]
npat = pat.iloc[0] if len(pat) else "n/a"
rows = []
for fr, fl in (("P", "P (at diagnosis)"), ("B", "B (at discharge)")):
    for co in ("EHR+ECG", "3-mod"):
        x = EV[(EV.frame == fl) & (EV.cohort == co)]
        if not len(x):
            continue
        n = int(x.patients.iloc[0])
        unit = "admissions" if fr == "P" else "discharges"
        rows.append([FRL[fr], co, "%s %s" % (cnt(n), unit)])
        rec(S, "%s %s: %s" % (FRL[fr], co, unit), n, d=0,
            src="optionB_outcome_events.csv",
            filt="frame %s, cohort %s" % (fl, co), text=cnt(n))
table(["Pipeline", "Cohort", "Size"], rows)
p("At diagnosis, the %s admissions come from %s patients (Table 1)."
  % (cnt(EV[(EV.frame == "P (at diagnosis)") &
            (EV.cohort == "EHR+ECG")].patients.iloc[0]), npat))
rec(S, "At diagnosis: patients", np.nan, src="table1.csv",
    filt="row 'Patients, n', first column", text=str(npat))
src("optionB_outcome_events.csv", "frame, cohort, horizon")
for fr, fl in (("P", "P (at diagnosis)"), ("B", "B (at discharge)")):
    for co in ("EHR+ECG", "3-mod"):
        x = EV[(EV.frame == fl) & (EV.cohort == co)].set_index("horizon")
        if not len(x):
            continue
        p("**%s, %s cohort: events by each horizon, n (%%)**" %
          (FRL[fr], co))
        p("")
        rows = []
        for oc, col in (("death", "death"), ("comp", "composite"),
                        ("cv", "cv_first")):
            r = [OCL[oc]]
            for hz in HZ:
                e = x.loc[hz, col]
                n = x.loc[hz, "patients"]
                t = "%s (%s)" % (cnt(e), pct(e / n)) if e >= 10 or e == 0 \
                    else cnt(e)
                r.append(rec(S, "%s %s %s events" % (FRL[fr], co, OCL[oc]),
                             e, d=0, src="optionB_outcome_events.csv",
                             filt="frame %s, cohort %s, column %s" %
                             (fl, co, col), horizon=hz, text=t))
            rows.append(r)
        table(["Outcome"] + HH, rows)
PS = rd("presentation_summary.csv")
if PS is not None:
    p("**Presentation data (at diagnosis)**")
    p("")
    rows = []
    for _, r in PS.iterrows():
        rows.append([r["block"], r["item"], r["value"]])
        rec(S, "%s: %s" % (r["block"], r["item"]), np.nan,
            src="presentation_summary.csv", filt="item '%s'" % r["item"],
            text=str(r["value"]))
    table(["Block", "Item", "Value"], rows)
    src("presentation_summary.csv", "all rows")
# events agree with Table 1
t30 = T1[T1.iloc[:, 0] == "Death by 30 days, n (%)"][first]
t365 = T1[T1.iloc[:, 0] == "Death by 365 days, n (%)"][first]
x = EV[(EV.frame == "P (at diagnosis)") & (EV.cohort == "EHR+ECG")]
x = x.set_index("horizon")
ok = (len(t30) and len(t365) and
      str(t30.iloc[0]).startswith(str(x.loc[30, "death"]) + " ") and
      str(t365.iloc[0]).startswith(str(x.loc[365, "death"]) + " "))
chk("deaths by 30 and 365 days: stage 4r agrees with Table 1",
    "%s / %s" % (x.loc[30, "death"], x.loc[365, "death"]),
    "PASS" if ok else "FAIL")
stop_if_fail()

# ================= 2. registered primary analysis =================
S = "2 Registered primary analysis"
h("2. Registered primary analysis (H1), Frame B, EHR+ECG cohort")
PR = rd("primary.csv")
p("Deaths are split at day 30 after discharge. The estimand is the ECG "
  "hazard ratio in days 0 to 30 over that in days 31 to 365.")
p("")
rows = []
for q in ("patients", "deaths days 0-30", "deaths days 31-365"):
    v = PR[PR.quantity == q].estimate
    if len(v):
        rows.append([q, rec(S, q, v.iloc[0], d=0, src="primary.csv",
                            filt="quantity '%s'" % q, text=cnt(v.iloc[0]))])
for q in ("HR per SD, ECG, days 0-30", "HR per SD, ECG, days 31-365",
          "HR per SD, EHR, days 0-30", "HR per SD, EHR, days 31-365"):
    x = PR[(PR.quantity == q) & (PR.method == "bootstrap percentile")]
    if len(x):
        r = x.iloc[0]
        rows.append([q, rec(S, q, r.estimate, r.lo, r.hi, 2, "primary.csv",
                            "quantity '%s', bootstrap percentile" % q)])
x = PR[PR.quantity.str.startswith("PRIMARY: ECG HR ratio")]
if len(x):
    r = x.iloc[0]
    rows.append(["ECG ratio, early / late (primary)",
                 rec(S, "ECG ratio, early / late", r.estimate, r.lo, r.hi,
                     2, "primary.csv", "quantity 'PRIMARY: ...'")])
    rows.append(["Decision", rec(S, "H1 decision", np.nan,
                                 src="primary.csv", filt="note",
                                 text=str(r.note))])
x = PR[PR.quantity.str.startswith("EHR HR ratio")]
if len(x):
    r = x.iloc[0]
    rows.append(["EHR ratio, early / late (no decision rule)",
                 rec(S, "EHR ratio, early / late", r.estimate, r.lo, r.hi,
                     2, "primary.csv", "quantity 'EHR HR ratio ...'")])
table(["Quantity", "Value (95% bootstrap interval)"], rows)
src("primary.csv", "bootstrap percentile rows")

p("**The ECG's early / late ratio across sensitivity analyses**")
p("")
rows = [["Registered primary analysis", "primary.csv",
         ci(PR[PR.quantity.str.startswith("PRIMARY")].estimate.iloc[0],
            PR[PR.quantity.str.startswith("PRIMARY")].lo.iloc[0],
            PR[PR.quantity.str.startswith("PRIMARY")].hi.iloc[0], 2)]]
SENS = []
RB = rd("primary_robust.csv")
if RB is not None:
    for v, lab in (("winsorised", "Scores winsorised"),
                   ("rank_normal", "Scores rank-normalised")):
        x = RB[(RB.variant == v) & RB.quantity.str.startswith("ECG HR")]
        SENS.append((lab, x, "primary_robust.csv", "variant %s" % v))
HP = rd("primary_hospice.csv")
if HP is not None:
    x = HP[HP.quantity.str.startswith("SENSITIVITY: ECG HR ratio")]
    SENS.append(("Hospice discharges excluded", x, "primary_hospice.csv",
                 "quantity 'SENSITIVITY: ECG HR ratio ...'"))
S5 = rd("s5_primary.csv")
if S5 is not None:
    for m, lab in (("(a)", "Three-modality cohort"),
                   ("(b)", "Three-modality cohort, adjusted for CTPA")):
        x = S5[S5.model.str.startswith(m) &
               S5.quantity.str.startswith("ECG HR ratio")]
        SENS.append((lab, x, "s5_primary.csv", "model %s" % m))
AB = rd("ab_ecg_ratio.csv")
if AB is not None:
    for v, lab in (("EHR 28", "Retuned EHR score, 28 features"),
                   ("EHR 79", "Larger EHR score, 79 features"),
                   ("EHR 79 + MI", "79 features plus missing indicators")):
        x = AB[(AB.ehr_version == v) & (AB.quantity == "ECG ratio")]
        SENS.append((lab, x, "ab_ecg_ratio.csv", "ehr_version %s" % v))
for lab, x, f, filt in SENS:
    if len(x) == 1:
        r = x.iloc[0]
        rows.append([lab, f, rec(S, "ECG ratio: " + lab, r.estimate, r.lo,
                                 r.hi, 2, f, filt)])
table(["Analysis", "Source", "ECG ratio (95% interval)"], rows)

# ================= 3. final models =================
S = "3 Final models"
h("3. Final models: discrimination and calibration (Tables 2 and S2)")
p("At diagnosis, death uses Super learner P. At discharge, death uses "
  "Super learner v7. The composite and CV readmission use the stage 18r "
  "super learners. ICI has no interval (stage 13c note).")
TL = rd("tables23_long.csv")
FM = rd("fin_metrics.csv")
OCT = {"Death": "death", "Death or CV readmission": "comp",
       "CV readmission, with death as a competing event": "cv"}
for tab, co in (("2", "EHR+ECG"), ("S2", "3-mod")):
    h("Table %s: %s cohort" % (tab, co), 3)
    for pl in ("At diagnosis", "At discharge"):
        for oc in OCT:
            x = TL[(TL.table == tab) & (TL.pipeline == pl) &
                   (TL.outcome == oc)]
            if not len(x):
                continue
            p("**%s: %s**" % (pl, oc))
            p("")
            rows = []
            for rw in x.row.unique():
                r = [rw]
                for hz in HZ:
                    y = x[(x.row == rw) & (x.horizon == hz)]
                    if not len(y):
                        r.append("")
                        continue
                    y = y.iloc[0]
                    d = 0 if rw == "Events" else (
                        2 if rw.startswith(("Observed", "Calibration"))
                        else 3)
                    t = cnt(y.estimate) if rw == "Events" else None
                    r.append(rec(S, "Table %s, %s, %s: %s" %
                                 (tab, pl, oc, rw), y.estimate, y.lo, y.hi,
                                 d, "tables23_long.csv",
                                 "table %s; %s" % (tab, y.source),
                                 horizon=hz, text=t))
                rows.append(r)
            fr = "P" if pl == "At diagnosis" else "B"
            r = ["Calibration intercept (95% CI)"]
            for hz in HZ:
                y = FM[(FM.frame == fr) & (FM.outcome == OCT[oc]) &
                       (FM.cohort == co) & (FM.horizon == hz) &
                       (FM.measure == "citl")]
                r.append(rec(S, "%s, %s, %s: calibration intercept" %
                             (pl, oc, co), y.estimate.iloc[0],
                             y.lo.iloc[0], y.hi.iloc[0], 2,
                             "fin_metrics.csv",
                             "frame %s, outcome %s, cohort %s, citl" %
                             (fr, OCT[oc], co), horizon=hz)
                         if len(y) else "")
            rows.append(r)
            table(["Measure"] + HH, rows)
    src("tables23_long.csv", "table %s; original stage named in the CSV. "
        "Calibration intercept: `results/fin_metrics.csv`, measure citl"
        % tab)

# ================= 4. added value =================
S = "4 Added value"
h("4. Added value over clinical scores and the EHR alone (Tables 3, S3)")
p("Each cell is the final model minus the comparator. At discharge, the "
  "death rows use v7-core (A33).")
for tab, co in (("3", "EHR+ECG"), ("S3", "3-mod")):
    h("Table %s: %s cohort" % (tab, co), 3)
    for pl in ("At diagnosis", "At discharge"):
        for oc in OCT:
            x = TL[(TL.table == tab) & (TL.pipeline == pl) &
                   (TL.outcome == oc)]
            if not len(x):
                continue
            p("**%s: %s**" % (pl, oc))
            p("")
            rows = []
            for rw in x.row.unique():
                r = [rw]
                for hz in HZ:
                    y = x[(x.row == rw) & (x.horizon == hz)].iloc[0]
                    t = ci(y.estimate, y.lo, y.hi, 3)
                    if y.excludes_zero is True or y.excludes_zero == "True":
                        t = "**" + t + "**"
                    r.append(rec(S, "Table %s, %s, %s: %s" %
                                 (tab, pl, oc, rw), y.estimate, y.lo, y.hi,
                                 3, "tables23_long.csv",
                                 "table %s; %s" % (tab, y.source),
                                 horizon=hz, text=t))
                rows.append(r)
            table(["Comparison"] + HH, rows)
    p("Bold: the 95% interval excludes zero.")
    src("tables23_long.csv", "table %s" % tab)
SM = rd("slP_summary.csv")
h("The scores' and the EHR-alone model's own AUROC at diagnosis (death)", 3)
for co in ("EHR+ECG", "3-mod"):
    x = SM[(SM.cohort == co) & SM["item"].str.startswith("AUROC")]
    p("**%s cohort**" % co)
    p("")
    rows = []
    for it in x["item"].unique():
        r = [it.replace("AUROC: ", "")]
        for hz in HZ:
            y = x[(x["item"] == it) & (x.horizon == hz)]
            r.append(rec(S, "%s: %s" % (co, it), y.estimate.iloc[0],
                         y.lo.iloc[0], y.hi.iloc[0], 3, "slP_summary.csv",
                         "cohort %s, item '%s'" % (co, it), horizon=hz)
                     if len(y) else "")
        rows.append(r)
    table(["Model or score (AUROC)"] + HH, rows)
    x = SM[(SM.cohort == co) &
           SM["item"].str.contains("Frame A v4, recalibrated - ")]
    if len(x):
        rows = []
        for it in x["item"].unique():
            r = [tidy(it)]
            for hz in (10, 30):
                y = x[(x["item"] == it) & (x.horizon == hz)]
                r.append(rec(S, "%s: %s" % (co, it), y.estimate.iloc[0],
                             y.lo.iloc[0], y.hi.iloc[0], 3,
                             "slP_summary.csv",
                             "cohort %s, item '%s'" % (co, it), horizon=hz)
                         if len(y) else "")
            rows.append(r)
        p("Frame A's mixed-time model (v4) minus the at-diagnosis model:")
        p("")
        table(["Contrast", "10 d", "30 d"], rows)
src("slP_summary.csv", "evaluation 'all'")

# ================= 5. decision curves =================
S = "5 Decision curves"
h("5. Decision curves (EHR+ECG cohort; exploratory)")
p("Net benefit counts correctly flagged patients per patient, minus a "
  "penalty for false alarms. Cells are the model's net benefit minus "
  "each strategy's, with 95% intervals.")
KY = rd("fin_dca_key.csv")
for fr in ("P", "B"):
    p("")
    p("**%s**" % FRL[fr])
    p("")
    rows = []
    stgs = ["Treat all"] + (["sPESI 1 or more", "PESI-10 class III to V"]
                            if fr == "P" else [])
    for oc in ("death", "comp", "cv"):
        for hz in (30, 365):
            for th in (0.05, 0.10, 0.20):
                x = KY[(KY.frame == fr) & (KY.outcome == oc) &
                       (KY.cohort == "EHR+ECG") & (KY.horizon == hz) &
                       (np.isclose(KY.threshold, th))]
                if not len(x):
                    continue
                r = [OCL[oc], "%d d" % hz, pct(th, 0)]
                m = x[(x.quantity == "net benefit") &
                      (x.strategy == "Model")]
                r.append(rec(S, "%s %s %dd NB model" % (FRL[fr], oc, hz),
                             m.estimate.iloc[0], m.lo.iloc[0],
                             m.hi.iloc[0], 3, "fin_dca_key.csv",
                             "frame %s, outcome %s, threshold %.2f, "
                             "net benefit, Model" % (fr, oc, th),
                             horizon=hz))
                for stg in stgs:
                    y = x[(x.quantity == "model minus strategy") &
                          (x.strategy == stg)]
                    r.append(rec(S, "%s %s %dd model minus %s" %
                                 (FRL[fr], oc, hz, stg),
                                 y.estimate.iloc[0], y.lo.iloc[0],
                                 y.hi.iloc[0], 3, "fin_dca_key.csv",
                                 "frame %s, outcome %s, threshold %.2f, "
                                 "model minus %s" % (fr, oc, th, stg),
                                 horizon=hz) if len(y) else "")
                rows.append(r)
    nm = {"Treat all": "Minus treat all",
          "sPESI 1 or more": "Minus sPESI 1+",
          "PESI-10 class III to V": "Minus PESI-10 III to V"}
    table(["Outcome", "Horizon", "Threshold", "Model net benefit"] +
          [nm[z] for z in stgs], rows)
src("fin_dca_key.csv", "cohort EHR+ECG; 10, 90 and 180 days are in the "
    "CSV")

# ================= 6. survival against binary =================
S = "6 Survival against binary"
BS = rd("binsurv_summary.csv")
if BS is not None:
    h("6. Survival super learner minus matched binary super learner "
      "(AUROC)")
    p("Positive values favour the survival model. CV excl: patients who "
      "died first are left out, as in the dissertation.")
    p("")
    TG = {"death": "Death", "comp": "Death or CV readmission",
          "cv": "CV readmission (death competing)",
          "cv_excl": "CV readmission (CV excl)"}
    rows = []
    for fr in ("P", "B"):
        for tg in TG:
            x = BS[(BS.frame == fr) & (BS.target == tg)]
            if not len(x):
                continue
            r = [FRL[fr], TG[tg]]
            for hz in HZ:
                y = x[x.horizon == hz]
                if not len(y):
                    r.append("")
                    continue
                y = y.iloc[0]
                t = ci(y["diff"], y.lo, y.hi, 3)
                if bool(y.excludes_zero):
                    t = "**" + t + "**"
                r.append(rec(S, "%s %s survival minus binary" %
                             (FRL[fr], TG[tg]), y["diff"], y.lo, y.hi, 3,
                             "binsurv_summary.csv",
                             "frame %s, target %s" % (fr, tg), horizon=hz,
                             text=t))
            rows.append(r)
    table(["Pipeline", "Outcome"] + HH, rows)
    n01 = int(BS.within_001.sum())
    p("Within 0.01 AUROC of each other: %d of %d cells (A33's rule). "
      "Bold: the interval excludes zero." % (n01, len(BS)))
    rec(S, "cells within 0.01 AUROC", n01, d=0, src="binsurv_summary.csv",
        filt="column within_001", text="%d of %d" % (n01, len(BS)))
    src("binsurv_summary.csv", "all rows")

# ================= 7. temporal validation =================
S = "7 Temporal validation"
h("7. Temporal validation (development 2010 to 2013, test 2014 onwards)")
p("The super learner is refitted on the development years and tested on "
  "later years. The internal CV model is the main cross-validated model "
  "scored on the same test admissions.")
for fr in ("P", "B"):
    DR = rd("tv2_%s_drift.csv" % fr)
    SU = rd("tv2_%s_summary.csv" % fr)
    CN = rd("tv2_%s_contrasts.csv" % fr)
    if SU is None:
        continue
    h(FRL[fr], 3)
    if DR is not None:
        x = DR[DR.cohort == "EHR+ECG"]
        rows = []
        for _, r in x.iterrows():
            rows.append([OCL[r.outcome], r.period, cnt(r.admissions),
                         cnt(r.events_30), cnt(r.events_365),
                         pct(r.km_risk_365)])
            rec(S, "%s drift %s %s" % (FRL[fr], r.outcome, r.period),
                r.events_365, d=0, src="tv2_%s_drift.csv" % fr,
                filt="cohort EHR+ECG, outcome %s, period %s" %
                (r.outcome, r.period), text=cnt(r.events_365))
        table(["Outcome", "Period", "Admissions", "Events by 30 d",
               "Events by 365 d", "Risk by 365 d"], rows)
        src("tv2_%s_drift.csv" % fr, "cohort EHR+ECG")
    rows = []
    for oc in ("death", "comp", "cv"):
        for mdl in ("Super learner, recalibrated", "Internal CV model"):
            for ms in ("auc", "oe", "slope"):
                x = SU[(SU.outcome == oc) & (SU.cohort == "EHR+ECG") &
                       (SU.period == "test, pooled") & (SU.model == mdl) &
                       (SU.measure == ms)]
                if not len(x):
                    continue
                r = [OCL[oc], mdl.replace(", recalibrated", ""),
                     {"auc": "AUROC", "oe": "O/E", "slope": "Slope"}[ms]]
                for hz in HZ:
                    y = x[x.horizon == hz]
                    r.append(rec(S, "%s %s %s %s test pooled" %
                                 (FRL[fr], oc, mdl, ms), y.estimate.iloc[0],
                                 y.lo.iloc[0], y.hi.iloc[0],
                                 3 if ms == "auc" else 2,
                                 "tv2_%s_summary.csv" % fr,
                                 "outcome %s, EHR+ECG, test pooled, %s, %s"
                                 % (oc, mdl, ms), horizon=hz)
                             if len(y) else "")
                rows.append(r)
        if CN is not None:
            x = CN[(CN.outcome == oc) & (CN.cohort == "EHR+ECG") &
                   (CN.period == "test, pooled") & (CN.measure == "auc")]
            for cn in x.contrast.unique():
                r = [OCL[oc], tidy(cn), "AUROC diff"]
                for hz in HZ:
                    y = x[(x.contrast == cn) & (x.horizon == hz)]
                    if not len(y):
                        r.append("")
                        continue
                    y = y.iloc[0]
                    t = ci(y.estimate, y.lo, y.hi, 3)
                    if bool(y.excludes_zero):
                        t = "**" + t + "**"
                    r.append(rec(S, "%s %s %s" % (FRL[fr], oc, cn),
                                 y.estimate, y.lo, y.hi, 3,
                                 "tv2_%s_contrasts.csv" % fr,
                                 "outcome %s, EHR+ECG, test pooled, %s"
                                 % (oc, cn), horizon=hz, text=t))
                rows.append(r)
    table(["Outcome", "Model or contrast", "Measure"] + HH, rows)
    src("tv2_%s_summary.csv and tv2_%s_contrasts.csv" % (fr, fr),
        "EHR+ECG cohort, test pooled; per-period rows are in the CSVs")

# ================= 8. model choices =================
S = "8 Model choices"
h("8. Model choices and hospice discharges")
rows = []
D7 = rd("sl7_decision.csv")
if D7 is not None:
    r = D7.iloc[0]
    rows.append(["Discharge death model (A30)", "%s; better at %s of %s "
                 "horizons, worse at %s" % (r.final_model_frame_B, r.better,
                                            r.horizons, r.worse),
                 "sl7_decision.csv"])
    rec(S, "final discharge death model", np.nan, src="sl7_decision.csv",
        filt="row 1", text=str(r.final_model_frame_B))
D2 = rd("optB2_decision.csv")
if D2 is not None:
    for _, r in D2.iterrows():
        rows.append(["Composite, %s (A34 SL+ rule)" % r.pipeline,
                     "%s (SL+ better at %s, worse at %s horizons)" %
                     (r.final_composite_model, r.better, r.worse),
                     "optB2_decision.csv"])
        rec(S, "composite model %s" % r.pipeline, np.nan,
            src="optB2_decision.csv", filt="pipeline %s" % r.pipeline,
            text=str(r.final_composite_model))
CC = rd("optB2_core_check.csv")
if CC is not None:
    mx = float(np.abs(CC.v7_minus_core).max())
    rows.append(["v7-core against v7 (A33)",
                 "largest AUROC gap %s; representative at every horizon: "
                 "%s" % (num(mx, 4), bool(CC.representative.all())),
                 "optB2_core_check.csv"])
    rec(S, "v7-core largest AUROC gap", mx, d=4,
        src="optB2_core_check.csv", filt="column v7_minus_core")
table(["Choice", "Result", "Source"], rows)
NH = rd("sl7nh_summary.csv")
if NH is not None:
    p("**The hospice indicator's contribution at discharge (death): v7 "
      "minus v7 without it, AUROC**")
    p("")
    rows = []
    for ev in ("all", "hospice excluded"):
        x = NH[(NH.evaluation == ev) & (NH.measure == "auc") &
               NH.contrast.str.contains("v7nh, recalibrated$")]
        r = [ev]
        for hz in HZ:
            y = x[x.horizon == hz]
            if not len(y):
                r.append("")
                continue
            y = y.iloc[0]
            t = ci(y.estimate, y.lo, y.hi, 3)
            if y.lo > 0 or y.hi < 0:
                t = "**" + t + "**"
            r.append(rec(S, "v7 minus v7nh AUROC (%s)" % ev, y.estimate,
                         y.lo, y.hi, 3, "sl7nh_summary.csv",
                         "evaluation %s, v7 - v7nh, auc" % ev, horizon=hz,
                         text=t))
        rows.append(r)
    table(["Evaluation"] + HH, rows)
    p("Bold: the 95% interval excludes zero.")
    src("sl7nh_summary.csv", "measure auc, contrast v7 minus v7nh")
OB = rd("optB2_sl_summary.csv")
if OB is not None:
    p("**Final discharge models with hospice discharges left out "
      "(AUROC, EHR+ECG cohort)**")
    p("")
    rows = []
    for oc, mdl in (("death", "Super learner v7, recalibrated"),
                    ("comp", "Super learner, recalibrated"),
                    ("cv", "Super learner, recalibrated")):
        r = [OCL[oc]]
        for hz in HZ:
            y = OB[(OB.frame == "B") & (OB.outcome == oc) &
                   (OB.cohort == "EHR+ECG") &
                   (OB.evaluation == "hospice excluded") &
                   (OB.model == mdl) & (OB.set == "EHR+ECG") &
                   (OB.horizon == hz) & (OB.measure == "auc")]
            r.append(rec(S, "B %s hospice excluded AUROC" % oc,
                         y.estimate.iloc[0], y.lo.iloc[0], y.hi.iloc[0], 3,
                         "optB2_sl_summary.csv",
                         "frame B, %s, EHR+ECG, hospice excluded, %s, auc"
                         % (oc, mdl), horizon=hz) if len(y) else "")
        rows.append(r)
    table(["Outcome"] + HH, rows)
    src("optB2_sl_summary.csv", "evaluation 'hospice excluded'")

# ================= 9. at-diagnosis sensitivity =================
S = "9 At-diagnosis sensitivity"
h("9. At-diagnosis sensitivity analyses (death, EHR+ECG cohort, AUROC)")
p("Strict history counts chronic conditions from earlier admissions "
  "only (A36). The second block keeps only admissions with presentation "
  "vital signs (D17).")
p("")
SPM = rd("slP_metrics.csv")
rows = []
for ev in ("all", "with presentation vital signs"):
    for mdl, st in (("Super learner P, recalibrated", "EHR+ECG"),
                    ("Super learner Psens, recalibrated", "EHR+ECG"),
                    ("sPESI", "score"), ("sPESI (sens)", "score"),
                    ("PESI-10", "score"), ("PESI-10 (sens)", "score")):
        r = [ev, mdl.replace(", recalibrated", "")]
        for hz in HZ:
            y = SPM[(SPM.cohort == "EHR+ECG") & (SPM.evaluation == ev) &
                    (SPM.model == mdl) & (SPM.set == st) &
                    (SPM.horizon == hz) & (SPM.measure == "auc")]
            r.append(rec(S, "%s %s AUROC" % (ev, mdl), y.estimate.iloc[0],
                         y.lo.iloc[0], y.hi.iloc[0], 3, "slP_metrics.csv",
                         "cohort EHR+ECG, evaluation %s, model %s, auc"
                         % (ev, mdl), horizon=hz) if len(y) else "")
        rows.append(r)
table(["Patients", "Model or score"] + HH, rows)
src("slP_metrics.csv", "cohort EHR+ECG, measure auc")

# ================= 10. interpretability =================
S = "10 Interpretability"
TM = rd("interpP_spline_terms.csv")
IM = rd("interpP_shap_importance.csv")
DM = rd("interpP_shap_domains.csv")
AG = rd("interpP_agreement.csv")
TW = rd("interpP_fig_twins.csv")
if TM is not None:
    h("10. Interpretability at diagnosis (death, EHR+ECG cohort)")
    x = TM[TM.cohort == "EHR+ECG"].sort_values("term_sd", ascending=False)
    rows = []
    for _, r in x.head(12).iterrows():
        lb = lab_of(r.feature, r["label"])
        rows.append([lb, r.domain, r.type,
                     rec(S, "spline term SD: %s" % lb, r.term_sd,
                         d=3, src="interpP_spline_terms.csv",
                         filt="cohort EHR+ECG, feature %s" % r.feature),
                     num(r.edf, 2)])
    p("**Penalised spline Cox: the 12 largest terms (spread of each "
      "term's log hazard ratio)**. Names follow the figures (stage 20q): "
      "Blood pH, ALT and CK-MB are 'pH', 'Alanine aminotransferase alt' "
      "and 'Creatine kinase mb isoenzyme' in the CSV.")
    p("")
    table(["Feature", "Domain", "Type", "Term SD", "edf"], rows)
    src("interpP_spline_terms.csv", "cohort EHR+ECG, sorted by term_sd")
if IM is not None:
    for mdl in IM.model.unique():
        x = IM[(IM.model == mdl) & (IM.horizon == 365)]
        x = x.sort_values("mean_abs_shap", ascending=False).head(10)
        rows = []
        for _, r in x.iterrows():
            lb = lab_of(r.feature, r["label"])
            rows.append([int(r["rank"]), lb, r.domain,
                         rec(S, "%s SHAP share 365 d: %s" % (mdl, lb),
                             r.share, d=3,
                             src="interpP_shap_importance.csv",
                             filt="model %s, horizon 365, feature %s" %
                             (mdl, r.feature), text=pct(r.share))])
        p("**%s: the 10 largest SurvSHAP(t) shares at 365 days**" % mdl)
        p("")
        table(["Rank", "Feature", "Domain", "Share"], rows)
    if TW is not None and len(TW):
        p("Charlson and Elixhauser twins are identical columns, so each "
          "model splits their credit: %s." %
          "; ".join(TW["label"].tolist()))
    src("interpP_shap_importance.csv", "horizon 365")
if DM is not None:
    p("**Share of importance by domain, SurvSHAP(t)**")
    p("")
    rows = []
    for mdl in DM.model.unique():
        for dmn in DM.domain.unique():
            r = [mdl, dmn]
            for hz in HZ:
                y = DM[(DM.model == mdl) & (DM.domain == dmn) &
                       (DM.horizon == hz)]
                r.append(rec(S, "%s domain share %s" % (mdl, dmn),
                             y.share.iloc[0], d=3,
                             src="interpP_shap_domains.csv",
                             filt="model %s, domain %s" % (mdl, dmn),
                             horizon=hz, text=pct(y.share.iloc[0]))
                         if len(y) else "")
            rows.append(r)
    table(["Model", "Domain"] + HH, rows)
    src("interpP_shap_domains.csv", "all rows")
if AG is not None:
    rows = []
    for _, r in AG.iterrows():
        rows.append([r.model, r.comparison,
                     rec(S, "Spearman %s" % r.model, r.spearman, d=2,
                         src="interpP_agreement.csv",
                         filt="model %s" % r.model)])
    p("**Agreement in feature ranking (Spearman)**")
    p("")
    table(["Model", "Comparison", "Spearman"], rows)
    src("interpP_agreement.csv", "all rows")

# ================= 11. CTPA blocks =================
S = "11 CTPA blocks"
CB = rd("ctpa_blocks_P_time.csv")
CBC = rd("ctpa_blocks_P_contrasts.csv")
if CB is not None:
    h("11. CTPA blocks at diagnosis (three-modality cohort, death)")
    QS = [("HR days 0-30 (per SD)", "HR 0 to 30 d", 2),
          ("HR days 31-365 (per SD)", "HR 31 to 365 d", 2),
          ("early / late ratio", "Early / late", 2),
          ("AUROC, death by day 30", "AUROC by 30 d", 3),
          ("AUROC, deaths in days 31-365 among day-30 survivors",
           "AUROC 31 to 365 d", 3)]
    for mdl in ("alone", "with EHR"):
        rows = []
        for blk in CB.block.unique():
            x = CB[(CB.block == blk) & (CB.model == mdl)]
            if not len(x):
                continue
            r = [blk, int(x.n_features.iloc[0])]
            for q, _, d in QS:
                y = x[x.quantity == q]
                r.append(rec(S, "%s (%s): %s" % (blk, mdl, q),
                             y.estimate.iloc[0], y.lo.iloc[0],
                             y.hi.iloc[0], d, "ctpa_blocks_P_time.csv",
                             "block %s, model %s, quantity %s" %
                             (blk, mdl, q)) if len(y) else "")
            rows.append(r)
        p("**Block score %s**" % ("on its own" if mdl == "alone" else
                                  "with the EHR score in the model"))
        p("")
        table(["Block", "Features"] + [q[1] for q in QS], rows)
    src("ctpa_blocks_P_time.csv", "all rows")
if CBC is not None:
    rows = []
    for _, r in CBC.iterrows():
        dd = 2 if r.quantity.startswith("ratio") else 3
        t = ci(r.estimate, r.lo, r.hi, dd)
        if bool(r.excludes_null):
            t = "**" + t + "**"
        rows.append([r.contrast, r.model, r.quantity,
                     rec(S, "%s (%s): %s" % (r.contrast, r.model,
                                             r.quantity), r.estimate, r.lo,
                         r.hi, dd, "ctpa_blocks_P_contrasts.csv",
                         "contrast %s, model %s" % (r.contrast, r.model),
                         text=t)])
    p("**Block contrasts** (bold: the interval excludes no difference)")
    p("")
    table(["Contrast", "Model", "Quantity", "Estimate"], rows)
    src("ctpa_blocks_P_contrasts.csv", "all rows")

# ================= 12. subgroups =================
S = "12 Subgroups"
SG = rd("subgroup_metrics.csv")
if SG is not None:
    h("12. Subgroups at discharge (Super learner v7, death, AUROC)")
    p("Cells with fewer than 10 deaths have no AUROC.")
    p("")
    x = SG[(SG.frame == "B") & (SG.measure == "auc")]
    rows = []
    for (v, g), y in x.groupby(["variable", "group"], sort=False):
        r = [v, g, cnt(y.patients.iloc[0])]
        for hz in (30, 365):
            z = y[y.horizon == hz]
            if not len(z):
                r += ["", ""]
                continue
            z = z.iloc[0]
            r.append(cnt(z.deaths))
            r.append(rec(S, "B v7 AUROC %s %s" % (v, g), z.estimate, z.lo,
                         z.hi, 3, "subgroup_metrics.csv",
                         "frame B, %s %s, auc" % (v, g), horizon=hz)
                     if np.isfinite(z.estimate) else "not estimated")
        rows.append(r)
    table(["Variable", "Group", "Discharges", "Deaths 30 d", "AUROC 30 d",
           "Deaths 365 d", "AUROC 365 d"], rows)
    src("subgroup_metrics.csv", "frame B, measure auc")

# ================= 13. checks =================
h("13. Checks across the final stages")
rows = []
nfail = 0
for f in CHECKS:
    x = rd(f)
    if x is None:
        rows.append([f, "not found", "", "", ""])
        continue
    np_, nw, nf = [(x.status == s).sum() for s in ("PASS", "WARN",
                                                   "FAIL")]
    nfail += nf
    rows.append([f, int(np_), int(nw), int(nf),
                 int((x.status == "INFO").sum())])
table(["Checks file", "PASS", "WARN", "FAIL", "INFO"], rows)
chk("FAIL rows across the final stages' checks files", nfail,
    "PASS" if nfail == 0 else "FAIL")

h("Source files", 2)
table(["File", "Last modified"],
      [[f, t.strftime("%Y-%m-%d %H:%M")] for f, t in sorted(USED.items())])

# ---------------- write ----------------
print("\n=== WRITE ===")
head = ["# P2 verified results", "",
        "Generated %s by `scripts/13e_verified_results.py` from the saved "
        "aggregate CSVs. Every number below is read by code from the "
        "file named under its table; `verified_results.csv` holds each "
        "one with its exact source rows. Intervals are 95%% "
        "patient-clustered bootstrap intervals unless stated. Bold marks "
        "an interval that excludes zero; a bound printed as 0.000 then "
        "lies just beyond zero. Counts of 1 to 9 are shown as <10."
        % datetime.now().strftime(
            "%Y-%m-%d %H:%M"), ""]
small = [r for r in ROWS if isinstance(r["text"], str) and
         r["text"].strip().isdigit() and 1 <= int(r["text"]) <= 9]
chk("counts from 1 to 9 written as numbers", len(small),
    "PASS" if not small else "FAIL")
stop_if_fail()
with open(os.path.join(RS, "verified_results.md"), "w",
          encoding="utf-8") as fh:
    fh.write("\n".join(head + MD) + "\n")
pd.DataFrame(ROWS).to_csv(os.path.join(RS, "verified_results.csv"),
                          index=False, encoding="utf-8-sig")
chk("numbers recorded", len(ROWS), "INFO")
chk("source files read", len(USED), "INFO")
ck = pd.DataFrame(CHK)
ck.to_csv(os.path.join(RS, "verified_checks.csv"), index=False)
print("\n  checks: %d PASS, %d WARN, %d FAIL" % (
    (ck.status == "PASS").sum(), (ck.status == "WARN").sum(),
    (ck.status == "FAIL").sum()))
print("\nOUTPUT FILES")
for f in ("verified_results.md", "verified_results.csv",
          "verified_checks.csv"):
    print("  ", os.path.join(RS, f))
print("  ", os.path.join(RS, "verified_log.txt"),
      " (terminal log)")
print("\nOutputs: verified_results.md, verified_checks.csv")