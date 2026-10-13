"""13c_results_tables.py
Tables 2 and 3: performance and added value of the final models
(PREREGISTRATION_ADDENDUM_2.md A38).

Reads aggregate results only; no model is refitted. Each table is built
for the EHR+ECG cohort and, as a supplementary table, for the
three-modality cohort with CTPA (Tables S2 and S3).

Table 2 covers each pipeline (at diagnosis, from admission; at discharge,
from discharge) and outcome (death; death or CV readmission; CV
readmission with death as a competing event) at 10, 30, 90, 180 and 365
days: events, AUROC, IPA and observed/expected (stages 18q and 18r, with
their bootstrap intervals), calibration slope (stage 13b, with its
interval), ICI (stage 13b, point estimate only, because resampling noise
inflates the bootstrap ICI and its percentile interval can sit above the
estimate) and, for CV readmission, the AUROC excluding patients who died
first. The final models are Super learner P (death) and the stage 18r
super learners at diagnosis, and Super learner v7 (death) and the stage
18r super learners at discharge.

Table 3 gives added value, final model minus comparator, with bootstrap
intervals. At diagnosis, death is compared with sPESI and PESI-10 (AUROC)
and the specificity gain at sPESI's sensitivity is given; every outcome is
compared with the EHR alone (AUROC and IPA). At discharge, the death rows
use v7-core (A33), whose EHR-alone version exists, and its distance from
v7 is given in the table note. Bold marks an interval that excludes zero.

Counts from 1 to 9 are shown as "<10". Checks (FAIL): inputs found, every
cell found, events agree with stage 13b, and stage 13b's AUROC point
estimates agree with the sources to within 1e-6.

Inputs:  results/slP_metrics.csv, slP_contrasts.csv, optB2_sl_summary.csv,
         optB2_sl_contrasts.csv, optB2_core_check.csv,
         optionB_outcome_events.csv, fin_metrics.csv (stage 13b)
Outputs: results/tables_2_3.xlsx (the formatted tables, four sheets),
         results/table2.csv, table3.csv, tableS2.csv, tableS3.csv,
         results/tables23_long.csv (every number, unrounded),
         results/tables23_checks.csv, results/tables23_log.txt (terminal
         log)
Usage:   python scripts/13c_results_tables.py
Notes:   Run 13b_calibration_dca_final.R first, since it writes
         fin_metrics.csv.
"""
import os
import sys

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

P2 = os.path.join(os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd())), "SurvivalProject")
RS = os.path.join(P2, "results")
HZ = [10, 30, 90, 180, 365]
TOL = 1e-6
NO_CI = {"ici"}
CHK = []


def chk(item, value, status):
    CHK.append({"item": item, "value": str(value), "status": status})
    print("  %-5s %-62s %s" % (status, item, value), flush=True)


def stop_if_fail():
    ck = pd.DataFrame(CHK)
    if (ck.status == "FAIL").any():
        ck.to_csv(os.path.join(RS, "tables23_checks.csv"), index=False)
        print("\nSTOP: a check failed. see: tables23_checks.csv")
        sys.exit(1)


# ---------------- inputs ----------------
print("=== INPUTS ===")
FILES = ["slP_metrics.csv", "slP_contrasts.csv", "optB2_sl_summary.csv",
         "optB2_sl_contrasts.csv", "optB2_core_check.csv",
         "optionB_outcome_events.csv", "fin_metrics.csv"]
miss = [f for f in FILES if not os.path.exists(os.path.join(RS, f))]
chk("input files found", miss[0] if miss else "all %d" % len(FILES),
    "FAIL" if miss else "PASS")
if "fin_metrics.csv" in miss:
    print("  (fin_metrics.csv comes from stage 13b: run it first)")
stop_if_fail()
rd = {f: pd.read_csv(os.path.join(RS, f)) for f in FILES}
SP = rd["slP_metrics.csv"]
SP = SP[SP.evaluation == "all"]
SC = rd["slP_contrasts.csv"]
SC = SC[SC.evaluation == "all"]
OB = rd["optB2_sl_summary.csv"]
OB = OB[OB.evaluation == "all"]
OC = rd["optB2_sl_contrasts.csv"]
OC = OC[OC.evaluation == "all"]
CORE = rd["optB2_core_check.csv"]
EV = rd["optionB_outcome_events.csv"]
FM = rd["fin_metrics.csv"]

SETN = {"EHR+ECG": "EHR+ECG", "3-mod": "EHR+ECG+CTPA"}
FRL = {"P": "At diagnosis", "B": "At discharge"}
EVF = {"P": "P (at diagnosis)", "B": "B (at discharge)"}
EVC = {"death": "death", "comp": "composite", "cv": "cv_first"}
OCL = {"death": "Death",
       "comp": "Death or CV readmission",
       "cv": "CV readmission, with death as a competing event"}
MODEL = {("P", "death"): "Super learner P, recalibrated",
         ("B", "death"): "Super learner v7, recalibrated"}
SL = "Super learner, recalibrated"
LONG = []


def model_of(fr, oc):
    return MODEL.get((fr, oc), SL)


def one(df, what):
    """exactly one row expected; returns it or None (and a FAIL)."""
    if len(df) == 1:
        return df.iloc[0]
    chk("cell found: %s" % what, "%d rows" % len(df), "FAIL")
    return None


def perf(fr, oc, co, h, m):
    st = SETN[co]
    mdl = model_of(fr, oc)
    if m in ("slope", "ici"):
        x = FM[(FM.frame == fr) & (FM.outcome == oc) & (FM.cohort == co) &
               (FM.horizon == h) & (FM.measure == m)]
        src = "stage 13b"
    elif fr == "P" and oc == "death":
        x = SP[(SP.cohort == co) & (SP.model == mdl) & (SP.set == st) &
               (SP.horizon == h) & (SP.measure == m)]
        src = "stage 18q (slP_metrics)"
    else:
        x = OB[(OB.frame == fr) & (OB.outcome == oc) & (OB.cohort == co) &
               (OB.model == mdl) & (OB.set == st) & (OB.horizon == h) &
               (OB.measure == m)]
        src = "stage 18r (optB2_sl_summary)"
    r = one(x, "%s %s %s %d %s" % (fr, oc, co, h, m))
    if r is None:
        return None
    lo, hi = (np.nan, np.nan) if m in NO_CI else (r.lo, r.hi)
    return {"estimate": r.estimate, "lo": lo, "hi": hi, "source": src}


CON_P = {"sPESI": ("Super learner P, recalibrated - sPESI", "auc"),
         "PESI-10": ("Super learner P, recalibrated - PESI-10", "auc"),
         "spec": ("Super learner P, recalibrated - sPESI",
                  "specificity at sPESI's sensitivity"),
         "EHR auc": ("Super learner P, recalibrated - Super learner P-EHR, "
                     "recalibrated", "auc"),
         "EHR ipa": ("Super learner P, recalibrated - Super learner P-EHR, "
                     "recalibrated", "ipa")}
CON_O = "Super learner, recalibrated - Super learner-EHR, recalibrated"


def gain(fr, oc, co, h, key):
    st = SETN[co]
    if fr == "P" and oc == "death":
        con, m = CON_P[key]
        x = SC[(SC.cohort == co) & (SC.set == st) & (SC.contrast == con) &
               (SC.horizon == h) & (SC.measure == m)]
        src = "stage 18q (slP_contrasts)"
    else:
        m = key.split()[-1]
        x = OC[(OC.frame == fr) & (OC.outcome == oc) & (OC.cohort == co) &
               (OC.set == st) & (OC.contrast == CON_O) & (OC.horizon == h) &
               (OC.measure == m)]
        src = "stage 18r (optB2_sl_contrasts)"
    r = one(x, "%s %s %s %d %s" % (fr, oc, co, h, key))
    if r is None:
        return None
    return {"estimate": r.estimate, "lo": r.lo, "hi": r.hi, "source": src,
            "excludes_zero": bool(r.lo > 0 or r.hi < 0)}


def num(v, d):
    s = ("{:.%df}" % d).format(v)
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    if s.startswith("-"):
        s = "\u2212" + s[1:]
    return s


def ci(c, d):
    if c is None:
        return "\u2014"
    if not np.isfinite(c["lo"]):
        return num(c["estimate"], d)
    return "%s (%s to %s)" % (num(c["estimate"], d), num(c["lo"], d),
                              num(c["hi"], d))


def events(fr, oc, co, h):
    x = EV[(EV.frame == EVF[fr]) & (EV.cohort == co) & (EV.horizon == h)]
    r = one(x, "events %s %s %d" % (fr, co, h))
    if r is None:
        return None, None
    return int(r[EVC[oc]]), int(r.patients)


def fmt_ev(e, n):
    if e is None:
        return "\u2014"
    if 0 < e < 10:
        return "<10"
    return "{:,} ({:.1f})".format(e, 100.0 * e / n)


# ---------------- build the rows ----------------
print("\n=== BUILD ===")
PERF_ROWS = [("AUROC (95% CI)", "auc", 3),
             ("IPA (95% CI)", "ipa", 3),
             ("Observed / expected (95% CI)", "oe", 2),
             ("Calibration slope (95% CI)", "slope", 2),
             ("ICI", "ici", 3)]
GAIN_P = [("AUROC gain over sPESI", "sPESI", 3),
          ("AUROC gain over PESI-10", "PESI-10", 3),
          ("Specificity gain at sPESI's sensitivity", "spec", 3),
          ("AUROC gain over the EHR alone", "EHR auc", 3),
          ("IPA gain over the EHR alone", "EHR ipa", 3)]
GAIN_O = [("AUROC gain over the EHR alone", "EHR auc", 3),
          ("IPA gain over the EHR alone", "EHR ipa", 3)]


def build(co, kind):
    rows = []
    for fr in ("P", "B"):
        for oc in ("death", "comp", "cv"):
            _, npat = events(fr, oc, co, HZ[0])
            unit = "admissions" if fr == "P" else "discharges"
            band = "%s (%s %s): %s" % (FRL[fr], "{:,}".format(npat or 0),
                                       unit, OCL[oc])
            if kind == "gain" and fr == "B" and oc == "death":
                band += " \u2020"
            rows.append(("band", band, None, None))
            if kind == "perf":
                vals = []
                for h in HZ:
                    e, n = events(fr, oc, co, h)
                    vals.append(fmt_ev(e, n))
                    LONG.append({"table": "2" if co == "EHR+ECG" else "S2",
                                 "pipeline": FRL[fr], "outcome": OCL[oc],
                                 "cohort": co, "row": "Events", "horizon": h,
                                 "estimate": e if (e is None or e == 0 or
                                                   e >= 10) else np.nan,
                                 "lo": np.nan, "hi": np.nan,
                                 "excludes_zero": np.nan,
                                 "source": "optionB_outcome_events"})
                rows.append(("row", "Events, n (%)", vals, [False] * 5))
                spec = list(PERF_ROWS)
                if oc == "cv":
                    spec.insert(1, ("AUROC excluding earlier deaths (95% CI)",
                                    "auc_excl", 3))
                for lab, m, d in spec:
                    vals = []
                    for h in HZ:
                        c = perf(fr, oc, co, h, m)
                        vals.append(ci(c, d))
                        if c is not None:
                            LONG.append({
                                "table": "2" if co == "EHR+ECG" else "S2",
                                "pipeline": FRL[fr], "outcome": OCL[oc],
                                "cohort": co, "row": lab, "horizon": h,
                                "estimate": c["estimate"], "lo": c["lo"],
                                "hi": c["hi"], "excludes_zero": np.nan,
                                "source": c["source"]})
                    rows.append(("row", lab, vals, [False] * 5))
            else:
                spec = GAIN_P if (fr == "P" and oc == "death") else GAIN_O
                for lab, key, d in spec:
                    vals = []
                    bold = []
                    for h in HZ:
                        c = gain(fr, oc, co, h, key)
                        vals.append(ci(c, d))
                        bold.append(bool(c is not None and
                                         c["excludes_zero"]))
                        if c is not None:
                            LONG.append({
                                "table": "3" if co == "EHR+ECG" else "S3",
                                "pipeline": FRL[fr], "outcome": OCL[oc],
                                "cohort": co, "row": lab, "horizon": h,
                                "estimate": c["estimate"], "lo": c["lo"],
                                "hi": c["hi"],
                                "excludes_zero": c["excludes_zero"],
                                "source": c["source"]})
                    rows.append(("row", lab, vals, bold))
    return rows


T = {"Table 2": build("EHR+ECG", "perf"), "Table S2": build("3-mod", "perf"),
     "Table 3": build("EHR+ECG", "gain"), "Table S3": build("3-mod", "gain")}
nmiss = sum(v == "\u2014" for rows in T.values() for k, _, vals, _ in rows
            if k == "row" for v in vals)
chk("every cell found (cells missing)", nmiss, "PASS" if nmiss == 0
    else "FAIL")
stop_if_fail()

# events and AUROC agree with stage 13b
bad_e, bad_a, n_a = 0, 0.0, 0
for fr in ("P", "B"):
    for oc in ("death", "comp", "cv"):
        for co in ("EHR+ECG", "3-mod"):
            for h in HZ:
                e, _ = events(fr, oc, co, h)
                x = FM[(FM.frame == fr) & (FM.outcome == oc) &
                       (FM.cohort == co) & (FM.horizon == h)]
                ev13 = str(x.events.iloc[0]) if len(x) else "missing"
                ev0 = "<10" if 0 < e < 10 else str(e)
                bad_e += int(ev13 != ev0)
                a13 = x[x.measure == "auc"].estimate
                a0 = perf(fr, oc, co, h, "auc")
                if len(a13) and a0 is not None:
                    bad_a = max(bad_a, abs(float(a13.iloc[0]) -
                                           a0["estimate"]))
                    n_a += 1
chk("events agree with stage 13b (cells that differ)", bad_e,
    "PASS" if bad_e == 0 else "FAIL")
chk("stage 13b AUROC agrees with the sources (n; max diff)",
    "%d; %.1e" % (n_a, bad_a), "PASS" if n_a == 60 and bad_a <= TOL
    else "FAIL")
core = float(np.abs(CORE.v7_minus_core).max())
chk("v7-core: largest AUROC distance from v7 (EHR+ECG cohort)",
    "%.4f" % core, "INFO")
stop_if_fail()

# ---------------- notes ----------------
HEADS = ["Measure"] + ["%d days" % h for h in HZ]
TITLES = {
    "Table 2": ("Table 2. Discrimination and calibration of the final "
                "models at each horizon (EHR+ECG cohort)"),
    "Table S2": ("Table S2. Discrimination and calibration of the final "
                 "models at each horizon (three-modality cohort, with "
                 "CTPA)"),
    "Table 3": ("Table 3. Added value of the final models over clinical "
                "scores and over the EHR alone (EHR+ECG cohort)"),
    "Table S3": ("Table S3. Added value of the final models over clinical "
                 "scores and over the EHR alone (three-modality cohort, "
                 "with CTPA)")}
N_PERF = (
    "Out-of-fold predictions of the final models, recalibrated per "
    "horizon. At diagnosis, death uses Super learner P; at discharge, "
    "Super learner v7. The composite and CV readmission use the stage 18r "
    "super learners. Time runs from admission at diagnosis and from "
    "discharge at discharge. Events count outcomes by each horizon. CV "
    "readmission counts a readmission before any death; the AUROC "
    "excluding earlier deaths leaves out patients who died first. "
    "Intervals are 95% patient-clustered bootstrap intervals from 2,000 "
    "samples. ICI has no interval, because resampling noise inflates the "
    "bootstrap ICI. AUROC, area under the ROC curve. IPA, index of "
    "prediction accuracy (Brier skill against the observed proportion). "
    "ICI, integrated calibration index (mean gap between the smoothed "
    "calibration curve and the prediction). CV, cardiovascular. Counts of "
    "1 to 9 are shown as <10.")
N_GAIN = (
    "Each cell is the final model minus the comparator, with a 95%% "
    "patient-clustered bootstrap interval (2,000 samples). Bold: the "
    "interval excludes zero (a bound shown as 0.000 lies just above "
    "zero). The EHR alone is a super learner over the "
    "same EHR candidates, without ECG or CTPA. The specificity gain "
    "compares specificities at sPESI's sensitivity (sPESI 1 or more as "
    "high risk). \u2020 At discharge, the death rows use v7-core, the "
    "super learner over the four EHR+ families (A33). Its AUROC is within "
    "%.3f of v7's at every horizon. sPESI, simplified PESI. PESI-10, PESI "
    "without altered mental status. IPA, index of prediction accuracy. "
    "CV, cardiovascular." % core)
NOTES = {"Table 2": N_PERF, "Table S2": N_PERF, "Table 3": N_GAIN,
         "Table S3": N_GAIN}

# ---------------- write ----------------
print("\n=== WRITE ===")
CSVN = {"Table 2": "table2.csv", "Table S2": "tableS2.csv",
        "Table 3": "table3.csv", "Table S3": "tableS3.csv"}
for name, rows in T.items():
    flat = []
    for kind, lab, vals, _ in rows:
        flat.append([lab] + ([""] * len(HZ) if kind == "band" else vals))
    pd.DataFrame(flat, columns=HEADS).to_csv(
        os.path.join(RS, CSVN[name]), index=False, encoding="utf-8-sig")
pd.DataFrame(LONG).to_csv(os.path.join(RS, "tables23_long.csv"),
                          index=False, encoding="utf-8-sig")

thin = Side(style="thin", color="FF000000")
grey = Side(style="thin", color="FFBFBFBF")
med = Side(style="medium", color="FF000000")
dbl = Side(style="double", color="FF000000")
cen = Alignment(horizontal="center", vertical="center", wrap_text=True)
f10 = Font(name="Calibri", size=10)
f10b = Font(name="Calibri", size=10, bold=True)
wb = Workbook()
wb.remove(wb.active)
for name in ("Table 2", "Table 3", "Table S2", "Table S3"):
    rows = T[name]
    ws = wb.create_sheet(name)
    ws.sheet_view.showGridLines = False
    ws.page_setup.orientation = "landscape"
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0
    nc = len(HEADS)
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=nc)
    c = ws.cell(row=1, column=1, value=TITLES[name])
    c.font = Font(name="Calibri", size=12, bold=True)
    c.alignment = Alignment(horizontal="left", vertical="center")
    for j, h in enumerate(HEADS, start=1):
        c = ws.cell(row=2, column=j, value=h)
        c.font = f10b
        c.fill = PatternFill("solid", fgColor="FFFFFFFF")
        c.alignment = cen
        c.border = Border(left=thin, right=thin, top=thin, bottom=dbl)
    ws.row_dimensions[2].height = 20
    r = 3
    first_band = True
    blocks = []
    for kind, lab, vals, bold in rows:
        if kind == "band":
            blocks.append([r, []])
            ws.merge_cells(start_row=r, start_column=1, end_row=r,
                           end_column=nc)
            for j in range(1, nc + 1):
                c = ws.cell(row=r, column=j)
                c.fill = PatternFill("solid", fgColor="FFF2F2F2")
                c.border = Border(left=thin, right=thin,
                                  top=dbl if first_band else None,
                                  bottom=med)
            c = ws.cell(row=r, column=1, value=lab)
            c.font = f10b
            c.alignment = cen
            ws.row_dimensions[r].height = 15.6
            first_band = False
        else:
            for j, v in enumerate([lab] + vals, start=1):
                c = ws.cell(row=r, column=j, value=v)
                c.font = f10b if (j > 1 and bold[j - 2]) else f10
                c.alignment = cen
            blocks[-1][1].append(r)
        r += 1
    for _, rr_list in blocks:
        for k, rr in enumerate(rr_list):
            top = med if k == 0 else None
            bot = med if k == len(rr_list) - 1 else grey
            for j in range(1, nc + 1):
                ws.cell(row=rr, column=j).border = Border(
                    left=thin, right=thin, top=top, bottom=bot)
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=nc)
    c = ws.cell(row=r, column=1, value=NOTES[name])
    c.font = Font(name="Calibri", size=9, italic=True)
    c.alignment = Alignment(horizontal="left", vertical="top",
                            wrap_text=True)
    ws.row_dimensions[r].height = 92
    ws.column_dimensions["A"].width = 40
    for j in range(2, nc + 1):
        ws.column_dimensions[chr(64 + j)].width = 21
xl = os.path.join(RS, "tables_2_3.xlsx")
wb.save(xl)

# ---------------- summary ----------------
print("\n=== TABLE 2 (EHR+ECG cohort), AUROC ===")
L = pd.DataFrame(LONG)
x = L[(L.table == "2") & (L.row == "AUROC (95% CI)")]
print(x.pivot_table(index=["pipeline", "outcome"], columns="horizon",
                    values="estimate").round(3).to_string())
print("\n=== TABLE 3 (EHR+ECG cohort), gains with intervals excluding 0 ===")
x = L[(L.table == "3")]
print(x.groupby(["pipeline", "outcome", "row"]).excludes_zero.sum()
      .astype(int).to_string())
ck = pd.DataFrame(CHK)
ck.to_csv(os.path.join(RS, "tables23_checks.csv"), index=False)
print("\n  checks: %d PASS, %d WARN, %d FAIL" % (
    (ck.status == "PASS").sum(), (ck.status == "WARN").sum(),
    (ck.status == "FAIL").sum()))
print("\nOUTPUT FILES")
for f in ("tables_2_3.xlsx", "table2.csv", "table3.csv", "tableS2.csv",
          "tableS3.csv", "tables23_long.csv", "tables23_checks.csv"):
    print("  ", os.path.join(RS, f))
print("  ", os.path.join(RS, "tables23_log.txt"),
      " (terminal log)")
print("\nOutputs: table2.csv, table3.csv, tableS2.csv, tableS3.csv,",
      "tables23_long.csv, tables23_checks.csv")