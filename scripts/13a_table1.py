"""13a_table1.py
Table 1: characteristics of the PE admissions at presentation and at
discharge (A36).

The columns are all at-diagnosis admissions (Frame P); died within 30 days
of admission; alive at 30 days; with a presentation CTPA (three-modality
cohort); discharged alive (Frame B, the discharge cohort); and data
available (% of all admissions, presentation measurements only). The rows
are demographics; first vital signs and labs at presentation (24 h before
to 24 h after admission, as in the at-diagnosis model); history as the
at-diagnosis model uses it (A36: chronic conditions from earlier or
index-admission codes) and, for cancer, heart failure and chronic lung
disease, from earlier admissions only (the strict sensitivity version);
sPESI and PESI-10 (A36 history); modalities at presentation (ECG within
24 h of admission, D19); the index admission (length of stay, in-hospital
death, discharge destination); and outcomes from admission (death, CV
readmission first, composite by 30 and 365 days) and, for the discharge
cohort, from discharge.

Counts from 1 to 9 are shown as "<10" (MIMIC data use). When one of the
two 30-day columns is "<10", the other is not shown either ("†"), since
the two add up to the all-admissions count. Decimals are fixed per row
and the same in every column. Race groups follow the subgroup analysis of
stage 22a (Portuguese counted as White); a mismatch with its counts is a
WARN, and the labels grouped as Other are listed in table1_checks.csv
(counts under 10 as "<10").
Version 3 (A36): fixed decimals per row, the complementary 30-day
suppression and the Portuguese grouping.

Inputs:  data/processed/presentation/frame_P_outcomes.csv, ehrP.csv,
         ehrP_sens.csv, scores_P.csv; data/processed/frame_B.csv;
         data/inputs/demographics.csv, followup_365.csv,
         mimic_labels_harmonised.csv, index_admission_times.csv
Outputs: results/table1.xlsx (the formatted table), results/table1.csv
         (same numbers, plain), results/table1_checks.csv,
         results/table1_log.txt (terminal log)
Usage:   python scripts/13a_table1.py
Notes:   index_admission_times.csv is taken from data/inputs, or else from
         Dissertation/figwork/data in the dissertation
         (pe-multimodal-prognosis-msc).
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
IN = os.path.join(P2, "data", "inputs")
PR = os.path.join(P2, "data", "processed")
PP = os.path.join(PR, "presentation")
RS = os.path.join(P2, "results")
CHK = []


def chk(item, value, status):
    CHK.append({"item": item, "value": str(value), "status": status})
    print("  %-5s %-56s %s" % (status, item, value), flush=True)


def num(s):
    return pd.to_numeric(s, errors="coerce").astype("float64")


def need(path):
    ok = os.path.exists(path)
    chk("found %s" % os.path.basename(path), "yes" if ok else "NO",
        "PASS" if ok else "FAIL")
    return ok


# ---------------- inputs ----------------
print("=== INPUTS ===")
files = [os.path.join(PP, f) for f in ("frame_P_outcomes.csv", "ehrP.csv",
                                       "ehrP_sens.csv", "scores_P.csv")]
files += [os.path.join(PR, "frame_B.csv")]
files += [os.path.join(IN, f) for f in ("demographics.csv",
                                        "followup_365.csv",
                                        "mimic_labels_harmonised.csv")]
ok = all([need(f) for f in files])
adm_f = [os.path.join(IN, "index_admission_times.csv"),
         os.path.join(DS, "Dissertation", "figwork", "data",
                      "index_admission_times.csv")]
adm_f = [f for f in adm_f if os.path.exists(f)]
chk("found index_admission_times.csv", adm_f[0] if adm_f else "NO",
    "PASS" if adm_f else "FAIL")
if not (ok and adm_f):
    pd.DataFrame(CHK).to_csv(os.path.join(RS, "table1_checks.csv"),
                             index=False)
    sys.exit("STOP: an input is missing. see: table1_checks.csv")


def rd(path):
    d = pd.read_csv(path)
    d = d.drop_duplicates("hadm_id")
    d["hadm_id"] = d["hadm_id"].astype("int64")
    return d.set_index("hadm_id")


F = rd(files[0])
ids = F.index
P = rd(files[1]).reindex(ids)
S = rd(files[2]).reindex(ids)
SC = rd(files[3]).reindex(ids)
B = rd(files[4])
DM = rd(files[5]).reindex(ids)
FU = rd(files[6]).reindex(ids)
LB = rd(files[7]).reindex(ids)
AT = rd(adm_f[0]).reindex(ids)

# ---------------- the five cohorts ----------------
died30 = (num(F["event_death"]) == 1) & (num(F["time_death"]) <= 30.5)
COLS = {
    "All admissions": pd.Series(True, index=ids),
    "Died within 30 days": died30,
    "Alive at 30 days": ~died30,
    "With a presentation CTPA": num(F["in_three_mod_p"]) == 1,
    "Discharged alive": pd.Series(ids.isin(B.index), index=ids)}
DIS = list(COLS)[-1]
for k, m in COLS.items():
    chk("column '%s': admissions" % k, int(m.sum()), "INFO")
chk("discharge cohort inside the at-diagnosis cohort",
    "%d of %d" % (int(B.index.isin(ids).sum()), len(B)),
    "PASS" if B.index.isin(ids).all() else "FAIL")


def fmt_n(n, d):
    if d == 0:
        return "—"
    if 0 < n < 10:
        return "<10"
    return "{:,} ({:.1f})".format(int(n), 100.0 * n / d)


def fmt_q(v, dig=None):
    v = v[np.isfinite(v)]
    if len(v) < 10:
        return "—"
    a, b, c = np.percentile(v, [25, 50, 75])
    if dig is None:
        dig = 0 if b >= 100 else 1 if b >= 1 else 2
    f = "{:.%df}" % dig
    return (f + " (" + f + "–" + f + ")").format(b, a, c)


ROWS = []


def band(title):
    ROWS.append(("band", title, None, None))


IDX30 = (1, 2)   # the died / alive at 30 days columns


def cat(label, x, cols=None, denom=None):
    """x: boolean Series (NaN = unknown, left out of the denominator)."""
    out = []
    for k, m in COLS.items():
        if cols is not None and k not in cols:
            out.append("—")
            continue
        mm = m.copy()
        if denom is not None:
            mm &= denom
        xx = x[mm]
        xx = xx[xx.notna()]
        out.append(fmt_n(int(xx.astype(bool).sum()), len(xx)))
    # the two 30-day columns add up to the first: if one is under 10,
    # the other would give it away, so it is not shown either
    a, b = IDX30
    if out[a] == "<10" and out[b] not in ("<10", "—"):
        out[b] = "†"
    elif out[b] == "<10" and out[a] not in ("<10", "—"):
        out[a] = "†"
    ROWS.append(("row", label, out, ""))


def cont(label, v, avail=False, cols=None, dig=None):
    v = num(v)
    out = []
    for k, m in COLS.items():
        if cols is not None and k not in cols:
            out.append("—")
            continue
        out.append(fmt_q(v[m].values, dig))
    av = "{:.1f}".format(100 * v.notna().mean()) if avail else ""
    ROWS.append(("row", label, out, av))


# ---------------- demographics ----------------
band("Demographics")
ROWS.append(("row", "Admissions, n",
              ["{:,}".format(int(m.sum())) for m in COLS.values()], ""))
pat = F["subject_id"]
ROWS.append(("row", "Patients, n",
              ["{:,}".format(pat[m].nunique()) for m in COLS.values()],
              ""))
cont("Age, years, median (IQR)", P["age"], avail=True, dig=0)
g = DM["gender"].astype(str).str.upper().str[0] if "gender" in DM else \
    None
fem = (g == "F") if g is not None else (num(P["sex_f"]) == 1)
cat("Female sex, n (%)", fem)
race_col = next((c for c in DM.columns if c.lower() in ("race",
                                                        "ethnicity")),
                None)
if race_col:
    r = DM[race_col].astype(str).str.upper()
    grp = np.select(
        [r.str.startswith("WHITE") | r.str.startswith("PORTUGUESE"),
         r.str.startswith("BLACK"),
         r.str.contains("HISPANIC|LATINO|SOUTH AMERICAN"),
         r.str.startswith("ASIAN"),
         r.str.contains("UNKNOWN|UNABLE|DECLINED|NAN")],
        ["White", "Black", "Hispanic or Latino", "Asian", "Unknown"],
        "Other")
    rg = pd.Series(grp, index=ids)
    for lab in ("White", "Black", "Hispanic or Latino", "Asian", "Other",
                "Unknown"):
        cat("Race and ethnicity: %s, n (%%)" % lab, rg == lab)
    exp = os.path.join(RS, "demographics_summary.csv")
    if os.path.exists(exp):
        e = pd.read_csv(exp)
        e = e[(e["frame"] == "A") & (e["cohort"] == "EHR+ECG") &
              (e["variable"] == "race_group")].set_index("group")["n"]
        got = rg.value_counts()
        same = all(int(got.get(k, 0)) == int(v) for k, v in e.items())
        chk("race groups match the subgroup analysis (stage 22a)",
            "; ".join("%s %d" % (k, int(got.get(k, 0))) for k in e.index),
            "PASS" if same else "WARN")
    for lab_, n_ in r[rg == "Other"].value_counts().items():
        chk("race label grouped as Other: %s" % lab_,
            "<10" if n_ < 10 else int(n_), "INFO")
else:
    chk("race column in demographics.csv", "not found", "WARN")

# ---------------- presentation ----------------
band("Vital signs at presentation, median (IQR)")
VIT = [("hr", "Heart rate, beats/min", 0),
       ("sbp", "Systolic BP, mmHg", 0),
       ("dbp", "Diastolic BP, mmHg", 0),
       ("rr", "Respiratory rate, /min", 0),
       ("spo2", "Oxygen saturation, %", 0),
       ("temp", "Temperature", 1)]
for k, lab, dg in VIT:
    c = k + "_first"
    if c not in P.columns:
        chk("vital sign column %s" % c, "not found", "WARN")
        continue
    if k == "temp":
        med = num(P[c]).median()
        lab = lab + (", °F" if med > 50 else ", °C")
    cont(lab, P[c], avail=True, dig=dg)
cat("No vital signs at presentation, n (%)",
    num(P["no_pres_vitals"]) == 1)
band("Laboratory tests at presentation, median (IQR)")
LAB = [("creatinine", None, "Creatinine, mg/dL", 1),
       ("bun", None, "Urea nitrogen, mg/dL", 0),
       ("sodium", None, "Sodium, mEq/L", 0),
       ("potassium", None, "Potassium, mEq/L", 1),
       ("bicarb", None, "Bicarbonate, mEq/L", 0),
       ("glucose", None, "Glucose, mg/dL", 0),
       ("hgb", None, "Haemoglobin, g/dL", 1),
       ("wbc", None, "White cell count, K/µL", 1),
       ("plt", None, "Platelets, K/µL", 0),
       (None, 50813, "Lactate, mmol/L", 1),
       (None, 51003, "Troponin T, ng/mL", 2),
       (None, 50963, "NT-proBNP, pg/mL", 0)]
for nm, iid, lab, dg in LAB:
    if nm is not None:
        c = nm + "_first" if nm + "_first" in P.columns else None
    else:
        hit = [x for x in P.columns
               if re.search(r"(^|_)%d(_|$)" % iid, x)
               and x.endswith("_first")]
        c = hit[0] if hit else None
    chk("lab column for %s" % lab.split(",")[0], c if c else "not found",
        "PASS" if c else "WARN")
    if c:
        cont(lab, P[c], avail=True, dig=dg)

# ---------------- history and scores ----------------
band("History as used at diagnosis (chronic conditions from earlier or "
     "index-admission codes)")
npr = num(P["n_prior"])
cat("Any earlier admission, n (%)", npr > 0)


def flag(D, *cols):
    have = [c for c in cols if c in D.columns]
    if not have:
        return None
    return (D[have].apply(num).fillna(0).max(axis=1) > 0)


for lab, cs in (("Cancer (any, incl. metastatic)", ("hx_ch_canc",
                                                    "hx_ch_metacanc")),
                ("Heart failure", ("hx_ch_chf",)),
                ("Chronic lung disease", ("hx_ch_cpd",)),
                ("Renal disease", ("hx_ch_rend",)),
                ("Diabetes", ("hx_ch_diab", "hx_ch_diabwc"))):
    f = flag(P, *cs)
    chk("history columns for %s" % lab, "found" if f is not None else
        "not found", "PASS" if f is not None else "WARN")
    if f is not None:
        cat(lab + ", n (%)", f)
if "hx_charlson_score" in P.columns:
    cont("Charlson score, median (IQR)", P["hx_charlson_score"], dig=0)
band("History from earlier admissions only (strict sensitivity "
     "version)")
for lab, cs in (("Cancer (any, incl. metastatic)", ("hx_ch_canc",
                                                    "hx_ch_metacanc")),
                ("Heart failure", ("hx_ch_chf",)),
                ("Chronic lung disease", ("hx_ch_cpd",))):
    f = flag(S, *cs)
    if f is not None:
        cat(lab + ", n (%)", f)
band("Risk scores at presentation")
sp = num(SC["spesi"])
cat("sPESI 0 (low risk), n (%)", sp == 0)
cat("sPESI 1 or more, n (%)", sp >= 1)
cont("PESI-10, median (IQR)", SC["pesi10"], dig=0)
cat("PESI-10 above 85 (class III–V), n (%)", num(SC["pesi10"]) > 85)
band("Modalities at presentation")
cat("Vital signs recorded, n (%)", num(F["has_pres_vitals"]) == 1)
cat("ECG within 24 h of admission, n (%)",
    num(F["has_pres_ecg"]) == 1)
cat("CTPA report at presentation, n (%)", num(F["in_three_mod_p"]) == 1)

# ---------------- index admission ----------------
band("Index admission")
a0 = pd.to_datetime(AT["admittime"], errors="coerce")
d0 = pd.to_datetime(AT["dischtime"], errors="coerce") if \
    "dischtime" in AT.columns else None
if d0 is not None:
    los = (d0 - a0).dt.total_seconds() / 86400
    cont("Length of stay, days, median (IQR)", los, dig=0)
    chk("length of stay available", "{:.1f}%".format(
        100 * los.notna().mean()), "INFO")
ih = LB["inhosp_days"].notna()
cat("Died in hospital, n (%)", ih)
loc = FU["discharge_location"].astype(str).str.upper().str.strip()
alive = ~ih
DEST = [("Home", loc == "HOME"),
        ("Home with home health care", loc == "HOME HEALTH CARE"),
        ("Skilled nursing facility", loc == "SKILLED NURSING FACILITY"),
        ("Rehabilitation", loc == "REHAB"),
        ("Hospice", loc == "HOSPICE")]
known = np.zeros(len(loc), dtype=bool)
for lab, m in DEST:
    known |= m.values
    cat("Discharged to: %s, n (%%)" % lab, m, denom=alive)
cat("Discharged to: other or unknown, n (%)",
    pd.Series(~known, index=ids), denom=alive)


# ---------------- outcomes ----------------
def ev_by(D, kind, h):
    cu = h + 0.5
    if kind == "death":
        return (num(D["event_death"]) == 1) & (num(D["time_death"]) <= cu)
    col = "event_cr" if kind == "cv" else "event_comp"
    val = 1
    return (num(D[col]) == val) & (num(D["time_first"]) <= cu)


band("Outcomes from admission (at-diagnosis frame)")
OUTS = (("death", "Death"), ("cv", "CV readmission before death"),
        ("comp", "Composite (CV readmission or death)"))
notB = [k for k in COLS if k != DIS]
for kind, lab in OUTS:
    for h in (30, 365):
        cat("%s by %d days, n (%%)" % (lab, h), ev_by(F, kind, h),
            cols=notB)
band("Outcomes from discharge (discharge frame)")
Bf = B.reindex(ids)
for kind, lab in OUTS:
    for h in (30, 365):
        x = ev_by(Bf, kind, h).where(COLS[DIS])
        cat("%s by %d days, n (%%)" % (lab, h), x, cols=[DIS])

# ---------------- write ----------------
print("\n=== WRITE ===")
heads = ["Characteristic"] + [
    "%s\n(n = {:,})".format(int(m.sum())) % k for k, m in COLS.items()]
heads += ["Data available, % of all admissions"]
flat = []
for kind, lab, vals, av in ROWS:
    if kind == "band":
        flat.append([lab] + [""] * (len(heads) - 1))
    else:
        flat.append([lab] + vals + [av])
pd.DataFrame(flat, columns=[h.replace("\n", " ") for h in heads]).to_csv(
    os.path.join(RS, "table1.csv"), index=False, encoding="utf-8-sig")
TITLE = ("Table 1. Characteristics of the pulmonary embolism admissions "
         "at presentation and at discharge")
NOTE = (
    "Values are n (%) or median (IQR). Presentation values are the first "
    "recorded from 24 h before to 24 h after admission; vital signs come "
    "from the ED or the ICU and are missing where neither recorded them. "
    "History counts chronic conditions coded on an earlier or the index "
    "admission, as the at-diagnosis model does; the strict rows count "
    "earlier admissions only. sPESI and PESI-10 use the at-diagnosis "
    "history, and PESI-10 is PESI without altered mental status. The ECG "
    "row counts an ECG from 24 h before to 24 h after admission. "
    "Discharge destination percentages are of patients discharged alive. "
    "Outcomes from admission count from the admission date; outcomes "
    "from discharge, for the discharge cohort only, count from the "
    "discharge date. CV, cardiovascular. Counts of 1 to 9 are shown as "
    "<10; † marks a count not shown because it would reveal a count under "
    "10 in the other 30-day column.")
wb = Workbook()
ws = wb.active
ws.title = "Table 1"
ws.sheet_view.showGridLines = False
nc = len(heads)
thin = Side(style="thin", color="FF000000")
grey = Side(style="thin", color="FFBFBFBF")
med = Side(style="medium", color="FF000000")
dbl = Side(style="double", color="FF000000")
cen = Alignment(horizontal="center", vertical="center", wrap_text=True)
f10 = Font(name="Calibri", size=10)
f10b = Font(name="Calibri", size=10, bold=True)
ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=nc)
c = ws.cell(row=1, column=1, value=TITLE)
c.font = Font(name="Calibri", size=12, bold=True)
c.alignment = Alignment(horizontal="left", vertical="center")
for j, h in enumerate(heads, start=1):
    c = ws.cell(row=2, column=j, value=h)
    c.font = f10b
    c.fill = PatternFill("solid", fgColor="FFFFFFFF")
    c.alignment = cen
    c.border = Border(left=thin, right=thin, top=thin, bottom=dbl)
ws.row_dimensions[2].height = 31.2
r = 3
first_band = True
blocks = []
for kind, lab, vals, av in ROWS:
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
        for j, v in enumerate([lab] + vals + [av], start=1):
            c = ws.cell(row=r, column=j, value=v)
            c.font = f10
            c.alignment = cen
        blocks[-1][1].append(r)
    r += 1
for _, rows in blocks:
    for k, rr in enumerate(rows):
        top = med if k == 0 else None
        bot = med if k == len(rows) - 1 else grey
        for j in range(1, nc + 1):
            ws.cell(row=rr, column=j).border = Border(
                left=thin, right=thin, top=top, bottom=bot)
ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=nc)
c = ws.cell(row=r, column=1, value=NOTE)
c.font = Font(name="Calibri", size=9, italic=True)
c.alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
ws.row_dimensions[r].height = 95
ws.column_dimensions["A"].width = 40
for j in range(2, nc + 1):
    ws.column_dimensions[chr(64 + j)].width = 17
xl = os.path.join(RS, "table1.xlsx")
wb.save(xl)
chk("rows written (bands)", "%d (%d)" % (
    sum(1 for x in ROWS if x[0] == "row"), len(blocks)), "INFO")
ck = pd.DataFrame(CHK)
ck.to_csv(os.path.join(RS, "table1_checks.csv"), index=False)
print("\n  checks: %d PASS, %d WARN, %d FAIL" % (
    (ck.status == "PASS").sum(), (ck.status == "WARN").sum(),
    (ck.status == "FAIL").sum()))
print("\nOUTPUT FILES")
for f in ("table1.xlsx", "table1.csv", "table1_checks.csv"):
    print("  ", os.path.join(RS, f))
print("  ", os.path.join(RS, "table1_log.txt"), " (terminal log)")
print("\nOutputs: table1.csv, table1_checks.csv")