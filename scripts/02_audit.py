"""02_audit.py
Event-time audit of the MIMIC and INSPECT data, before any frame is built.

Descriptive only; no model is fitted. It answers four questions.
  A. Are the MIMIC event times coherent? Negative times, in-hospital
     deaths without a dod, cardiovascular readmissions dated inside the
     index stay or after death, same-day ties, ECG timing.
  B. Does the event coding reproduce the dissertation label file at day
     30, row by row? (FAIL if not.)
  C. How many events are there at each horizon (10, 30, 90, 180 and 365
     days) in three frames? A_admission_all: clock at admission, every
     index admission. A_admission_first: clock at admission, first
     admission per patient. B_discharge_first: clock at discharge, first
     admission per patient, alive at discharge. Cardiovascular
     readmission is counted only to day 30, where the label file stops
     (03_bq_followup.py extends it).
  D. Does INSPECT's native time-to-event mortality (tte in minutes, with a
     censoring flag) support survival analysis? Censoring pattern,
     Kaplan-Meier mortality per horizon beside MIMIC's, whether early
     censored patients differ, and whether its deaths agree with the
     dissertation labels.

Inputs:  data/inputs/{mimic_labels_harmonised,index_admission_times,
         cohorts,ecg_timing,inspect_ehr}.csv; copied into
         data/inputs/inspect/ (MD5 in manifest_inspect.csv): side project
         data/raw/inspect/{labels,study_mapping}_*.tsv and
         data/processed/inspect_imp_features_v2.csv (only to check the
         index-scan choice), dissertation figwork/data/
         inspect_labels_final.csv
Outputs: results/audit_checks.csv, results/audit_horizons.csv,
         results/audit_log.txt; data/inputs/inspect/ (not in the repo)
Usage:   python scripts/02_audit.py
Notes:   needs scikit-survival (Kaplan-Meier estimator).
"""
import os
import glob
import shutil
import hashlib
import numpy as np
import pandas as pd
from sksurv.nonparametric import (
    kaplan_meier_estimator)

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
SIDE = os.path.join(DS, "SideProject")
FIG = os.path.join(DS, "Dissertation", "figwork",
                   "data")
INP = os.path.join(P2, "data", "inputs")
INS = os.path.join(INP, "inspect")
RES = os.path.join(P2, "results")
RAWI = os.path.join(SIDE, "data", "raw",
                    "inspect")
for p_ in (INS, RES):
    os.makedirs(p_, exist_ok=True)

KEY = ["subject_id", "hadm_id"]
H = [10, 30, 90, 180, 365]
GAP_D = 30   # index scan within 30 d of PE date
CHK, HOR, MAN = [], [], []
MIM, KMI = {}, {}


def banner(t):
    print("", flush=True)
    print("=" * 70)
    print(t)
    print("=" * 70, flush=True)


def rec(sec, coh, item, value, status="INFO"):
    CHK.append({"section": sec, "cohort": coh,
                "item": item, "value": value,
                "status": status})
    print("  %-4s %-9s %-52s %s"
          % (status, coh, item, value),
          flush=True)


def fl(sec, coh, item, mask, bad=True):
    v = int(mask.sum())
    rec(sec, coh, item, v,
        "WARN" if (bad and v) else "INFO")


def hor(ds, frame, coh, endpoint, h, n, ev,
        risk, method):
    HOR.append({"dataset": ds, "frame": frame,
                "cohort": coh,
                "endpoint": endpoint,
                "horizon_days": h, "n": n,
                "events": ev,
                "risk": round(float(risk), 4),
                "method": method})


def flt(s):
    return pd.to_numeric(
        s, errors="coerce").astype("float64")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20),
                      b""):
            h.update(b)
    return h.hexdigest()


def latest(folder, prefix):
    hits = sorted(glob.glob(os.path.join(
        folder, prefix + "*.tsv")))
    return hits[-1] if hits else None


def counts(ds, frame, coh, n, t_d, t_cv, hs,
           with_cv):
    """Events per horizon. A CV readmission and a
    death on the same day count as CV, which is
    the DISSERTATION rule."""
    for h in hs:
        de = (t_d >= 0) & (t_d <= h)
        hor(ds, frame, coh, "death", h, n,
            int(de.sum()), de.sum() / n,
            "empirical")
        MIM[(coh, frame, h)] = de.sum() / n
        if not with_cv:
            continue
        cv = (t_cv >= 0) & (t_cv <= h)
        dfs = de & (~cv | (t_d < t_cv))
        cfs = cv & ~dfs
        for nm, m_ in (("composite", de | cv),
                       ("cv_first_competing", cfs),
                       ("death_first", dfs)):
            hor(ds, frame, coh, nm, h, n,
                int(m_.sum()), m_.sum() / n,
                "empirical")


def first_adm(df):
    return (df.sort_values(
        ["subject_id", "admittime", "hadm_id"],
        kind="mergesort")
        .drop_duplicates("subject_id"))


# ================= A. MIMIC =================
banner("A. MIMIC EVENT TIMES (checks on the"
       " EHR+ECG cohort)")
lab = pd.read_csv(os.path.join(
    INP, "mimic_labels_harmonised.csv"))
for c in ("dod_days", "inhosp_days", "cv_days",
          "death_days"):
    lab[c] = flt(lab[c])
adm = pd.read_csv(os.path.join(
    INP, "index_admission_times.csv"))
adm = adm.drop_duplicates("hadm_id").copy()
at = pd.to_datetime(adm["admittime"],
                    errors="coerce")
dt = pd.to_datetime(adm["dischtime"],
                    errors="coerce")
adm["admittime"] = at
adm["dischtime"] = dt
adm["los_d"] = flt((dt.dt.normalize()
                    - at.dt.normalize()).dt.days)
adm["los_h"] = ((dt - at).dt.total_seconds()
                / 3600.0)
coh = pd.read_csv(os.path.join(INP, "cohorts.csv"))
ecgt = pd.read_csv(os.path.join(
    INP, "ecg_timing.csv"))

d = coh.merge(lab, on=KEY, how="left")
d = d.merge(adm[["hadm_id", "admittime",
                 "dischtime", "los_d", "los_h"]],
            on="hadm_id", how="left")
d = d.merge(ecgt[["hadm_id", "ecg_h_last"]],
            on="hadm_id", how="left")
d["t_death"] = np.fmin(d["dod_days"],
                       d["inhosp_days"])
E = d[d["in_ehr_ecg"] == 1].copy()
T = d[d["in_three_mod"] == 1].copy()
C = "EHR+ECG"

rec("A", C, "admissions", len(E))
rec("A", C, "patients", E["subject_id"].nunique())
rec("A", "3-mod", "admissions", len(T))
rec("A", "3-mod", "patients",
    T["subject_id"].nunique())
nl = int(E["death_30d"].isna().sum())
rec("A", C, "admissions with no label row", nl,
    "PASS" if nl == 0 else "FAIL")
na_ = int(E["admittime"].isna().sum())
nd_ = int(E["dischtime"].isna().sum())
rec("A", C, "missing admittime / dischtime",
    "%d / %d" % (na_, nd_),
    "PASS" if na_ + nd_ == 0 else "WARN")
rec("A", C, "median length of stay (days)",
    float(E["los_d"].median()))

oh = E["inhosp_days"].isna()
fl("A", C, "negative length of stay",
   E["los_d"] < 0)
fl("A", C, "dod_days < 0", E["dod_days"] < 0)
fl("A", C, "inhosp_days < 0",
   E["inhosp_days"] < 0)
fl("A", C, "cv_days < 0", E["cv_days"] < 0)
fl("A", C, "in-hospital death with no dod",
   ~oh & E["dod_days"].isna())
fl("A", C, "dod and in-hospital death >1 d apart",
   (E["dod_days"] - E["inhosp_days"]).abs() > 1)
fl("A", C, "in-hospital death >1 d from discharge",
   (E["inhosp_days"] - E["los_d"]).abs() > 1)
fl("A", C, "dod BEFORE discharge, not in-hospital",
   oh & (E["dod_days"] < E["los_d"]))
fl("A", C, "dod ON discharge date, not in-hosp",
   oh & (E["dod_days"] == E["los_d"]), bad=False)
fl("A", C, "CV readmission on admission day",
   E["cv_days"] == 0)
fl("A", C, "CV readmission dated inside stay",
   E["cv_days"] < E["los_d"])
fl("A", C, "CV readmission after death",
   E["cv_days"] > E["t_death"])
fl("A", C, "CV readmission and death same day",
   E["cv_days"] == E["t_death"], bad=False)
fl("A", C, "deaths on day 0 or 1",
   E["t_death"] <= 1, bad=False)
rec("A", C, "any dod recorded",
    int(E["dod_days"].notna().sum()))
rec("A", C, "latest dod (days after admission)",
    float(E["dod_days"].max()))
fl("A", C, "deaths >365 d after admission",
   E["dod_days"] > 365, bad=False)
fl("A", C, "deaths >365 d after index discharge",
   E["dod_days"] > E["los_d"] + 365, bad=False)
fl("A", C, "last in-window ECG after admission",
   E["ecg_h_last"] > 0, bad=False)
fl("A", C, "last in-window ECG >24 h in",
   E["ecg_h_last"] > 24, bad=False)
fl("A", C, "last in-window ECG after discharge",
   E["ecg_h_last"] > E["los_h"], bad=False)

# ================= B. CODING =================
banner("B. P2 EVENT CODING vs DISSERTATION LABELS"
       " AT DAY 30 (row by row)")
td, cv = E["t_death"], E["cv_days"]
de30 = (td >= 0) & (td <= 30)
cv30 = (cv >= 0) & (cv <= 30)
df30 = de30 & (~cv30 | (td < cv))
cf30 = cv30 & ~df30
for nm, mine in (("death_30d", de30),
                 ("death_first", df30),
                 ("cv_first", cf30),
                 ("composite_30d", de30 | cv30)):
    ref = flt(E[nm]).fillna(0).astype(int).values
    nb = int((mine.astype(int).values
              != ref).sum())
    rec("B", C, "%s: rows that disagree" % nm, nb,
        "PASS" if nb == 0 else "FAIL")

# ================= C. HORIZONS =================
banner("C. EVENTS PER HORIZON")
for cn, df_ in (("EHR+ECG", E), ("3-mod", T)):
    n = len(df_)
    counts("MIMIC", "A_admission_all", cn, n,
           df_["t_death"], df_["cv_days"],
           [10, 30], True)
    counts("MIMIC", "A_admission_all", cn, n,
           df_["t_death"], None,
           [90, 180, 365], False)
    f1 = first_adm(df_)
    n1 = len(f1)
    counts("MIMIC", "A_admission_first", cn, n1,
           f1["t_death"], f1["cv_days"],
           [10, 30], True)
    counts("MIMIC", "A_admission_first", cn, n1,
           f1["t_death"], None,
           [90, 180, 365], False)
    died = f1["inhosp_days"].notna()
    nodis = f1["dischtime"].isna()
    early = ~died & (f1["dod_days"] < f1["los_d"])
    b = f1[~died & ~nodis & ~early]
    rec("C", cn, "first admissions", n1)
    rec("C", cn, "Frame B excl: died in hospital",
        int(died.sum()))
    rec("C", cn, "Frame B excl: no discharge time",
        int(nodis.sum()))
    rec("C", cn, "Frame B excl: dod before discharge",
        int(early.sum()),
        "WARN" if early.sum() else "INFO")
    rec("C", cn, "Frame B cohort (alive at discharge)",
        len(b))
    tb = b["dod_days"] - b["los_d"]
    counts("MIMIC", "B_discharge_first", cn, len(b),
           tb, None, H, False)
    rec("C", cn, "conditional cv_first cohort, day 30",
        int((df_["death_first"] == 0).sum()))

s = E.sort_values(["subject_id", "admittime",
                   "hadm_id"],
                  kind="mergesort").copy()
s["nxt"] = s.groupby("subject_id")[
    "admittime"].shift(-1)
gp = flt((s["nxt"].dt.normalize()
          - s["dischtime"].dt.normalize()).dt.days)
k_ = E.groupby("subject_id").size()
rec("C", C, "patients with 2+ index admissions",
    int((k_ >= 2).sum()))
rec("C", C, "most index admissions, one patient",
    int(k_.max()))
fl("C", C, "next PE admission <=30 d after disch",
   gp <= 30, bad=False)
fl("C", C, "next PE admission <=365 d after disch",
   gp <= 365, bad=False)


# ================= D. INSPECT =================
def inspect_section():
    banner("D. INSPECT NATIVE TIME-TO-EVENT"
           " MORTALITY")
    src = {
        "labels": latest(RAWI, "labels_"),
        "study_mapping": latest(
            RAWI, "study_mapping_"),
        "inspect_labels_final": os.path.join(
            FIG, "inspect_labels_final.csv"),
        "inspect_imp_features_v2": os.path.join(
            SIDE, "data", "processed",
            "inspect_imp_features_v2.csv")}
    need = ("labels", "study_mapping",
            "inspect_labels_final")
    loc = {}
    for k, p in src.items():
        if not p or not os.path.exists(p):
            rec("D", "INSPECT",
                "input not found: " + k, str(p),
                "FAIL" if k in need else "WARN")
            continue
        q = os.path.join(INS, os.path.basename(p))
        shutil.copy2(p, q)
        ok = md5(p) == md5(q)
        MAN.append({"file": os.path.basename(p),
                    "source": p, "md5": md5(q),
                    "identical": int(ok)})
        rec("D", "INSPECT",
            "copied " + os.path.basename(p),
            "identical" if ok else "MISMATCH",
            "PASS" if ok else "FAIL")
        loc[k] = q
    if not all(k in loc for k in need):
        return

    il = pd.read_csv(loc["labels"], sep="\t")
    if "person_id" in il.columns:
        il = il.drop(columns=["person_id"])
    mp = pd.read_csv(loc["study_mapping"], sep="\t")
    mp = mp[["impression_id", "person_id",
             "procedure_DATETIME"]]
    hl = pd.read_csv(loc["inspect_labels_final"])
    hl = hl[[c for c in ("person_id", "pe_date",
                         "days_to_death")
             if c in hl.columns]].copy()
    ie = pd.read_csv(os.path.join(
        INP, "inspect_ehr.csv"))
    rec("D", "INSPECT", "label rows / columns",
        "%d / %d" % il.shape)

    tte = [c for c in il.columns
           if c.lower().startswith("tte")
           and "mort" in c.lower()]
    cen = [c for c in il.columns
           if "censor" in c.lower()
           and "mort" in c.lower()]
    rec("D", "INSPECT", "tte mortality column",
        ",".join(tte) or "none")
    rec("D", "INSPECT", "censoring mortality column",
        ",".join(cen) or "none")
    if len(tte) != 1 or len(cen) != 1:
        rec("D", "INSPECT",
            "exactly one tte and one censor column",
            "%d / %d" % (len(tte), len(cen)),
            "FAIL")
        print("  label columns:", list(il.columns))
        return
    tte, cen = tte[0], cen[0]

    # index scan: closest to the PE date, as in
    # the side project; checked against its file
    g = il.merge(mp, on="impression_id",
                 how="inner")
    g["pdt"] = pd.to_datetime(
        g["procedure_DATETIME"], errors="coerce")
    hl["pdate"] = pd.to_datetime(
        hl["pe_date"], errors="coerce")
    g = g[g["person_id"].isin(
        set(ie["person_id"]))]
    g = g.merge(hl, on="person_id", how="left")
    g["gap"] = (g["pdt"] - g["pdate"]).abs()
    g = g[g["gap"] <= pd.Timedelta(days=GAP_D)]
    g = (g.sort_values(["person_id", "gap"],
                       kind="mergesort")
         .drop_duplicates("person_id"))
    rec("D", "INSPECT",
        "persons with an index scan (of %d)"
        % len(ie), len(g),
        "PASS" if len(g) == len(ie) else "WARN")

    if ("inspect_imp_features_v2" in loc
            and "1_month_mortality" in g.columns):
        sp = pd.read_csv(
            loc["inspect_imp_features_v2"],
            usecols=["person_id",
                     "1_month_mortality"])
        m_ = g[["person_id",
                "1_month_mortality"]].merge(
            sp, on="person_id",
            suffixes=("", "_sp"))
        a_ = m_["1_month_mortality"].astype(
            str).str.lower()
        b_ = m_["1_month_mortality_sp"].astype(
            str).str.lower()
        nd = int((a_ != b_).sum())
        rec("D", "INSPECT",
            "matched to side-project index file",
            len(m_))
        rec("D", "INSPECT",
            "1-month label disagreements vs side",
            nd, "PASS" if nd == 0 else "WARN")

    g["tte_d"] = flt(g[tte]) / 1440.0
    cz = g[cen].astype(str).str.strip().str.lower()
    g["cens"] = cz.isin(["true", "1", "1.0", "yes"])
    bad = g["tte_d"].isna() | (g["tte_d"] < 0)
    rec("D", "INSPECT", "tte missing or negative",
        int(bad.sum()),
        "WARN" if bad.sum() else "INFO")
    k = g[~bad].copy()
    ev = (~k["cens"]).to_numpy(dtype=bool)
    td = k["tte_d"].to_numpy(dtype=float)
    n = len(k)
    rec("D", "INSPECT", "persons analysed", n)
    rec("D", "INSPECT", "deaths observed",
        int(ev.sum()))
    rec("D", "INSPECT", "censored",
        int((~ev).sum()))
    if ev.any():
        rec("D", "INSPECT",
            "median days to death (deaths)",
            round(float(np.median(td[ev])), 1))
    if (~ev).any():
        rec("D", "INSPECT",
            "median days to censoring (censored)",
            round(float(np.median(td[~ev])), 1))
    xr, sr = kaplan_meier_estimator(
        ev, td, reverse=True)
    lo = np.where(sr <= 0.5)[0]
    mfu = float(xr[lo[0]]) if len(lo) else np.nan
    rec("D", "INSPECT",
        "median follow-up, reverse KM (days)",
        round(mfu, 1))

    # how the native labels relate to tte
    for lb, cuts in (
            ("1_month_mortality", (30, 31)),
            ("6_month_mortality", (180, 182, 183)),
            ("12_month_mortality", (365, 366))):
        if lb not in k.columns:
            continue
        nat = (k[lb].astype(str).str.strip()
               .str.lower().to_numpy())
        vc = pd.Series(nat).value_counts().to_dict()
        rec("D", "INSPECT", "native %s" % lb, vc)
        for cut in cuts:
            dv = np.where(
                ev & (td <= cut), "true",
                np.where(td > cut, "false",
                         "censored"))
            rec("D", "INSPECT",
                "  agreement with tte cut at %d d"
                % cut,
                round(float((dv == nat).mean()), 4))

    x, sv = kaplan_meier_estimator(ev, td)
    for h in H:
        i = np.searchsorted(x, h, side="right") - 1
        s_h = 1.0 if i < 0 else float(sv[i])
        evh = int((ev & (td <= h)).sum())
        ceh = int((~ev & (td <= h)).sum())
        alive = int((td > h).sum())
        KMI[h] = 1.0 - s_h
        hor("INSPECT", "scan_origin", "PE persons",
            "death", h, n, evh, 1.0 - s_h,
            "kaplan_meier")
        if evh + alive:
            hor("INSPECT", "scan_origin",
                "PE persons", "death", h,
                evh + alive, evh,
                evh / (evh + alive),
                "binary_censored_dropped")
        rec("D", "INSPECT",
            "censored on or before day %d" % h,
            "%d (%.1f%%)" % (ceh, 100.0 * ceh / n))

    # are early-censored patients different?
    kk = k.merge(ie[["person_id", "age", "cancer",
                     "heart_failure"]],
                 on="person_id", how="left")
    ce = (kk["cens"]
          & (kk["tte_d"] <= 365)).to_numpy()
    rec("D", "INSPECT", "censored within a year",
        int(ce.sum()))
    for v in ("age", "cancer", "heart_failure"):
        a = flt(kk.loc[ce, v]).mean()
        b2 = flt(kk.loc[~ce, v]).mean()
        rec("D", "INSPECT",
            "mean %s: censored <1y vs rest" % v,
            "%.3f vs %.3f" % (a, b2))

    # same deaths as the dissertation labels?
    if "days_to_death" in kk.columns:
        hd = flt(kk["days_to_death"]).to_numpy()
        for h in (30, 365):
            a_ = (hd >= 0) & (hd <= h)
            b_ = ev & (td <= h)
            rec("D", "INSPECT",
                "deaths <=%d d: both/DISS only/"
                "native only" % h,
                "%d / %d / %d"
                % ((a_ & b_).sum(), (a_ & ~b_).sum(),
                   (~a_ & b_).sum()))

    m30 = MIM.get(("EHR+ECG", "A_admission_first",
                   30))
    m365 = MIM.get(("EHR+ECG", "A_admission_first",
                    365))
    if m30 and m365 and KMI.get(30):
        rec("E", "compare",
            "1-year / 30-day mortality, MIMIC vs"
            " INSPECT",
            "%.2f vs %.2f" % (m365 / m30,
                              KMI[365] / KMI[30]))


inspect_section()

# ================= SUMMARY =================
banner("SUMMARY")
ck = pd.DataFrame(CHK)
hz = pd.DataFrame(HOR)
ck.to_csv(os.path.join(RES, "audit_checks.csv"),
          index=False)
hz.to_csv(os.path.join(RES, "audit_horizons.csv"),
          index=False)
if MAN:
    pd.DataFrame(MAN).to_csv(os.path.join(
        INS, "manifest_inspect.csv"), index=False)

print(hz.pivot_table(
    index=["dataset", "frame", "cohort",
           "endpoint", "method"],
    columns="horizon_days", values="risk")
    .round(4).to_string())
nf = int((ck["status"] == "FAIL").sum())
nw = int((ck["status"] == "WARN").sum())
print("")
print("FAIL %d   WARN %d" % (nf, nw))
if nf + nw:
    print(ck[ck["status"].isin(["FAIL", "WARN"])]
          [["section", "cohort", "item", "value",
            "status"]].to_string(index=False))
print("")
print("OUTPUT FILES")
print("  " + os.path.join(RES, "audit_checks.csv"))
print("  " + os.path.join(RES, "audit_horizons.csv"))
print("  " + os.path.join(RES, "audit_log.txt"),
      "(terminal log)")
print("  " + INS + "  (patient-level, local only)")
print("")
print("Outputs: audit_checks.csv and"
      " audit_horizons.csv")