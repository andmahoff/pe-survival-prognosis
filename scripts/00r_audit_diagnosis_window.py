"""00r_audit_diagnosis_window.py
Audit of the at-diagnosis window (A32, D17, D18).

Checks, block by block, whether the information left out of the
at-diagnosis model lies outside a reasonable diagnostic period, and
whether anything it kept does not. Results are printed and saved as
aggregates (shares, medians, correlations; counts of 1 to 9 are shown as
<10); the patient-level rows pulled from BigQuery stay in
data/processed/audit.
1. Admission times: real admittime and ED arrival (edregtime) from
   MIMIC-IV; deaths on the day of admission or the day after.
2. ECG: real acquisition times from the first source with clock times
   (BigQuery MIMIC-IV-ECG record_list, then raw WFDB .hea headers); hours
   from admission; how many ECGs dropped by D18 (dated the day after
   admission or later) were within 24 h; first-ECG coverage under clock
   windows, including the side project's -12 h to +48 h.
3. Labs: for each lab the models use, coverage at presentation (24 h
   before to 24 h after admission), within 48 h and over the stay; share
   of measurements in the first 24 and 48 h; rank correlation of the first
   value with the stay mean and with the last value.
4. Vital signs: when the first heart rate from ED triage, ED vital signs
   and ICU charts falls, and coverage by window; which definition
   reproduces the side project's ehr28 heart rate.
5. History: each Charlson and Elixhauser category among discharge-cohort
   admissions, from earlier admissions only (as in the model) and from
   the index admission's codes; chronic, or possibly acquired in the stay.
6. Kept features: every at-diagnosis EHR column belongs to a known
   presentation group; n_prior and n_codes count earlier admissions.

Inputs:  data/processed/presentation/{frame_P,ehrP,ehrP_sens}.csv,
         data/processed/discharge/ehrplusB.csv, data/inputs/{ehr28,
         ehr_extra,pres_labs_first,ecg_timing}.csv, results/
         ehr_extra_dictionary.csv (if present), BigQuery (MIMIC-IV);
         side project: data/raw/ecg/ecg_record_index.csv and raw ECG
         headers (if present)
Outputs: results/diagwin_summary.csv, diagwin_labs.csv,
         diagwin_history.csv, diagwin_ecg.csv, diagwin_log.txt;
         data/processed/audit/*.csv (patient level, not in the repo)
Usage:   python scripts/00r_audit_diagnosis_window.py
Notes:   set GCP_PROJECT to the BigQuery billing project. The script
         stops if an input file or a MIMIC-IV dataset is missing.
"""
import os
import re
import sys
import time

import numpy as np
import pandas as pd

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
SP = os.path.join(DS, "SideProject")
IN = os.path.join(P2, "data", "inputs")
PR = os.path.join(P2, "data", "processed")
PP = os.path.join(PR, "presentation")
AU = os.path.join(PR, "audit")
RS = os.path.join(P2, "results")
os.makedirs(AU, exist_ok=True)
PROJECT = os.environ["GCP_PROJECT"]
HOSP = ["physionet-data.mimiciv_hosp",
        "physionet-data.mimiciv_3_1_hosp"]
ICU = ["physionet-data.mimiciv_icu",
       "physionet-data.mimiciv_3_1_icu"]
ED = ["physionet-data.mimiciv_ed"]
ECGDS = ["physionet-data.mimiciv_ecg",
         "physionet-data.mimic_iv_ecg"]
# MIMIC-IV item IDs of the 18 registered labs (ehr28's labs)
REG = {"creatinine": 50912, "sodium": 50983, "potassium": 50971,
       "bun": 51006, "glucose": 50931, "calcium": 50893,
       "bicarb": 50882, "chloride": 50902, "hct": 51221,
       "plt": 51265, "wbc": 51301, "hgb": 51222, "aniongap": 50868,
       "rbc": 51279, "mchc": 51249, "mch": 51248, "mcv": 51250,
       "rdw": 51277}
# categories that can be coded for an event during the index stay
ACUTE = {"mi", "cevd", "carit", "pcd", "fed", "coag", "blane",
         "wloss", "para", "ond"}
t_start = time.time()
OUT = []


def put(section, item, value, status="INFO"):
    OUT.append({"section": section, "item": item,
                "value": str(value), "status": status})
    print("  %-5s %-58s %s" % (status, item, value), flush=True)


def head(t):
    print("\n" + "=" * 70 + "\n" + t + "\n" + "=" * 70, flush=True)


def num(s):
    return pd.to_numeric(s, errors="coerce").astype("float64")


def pct(x):
    return "%.1f%%" % (100 * x) if np.isfinite(x) else "n/a"


def cnt(n):
    """Counts of 1 to 9 shown as <10 (MIMIC small-cell rule)."""
    n = int(n)
    return "<10" if 1 <= n <= 9 else "%d" % n


def q3(v):
    v = np.asarray(v, dtype="float64")
    v = v[np.isfinite(v)]
    if not len(v):
        return "n/a"
    a, b, c = np.percentile(v, [25, 50, 75])
    return "%.1f (%.1f to %.1f)" % (b, a, c)


def save_all():
    pd.DataFrame(OUT).to_csv(os.path.join(RS, "diagwin_summary.csv"),
                             index=False)


# ---------------- 0. inputs ----------------
head("0. INPUTS")
need = [os.path.join(PP, f) for f in ("frame_P.csv", "ehrP.csv",
                                      "ehrP_sens.csv")]
need += [os.path.join(IN, f) for f in ("ehr28.csv", "ehr_extra.csv",
                                       "pres_labs_first.csv")]
need += [os.path.join(SP, "data", "raw", "ecg", "ecg_record_index.csv")]
miss = [f for f in need if not os.path.exists(f)]
put("inputs", "input files found",
    "all" if not miss else os.path.basename(miss[0]),
    "PASS" if not miss else "FAIL")
if miss:
    save_all()
    sys.exit(1)
fp = pd.read_csv(os.path.join(PP, "frame_P.csv"))
fp["hadm_id"] = fp["hadm_id"].astype("int64")
fp["subject_id"] = fp["subject_id"].astype("int64")
HADM = sorted(fp["hadm_id"].unique().tolist())
SUBJ = sorted(fp["subject_id"].unique().tolist())
put("inputs", "at-diagnosis admissions / patients",
    "%d / %d" % (len(HADM), len(SUBJ)))
E = pd.read_csv(os.path.join(SP, "data", "raw", "ecg",
                             "ecg_record_index.csv"))
E = E[E["hadm_id"].isin(HADM)].drop_duplicates(["rec", "hadm_id"])
E["rec"] = E["rec"].astype(str).str.replace("\\", "/", regex=False)
E["rec"] = E["rec"].str.strip("/")
E["rec"] = np.where(E["rec"].str.startswith("files/"), E["rec"],
                    "files/" + E["rec"])
put("inputs", "ECG records in the cohort", len(E))
pl = pd.read_csv(os.path.join(IN, "pres_labs_first.csv"), nrows=1)
items = dict((k, v) for k, v in REG.items())
for c in pl.columns:
    m = re.match(r"lab_(\d+)_(.+)_first$", c)
    if m:
        items[m.group(2)] = int(m.group(1))
put("inputs", "labs audited (registered + extra)", len(items))

# ---------------- BigQuery: datasets, probe ----------------
head("1. BIGQUERY: DATASETS AND DRY RUNS")
from google.cloud import bigquery  # noqa: E402

client = bigquery.Client(project=PROJECT)


def cfg(dry=False, **arrays):
    p = [bigquery.ArrayQueryParameter(k, "INT64", v)
         for k, v in arrays.items()]
    return bigquery.QueryJobConfig(dry_run=dry, use_query_cache=not dry,
                                   query_parameters=p)


def pick(cands, table):
    for ds in cands:
        sql = "SELECT 1 FROM `%s.%s` LIMIT 1" % (ds, table)
        try:
            client.query(sql, job_config=cfg(True))
            return ds
        except Exception as e:
            msg = str(e).splitlines()[0][:90]
            print("    %s.%s not usable: %s" % (ds, table, msg))
    return None


H = pick(HOSP, "admissions")
I = pick(ICU, "chartevents")
D = pick(ED, "edstays")
G = pick(ECGDS, "record_list")
put("bigquery", "hosp / icu / ed datasets", "%s / %s / %s" % (H, I, D),
    "PASS" if H and I and D else "FAIL")
put("bigquery", "MIMIC-IV-ECG record_list on BigQuery",
    G if G else "not available (headers on this computer tried next)")
if not (H and I and D):
    save_all()
    sys.exit(1)
STUDY = sorted(num(E["ecg_study_id"]).dropna().astype("int64")
               .unique().tolist())
ITEMS = sorted(set(items.values()))

C_ADM = """c AS (
  SELECT subject_id, hadm_id, admittime, dischtime
  FROM `{H}.admissions` WHERE hadm_id IN UNNEST(@hadm))"""
Q = {}
Q["admissions"] = ("""
SELECT subject_id, hadm_id, admittime, edregtime, dischtime
FROM `{H}.admissions` WHERE hadm_id IN UNNEST(@hadm)
""".format(H=H), dict(hadm=HADM))
if G:
    Q["ecg"] = ("""
SELECT study_id, ecg_time FROM `{G}.record_list`
WHERE study_id IN UNNEST(@study)
""".format(G=G), dict(study=STUDY))
Q["labs"] = ("""
WITH {C},
l AS (
  SELECT c.hadm_id, le.itemid, le.valuenum AS v,
    DATETIME_DIFF(le.charttime, c.admittime, MINUTE) / 60.0 AS h
  FROM `{H}.labevents` le
  JOIN c ON le.subject_id = c.subject_id
  WHERE le.itemid IN UNNEST(@items) AND le.valuenum IS NOT NULL
    AND le.charttime >= DATETIME_SUB(c.admittime, INTERVAL 24 HOUR)
    AND le.charttime <= c.dischtime),
a AS (
  SELECT hadm_id, itemid, COUNT(*) AS n_all,
    COUNTIF(h <= 24) AS n24, COUNTIF(h <= 48) AS n48,
    ARRAY_AGG(v ORDER BY h LIMIT 1)[OFFSET(0)] AS v_first,
    ARRAY_AGG(v ORDER BY h DESC LIMIT 1)[OFFSET(0)] AS v_last,
    AVG(v) AS v_mean, MIN(h) AS h_first
  FROM l GROUP BY hadm_id, itemid),
r AS (
  SELECT *,
    PERCENT_RANK() OVER (PARTITION BY itemid ORDER BY v_first) AS rf,
    PERCENT_RANK() OVER (PARTITION BY itemid ORDER BY v_mean) AS rm,
    PERCENT_RANK() OVER (PARTITION BY itemid ORDER BY v_last) AS rl
  FROM a)
SELECT itemid, COUNT(*) AS adm_stay,
  COUNTIF(h_first <= 24) AS adm_24, COUNTIF(h_first <= 48) AS adm_48,
  SUM(n_all) AS meas_all, SUM(n24) AS meas_24, SUM(n48) AS meas_48,
  APPROX_QUANTILES(n_all, 2)[OFFSET(1)] AS med_n,
  APPROX_QUANTILES(h_first, 4)[OFFSET(1)] AS h_first_q25,
  APPROX_QUANTILES(h_first, 4)[OFFSET(2)] AS h_first_q50,
  APPROX_QUANTILES(h_first, 4)[OFFSET(3)] AS h_first_q75,
  CORR(rf, rm) AS rho_first_mean, CORR(rf, rl) AS rho_first_last,
  COUNTIF(n_all > n24) AS adm_values_after_24h
FROM r GROUP BY itemid
""".format(C=C_ADM.format(H=H), H=H), dict(hadm=HADM, items=ITEMS))
Q["heart_rate"] = ("""
WITH {C},
hr AS (
  SELECT e.subject_id, e.intime AS t,
    CAST(t.heartrate AS FLOAT64) AS v, 'triage' AS src
  FROM `{D}.triage` t JOIN `{D}.edstays` e ON t.stay_id = e.stay_id
  WHERE e.subject_id IN (SELECT subject_id FROM c)
  UNION ALL
  SELECT e.subject_id, vs.charttime,
    CAST(vs.heartrate AS FLOAT64), 'ed'
  FROM `{D}.vitalsign` vs JOIN `{D}.edstays` e
    ON vs.stay_id = e.stay_id
  WHERE e.subject_id IN (SELECT subject_id FROM c)
  UNION ALL
  SELECT ce.subject_id, ce.charttime,
    CAST(ce.valuenum AS FLOAT64), 'icu'
  FROM `{I}.chartevents` ce
  WHERE ce.itemid = 220045
    AND ce.subject_id IN (SELECT subject_id FROM c)),
j AS (
  SELECT c.hadm_id, hr.src, hr.v,
    DATETIME_DIFF(hr.t, c.admittime, MINUTE) / 60.0 AS h,
    (hr.t >= DATETIME_SUB(c.admittime, INTERVAL 24 HOUR)
     AND hr.t <= c.dischtime) AS in_idx
  FROM c JOIN hr ON hr.subject_id = c.subject_id
  WHERE hr.v BETWEEN 20 AND 300)
SELECT hadm_id,
  MIN(IF(in_idx, h, NULL)) AS h_first,
  MIN(IF(in_idx AND src = 'triage', h, NULL)) AS h_first_triage,
  MIN(IF(in_idx AND src = 'ed', h, NULL)) AS h_first_ed,
  MIN(IF(in_idx AND src = 'icu', h, NULL)) AS h_first_icu,
  AVG(IF(in_idx, v, NULL)) AS hr_index,
  AVG(v) AS hr_patient,
  AVG(IF(in_idx AND src != 'ed', v, NULL)) AS hr_index_no_edv,
  AVG(IF(src != 'ed', v, NULL)) AS hr_patient_no_edv,
  COUNTIF(in_idx) AS n_index, COUNTIF(NOT in_idx) AS n_outside
FROM j GROUP BY hadm_id
""".format(C=C_ADM.format(H=H), D=D, I=I), dict(hadm=HADM))

gb = {}
for k, (sql, arr) in Q.items():
    try:
        b = client.query(sql, job_config=cfg(True, **arr))
        gb[k] = b.total_bytes_processed / 1e9
        put("bigquery", "query '%s': GB to scan" % k, "%.2f" % gb[k],
            "PASS")
    except Exception as e:
        put("bigquery", "query '%s': dry run" % k,
            str(e).splitlines()[0][:120], "FAIL")
save_all()
if any(r["status"] == "FAIL" for r in OUT):
    print("\nSTOP: a dry run failed. see: diagwin_summary.csv")
    sys.exit(1)
put("bigquery", "total GB to scan (free tier 1,000 GB / month)",
    "%.1f" % sum(gb.values()))
print("  starting in 5 s (Ctrl+C to stop)...", flush=True)
time.sleep(5)
R = {}
for k, (sql, arr) in Q.items():
    t0 = time.time()
    R[k] = client.query(sql, job_config=cfg(False, **arr)).to_dataframe()
    put("bigquery", "query '%s': rows" % k,
        "%d (%.0f s)" % (len(R[k]), time.time() - t0))
    R[k].to_csv(os.path.join(AU, "bq_%s.csv" % k), index=False)

# ---------------- 2. admission times and early deaths ----------------
head("2. ADMISSION TIMES AND EARLY DEATHS")
A = R["admissions"].copy()
for c in ("admittime", "edregtime", "dischtime"):
    A[c] = pd.to_datetime(A[c], errors="coerce")
A["hadm_id"] = A["hadm_id"].astype("int64")
cov = np.mean(np.isin(HADM, A["hadm_id"]))
put("admission", "cohort admissions found in MIMIC-IV", pct(cov),
    "PASS" if cov >= 0.99 else "FAIL")
mid = (A["admittime"].dt.hour == 0) & (A["admittime"].dt.minute == 0)
put("admission", "real admission times at exactly midnight", pct(
    mid.mean()))
ed_h = (A["admittime"] - A["edregtime"]).dt.total_seconds() / 3600
put("admission", "admissions through the ED (edregtime present)",
    pct(A["edregtime"].notna().mean()))
put("admission", "hours from ED arrival to admission (median, IQR)",
    q3(ed_h))
d_day = num(fp["time_death"]) - 0.5
ev = num(fp["event_death"]) == 1
put("admission", "deaths on the day of admission / by the day after",
    "%s / %s of %s within a year" % (
        cnt((ev & (d_day <= 0)).sum()), cnt((ev & (d_day <= 1)).sum()),
        cnt(ev.sum())))
AD = A.set_index("hadm_id")

# ---------------- 3. ECG ----------------
head("3. ECG: REAL ACQUISITION TIMES")
E["t_date"] = pd.to_datetime(E["ecg_charttime"], errors="coerce")
E["hadm_id"] = E["hadm_id"].astype("int64")
src, T = None, None
if G is not None and "ecg" in R and len(R["ecg"]):
    g = R["ecg"].copy()
    g["study_id"] = g["study_id"].astype("int64")
    g["t"] = pd.to_datetime(g["ecg_time"], errors="coerce")
    e2 = E.assign(study_id=num(E["ecg_study_id"]))
    e2 = e2.merge(g[["study_id", "t"]], on="study_id", how="left")
    if e2["t"].notna().mean() >= 0.9:
        src, T = "BigQuery record_list", e2


def find_root():
    """Folder holding MIMIC-IV-ECG's files/pNNNN/... tree, if any."""
    probe = E["rec"].head(50).tolist()
    starts = [os.path.join(SP, "data", "raw", "ecg"),
              os.path.join(SP, "data", "raw"), os.path.join(SP, "data"),
              SP, DS]
    for s in starts:
        for r in probe:
            if os.path.exists(os.path.join(s, *r.split("/")) + ".hea"):
                return s
    for top in (SP, os.path.join(DS, "Dissertation")):
        for dirpath, dirnames, _ in os.walk(top):
            if "files" in dirnames:
                for r in probe:
                    f = os.path.join(dirpath, *r.split("/")) + ".hea"
                    if os.path.exists(f):
                        return dirpath
            dirnames[:] = [d for d in dirnames
                           if not re.match(r"^p\d{2,8}$", d)
                           and not d.startswith(".")
                           and d not in SKIP]
    return None


SKIP = {"venv", "venv_surv", "__pycache__", "node_modules", "Lib",
        "site-packages"}


def read_hea(path):
    try:
        with open(path, "r", errors="ignore") as fh:
            tok = fh.readline().split()
    except OSError:
        return None
    tm = next((x for x in tok if ":" in x), None)
    dt = next((x for x in tok if "/" in x), None)
    if tm is None or dt is None:
        return None
    return pd.to_datetime(dt + " " + tm.split(".")[0],
                          format="%d/%m/%Y %H:%M:%S", errors="coerce")


if src is None:
    root = find_root()
    put("ecg", "raw ECG headers (.hea) on this computer",
        root if root else "not found")
    if root:
        ts = [read_hea(os.path.join(root, *r.split("/")) + ".hea")
              for r in E["rec"]]
        T = E.assign(t=pd.to_datetime(pd.Series(ts, index=E.index)))
        if T["t"].notna().mean() >= 0.9:
            src = "raw WFDB headers (.hea)"
put("ecg", "source of real ECG times", src if src else "none found",
    "PASS" if src else "WARN")
ecg_rows = []
ET = os.path.join(IN, "ecg_timing.csv")
if os.path.exists(ET):
    et = pd.read_csv(ET)
    h = num(et["ecg_h_first"]).dropna()
    put("ecg", "ecg_timing.csv: first-ECG hours that are whole days",
        pct(np.mean(np.isclose(h / 24, np.round(h / 24)))))
if src:
    T = T[T["t"].notna()].copy()
    put("ecg", "records with a real time", "%d of %d (%s)" % (
        len(T), len(E), pct(len(T) / len(E))))
    same = (T["t"].dt.normalize() == T["t_date"].dt.normalize()).mean()
    put("ecg", "real date equals the date in ecg_record_index", pct(same),
        "PASS" if same >= 0.99 else "WARN")
    put("ecg", "real times at exactly midnight", pct(
        ((T["t"].dt.hour == 0) & (T["t"].dt.minute == 0)).mean()))
    T["admit"] = T["hadm_id"].map(AD["admittime"])
    T["edreg"] = T["hadm_id"].map(AD["edregtime"])
    T["h"] = (T["t"] - T["admit"]).dt.total_seconds() / 3600
    T["h_ed"] = (T["t"] - T["edreg"]).dt.total_seconds() / 3600
    T["day"] = (T["t_date"].dt.normalize()
                - T["admit"].dt.normalize()).dt.days
    T[["rec", "subject_id", "hadm_id", "t", "h", "h_ed", "day"]].to_csv(
        os.path.join(AU, "ecg_times_real.csv"), index=False)
    put("ecg", "hours from admission, all records (median, IQR)",
        q3(T["h"]))
    put("ecg", "records within the side project's -12 h to +48 h",
        pct(T["h"].between(-12, 48).mean()))
    kept = T["day"].between(-1, 0)
    put("ecg", "records D18 kept (dated day before / day of admission)",
        "%d (%s)" % (kept.sum(), pct(kept.mean())))
    drop = T[~kept & (T["day"] >= 1)]
    put("ecg", "records D18 dropped as dated after admission",
        cnt(len(drop)))
    put("ecg", "  of those, within 24 h of admission by the clock",
        pct(drop["h"].le(24).mean()) if len(drop) else "n/a")
    put("ecg", "  of those, within 12 h of admission by the clock",
        pct(drop["h"].le(12).mean()) if len(drop) else "n/a")
    put("ecg", "records D18 kept that were after 24 h by the clock",
        cnt((kept & T["h"].gt(24)).sum()))
    wins = {"D18 (dates: day before or day of admission)": None,
            "-24 h to +24 h": (-24, 24), "-12 h to +24 h": (-12, 24),
            "-24 h to +12 h": (-24, 12), "-24 h to +48 h": (-24, 48),
            "-12 h to +48 h (side project)": (-12, 48)}
    for nm, w in wins.items():
        sel = kept if w is None else T["h"].between(w[0], w[1])
        hz = T[sel].groupby("hadm_id")["h"].min()
        ecg_rows.append({
            "window": nm,
            "admissions_with_ecg": len(hz),
            "share_of_cohort": round(len(hz) / len(HADM), 4),
            "first_ecg_hours_median": round(float(hz.median()), 2),
            "first_ecg_hours_q75": round(float(hz.quantile(.75)), 2),
            "records_used": int(sel.sum())})
    fc = T[T["h"].between(-24, 24)].sort_values("h").groupby(
        "hadm_id").head(1).set_index("hadm_id")["rec"]
    fd = T[kept].sort_values(["day", "h"]).groupby("hadm_id").head(1)
    fd = fd.set_index("hadm_id")["rec"]
    both = fc.index.intersection(fd.index)
    put("ecg", "admissions whose clock-first ECG (-24 h to +24 h) "
        "differs from D18's", "%s of %s" % (
            cnt((fc[both] != fd[both]).sum()), cnt(len(both))))
    print(pd.DataFrame(ecg_rows).to_string(index=False), flush=True)
else:
    put("ecg", "next step", "copy mimic_pe_mace_cohort.csv (real "
        "ecg_charttime) from the cluster, or supply this file", "WARN")
pd.DataFrame(ecg_rows).to_csv(os.path.join(RS, "diagwin_ecg.csv"),
                              index=False)

# ---------------- 4. labs ----------------
head("4. LABS: WHEN THE INFORMATION IS RECORDED")
L = R["labs"].copy()
L["itemid"] = L["itemid"].astype("int64")
inv = {v: k for k, v in items.items()}
L.insert(1, "lab", L["itemid"].map(inv))
n = float(len(HADM))
for c in L.columns[2:]:
    L[c] = num(L[c])
L["cov_presentation_24h"] = L["adm_24"] / n
L["cov_48h"] = L["adm_48"] / n
L["cov_stay"] = L["adm_stay"] / n
L["share_meas_first_24h"] = L["meas_24"] / L["meas_all"]
L["share_meas_first_48h"] = L["meas_48"] / L["meas_all"]
L["share_adm_values_after_24h"] = L["adm_values_after_24h"] / L[
    "adm_stay"]
L["registered"] = L["lab"].isin(REG.keys())
cols = ["itemid", "lab", "registered", "cov_presentation_24h", "cov_48h",
        "cov_stay", "med_n", "h_first_q50", "share_meas_first_24h",
        "share_meas_first_48h", "share_adm_values_after_24h",
        "rho_first_mean", "rho_first_last"]
L = L[cols].sort_values(["registered", "cov_stay"],
                        ascending=[False, False])
L.round(4).to_csv(os.path.join(RS, "diagwin_labs.csv"), index=False)
for grp, d in (("registered 18", L[L["registered"]]),
               ("extra", L[~L["registered"]])):
    put("labs", "%s: coverage at presentation / 48 h / stay (median)"
        % grp, "%s / %s / %s" % (
            pct(d["cov_presentation_24h"].median()),
            pct(d["cov_48h"].median()), pct(d["cov_stay"].median())))
    put("labs", "%s: share of stay measurements in first 24 / 48 h"
        % grp, "%s / %s" % (pct(d["share_meas_first_24h"].median()),
                           pct(d["share_meas_first_48h"].median())))
    put("labs", "%s: first value vs stay mean / vs last (rank r)" % grp,
        "%.2f / %.2f" % (d["rho_first_mean"].median(),
                         d["rho_first_last"].median()))
gain = L["cov_48h"] - L["cov_presentation_24h"]
put("labs", "labs whose coverage rises by 10+ points from 24 h to 48 h",
    ", ".join(L.loc[gain >= 0.10, "lab"].astype(str)) or "none")

# ---------------- 5. vital signs ----------------
head("5. VITAL SIGNS (HEART RATE)")
V = R["heart_rate"].copy()
V["hadm_id"] = V["hadm_id"].astype("int64")
for c in V.columns[1:]:
    V[c] = num(V[c])
V = pd.DataFrame({"hadm_id": HADM}).merge(V, on="hadm_id", how="left")
hf = V["h_first"]
for nm, lo, hi in (("presentation (-24 h to +24 h)", -24, 24),
                   ("within 48 h", -24, 48),
                   ("any time in the stay", -24, 1e9)):
    put("vitals", "heart rate recorded: %s" % nm,
        pct(hf.between(lo, hi).mean()))
for s in ("triage", "ed", "icu"):
    put("vitals", "first heart rate within 24 h from %s" % s,
        pct(V["h_first_%s" % s].between(-24, 24).mean()))
put("vitals", "admissions whose patient has heart rates outside the "
    "index stay", pct((V["n_outside"] > 0).mean()))
e28 = pd.read_csv(os.path.join(IN, "ehr28.csv"))
e28 = e28.drop_duplicates("hadm_id").set_index("hadm_id")
if "hr" in e28.columns:
    x = num(e28["hr"].reindex(V["hadm_id"])).values
    for c in ("hr_index", "hr_patient", "hr_index_no_edv",
              "hr_patient_no_edv"):
        y = V[c].values
        ok = np.isfinite(x) & np.isfinite(y)
        put("vitals", "ehr28 heart rate reproduced by %s (within 0.5)"
            % c, "%s of %d" % (pct(np.mean(np.abs(x[ok] - y[ok]) <= .5)),
                              ok.sum()))
    ok = np.isfinite(x)
    put("vitals", "ehr28 heart rate present / present without any "
        "index-stay heart rate", "%s / %s" % (
            pct(ok.mean()), cnt((ok & ~np.isfinite(V["hr_index"])).sum())))

# ---------------- 6. history ----------------
head("6. HISTORY: EARLIER ADMISSIONS vs THE INDEX ADMISSION")
P = pd.read_csv(os.path.join(PP, "ehrP.csv")).set_index("hadm_id")
S = pd.read_csv(os.path.join(PP, "ehrP_sens.csv")).set_index("hadm_id")
bf = os.path.join(PR, "discharge", "ehrplusB.csv")
hist = []
if os.path.exists(bf):
    B = pd.read_csv(bf).set_index("hadm_id")
    ids = B.index.intersection(P.index)
    put("history", "discharge-cohort admissions compared", len(ids))
    for c in P.columns:
        m = re.match(r"hx_(ch|el)_(.+)$", c)
        if not m or (m.group(1) + "_" + m.group(2)) not in B.columns:
            continue
        k = m.group(2)
        e_ = num(P.loc[ids, c]).mean()
        i_ = num(B.loc[ids, m.group(1) + "_" + k]).mean()
        hist.append({
            "index": "Charlson" if m.group(1) == "ch" else "Elixhauser",
            "category": k,
            "type": "possibly acquired during the stay" if k in ACUTE
            else "chronic history",
            "prev_earlier_only": round(e_, 4),
            "prev_index_codes": round(i_, 4),
            "index_minus_earlier": round(i_ - e_, 4)})
H_ = pd.DataFrame(hist)
H_.to_csv(os.path.join(RS, "diagwin_history.csv"), index=False)
if len(H_):
    for ty, d in H_.groupby("type"):
        put("history", "%s: categories / median prevalence earlier vs "
            "index" % ty, "%d / %s vs %s" % (
                len(d), pct(d["prev_earlier_only"].median()),
                pct(d["prev_index_codes"].median())))
ex = pd.read_csv(os.path.join(IN, "ehr_extra.csv"))
ex = ex.drop_duplicates("hadm_id").set_index("hadm_id")
np0 = num(ex["n_prior"].reindex(HADM))
put("history", "admissions with no earlier admission (history unseen)",
    pct((np0 == 0).mean()))
for k, nm in (("cancer", "canc"), ("heart_failure", "chf"),
              ("copd", "cpd")):
    if k in e28.columns:
        ehr28v = num(e28[k].reindex(HADM)).mean()
        ear = num(P.get("hx_ch_" + nm, pd.Series(dtype=float)).reindex(
            HADM)).mean()
        sen = num(S.get("hx_ch_" + nm, pd.Series(dtype=float)).reindex(
            HADM)).mean()
        put("history", "%s: ehr28 / earlier only / earlier or index" % k,
            "%s / %s / %s" % (pct(ehr28v), pct(ear), pct(sen)))

# ---------------- 7. kept features ----------------
head("7. FEATURES THE AT-DIAGNOSIS MODEL KEPT")
grp = {"vital (first in window)": r"^(hr|sbp|dbp|rr|spo2|temp)_first$",
       "lab (first in window)": r"_first$",
       "history (earlier admissions)": r"^hx_",
       "demographics": r"^(age|sex_f)$",
       "missing-vitals indicator": r"^no_pres_vitals$",
       "earlier admissions / codes": r"^(n_prior|n_codes)$"}
left = list(P.columns)
for g, rx in grp.items():
    hit = [c for c in left if re.search(rx, c)]
    left = [c for c in left if c not in hit]
    put("kept", g, len(hit))
put("kept", "columns in no known presentation group",
    ", ".join(left) if left else "none", "PASS" if not left else "WARN")
dic = os.path.join(RS, "ehr_extra_dictionary.csv")
if os.path.exists(dic):
    dd = pd.read_csv(dic).set_index("column")
    for c in ("n_prior", "n_codes"):
        if c in dd.index:
            put("kept", "%s is defined as" % c, dd.loc[c, "label"],
                "PASS" if "prior" in str(dd.loc[c, "label"]).lower() or
                "before" in str(dd.loc[c, "label"]).lower() else "WARN")
save_all()
print("\n  elapsed %.1f min" % ((time.time() - t_start) / 60))
print("\nOUTPUT FILES")
for f in ("diagwin_summary.csv", "diagwin_labs.csv",
          "diagwin_history.csv", "diagwin_ecg.csv"):
    print("  ", os.path.join(RS, f))
print("  ", os.path.join(RS, "diagwin_log.txt"),
      " (terminal log)")
print("  ", AU, " (patient-level BigQuery pulls, local only)")
print("\nOutputs: diagwin_summary.csv, diagwin_labs.csv,",
      "diagwin_history.csv, diagwin_ecg.csv")