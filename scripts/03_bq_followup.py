"""03_bq_followup.py
One-year follow-up pulled from MIMIC-IV on BigQuery.

Extends the dissertation's 30-day labels. For every cohort admission it
pulls the days to the first later admission with a cardiovascular code
(cv_days_all, no upper bound; MI I21/I22/410, heart failure I50/428,
cardiac arrest I46/4275, AF I48/4273, ischaemic stroke I63/433/434), the
discharge location (for the hospice sensitivity analysis), the anchor year
and an approximate calendar year (for the end-of-data truncation check),
and the days to the patient's last MIMIC discharge (deaths are recorded
for up to one year after it). The codes and the readmission rule are those
of the dissertation's script 80: any different hadm_id admitted on or
after the index admission date. dod_days is pulled again as a cross-check.

Checks: one row per cohort admission; discharge dates match
index_admission_times.csv; dod_days equals the label file's; cv_days_all
capped at day 30 reproduces the label file's cv_days row by row. If any
fails, no patient-level file is written.

Inputs:  data/inputs/cohorts.csv, mimic_labels_harmonised.csv,
         index_admission_times.csv; BigQuery (MIMIC-IV 3.1)
Outputs: data/inputs/followup_365.csv (patient level, not in the repo),
         results/bq_followup_checks.csv, results/bq_followup_horizons.csv,
         results/bq_followup_log.txt
Usage:   GCP_PROJECT=<billing project> python scripts/03_bq_followup.py
Notes:   one BigQuery query; run with google-cloud-bigquery 3.25.0 and
         db-dtypes 1.2.0.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from google.cloud import bigquery

T0 = time.time()
pd.set_option("display.width", 200)

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
RES = os.path.join(P2, "results")
os.makedirs(RES, exist_ok=True)
PROJECT = os.environ["GCP_PROJECT"]
HOSP = "physionet-data.mimiciv_3_1_hosp"
H = [10, 30, 90, 180, 365]
CHK, HOR = [], []


def banner(t):
    print("", flush=True)
    print("=" * 70)
    print(t)
    print("=" * 70, flush=True)


def rec(sec, item, value, status="INFO"):
    CHK.append({"section": sec, "item": item,
                "value": value, "status": status})
    print("  %-4s %-54s %s"
          % (status, item, value), flush=True)


def hor(coh, endpoint, h, n, ev):
    HOR.append({"dataset": "MIMIC",
                "frame": "B_discharge_first",
                "cohort": coh,
                "endpoint": endpoint,
                "horizon_days": h, "n": n,
                "events": int(ev),
                "risk": round(float(ev) / n, 4)
                if n else np.nan})


def flt(s):
    return pd.to_numeric(
        s, errors="coerce").astype("float64")


def write_checks():
    pd.DataFrame(CHK).to_csv(os.path.join(
        RES, "bq_followup_checks.csv"), index=False)
    pd.DataFrame(HOR).to_csv(os.path.join(
        RES, "bq_followup_horizons.csv"),
        index=False)


# ============ 0. ACCESS ============
banner("0. BIGQUERY ACCESS")
try:
    client = bigquery.Client(project=PROJECT)
    list(client.query("SELECT 1 AS ok").result())
    print("  access OK, project " + PROJECT)
except Exception as e:
    print("  BigQuery access FAILED:")
    print("  " + repr(e)[:300])
    print("  Fix, in a plain PowerShell window:")
    print("    gcloud auth application-default login")
    print("    (sign in with the account linked to PhysioNet)")
    print("    gcloud auth application-default"
          " set-quota-project " + PROJECT)
    sys.exit(1)

# ============ 1. PULL ============
banner("1. PULL")
coh = pd.read_csv(os.path.join(INP, "cohorts.csv"))
lab = pd.read_csv(os.path.join(
    INP, "mimic_labels_harmonised.csv"))
adm = pd.read_csv(os.path.join(
    INP, "index_admission_times.csv"))
hadms = ",".join(str(int(x)) for x in
                 sorted(coh["hadm_id"].unique()))
print("  cohort admissions:",
      coh["hadm_id"].nunique(), flush=True)

SQL = """
WITH idx AS (
  SELECT subject_id, hadm_id, admittime,
         dischtime, discharge_location,
         admission_type
  FROM `{hosp}.admissions`
  WHERE hadm_id IN ({h})
),
mace_adm AS (
  SELECT DISTINCT a.subject_id, a.hadm_id,
         a.admittime
  FROM `{hosp}.admissions` a
  JOIN `{hosp}.diagnoses_icd` d
    ON a.hadm_id = d.hadm_id
  WHERE (d.icd_version = 10
         AND SUBSTR(d.icd_code, 1, 3) IN
             ('I21','I22','I50','I46',
              'I48','I63'))
     OR (d.icd_version = 9
         AND (SUBSTR(d.icd_code, 1, 3) IN
              ('410','428','433','434')
           OR SUBSTR(d.icd_code, 1, 4) IN
              ('4275','4273')))
),
cv AS (
  SELECT i.hadm_id AS ih,
         MIN(DATE_DIFF(DATE(m.admittime),
             DATE(i.admittime), DAY)) AS cvd
  FROM idx i
  JOIN mace_adm m
    ON m.subject_id = i.subject_id
   AND m.hadm_id != i.hadm_id
   AND DATE_DIFF(DATE(m.admittime),
       DATE(i.admittime), DAY) >= 0
  GROUP BY 1
),
lastd AS (
  SELECT subject_id,
         MAX(dischtime) AS last_disch
  FROM `{hosp}.admissions`
  WHERE subject_id IN
        (SELECT subject_id FROM idx)
  GROUP BY 1
)
SELECT
  i.subject_id, i.hadm_id, i.dischtime,
  i.discharge_location, i.admission_type,
  DATE_DIFF(DATE(p.dod), DATE(i.admittime),
            DAY) AS dod_days_bq,
  p.anchor_year, p.anchor_year_group,
  EXTRACT(YEAR FROM i.admittime)
      AS admit_year_shifted,
  cv.cvd AS cv_days_all,
  DATE_DIFF(DATE(l.last_disch),
            DATE(i.admittime), DAY)
      AS last_disch_days
FROM idx i
LEFT JOIN `{hosp}.patients` p
  ON p.subject_id = i.subject_id
LEFT JOIN cv ON cv.ih = i.hadm_id
LEFT JOIN lastd l
  ON l.subject_id = i.subject_id
""".format(hosp=HOSP, h=hadms)

print("  querying ...", flush=True)
t = client.query(SQL).to_dataframe()
print("  rows returned:", len(t), flush=True)
for c in ("subject_id", "hadm_id"):
    t[c] = pd.to_numeric(t[c]).astype("int64")
for c in ("dod_days_bq", "cv_days_all",
          "last_disch_days", "anchor_year",
          "admit_year_shifted"):
    t[c] = flt(t[c])
t["dischtime"] = pd.to_datetime(
    t["dischtime"], errors="coerce")

# ============ 2. HARD CHECKS ============
banner("2. HARD CHECKS")
nu = t["hadm_id"].nunique()
nc = coh["hadm_id"].nunique()
rec("check", "rows / unique admissions / cohort",
    "%d / %d / %d" % (len(t), nu, nc),
    "PASS" if len(t) == nu == nc else "FAIL")

a2 = adm.drop_duplicates("hadm_id").copy()
a2["at_local"] = pd.to_datetime(
    a2["admittime"], errors="coerce")
a2["dt_local"] = pd.to_datetime(
    a2["dischtime"], errors="coerce")
m = t.merge(a2[["hadm_id", "at_local",
                "dt_local"]],
            on="hadm_id", how="left")
bd = int((m["dischtime"].dt.normalize()
          != m["dt_local"].dt.normalize()).sum())
rec("check", "discharge dates that differ", bd,
    "PASS" if bd == 0 else "FAIL")

lb = lab[["hadm_id", "dod_days", "cv_days",
          "inhosp_days"]].drop_duplicates(
    "hadm_id")
m = m.merge(lb, on="hadm_id", how="left")
for c in ("dod_days", "cv_days", "inhosp_days"):
    m[c] = flt(m[c])
same = ((m["dod_days"].isna()
         & m["dod_days_bq"].isna())
        | (m["dod_days"] == m["dod_days_bq"]))
nd = int((~same).sum())
rec("check", "dod_days rows that disagree", nd,
    "PASS" if nd == 0 else "FAIL")
cap = m["cv_days_all"].where(
    m["cv_days_all"] <= 30)
same = ((cap.isna() & m["cv_days"].isna())
        | (cap == m["cv_days"]))
nv = int((~same).sum())
rec("check", "cv_days capped at 30: rows disagree",
    nv, "PASS" if nv == 0 else "FAIL")

nf = sum(1 for r in CHK if r["status"] == "FAIL")
if nf:
    write_checks()
    banner("STOPPED: %d hard check(s) failed, no"
           " patient data written" % nf)
    sys.exit(1)

# ============ 3. DESCRIPTION ============
banner("3. FOLLOW-UP DESCRIPTION")
m = m.merge(coh[["hadm_id", "in_three_mod"]],
            on="hadm_id", how="left")
m["los_d"] = flt((m["dt_local"].dt.normalize()
                  - m["at_local"].dt.normalize())
                 .dt.days)
g = m["anchor_year_group"].astype(str)
off = m["admit_year_shifted"] - m["anchor_year"]
m["approx_year_lo"] = flt(g.str[:4]) + off
m["approx_year_hi"] = flt(g.str[-4:]) + off
data_end = int(np.nanmax(flt(g.str[-4:])))
rec("info", "latest anchor year group ends",
    data_end)
for k_, v_ in (m["anchor_year_group"]
               .value_counts().sort_index()
               .items()):
    rec("info", "anchor_year_group " + str(k_),
        int(v_))

rec("info", "any later MACE admission (any time)",
    int(m["cv_days_all"].notna().sum()))
x = int((m["cv_days_all"] < m["los_d"]).sum())
rec("info", "MACE admission dated inside stay", x,
    "WARN" if x else "INFO")
x = int((m["last_disch_days"] < m["los_d"]).sum())
rec("info", "last MIMIC discharge before index",
    x, "WARN" if x else "INFO")
died = m["discharge_location"].astype(
    str).str.upper().eq("DIED")
x = int((died != m["inhosp_days"].notna()).sum())
rec("info", "'DIED' location vs in-hospital death"
    " mismatch", x, "WARN" if x else "INFO")

# Frame B: first admission, alive at discharge
m = m.sort_values(["subject_id", "at_local",
                   "hadm_id"], kind="mergesort")
f1 = m.drop_duplicates("subject_id")
ok = (f1["inhosp_days"].isna()
      & ~(f1["dod_days"] < f1["los_d"]))
B = f1[ok].copy()
B["tB_death"] = B["dod_days"] - B["los_d"]
B["tB_cv"] = B["cv_days_all"] - B["los_d"]
x = int((B["tB_cv"] < 0).sum())
rec("frameB", "CV readmission before discharge", x,
    "WARN" if x else "INFO")

for cn, df_ in (("EHR+ECG", B),
                ("3-mod", B[B["in_three_mod"]
                            == 1])):
    n = len(df_)
    rec("frameB", cn + " cohort", n)
    hs = df_["discharge_location"].astype(
        str).str.upper().str.contains("HOSPICE")
    rec("frameB", cn + " discharged to hospice",
        int(hs.sum()))
    td, tc = df_["tB_death"], df_["tB_cv"]
    for h in H:
        de = (td >= 0) & (td <= h)
        cv = (tc >= 0) & (tc <= h)
        dfs = de & (~cv | (td < tc))
        cfs = cv & ~dfs
        for nm, mk in (("death", de),
                       ("composite", de | cv),
                       ("cv_first_competing", cfs),
                       ("death_first", dfs)):
            hor(cn, nm, h, n, mk.sum())

print("")
print("  discharge locations, EHR+ECG Frame B:")
for k_, v_ in (B["discharge_location"]
               .fillna("missing")
               .value_counts().head(12).items()):
    rec("frameB", "  location " + str(k_), int(v_))

# end-of-data truncation check
late = B["approx_year_hi"] >= data_end
early = B["approx_year_lo"] >= data_end
rec("trunc", "Frame B, possibly in final year",
    int(late.sum()))
rec("trunc", "Frame B, certainly in final year",
    int(early.sum()))
for lbl, sel in (("possibly final year", late),
                 ("earlier", ~late)):
    s = B[sel]
    n = len(s)
    if n == 0:
        continue
    d365 = ((s["tB_death"] >= 0)
            & (s["tB_death"] <= 365)).sum()
    c365 = ((s["tB_cv"] >= 0)
            & (s["tB_cv"] <= 365)).sum()
    rec("trunc", "%s: n / 1y death / 1y CV" % lbl,
        "%d / %.4f / %.4f"
        % (n, d365 / n, c365 / n))

# ============ 4. WRITE ============
banner("4. WRITING")
out = m[["subject_id", "hadm_id",
         "discharge_location", "admission_type",
         "anchor_year", "anchor_year_group",
         "admit_year_shifted", "approx_year_lo",
         "approx_year_hi", "cv_days_all",
         "last_disch_days", "dod_days_bq"]]
fp = os.path.join(INP, "followup_365.csv")
out.to_csv(fp, index=False)
print("  wrote", fp, out.shape)
write_checks()

hz = pd.DataFrame(HOR)
print("")
print(hz.pivot_table(
    index=["cohort", "endpoint"],
    columns="horizon_days", values="events")
    .astype(int).to_string())
print("")
print("elapsed %.1f s" % (time.time() - T0))
print("")
print("OUTPUT FILES")
print("  " + fp + "  (patient-level, local only)")
print("  " + os.path.join(
    RES, "bq_followup_checks.csv"))
print("  " + os.path.join(
    RES, "bq_followup_horizons.csv"))
print("  " + os.path.join(
    RES, "bq_followup_log.txt"),
    "(terminal log)")
print("")
print("Outputs: bq_followup_checks.csv and"
      " bq_followup_horizons.csv")