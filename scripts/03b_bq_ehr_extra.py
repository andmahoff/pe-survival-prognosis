"""03b_bq_ehr_extra.py
Extra MIMIC EHR features from BigQuery for every admission in ehr28.csv.

Pulls:
1. SpO2, built like ehr28's vitals: ICU chartevents (itemid 220277), then
   ED triage (o2sat) where still missing.
2. Presentation vitals for sPESI (hr_first, sbp_first, spo2_first):
   triage of the ED visit that led to the index admission, else that
   visit's first ED vital signs, else the first ICU value within 24 h of
   admission. Values outside hr 20-220, sbp 40-260 and spo2 50-100 are
   treated as missing and the next source is tried.
3. The dissertation's MIMIC-only EHR features (script 135's ceiling set)
   that are not in ehr28: every numeric lab with coverage of at least 0.15
   in these admissions, minus labs that duplicate an ehr28 column (same
   itemid, or Spearman >= 0.98); neutrophil/lymphocyte ratio; sex; number
   of prior admissions; number of diagnosis codes on prior admissions.

Validation comes first. ehr28's 5 vitals and 18 labs are re-derived under
every candidate definition, and the one that reproduces ehr28 is used for
the new columns. Vitals vary the window (patient history to index
discharge, index admission only, or 7 days before admission to index
discharge), the clipping (raw values or means) and the combination of two
itemids (pooled or mean of itemid means), with ICU values first, then ward
omr blood pressure where sbp is missing, then ED triage. Labs vary the
window, and each ehr28 lab takes the itemid (or pair) that matches best.
A value matches when |new - ehr28| <= 0.001 x max(1, |ehr28|). With
nothing patient level written, the script stops unless the chosen vitals
definition reproduces hr and sbp for at least 95% of admissions and the
chosen lab window reproduces the ehr28 labs for at least 90% of
admissions (median over the 18 labs).

Inputs:  data/inputs/ehr28.csv, BigQuery (MIMIC-IV 3.1)
Outputs: data/inputs/ehr_extra.csv (patient level, not in the repo),
         results/ehr_extra_checks.csv, results/ehr_extra_validation.csv,
         results/ehr_extra_dictionary.csv, results/ehr_extra_log.txt
Usage:   GCP_PROJECT=<billing project> python scripts/03b_bq_ehr_extra.py
Notes:   every query is dry-run first, the GB to be scanned are printed
         and the script pauses 5 s. Results are converted to float64,
         int64 or bool on arrival, because pandas 2.0.3 passes BigQuery's
         nullable types to numpy as generic objects.
"""
import os
import re
import sys
import time
import warnings
import numpy as np
import pandas as pd
from google.cloud import bigquery

warnings.filterwarnings("ignore", message=".*BigQuery Storage.*")
T0 = time.time()
pd.set_option("display.width", 200)
DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
RES = os.path.join(P2, "results")
os.makedirs(RES, exist_ok=True)
PROJECT = os.environ["GCP_PROJECT"]
HOSP = "physionet-data.mimiciv_3_1_hosp"
ICU = "physionet-data.mimiciv_3_1_icu"
ED = "physionet-data.mimiciv_ed"
COV_MIN = 0.15
TOL = 1e-3
VIT5 = ["hr", "sbp", "dbp", "rr", "temp"]
LAB18 = ["creatinine", "sodium", "potassium", "bun", "glucose",
         "calcium", "bicarb", "chloride", "hct", "plt", "wbc", "hgb",
         "aniongap", "rbc", "mchc", "mch", "mcv", "rdw"]
VAR_OF = {220045: "hr", 220050: "sbp", 220179: "sbp", 220051: "dbp",
          220180: "dbp", 220210: "rr", 220277: "spo2", 223761: "temp",
          223762: "temp"}
BOUNDS = {"hr": (20, 220), "sbp": (40, 260), "dbp": (20, 180),
          "rr": (4, 60), "spo2": (0, 100), "temp": (28, 43)}
LAB_CLIP = {"wbc": (0, 100), "creatinine": (0, 20)}
PRES_OK = {"hr": (20, 220), "sbp": (40, 260), "spo2": (50, 100)}
WINDOWS = {"subj": "patient history to index discharge",
           "hadm": "index admission only",
           "7d": "7 days before admission to index discharge"}
KEYS = ("hadm_id", "itemid", "stay_id")
CHK, VAL, DIC = [], [], []


def banner(t):
    print("", flush=True)
    print("=" * 70)
    print(t)
    print("=" * 70, flush=True)


def rec(sec, item, value, status="INFO"):
    CHK.append({"section": sec, "item": item, "value": value,
                "status": status})
    print("  %-4s %-58s %s" % (status, item, value), flush=True)


def write_small():
    pd.DataFrame(CHK).to_csv(os.path.join(
        RES, "ehr_extra_checks.csv"), index=False)
    pd.DataFrame(VAL).to_csv(os.path.join(
        RES, "ehr_extra_validation.csv"), index=False)
    pd.DataFrame(DIC).to_csv(os.path.join(
        RES, "ehr_extra_dictionary.csv"), index=False)


def stop(msg):
    write_small()
    banner("STOPPED: " + msg + " -- no patient data written")
    print("Outputs: ehr_extra_checks.csv, ehr_extra_validation.csv,"
          " ehr_extra_log.txt")
    sys.exit(1)


def clean(df):
    """Nullable BigQuery types -> plain numpy types."""
    df = df.copy()
    for c in df.columns:
        dt = str(df[c].dtype)
        if c in KEYS:
            df[c] = pd.to_numeric(df[c]).astype("int64")
        elif dt == "boolean":
            df[c] = df[c].fillna(False).astype(bool)
        elif dt in ("Int8", "Int16", "Int32", "Int64", "Float32",
                    "Float64"):
            df[c] = df[c].astype("float64")
    return df


def match(a, b):
    a = pd.to_numeric(a, errors="coerce").astype("float64")
    b = pd.to_numeric(b, errors="coerce").astype("float64")
    both = a.notna() & b.notna()
    n = int(both.sum())
    if n == 0:
        return {"n_both": 0, "share_match": np.nan, "corr": np.nan,
                "ehr28_only": int(b.notna().sum()),
                "new_only": int(a.notna().sum())}
    x = a[both].values
    y = b[both].values
    ok = np.abs(x - y) <= TOL * np.maximum(1.0, np.abs(y))
    c = float(np.corrcoef(x, y)[0, 1]) if n > 2 else np.nan
    return {"n_both": n, "share_match": round(float(ok.mean()), 4),
            "corr": round(c, 4),
            "ehr28_only": int((b.notna() & a.isna()).sum()),
            "new_only": int((a.notna() & b.isna()).sum())}


def spearman(x, y, min_n=200):
    both = x.notna() & y.notna()
    if both.sum() < min_n:
        return np.nan
    return float(np.corrcoef(x[both].rank().values,
                             y[both].rank().values)[0, 1])


# ============ 0. ACCESS AND INPUTS ============
banner("0. ACCESS AND INPUTS")
try:
    client = bigquery.Client(project=PROJECT)
    list(client.query("SELECT 1 AS ok").result())
    print("  access OK, project " + PROJECT)
except Exception as e:
    print("  BigQuery access FAILED: " + repr(e)[:300])
    print("  Fix, in a plain PowerShell window:")
    print("    gcloud auth application-default login")
    print("    (sign in with the account linked to PhysioNet)")
    print("    gcloud auth application-default set-quota-project "
          + PROJECT)
    sys.exit(1)
e28 = pd.read_csv(os.path.join(INP, "ehr28.csv"))
rec("input", "ehr28 rows / unique admissions",
    "%d / %d" % (len(e28), e28["hadm_id"].nunique()))
E = e28.drop_duplicates("hadm_id").copy()
E["hadm_id"] = E["hadm_id"].astype("int64")
E = E.set_index("hadm_id")
ids = E.index.values
for c in VIT5 + LAB18:
    if c not in E.columns:
        stop("ehr28 has no column '%s'" % c)
    E[c] = pd.to_numeric(E[c], errors="coerce").astype("float64")
H = ",".join(str(int(x)) for x in ids)
IDX = """idx AS (
  SELECT subject_id, hadm_id, admittime, dischtime
  FROM `{hosp}.admissions` WHERE hadm_id IN ({h}))""".format(
    hosp=HOSP, h=H)

# ============ 1. QUERIES ============
Q = {}
Q["labs"] = """WITH {idx},
lv AS (
  SELECT i.hadm_id AS ih, l.itemid, l.valuenum AS v,
    IFNULL(l.hadm_id = i.hadm_id, FALSE) AS wh,
    IFNULL(l.charttime >= DATETIME_SUB(i.admittime, INTERVAL 7 DAY),
           FALSE) AS w7
  FROM idx i JOIN `{hosp}.labevents` l ON l.subject_id = i.subject_id
  WHERE l.valuenum IS NOT NULL AND l.charttime <= i.dischtime)
SELECT ih AS hadm_id, itemid,
  SUM(v) AS s_subj, COUNT(*) AS n_subj,
  SUM(IF(wh, v, 0)) AS s_hadm, COUNTIF(wh) AS n_hadm,
  SUM(IF(w7, v, 0)) AS s_7d, COUNTIF(w7) AS n_7d
FROM lv GROUP BY 1, 2""".format(idx=IDX, hosp=HOSP)
Q["icu"] = """WITH {idx},
cv AS (
  SELECT i.hadm_id AS ih, c.itemid,
    CASE WHEN c.itemid = 223761 THEN (c.valuenum - 32) * 5 / 9
         ELSE c.valuenum END AS v,
    IFNULL(c.hadm_id = i.hadm_id, FALSE) AS wh,
    IFNULL(c.charttime >= DATETIME_SUB(i.admittime, INTERVAL 7 DAY),
           FALSE) AS w7
  FROM idx i JOIN `{icu}.chartevents` c ON c.subject_id = i.subject_id
  WHERE c.itemid IN (220045, 220050, 220179, 220051, 220180, 220210,
                     220277, 223761, 223762)
    AND c.valuenum IS NOT NULL AND c.charttime <= i.dischtime),
cc AS (
  SELECT ih, itemid, v, wh, w7,
    LEAST(GREATEST(v,
      CASE WHEN itemid = 220045 THEN 20
           WHEN itemid IN (220050, 220179) THEN 40
           WHEN itemid IN (220051, 220180) THEN 20
           WHEN itemid = 220210 THEN 4
           WHEN itemid = 220277 THEN 0 ELSE 28 END),
      CASE WHEN itemid = 220045 THEN 220
           WHEN itemid IN (220050, 220179) THEN 260
           WHEN itemid IN (220051, 220180) THEN 180
           WHEN itemid = 220210 THEN 60
           WHEN itemid = 220277 THEN 100 ELSE 43 END) AS vc
  FROM cv)
SELECT ih AS hadm_id, itemid,
  SUM(v) AS s_subj, SUM(vc) AS sc_subj, COUNT(*) AS n_subj,
  SUM(IF(wh, v, 0)) AS s_hadm, SUM(IF(wh, vc, 0)) AS sc_hadm,
  COUNTIF(wh) AS n_hadm,
  SUM(IF(w7, v, 0)) AS s_7d, SUM(IF(w7, vc, 0)) AS sc_7d,
  COUNTIF(w7) AS n_7d
FROM cc GROUP BY 1, 2""".format(idx=IDX, icu=ICU)
Q["omr"] = """WITH {idx}
SELECT i.hadm_id, o.result_value,
  o.chartdate >= DATE(i.admittime) AS wh,
  o.chartdate >= DATE_SUB(DATE(i.admittime), INTERVAL 7 DAY) AS w7
FROM idx i JOIN `{hosp}.omr` o ON o.subject_id = i.subject_id
WHERE o.result_name = 'Blood Pressure' AND o.result_value LIKE '%/%'
  AND o.chartdate <= DATE(i.dischtime)""".format(idx=IDX, hosp=HOSP)
Q["triage"] = """WITH {idx}
SELECT i.hadm_id, e.stay_id, e.intime,
  IFNULL(e.hadm_id = i.hadm_id, FALSE) AS wh,
  IFNULL(e.intime >= DATETIME_SUB(i.admittime, INTERVAL 7 DAY),
         FALSE) AS w7,
  t.heartrate, t.resprate, t.o2sat, t.sbp, t.dbp, t.temperature
FROM idx i JOIN `{ed}.edstays` e ON e.subject_id = i.subject_id
JOIN `{ed}.triage` t ON t.stay_id = e.stay_id
WHERE e.intime <= i.dischtime""".format(idx=IDX, ed=ED)
Q["edvitals"] = """WITH {idx}
SELECT i.hadm_id, v.charttime, v.heartrate, v.sbp, v.o2sat
FROM idx i JOIN `{ed}.edstays` e
  ON e.subject_id = i.subject_id AND e.hadm_id = i.hadm_id
JOIN `{ed}.vitalsign` v ON v.stay_id = e.stay_id""".format(
    idx=IDX, ed=ED)
Q["icu_first"] = """WITH {idx}
SELECT i.hadm_id, c.itemid, c.charttime, c.valuenum
FROM idx i JOIN `{icu}.chartevents` c ON c.hadm_id = i.hadm_id
WHERE c.itemid IN (220045, 220050, 220179, 220277)
  AND c.valuenum IS NOT NULL AND c.charttime >= i.admittime
  AND c.charttime <= DATETIME_ADD(i.admittime, INTERVAL 24 HOUR)""".format(
    idx=IDX, icu=ICU)
Q["demo"] = """WITH {idx},
pa AS (
  SELECT i.hadm_id AS ih, a.hadm_id AS ph
  FROM idx i JOIN `{hosp}.admissions` a
    ON a.subject_id = i.subject_id AND a.dischtime < i.admittime),
np AS (SELECT ih, COUNT(*) AS n_prior FROM pa GROUP BY ih),
nc AS (SELECT pa.ih, COUNT(*) AS n_codes
       FROM pa JOIN `{hosp}.diagnoses_icd` d ON d.hadm_id = pa.ph
       GROUP BY pa.ih)
SELECT i.hadm_id, p.gender, IFNULL(np.n_prior, 0) AS n_prior,
  IFNULL(nc.n_codes, 0) AS n_codes
FROM idx i
LEFT JOIN `{hosp}.patients` p ON p.subject_id = i.subject_id
LEFT JOIN np ON np.ih = i.hadm_id
LEFT JOIN nc ON nc.ih = i.hadm_id""".format(idx=IDX, hosp=HOSP)
Q["d_labitems"] = ("SELECT itemid, label, fluid, category FROM `%s.d_labitems`"
                   % HOSP)

banner("1. PROBE: DRY RUNS")
tot = 0.0
for k, sql in Q.items():
    try:
        job = client.query(sql, job_config=bigquery.QueryJobConfig(
            dry_run=True, use_query_cache=False))
    except Exception as e:
        rec("probe", "query '%s' dry run" % k, repr(e)[:160], "FAIL")
        stop("query '%s' did not compile" % k)
    gb = job.total_bytes_processed / 1e9
    tot += gb
    rec("probe", "query '%s': GB to scan" % k, "%.2f" % gb)
rec("probe", "total GB to scan (free tier: 1,000 GB/month)",
    "%.1f" % tot)
print("  starting in 5 s (Ctrl+C to stop)...", flush=True)
time.sleep(5)

banner("2. PULL")
R = {}
for k, sql in Q.items():
    t = time.time()
    R[k] = clean(client.query(sql).to_dataframe())
    rec("pull", "query '%s': rows" % k, "%d (%.0f s)"
        % (len(R[k]), time.time() - t))

# ============ 3. VITALS VALIDATION ============
banner("3. VITALS: REPRODUCE ehr28")
IC = R["icu"].copy()
IC["var"] = IC["itemid"].map(VAR_OF)


def icu_means(w, clip, comb):
    d = IC[IC["n_" + w] > 0].copy()
    d["s"] = d["sc_" + w] if clip == "raw" else d["s_" + w]
    d["n"] = d["n_" + w]
    if comb == "pooled":
        g = d.groupby(["hadm_id", "var"])[["s", "n"]].sum()
        m = g["s"] / g["n"]
    else:
        d["m"] = d["s"] / d["n"]
        m = d.groupby(["hadm_id", "var"])["m"].mean()
    m = m.unstack("var").reindex(index=ids,
                                 columns=list(BOUNDS.keys()))
    m = m.astype("float64")
    if clip == "mean":
        for v, (lo, hi) in BOUNDS.items():
            m[v] = m[v].clip(lo, hi)
    return m


OM = R["omr"].copy()
pp = OM["result_value"].astype(str).str.split("/", n=1, expand=True)
OM["sbp"] = pd.to_numeric(pp[0], errors="coerce")
OM["dbp"] = pd.to_numeric(pp[1], errors="coerce")
OM = OM.dropna(subset=["sbp", "dbp"])
OM["sbp"] = OM["sbp"].clip(*BOUNDS["sbp"])
OM["dbp"] = OM["dbp"].clip(*BOUNDS["dbp"])
TR = R["triage"].copy()
TRV = pd.DataFrame({"hadm_id": TR["hadm_id"]})
for src, v in [("heartrate", "hr"), ("resprate", "rr"),
               ("o2sat", "spo2"), ("sbp", "sbp"), ("dbp", "dbp"),
               ("temperature", "temp")]:
    x = pd.to_numeric(TR[src], errors="coerce").astype("float64")
    if v == "temp":
        x = x.where(x < 50, (x - 32) * 5 / 9)
    TRV[v] = x.clip(*BOUNDS[v])


def sub(df, w):
    if w == "subj":
        return df
    return df[df["wh" if w == "hadm" else "w7"].astype(bool)]


def build_vitals(w, clip, comb):
    V = icu_means(w, clip, comb)
    need = V["sbp"].isna()
    o = sub(OM, w).groupby("hadm_id")[["sbp", "dbp"]].mean().reindex(ids)
    for c in ("sbp", "dbp"):
        V.loc[need, c] = V.loc[need, c].fillna(o.loc[need, c])
    tv = TRV.loc[sub(TR, w).index]
    t = tv.groupby("hadm_id")[list(BOUNDS.keys())].mean().reindex(ids)
    for c in BOUNDS:
        V[c] = V[c].fillna(t[c])
    return V


best_v, best_s = None, -1.0
for w in WINDOWS:
    for clip in ("raw", "mean"):
        for comb in ("pooled", "mean_of_itemids"):
            V = build_vitals(w, clip, comb)
            ss = []
            for c in VIT5:
                m = match(V[c], E[c])
                ss.append(m["share_match"])
                VAL.append(dict(domain="vitals", window=w, clip=clip,
                                combine=comb, feature=c, itemids="",
                                label="", **m))
            s = float(np.nanmean(ss))
            print("  %-5s clip %-4s %-15s mean share %.4f  (hr %.4f,"
                  " sbp %.4f)" % (w, clip, comb, s, ss[0], ss[1]))
            if s > best_s:
                best_v, best_s = (w, clip, comb), s
V = build_vitals(*best_v)
mh = match(V["hr"], E["hr"])["share_match"]
ms = match(V["sbp"], E["sbp"])["share_match"]
rec("vitals", "chosen definition (window, clip, combine)",
    " / ".join(best_v))
rec("vitals", "share of ehr28 hr reproduced", mh,
    "PASS" if mh >= 0.95 else "FAIL")
rec("vitals", "share of ehr28 sbp reproduced", ms,
    "PASS" if ms >= 0.95 else "FAIL")
if mh < 0.95 or ms < 0.95:
    stop("ehr28 vitals not reproduced; check 03_mimic_vitals.py")

# ============ 4. LABS VALIDATION ============
banner("4. LABS: REPRODUCE ehr28")
LB = R["labs"].copy()
DL = R["d_labitems"].copy().set_index("itemid")
cov0 = LB[LB["n_subj"] > 0].groupby("itemid")["hadm_id"].nunique()
items = cov0[cov0 >= 0.005 * len(ids)].index.values.astype("int64")
LB = LB[LB["itemid"].isin(items)]
SN = {}
for w in WINDOWS:
    s = LB.pivot(index="hadm_id", columns="itemid",
                 values="s_" + w).reindex(index=ids, columns=items)
    n = LB.pivot(index="hadm_id", columns="itemid",
                 values="n_" + w).reindex(index=ids, columns=items)
    SN[w] = (s.fillna(0).values.astype("float64"),
             n.fillna(0).values.astype("float64"))


def share_vec(M, y):
    Y = y[:, None]
    with np.errstate(invalid="ignore"):
        ok = np.abs(M - Y) <= TOL * np.maximum(1.0, np.abs(Y))
    both = ~np.isnan(M) & ~np.isnan(Y)
    nb = both.sum(0)
    return np.where(nb >= 50, (ok & both).sum(0) / np.maximum(nb, 1),
                    0.0), nb


def lab_label(j):
    return (str(DL["label"].get(j, "?")) + " | "
            + str(DL["fluid"].get(j, "?")))


LABMAP = {}
wscore = {}
for w in WINDOWS:
    S, N = SN[w]
    with np.errstate(invalid="ignore", divide="ignore"):
        M = np.where(N > 0, S / np.where(N > 0, N, 1), np.nan)
    shares = []
    LABMAP[w] = {}
    for c in LAB18:
        y = E[c].reindex(ids).values.astype(float)
        Mc = M.copy()
        if c in LAB_CLIP:
            Mc = np.clip(Mc, *LAB_CLIP[c])
        sh, nb = share_vec(Mc, y)
        j = int(np.argmax(sh))
        best_set, best = [j], sh[j]
        for k in np.argsort(-sh)[:8]:
            if k == j:
                continue
            num = S[:, j] + S[:, k]
            den = N[:, j] + N[:, k]
            with np.errstate(invalid="ignore", divide="ignore"):
                pm = np.where(den > 0, num / np.where(den > 0, den, 1),
                              np.nan)
            if c in LAB_CLIP:
                pm = np.clip(pm, *LAB_CLIP[c])
            s2, _ = share_vec(pm[:, None], y)
            if s2[0] > best + 0.01:
                best_set, best = [j, int(k)], s2[0]
        LABMAP[w][c] = [int(items[q]) for q in best_set]
        shares.append(best)
        VAL.append(dict(domain="labs", window=w, clip="", combine="",
                        feature=c,
                        itemids="+".join(str(items[q]) for q in best_set),
                        label="; ".join(lab_label(items[q])
                                        for q in best_set),
                        n_both=int(nb[j]), share_match=round(best, 4),
                        corr=np.nan, ehr28_only=np.nan, new_only=np.nan))
    wscore[w] = float(np.median(shares))
    print("  window %-5s median share over 18 labs %.4f (min %.4f)"
          % (w, wscore[w], min(shares)))
wl = max(wscore, key=wscore.get)
rec("labs", "chosen window", wl + " (" + WINDOWS[wl] + ")")
rec("labs", "median share of ehr28 labs reproduced",
    round(wscore[wl], 4), "PASS" if wscore[wl] >= 0.90 else "FAIL")
if wscore[wl] < 0.90:
    stop("ehr28 labs not reproduced; check the SIDE PROJECT lab script")

# ============ 5. NEW FEATURES ============
banner("5. NEW FEATURES")
S, N = SN[wl]
with np.errstate(invalid="ignore", divide="ignore"):
    M = pd.DataFrame(np.where(N > 0, S / np.where(N > 0, N, 1), np.nan),
                     index=ids, columns=items)
used = {i for c in LAB18 for i in LABMAP[wl][c]}
new = pd.DataFrame(index=ids)
new.index.name = "hadm_id"
for j in items:
    lab_ = str(DL["label"].get(j, ""))
    cov = float(M[j].notna().mean())
    row = {"column": "", "itemids": str(j), "label": lab_,
           "fluid": str(DL["fluid"].get(j, "")),
           "category": str(DL["category"].get(j, "")),
           "coverage": round(cov, 4)}
    if cov < COV_MIN:
        continue
    if j in used:
        row["status"] = "excluded: the ehr28 source itemid"
    elif len(lab_.strip()) <= 1:
        row["status"] = "excluded: flag item"
    elif M[j].nunique() <= 1:
        row["status"] = "excluded: constant"
    else:
        rs = {c: spearman(M[j], E[c].reindex(ids)) for c in LAB18}
        rs = {c: abs(r) for c, r in rs.items() if not np.isnan(r)}
        top = max(rs, key=rs.get) if rs else None
        if top is not None and rs[top] >= 0.98:
            row["status"] = "excluded: duplicates ehr28 %s (rho %.3f)" % (
                top, rs[top])
        else:
            slug = re.sub("[^a-z0-9]+", "_", lab_.lower()).strip("_")[:28]
            col = "lab_%d_%s" % (j, slug)
            new[col] = M[j].values
            row["column"] = col
            row["status"] = "kept"
            q = M[j].quantile([0.25, 0.5, 0.75])
            row.update(p25=q.iloc[0], median=q.iloc[1], p75=q.iloc[2])
    DIC.append(row)
nk = sum(1 for r in DIC if r.get("status") == "kept")
rec("new", "labs with coverage >= %.2f: kept / excluded" % COV_MIN,
    "%d / %d" % (nk, len(DIC) - nk))

new["spo2"] = V["spo2"].reindex(ids).values
cov = float(new["spo2"].notna().mean())
rec("new", "SpO2 (averaged like ehr28): coverage / median",
    "%.3f / %.1f" % (cov, new["spo2"].median()),
    "PASS" if cov >= 0.5 else "WARN")
DIC.append({"column": "spo2", "itemids": "220277 + ED triage o2sat",
            "label": "SpO2, averaged as ehr28 vitals", "coverage":
            round(cov, 4), "median": new["spo2"].median(),
            "status": "kept"})
if 51256 in M.columns and 51244 in M.columns:
    ly = M[51244].where(M[51244] > 0)
    new["nlr"] = (M[51256] / ly).values
    DIC.append({"column": "nlr", "itemids": "51256/51244",
                "label": "neutrophil % / lymphocyte %",
                "coverage": round(float(new["nlr"].notna().mean()), 4),
                "median": new["nlr"].median(), "status": "kept"})

# presentation vitals
t = TR[TR["wh"]].sort_values("intime")
tri = pd.DataFrame({"hadm_id": t["hadm_id"],
                    "hr": pd.to_numeric(t["heartrate"], errors="coerce"),
                    "sbp": pd.to_numeric(t["sbp"], errors="coerce"),
                    "spo2": pd.to_numeric(t["o2sat"], errors="coerce")})
ev = R["edvitals"].sort_values("charttime")
edv = pd.DataFrame({"hadm_id": ev["hadm_id"],
                    "hr": pd.to_numeric(ev["heartrate"], errors="coerce"),
                    "sbp": pd.to_numeric(ev["sbp"], errors="coerce"),
                    "spo2": pd.to_numeric(ev["o2sat"], errors="coerce")})
fi = R["icu_first"].sort_values("charttime").copy()
fi["var"] = fi["itemid"].map({220045: "hr", 220050: "sbp", 220179: "sbp",
                              220277: "spo2"})
icf = fi.pivot_table(index="hadm_id", columns="var", values="valuenum",
                     aggfunc="first")


def first_valid(df, v):
    x = df[["hadm_id", v]].copy()
    lo, hi = PRES_OK[v]
    x = x[(x[v] >= lo) & (x[v] <= hi)]
    return x.groupby("hadm_id")[v].first().reindex(ids)


for v in PRES_OK:
    a = first_valid(tri, v)
    b = first_valid(edv, v)
    lo, hi = PRES_OK[v]
    if v in icf.columns:
        c = icf[v].where((icf[v] >= lo) & (icf[v] <= hi)).reindex(ids)
    else:
        c = pd.Series(np.nan, index=ids)
    val = a.fillna(b).fillna(c)
    src = np.where(a.notna(), "triage", np.where(
        b.notna(), "ed_vitals", np.where(c.notna(), "icu_24h", "")))
    new[v + "_first"] = val.values
    new[v + "_first_src"] = src
    rec("presentation", "%s_first: coverage (triage/ED/ICU shares)" % v,
        "%.3f (%.3f/%.3f/%.3f)" % (val.notna().mean(),
                                   (src == "triage").mean(),
                                   (src == "ed_vitals").mean(),
                                   (src == "icu_24h").mean()))
    DIC.append({"column": v + "_first", "itemids": "triage > ED vitals"
                " > ICU 24h", "label": "%s at presentation" % v,
                "coverage": round(float(val.notna().mean()), 4),
                "median": val.median(), "status": "kept"})

d = R["demo"].set_index("hadm_id").reindex(ids)
new["sex_f"] = (d["gender"] == "F").astype(float).where(
    d["gender"].notna()).values
new["n_prior"] = d["n_prior"].values
new["n_codes"] = d["n_codes"].values
for c in ("sex_f", "n_prior", "n_codes"):
    DIC.append({"column": c, "itemids": "", "label": {
        "sex_f": "female sex", "n_prior": "admissions before the index",
        "n_codes": "diagnosis codes on prior admissions"}[c],
        "coverage": round(float(new[c].notna().mean()), 4),
        "median": float(np.nanmedian(new[c])), "status": "kept"})
rec("new", "female share / median prior admissions",
    "%.3f / %.0f" % (np.nanmean(new["sex_f"]),
                     np.nanmedian(new["n_prior"])))

# ============ 6. WRITE ============
banner("6. WRITE")
out = new.reset_index()
out.to_csv(os.path.join(INP, "ehr_extra.csv"), index=False)
rec("write", "ehr_extra.csv rows / columns",
    "%d / %d" % (out.shape[0], out.shape[1] - 1))
write_small()
nf = sum(1 for r in CHK if r["status"] == "FAIL")
print("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min"
      % (sum(1 for r in CHK if r["status"] == "PASS"),
         sum(1 for r in CHK if r["status"] == "WARN"), nf,
         (time.time() - T0) / 60))
print("\nOUTPUT FILES")
print("  " + os.path.join(INP, "ehr_extra.csv") + "   (local only)")
for f in ("ehr_extra_checks.csv", "ehr_extra_validation.csv",
          "ehr_extra_dictionary.csv"):
    print("  " + os.path.join(RES, f))
print("  " + os.path.join(RES, "ehr_extra_log.txt")
      + "   (terminal log)")
print("\nOutputs: ehr_extra_checks.csv, ehr_extra_validation.csv,"
      " ehr_extra_dictionary.csv")