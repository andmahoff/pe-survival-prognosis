"""03e_bq_presentation.py
Information available at presentation for every Frame A admission (A32).

Pulled from MIMIC-IV on BigQuery for the at-diagnosis model. The window
runs from 24 h before to 24 h after the admission time, so the emergency
department stay that led to the admission falls inside it.
  labs     the first value in the window of every lab analyte of the
           79-feature EHR: the 18 registered labs (item IDs as matched by
           03b_bq_ehr_extra.py, read from ehr_extra_validation.csv) and
           the extra labs (item IDs from the lab_<itemid>_ column names of
           ehr_extra.csv)
  vitals   the first plausible value in the window of heart rate, systolic
           and diastolic blood pressure, respiratory rate, oxygen
           saturation and temperature, from ED triage (timed at ED
           arrival), ED vital signs and ICU charting. Values outside the
           plausible ranges of 03b are discarded, not clipped, and
           Fahrenheit temperatures are converted.
  history  the diagnosis codes of every earlier admission (discharged
           before the index admission began), for comorbidity history
  index    the index admission's own codes (assigned at discharge), used
           only in the pre-specified sensitivity analysis of A32

Checks: item IDs found for all 18 registered labs (FAIL); admissions with
any lab and with any vital sign in the window (FAIL below 90%).

Inputs:  data/processed/frame_A.csv, data/inputs/ehr_extra.csv,
         results/ehr_extra_validation.csv; BigQuery (MIMIC-IV 3.1)
Outputs: data/inputs/pres_labs_first.csv, pres_vitals_first.csv,
         dx_prior.csv, dx_index_A.csv (patient level, not in the repo),
         results/pres_pull_summary.csv (coverage only),
         results/pres_pull_log.txt
Usage:   GCP_PROJECT=<billing project> python scripts/03e_bq_presentation.py
Notes:   every query is dry-run first; the GB to be scanned are printed
         and the script pauses 5 s. Item IDs from 03b are accepted
         whether saved as numbers ("50912.0") or as text. Run with
         google-cloud-bigquery and db-dtypes.
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
DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
PRO = os.path.join(P2, "data", "processed")
INP = os.path.join(P2, "data", "inputs")
RES = os.path.join(P2, "results")
PROJECT = os.environ["GCP_PROJECT"]
HOSP = "physionet-data.mimiciv_3_1_hosp"
ICU = "physionet-data.mimiciv_3_1_icu"
ED = "physionet-data.mimiciv_ed"
HOURS = 24
LAB18 = ["creatinine", "sodium", "potassium", "bun", "glucose",
         "calcium", "bicarb", "chloride", "hct", "plt", "wbc", "hgb",
         "aniongap", "rbc", "mchc", "mch", "mcv", "rdw"]
BOUNDS = {"hr": (20, 220), "sbp": (40, 260), "dbp": (20, 180),
          "rr": (4, 60), "spo2": (50, 100), "temp": (28, 43)}
ICU_ITEMS = {220045: "hr", 220050: "sbp", 220179: "sbp", 220051: "dbp",
             220180: "dbp", 220210: "rr", 220277: "spo2", 223761: "temp",
             223762: "temp"}
SUM = []


def rec(item, value, status="INFO"):
    SUM.append({"item": item, "value": value, "status": status})
    print("  %-5s %-58s %s" % (status, item, value), flush=True)


def stop(msg):
    pd.DataFrame(SUM).to_csv(os.path.join(RES, "pres_pull_summary.csv"),
                             index=False)
    print("\nSTOPPED: %s. see: pres_pull_summary.csv" % msg)
    sys.exit(1)


def lab_itemids():
    """18 registered labs from stage 3b's validation, extras from names."""
    v = pd.read_csv(os.path.join(RES, "ehr_extra_validation.csv"))
    v = v[v["domain"] == "labs"]
    m18 = {}
    for c in LAB18:
        x = v[v["feature"] == c].sort_values("share_match",
                                             ascending=False)
        if len(x):
            # item IDs may have been read as numbers ("50912.0") when a
            # lab matched a single item, or as text ("50912+50809")
            m18[c] = [int(float(i)) for i in
                      str(x["itemids"].iloc[0]).split("+")]
    hdr = pd.read_csv(os.path.join(INP, "ehr_extra.csv"), nrows=1).columns
    extra = {}
    for c in hdr:
        mm = re.match(r"^lab_(\d+)_", c)
        if mm and not re.search("vancomycin|_50825_", c):
            extra[c] = [int(mm.group(1))]
    return m18, extra


def first_per(df, keys, tcol):
    """earliest row per key group (ties broken by row order)."""
    df = df.sort_values(keys + [tcol], kind="mergesort")
    return df.drop_duplicates(keys, keep="first")


def build_labs(long, m18, extra):
    out = {}
    for name, ids in list(m18.items()) + list(extra.items()):
        x = long[long["itemid"].isin(ids)]
        x = first_per(x, ["hadm_id"], "h_rel")
        out[name + "_first"] = x.set_index("hadm_id")["valuenum"]
    return pd.DataFrame(out)


def build_vitals(tri, edv, icu):
    rows = []
    for v in ["hr", "sbp", "dbp", "rr", "spo2", "temp"]:
        for src, d, col in (("triage", tri, v), ("ed", edv, v),
                            ("icu", icu, v)):
            if col not in d.columns:
                continue
            x = d[["hadm_id", "h_rel", col]].rename(columns={col: "val"})
            x = x[x["val"].notna()].copy()
            x["vital"] = v
            x["src"] = src
            rows.append(x)
    a = pd.concat(rows, ignore_index=True)
    a["val"] = a["val"].astype(float)
    hot = (a["vital"] == "temp") & (a["val"] > 50)
    a.loc[hot, "val"] = (a.loc[hot, "val"] - 32) * 5 / 9
    lo = a["vital"].map(lambda k: BOUNDS[k][0])
    hi = a["vital"].map(lambda k: BOUNDS[k][1])
    a = a[(a["val"] >= lo) & (a["val"] <= hi)]
    a = first_per(a, ["hadm_id", "vital"], "h_rel")
    w = a.pivot(index="hadm_id", columns="vital", values="val")
    w.columns = [c + "_first" for c in w.columns]
    s = a.pivot(index="hadm_id", columns="vital", values="src")
    s.columns = [c + "_first_src" for c in s.columns]
    return w.join(s)


def main():
    t0 = time.time()
    fa = pd.read_csv(os.path.join(PRO, "frame_A.csv"),
                     usecols=["subject_id", "hadm_id"])
    ids = sorted(int(x) for x in fa["hadm_id"].unique())
    rec("Frame A index admissions", len(ids))
    m18, extra = lab_itemids()
    rec("registered labs with item IDs (of 18)", len(m18),
        "PASS" if len(m18) == 18 else "FAIL")
    rec("extra labs with item IDs", len(extra))
    if len(m18) < 18:
        stop("not every registered lab has item IDs")
    items = sorted({i for v in list(m18.values()) + list(extra.values())
                    for i in v})
    H = ",".join(str(x) for x in ids)
    I = ",".join(str(x) for x in items)
    idx = ("idx AS (SELECT subject_id, hadm_id, admittime FROM "
           "`%s.admissions` WHERE hadm_id IN (%s))" % (HOSP, H))
    win = ("BETWEEN DATETIME_SUB(i.admittime, INTERVAL %d HOUR) "
           "AND DATETIME_ADD(i.admittime, INTERVAL %d HOUR)"
           % (HOURS, HOURS))
    Q = {}
    Q["labs"] = (
        "WITH %s SELECT i.hadm_id, l.itemid, l.valuenum, "
        "DATETIME_DIFF(l.charttime, i.admittime, MINUTE) / 60 AS h_rel "
        "FROM idx i JOIN `%s.labevents` l ON l.subject_id = i.subject_id "
        "WHERE l.itemid IN (%s) AND l.valuenum IS NOT NULL "
        "AND l.charttime %s" % (idx, HOSP, I, win))
    Q["triage"] = (
        "WITH %s SELECT i.hadm_id, "
        "DATETIME_DIFF(e.intime, i.admittime, MINUTE) / 60 AS h_rel, "
        "t.heartrate AS hr, t.sbp, t.dbp, t.resprate AS rr, "
        "t.o2sat AS spo2, t.temperature AS temp "
        "FROM idx i JOIN `%s.edstays` e ON e.subject_id = i.subject_id "
        "JOIN `%s.triage` t ON t.stay_id = e.stay_id "
        "WHERE e.intime %s" % (idx, ED, ED, win))
    Q["edvitals"] = (
        "WITH %s SELECT i.hadm_id, "
        "DATETIME_DIFF(v.charttime, i.admittime, MINUTE) / 60 AS h_rel, "
        "v.heartrate AS hr, v.sbp, v.dbp, v.resprate AS rr, "
        "v.o2sat AS spo2, v.temperature AS temp "
        "FROM idx i JOIN `%s.edstays` e ON e.subject_id = i.subject_id "
        "JOIN `%s.vitalsign` v ON v.stay_id = e.stay_id "
        "WHERE v.charttime %s" % (idx, ED, ED, win))
    Q["icu"] = (
        "WITH %s SELECT i.hadm_id, c.itemid, c.valuenum, "
        "DATETIME_DIFF(c.charttime, i.admittime, MINUTE) / 60 AS h_rel "
        "FROM idx i JOIN `%s.chartevents` c ON c.subject_id = i.subject_id "
        "WHERE c.itemid IN (%s) AND c.valuenum IS NOT NULL "
        "AND c.charttime %s"
        % (idx, ICU, ",".join(str(k) for k in ICU_ITEMS), win))
    Q["dx_prior"] = (
        "WITH %s SELECT i.hadm_id, d.icd_code, d.icd_version "
        "FROM idx i JOIN `%s.admissions` a ON a.subject_id = i.subject_id "
        "AND a.dischtime < i.admittime "
        "JOIN `%s.diagnoses_icd` d ON d.hadm_id = a.hadm_id"
        % (idx, HOSP, HOSP))
    Q["dx_index"] = (
        "SELECT hadm_id, seq_num, icd_code, icd_version FROM "
        "`%s.diagnoses_icd` WHERE hadm_id IN (%s)" % (HOSP, H))
    client = bigquery.Client(project=PROJECT)
    print("\n=== PROBE: DRY RUNS ===")
    tot = 0.0
    for k, sql in Q.items():
        try:
            job = client.query(sql, job_config=bigquery.QueryJobConfig(
                dry_run=True, use_query_cache=False))
        except Exception as e:
            rec("dry run: " + k, str(e)[:150], "FAIL")
            stop("query %s did not compile (sign in again with gcloud "
                 "auth application-default login if it is an access "
                 "error)" % k)
        tot += job.total_bytes_processed / 1e9
        print("  %-9s %.3f GB" % (k, job.total_bytes_processed / 1e9))
    print("  total %.2f GB; starting in 5 s (Ctrl+C to stop)..." % tot,
          flush=True)
    time.sleep(5)
    R = {}
    for k, sql in Q.items():
        R[k] = client.query(sql).to_dataframe()
        for c in R[k].columns:
            if c in ("hadm_id", "itemid", "icd_version", "seq_num"):
                R[k][c] = pd.to_numeric(R[k][c]).astype("int64")
            elif c == "icd_code":
                R[k][c] = R[k][c].astype(str).str.strip()
            else:
                R[k][c] = pd.to_numeric(R[k][c], errors="coerce") \
                    .astype("float64")
        print("  %-9s %d rows" % (k, len(R[k])), flush=True)
    labs = build_labs(R["labs"], m18, extra)
    icu = R["icu"].copy()
    icu["vital"] = icu["itemid"].map(ICU_ITEMS)
    icu = icu.pivot_table(index=["hadm_id", "h_rel"], columns="vital",
                          values="valuenum", aggfunc="first").reset_index()
    vit = build_vitals(R["triage"], R["edvitals"], icu)
    labs = labs.reindex(ids)
    vit = vit.reindex(ids)
    labs.index.name = vit.index.name = "hadm_id"
    labs.reset_index().to_csv(os.path.join(INP, "pres_labs_first.csv"),
                              index=False)
    vit.reset_index().to_csv(os.path.join(INP, "pres_vitals_first.csv"),
                             index=False)
    R["dx_prior"].to_csv(os.path.join(INP, "dx_prior.csv"), index=False)
    R["dx_index"].to_csv(os.path.join(INP, "dx_index_A.csv"), index=False)
    anylab = labs.notna().any(axis=1).mean()
    anyvit = vit.filter(regex="_first$").notna().any(axis=1).mean()
    rec("admissions with any first lab in the window",
        "%.1f%%" % (100 * anylab), "PASS" if anylab >= 0.9 else "FAIL")
    rec("admissions with any first vital in the window",
        "%.1f%%" % (100 * anyvit), "PASS" if anyvit >= 0.9 else "FAIL")
    for c in labs.columns[:18]:
        rec("coverage " + c, "%.1f%%" % (100 * labs[c].notna().mean()))
    for c in [c for c in vit.columns if c.endswith("_first")]:
        rec("coverage " + c, "%.1f%%" % (100 * vit[c].notna().mean()))
    rec("extra labs: median coverage", "%.1f%%" % (
        100 * labs.iloc[:, 18:].notna().mean().median()))
    hp = R["dx_prior"]["hadm_id"].nunique() / len(ids)
    rec("admissions with any earlier diagnosis codes", "%.1f%%" %
        (100 * hp))
    rec("admissions with index codes", "%.1f%%" % (
        100 * R["dx_index"]["hadm_id"].nunique() / len(ids)))
    pd.DataFrame(SUM).to_csv(os.path.join(RES, "pres_pull_summary.csv"),
                             index=False)
    print("\n  elapsed %.1f min" % ((time.time() - t0) / 60))
    print("\nOUTPUT FILES")
    for f in ("pres_labs_first.csv", "pres_vitals_first.csv",
              "dx_prior.csv", "dx_index_A.csv"):
        print("  " + os.path.join(INP, f), " (local only)")
    print("  " + os.path.join(RES, "pres_pull_summary.csv"))
    print("  " + os.path.join(RES, "pres_pull_log.txt"),
          " (terminal log)")
    print("\nOutputs: pres_pull_summary.csv")
    if any(r["status"] == "FAIL" for r in SUM):
        sys.exit(1)


if __name__ == "__main__":
    main()