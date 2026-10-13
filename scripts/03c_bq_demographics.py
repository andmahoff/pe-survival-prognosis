"""03c_bq_demographics.py
Sex, recorded race and anchor age for every Frame A and B admission (A29).

Pulled from MIMIC-IV on BigQuery for the subgroup evaluation
(22a_subgroups.R). Race comes from admissions.race and is grouped as White
(WHITE..., PORTUGUESE), Black (BLACK...), Hispanic or Latino (HISPANIC...,
SOUTH AMERICAN), Asian (ASIAN...), Unknown (UNKNOWN, UNABLE TO OBTAIN,
PATIENT DECLINED TO ANSWER, missing) and Other (everything else). Race is
not a model input; it is used only to evaluate the models.

Checks: every frame admission is found (FAIL below 99%); group counts.

Inputs:  data/processed/frame_A.csv, frame_B.csv; BigQuery (MIMIC-IV 3.1)
Outputs: data/inputs/demographics.csv (patient level, not in the repo),
         results/demographics_summary.csv (counts only),
         results/demographics_log.txt
Usage:   GCP_PROJECT=<billing project> python scripts/03c_bq_demographics.py
Notes:   one query, dry-run first; the GB to be scanned are printed and
         the script pauses 5 s. Run with google-cloud-bigquery and
         db-dtypes.
"""
import os
import sys
import time
import warnings
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


def group(r):
    r = str(r).upper().strip()
    if r.startswith("WHITE") or r == "PORTUGUESE":
        return "White"
    if r.startswith("BLACK"):
        return "Black"
    if r.startswith("HISPANIC") or r.startswith("SOUTH AMERICAN"):
        return "Hispanic or Latino"
    if r.startswith("ASIAN"):
        return "Asian"
    if r in ("UNKNOWN", "UNABLE TO OBTAIN", "PATIENT DECLINED TO ANSWER",
             "NAN", "NONE", ""):
        return "Unknown"
    return "Other"


def main():
    fr = pd.concat([pd.read_csv(os.path.join(PRO, "frame_%s.csv" % f),
                                usecols=["subject_id", "hadm_id"])
                    for f in ("A", "B")]).drop_duplicates("hadm_id")
    ids = [int(x) for x in fr["hadm_id"]]
    print("admissions in Frames A and B: %d" % len(ids))
    sql = ("SELECT a.subject_id, a.hadm_id, a.race, p.gender, "
           "p.anchor_age FROM `%s.admissions` a JOIN `%s.patients` p "
           "USING (subject_id) WHERE a.hadm_id IN UNNEST(@ids)"
           % (HOSP, HOSP))
    par = [bigquery.ArrayQueryParameter("ids", "INT64", ids)]
    client = bigquery.Client(project=PROJECT)
    try:
        dry = client.query(sql, job_config=bigquery.QueryJobConfig(
            query_parameters=par, dry_run=True, use_query_cache=False))
    except Exception as e:
        print("FAIL: BigQuery access (%s). Sign in again with: gcloud auth "
              "application-default login (the account linked to PhysioNet)"
              % str(e)[:200])
        sys.exit(1)
    print("dry run: %.3f GB" % (dry.total_bytes_processed / 1e9))
    print("starting in 5 s (Ctrl+C to stop)...", flush=True)
    time.sleep(5)
    d = client.query(sql, job_config=bigquery.QueryJobConfig(
        query_parameters=par)).to_dataframe()
    d = pd.DataFrame({"subject_id": d["subject_id"].astype("int64"),
                      "hadm_id": d["hadm_id"].astype("int64"),
                      "race": d["race"].astype(str),
                      "gender": d["gender"].astype(str),
                      "anchor_age": pd.to_numeric(d["anchor_age"],
                                                  errors="coerce")
                      .astype(float)})
    d["race_group"] = d["race"].map(group)
    d = d.drop_duplicates("hadm_id")
    found = fr["hadm_id"].isin(d["hadm_id"]).mean()
    status = "PASS" if found >= 0.99 else "FAIL"
    print("  %-5s admissions found: %.1f%%" % (status, 100 * found))
    d.to_csv(os.path.join(INP, "demographics.csv"), index=False)
    rows = []
    for f in ("A", "B"):
        fx = pd.read_csv(os.path.join(PRO, "frame_%s.csv" % f),
                         usecols=["hadm_id", "in_three_mod"])
        m = fx.merge(d, on="hadm_id", how="left")
        for co, sub in (("EHR+ECG", m), ("3-mod", m[m["in_three_mod"] ==
                                                     1])):
            for var in ("race_group", "gender"):
                for k, v in sub[var].fillna("missing").value_counts() \
                        .items():
                    rows.append({"frame": f, "cohort": co, "variable": var,
                                 "group": k, "n": int(v)})
    s = pd.DataFrame(rows)
    s["check"] = status
    s.to_csv(os.path.join(RES, "demographics_summary.csv"), index=False)
    print(s[(s["frame"] == "B") & (s["cohort"] == "EHR+ECG")]
          .to_string(index=False))
    print("\nOUTPUT FILES")
    print("  " + os.path.join(INP, "demographics.csv"), " (local only)")
    print("  " + os.path.join(RES, "demographics_summary.csv"))
    print("  " + os.path.join(RES, "demographics_log.txt"),
          " (terminal log)")
    print("\nOutputs: demographics_summary.csv")
    if status == "FAIL":
        sys.exit(1)


if __name__ == "__main__":
    main()