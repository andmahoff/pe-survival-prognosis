"""03d_bq_diagnoses.py
Diagnosis codes of every Frame B index admission (A30).

Pulled from MIMIC-IV on BigQuery (diagnoses_icd) for the comorbidity
scores of 06j_level1_ehrplus.R. The codes are assigned at discharge, so
they are used in Frame B only (prediction at discharge); in Frame A they
would leak the outcome.

Checks: admissions with at least one code (FAIL below 99%); ICD-9 and
ICD-10 split.

Inputs:  data/processed/frame_B.csv; BigQuery (MIMIC-IV 3.1)
Outputs: data/inputs/dx_index.csv (patient level, not in the repo),
         results/dx_index_summary.csv (counts only),
         results/dx_index_log.txt
Usage:   GCP_PROJECT=<billing project> python scripts/03d_bq_diagnoses.py
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


def main():
    fb = pd.read_csv(os.path.join(PRO, "frame_B.csv"), usecols=["hadm_id"])
    ids = [int(x) for x in fb["hadm_id"].unique()]
    print("Frame B index admissions: %d" % len(ids))
    sql = ("SELECT hadm_id, seq_num, icd_code, icd_version FROM "
           "`%s.diagnoses_icd` WHERE hadm_id IN UNNEST(@ids)" % HOSP)
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
    d = pd.DataFrame({"hadm_id": d["hadm_id"].astype("int64"),
                      "seq_num": pd.to_numeric(d["seq_num"],
                                               errors="coerce")
                      .astype(float),
                      "icd_code": d["icd_code"].astype(str).str.strip(),
                      "icd_version": d["icd_version"].astype("int64")})
    d = d.sort_values(["hadm_id", "seq_num"])
    found = fb["hadm_id"].isin(d["hadm_id"]).mean()
    status = "PASS" if found >= 0.99 else "FAIL"
    print("  %-5s admissions with codes: %.1f%%" % (status, 100 * found))
    d.to_csv(os.path.join(INP, "dx_index.csv"), index=False)
    per = d.groupby("hadm_id").size()
    ver = d.drop_duplicates("hadm_id")["icd_version"].value_counts()
    s = pd.DataFrame([
        {"item": "admissions with codes", "value": "%.1f%%" % (100 * found),
         "status": status},
        {"item": "codes per admission (median)", "value": "%g" %
         per.median(), "status": "INFO"},
        {"item": "admissions coded in ICD-9 / ICD-10",
         "value": "%d / %d" % (ver.get(9, 0), ver.get(10, 0)),
         "status": "INFO"}])
    s.to_csv(os.path.join(RES, "dx_index_summary.csv"), index=False)
    print(s.to_string(index=False))
    print("\nOUTPUT FILES")
    print("  " + os.path.join(INP, "dx_index.csv"), " (local only)")
    print("  " + os.path.join(RES, "dx_index_summary.csv"))
    print("  " + os.path.join(RES, "dx_index_log.txt"),
          " (terminal log)")
    print("\nOutputs: dx_index_summary.csv")
    if status == "FAIL":
        sys.exit(1)


if __name__ == "__main__":
    main()