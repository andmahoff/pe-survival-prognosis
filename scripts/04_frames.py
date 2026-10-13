"""04_frames.py
Analysis frames and cross-validation folds for MIMIC (D1).

Builds the two MIMIC analysis frames from the frozen inputs and the
03_bq_followup.py pull, codes every endpoint as (time, event), fixes the
feature whitelist and draws the folds shared by all later scripts. No
model is fitted.

Time convention (D1). Dates are whole days and day d runs from d to d + 1.
An event on day d is placed at d + 0.5 and event-free patients are
censored at tau + 1, the end of day tau. Horizon h is read as time <=
h + 0.5, so a death on day h counts, as in death_30d, and the C-index is
truncated at h + 1. A cardiovascular readmission and a death on the same
day count as cardiovascular (the dissertation's rule).

Frame A: clock at the index admission date; all EHR+ECG admissions;
horizon 30 days. Endpoints: time_death/event_death, time_first/event_comp
(first of readmission or death) and event_cr (1 = readmission first,
2 = death first, 0 = event-free).
Frame B: clock at the index discharge date; first admission per patient;
horizon 365 days. Excluded: in-hospital deaths (deathtime or discharge
location DIED), dod before discharge, and an in-window ECG after
discharge (every feature must be known when the clock starts). The
three-modality subset is the patients whose first admission has a CTPA,
nested in the EHR+ECG frame. Hospice and "index possibly in the final
data year" are kept as flags for sensitivity analyses.
Folds: StratifiedGroupKFold, 5 folds for each of the six side-project
seeds (42, 7, 13, 1, 2, 3), grouped on subject_id and stratified on
event_cr, drawn separately for each frame and cohort.

Checks (no frame is written if one fails): Frame A matches the
dissertation labels at day 30 and Frame B the event counts of
03_bq_followup.py; times obey D1; whitelisted features exist and none is
a label or timing column; each patient is in one fold per seed.

Inputs:  data/inputs/*.csv, results/bq_followup_horizons.csv
Outputs: results/frames_summary.csv, results/frames_events.csv,
         results/frames_log.txt; data/processed/frame_A.csv, frame_B.csv,
         folds.csv, feature_whitelist.csv (patient level, not in the repo)
Usage:   python scripts/04_frames.py
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from sklearn.model_selection import (
    StratifiedGroupKFold)

T0 = time.time()
pd.set_option("display.width", 200)

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")
for p_ in (PROC, RES):
    os.makedirs(p_, exist_ok=True)

KEY = ["subject_id", "hadm_id"]
TAU_A, TAU_B = 30, 365
HA, HB = [10, 30], [10, 30, 90, 180, 365]
NFOLD = 5
ENDS = ("death", "composite",
        "cv_first_competing", "death_first")
SUM, EVT, FOLDS = [], [], []


def banner(t):
    print("", flush=True)
    print("=" * 70)
    print(t)
    print("=" * 70, flush=True)


def rec(sec, fr, coh, item, value,
        status="INFO"):
    SUM.append({"section": sec, "frame": fr,
                "cohort": coh, "item": item,
                "value": value, "status": status})
    print("  %-4s %-2s %-8s %-48s %s"
          % (status, fr, coh, item, value),
          flush=True)


def flt(s):
    return pd.to_numeric(
        s, errors="coerce").astype("float64")


def code(t_d, t_cv, tau):
    """(time, event) coding. t_d, t_cv are days
    from the clock start (NaN = no event).
    Events at day + 0.5; event-free censored at
    the end of day tau, i.e. tau + 1."""
    de = (t_d >= 0) & (t_d <= tau)
    cv = (t_cv >= 0) & (t_cv <= tau)
    cvf = cv & (~de | (t_cv <= t_d))
    dthf = de & ~cvf
    t_death = np.where(de, t_d + 0.5, tau + 1.0)
    t_first = np.where(
        cvf, t_cv + 0.5,
        np.where(dthf, t_d + 0.5, tau + 1.0))
    e_cr = np.where(cvf, 1, np.where(dthf, 2, 0))
    return (t_death, de.astype(int), t_first,
            e_cr.astype(int),
            (e_cr > 0).astype(int))


def put(df, t_d, t_cv, tau):
    a, b, c, d_, e = code(
        np.asarray(t_d, float),
        np.asarray(t_cv, float), tau)
    df["time_death"], df["event_death"] = a, b
    df["time_first"], df["event_cr"] = c, d_
    df["event_comp"] = e
    return df


def check_times(df, fr, coh, tau):
    ev_t = pd.concat([
        df.loc[df["event_death"] == 1,
               "time_death"],
        df.loc[df["event_comp"] == 1,
               "time_first"]])
    ce_t = pd.concat([
        df.loc[df["event_death"] == 0,
               "time_death"],
        df.loc[df["event_comp"] == 0,
               "time_first"]])
    ok = bool((ev_t <= tau + 0.5).all()
              and (ce_t == tau + 1.0).all())
    rec("check", fr, coh,
        "events <= %.1f and censoring at %.1f"
        % (tau + 0.5, tau + 1.0),
        "yes" if ok else "NO",
        "PASS" if ok else "FAIL")


def events(df, fr, coh, hs):
    out, n = [], len(df)
    for h in hs:
        cut = h + 0.5
        tf = df["time_first"] <= cut
        mk = {"death": (df["event_death"] == 1)
              & (df["time_death"] <= cut),
              "composite": (df["event_comp"] == 1)
              & tf,
              "cv_first_competing":
              (df["event_cr"] == 1) & tf,
              "death_first":
              (df["event_cr"] == 2) & tf}
        for nm in ENDS:
            v = int(mk[nm].sum())
            out.append({"frame": fr, "cohort": coh,
                        "endpoint": nm,
                        "horizon_days": h, "n": n,
                        "events": v,
                        "risk": round(v / n, 4)})
    return out


def make_folds(df, fr, coh):
    y = df["event_cr"].to_numpy()
    grp = df["subject_id"].to_numpy()
    for s in SEEDS:
        sgk = StratifiedGroupKFold(
            n_splits=NFOLD, shuffle=True,
            random_state=s)
        f_ = np.full(len(df), -1)
        for k, (_, te) in enumerate(sgk.split(
                np.zeros(len(df)), y, groups=grp)):
            f_[te] = k
        FOLDS.append(pd.DataFrame({
            "frame": fr, "cohort": coh, "seed": s,
            "subject_id": df["subject_id"].values,
            "hadm_id": df["hadm_id"].values,
            "fold": f_}))
        t_ = pd.DataFrame({"g": grp, "f": f_,
                           "y": y})
        ok = bool((f_ >= 0).all() and t_.groupby(
            "g")["f"].nunique().max() == 1)
        e1 = (t_[t_["y"] == 1].groupby("f").size()
              .reindex(range(NFOLD), fill_value=0))
        e2 = (t_[t_["y"] == 2].groupby("f").size()
              .reindex(range(NFOLD), fill_value=0))
        rec("folds", fr, coh,
            "seed %d: one fold/patient; CV, death"
            " per fold" % s,
            "%s; %d-%d, %d-%d"
            % ("yes" if ok else "NO", e1.min(),
               e1.max(), e2.min(), e2.max()),
            "PASS" if ok else "FAIL")


# ============ LOAD ============
banner("LOAD")
lab = pd.read_csv(os.path.join(
    INP, "mimic_labels_harmonised.csv"))
for c in ("dod_days", "inhosp_days", "cv_days"):
    lab[c] = flt(lab[c])
adm = pd.read_csv(os.path.join(
    INP, "index_admission_times.csv"))
adm = adm.drop_duplicates("hadm_id").copy()
adm["at"] = pd.to_datetime(adm["admittime"],
                           errors="coerce")
adm["dt"] = pd.to_datetime(adm["dischtime"],
                           errors="coerce")
adm["los_d"] = flt((adm["dt"].dt.normalize()
                    - adm["at"].dt.normalize())
                   .dt.days)
adm["los_h"] = ((adm["dt"] - adm["at"])
                .dt.total_seconds() / 3600.0)
coh = pd.read_csv(os.path.join(INP, "cohorts.csv"))
ecgt = pd.read_csv(os.path.join(
    INP, "ecg_timing.csv"))
fu = pd.read_csv(os.path.join(
    INP, "followup_365.csv"))
fu["cv_days_all"] = flt(fu["cv_days_all"])
fu["approx_year_hi"] = flt(fu["approx_year_hi"])
cfg = pd.read_csv(os.path.join(
    INP, "side_config.csv")).set_index("key")
SEEDS = [int(s) for s in
         str(cfg.loc["seeds", "value"]).split(",")]
print("  seeds:", SEEDS)

d = coh.merge(lab, on=KEY, how="left")
d = d.merge(adm[["hadm_id", "at", "los_d",
                 "los_h"]], on="hadm_id", how="left")
d = d.merge(ecgt[["hadm_id", "ecg_h_last"]],
            on="hadm_id", how="left")
d = d.merge(fu[["hadm_id", "discharge_location",
                "cv_days_all", "approx_year_hi"]],
            on="hadm_id", how="left")
d["t_death"] = np.fmin(d["dod_days"],
                       d["inhosp_days"])
print("  rows:", len(d), flush=True)

# ============ FRAME A ============
banner("FRAME A: clock at admission, 30 days")
A = d[d["in_ehr_ecg"] == 1].sort_values(
    ["subject_id", "at", "hadm_id"],
    kind="mergesort").copy()
A["first_adm"] = (~A.duplicated(
    "subject_id")).astype(int)
A = put(A, A["t_death"], A["cv_days"], TAU_A)
for lc, mine in (
        ("death_30d", A["event_death"]),
        ("cv_first", (A["event_cr"] == 1)),
        ("death_first", (A["event_cr"] == 2)),
        ("composite_30d", A["event_comp"])):
    ref = flt(A[lc]).fillna(0).astype(int).values
    nb = int((mine.astype(int).values
              != ref).sum())
    rec("check", "A", "EHR+ECG",
        "%s rows that disagree" % lc, nb,
        "PASS" if nb == 0 else "FAIL")
nb = int((A["time_first"]
          > A["time_death"]).sum())
rec("check", "A", "EHR+ECG",
    "first-event time after death time", nb,
    "PASS" if nb == 0 else "FAIL")
check_times(A, "A", "EHR+ECG", TAU_A)
rec("info", "A", "EHR+ECG", "admissions", len(A))
rec("info", "A", "3-mod", "admissions",
    int((A["in_three_mod"] == 1).sum()))

# ============ FRAME B ============
banner("FRAME B: clock at discharge, 365 days")
f1 = A[A["first_adm"] == 1].copy()
r_ih = f1["inhosp_days"].notna()
r_early = f1["dod_days"] < f1["los_d"]
pre = f1[~r_ih & ~r_early].copy()
pre = put(pre, pre["dod_days"] - pre["los_d"],
          pre["cv_days_all"] - pre["los_d"],
          TAU_B)

# check 2: reproduce the stage 3 counts
rp = os.path.join(RES, "bq_followup_horizons.csv")
if os.path.exists(rp):
    ref = pd.read_csv(rp)
    mine = pd.DataFrame(
        events(pre, "B", "EHR+ECG", HB)
        + events(pre[pre["in_three_mod"] == 1],
                 "B", "3-mod", HB))
    j = mine.merge(ref, on=["cohort", "endpoint",
                            "horizon_days"],
                   suffixes=("", "_ref"))
    nb = int((j["events"]
              != j["events_ref"]).sum())
    rec("check", "B", "both",
        "cells differing from stage 3 (of %d)"
        % len(j), nb,
        "PASS" if nb == 0 and len(j) == 40
        else "FAIL")
else:
    rec("check", "B", "both",
        "stage 3 horizons file", "not found",
        "FAIL")

r_died = pre["discharge_location"].astype(
    str).str.upper().eq("DIED")
r_ecg = pre["ecg_h_last"] > pre["los_h"]
B = pre[~r_died & ~r_ecg].copy()
check_times(B, "B", "EHR+ECG", TAU_B)
rec("excl", "B", "EHR+ECG", "first admissions",
    len(f1))
rec("excl", "B", "EHR+ECG",
    "excluded: in-hospital death (deathtime)",
    int(r_ih.sum()))
rec("excl", "B", "EHR+ECG",
    "excluded: dod before discharge",
    int(r_early.sum()))
rec("excl", "B", "EHR+ECG",
    "excluded: discharge location DIED",
    int(r_died.sum()))
rec("excl", "B", "EHR+ECG",
    "excluded: ECG after discharge",
    int((r_ecg & ~r_died).sum()))
rec("excl", "B", "EHR+ECG", "Frame B cohort",
    len(B))
rec("excl", "B", "3-mod", "Frame B cohort",
    int((B["in_three_mod"] == 1).sum()))

B["hospice"] = B["discharge_location"].astype(
    str).str.upper().str.contains(
    "HOSPICE").astype(int)
grp_end = flt(fu["anchor_year_group"].astype(
    str).str[-4:])
data_end = int(np.nanmax(grp_end))
B["final_year_possible"] = (
    B["approx_year_hi"] >= data_end).astype(int)
for cn, sub in (("EHR+ECG", B),
                ("3-mod", B[B["in_three_mod"]
                            == 1])):
    rec("flags", "B", cn, "hospice discharges",
        int(sub["hospice"].sum()))
    rec("flags", "B", cn,
        "index possibly in final year (>=%d)"
        % data_end,
        int(sub["final_year_possible"].sum()))
td_ = B["dod_days"] - B["los_d"]
tc_ = B["cv_days_all"] - B["los_d"]
rec("info", "B", "EHR+ECG",
    "same-day CV and death within a year (CV)",
    int(((td_ == tc_) & (td_ <= 365)).sum()))

# ============ EVENTS ============
EVT += events(A, "A", "EHR+ECG", HA)
EVT += events(A[A["in_three_mod"] == 1], "A",
              "3-mod", HA)
EVT += events(B, "B", "EHR+ECG", HB)
EVT += events(B[B["in_three_mod"] == 1], "B",
              "3-mod", HB)

# ============ WHITELIST ============
banner("FEATURE WHITELIST AND COVERAGE")
fc = pd.read_csv(os.path.join(
    INP, "feature_columns.csv"))
files = {"ehr": "ehr28.csv", "ecg": "ecg_opt.csv",
         "ctpa": "ctpa46.csv"}
mats = {m: pd.read_csv(os.path.join(INP, f_))
        .drop_duplicates("hadm_id")
        for m, f_ in files.items()}
absent = [c for m, c in zip(fc["modality"],
                            fc["column"])
          if c not in mats[m].columns]
rec("check", "-", "-",
    "whitelisted columns missing from matrix",
    len(absent), "PASS" if not absent else "FAIL")
banned = (set(lab.columns) - set(KEY)) | {
    "time_death", "event_death", "time_first",
    "event_cr", "event_comp", "hospice",
    "final_year_possible", "first_adm"}
words = ("days", "death", "_first", "composite",
         "_30d", "inhosp", "mace", "label",
         "time_", "event_")
hits = [c for c in fc["column"]
        if c in banned
        or any(w in c.lower() for w in words)]
rec("check", "-", "-",
    "label or timing names among %d features"
    % len(fc), len(hits),
    "PASS" if not hits else "FAIL")
if hits:
    print("    offending:", hits[:10])
for m in ("ehr", "ecg", "ctpa"):
    rec("info", "-", "-", "%s feature columns" % m,
        int((fc["modality"] == m).sum()))

for fr, df_ in (("A", A), ("B", B)):
    for cn, sub in (("EHR+ECG", df_),
                    ("3-mod", df_[df_["in_three_mod"]
                                  == 1])):
        for m in ("ehr", "ecg", "ctpa"):
            if m == "ctpa" and cn == "EHR+ECG":
                continue
            cols = fc.loc[fc["modality"] == m,
                          "column"].tolist()
            jj = sub[["hadm_id"]].merge(
                mats[m][["hadm_id"] + cols],
                on="hadm_id", how="left",
                indicator=True)
            has = jj["_merge"] == "both"
            cov = float(has.mean())
            miss = float(jj.loc[has, cols]
                         .isna().mean().mean())
            rec("coverage", fr, cn,
                "%s: rows present / cells missing"
                % m, "%.4f / %.4f" % (cov, miss),
                "PASS" if cov == 1.0 else "WARN")

# ============ FOLDS ============
banner("FOLDS")
make_folds(A, "A", "EHR+ECG")
make_folds(A[A["in_three_mod"] == 1].copy(), "A",
           "3-mod")
make_folds(B, "B", "EHR+ECG")
make_folds(B[B["in_three_mod"] == 1].copy(), "B",
           "3-mod")

# ============ WRITE ============
banner("WRITE")
sm = pd.DataFrame(SUM)
ev = pd.DataFrame(EVT)
sm.to_csv(os.path.join(RES, "frames_summary.csv"),
          index=False)
ev.to_csv(os.path.join(RES, "frames_events.csv"),
          index=False)
nf = int((sm["status"] == "FAIL").sum())
if nf:
    banner("STOPPED: %d hard check(s) failed,"
           " no frame written" % nf)
    print(sm[sm["status"] == "FAIL"][
        ["frame", "cohort", "item", "value"]]
        .to_string(index=False))
    sys.exit(1)

cA = ["subject_id", "hadm_id", "in_three_mod",
      "first_adm", "time_death", "event_death",
      "time_first", "event_cr", "event_comp"]
cB = ["subject_id", "hadm_id", "in_three_mod",
      "hospice", "final_year_possible",
      "time_death", "event_death", "time_first",
      "event_cr", "event_comp"]
out = {"frame_A.csv": A[cA], "frame_B.csv": B[cB],
       "folds.csv": pd.concat(FOLDS,
                              ignore_index=True)}
fc["file"] = fc["modality"].map(files)
out["feature_whitelist.csv"] = fc
for nm, df_ in out.items():
    p = os.path.join(PROC, nm)
    df_.to_csv(p, index=False)
    print("  wrote %-24s %7d x %d"
          % (nm, len(df_), df_.shape[1]))

print("")
print(ev.pivot_table(
    index=["frame", "cohort", "endpoint"],
    columns="horizon_days", values="events")
    .to_string())
print("")
print("elapsed %.1f s" % (time.time() - T0))
print("")
print("OUTPUT FILES")
for nm in out:
    print("  " + os.path.join(PROC, nm)
          + "  (patient-level, local only)")
print("  " + os.path.join(RES, "frames_summary.csv"))
print("  " + os.path.join(RES, "frames_events.csv"))
print("  " + os.path.join(RES, "frames_log.txt"),
      "(terminal log)")
print("")
print("Outputs: frames_summary.csv and"
      " frames_events.csv")