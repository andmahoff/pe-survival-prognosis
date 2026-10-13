"""04b_inspect_frame.py
INSPECT frame and its MIMIC partner for the secondary analysis (D1).

Frame I: INSPECT, 3,300 PE-positive persons, clock at the index CTPA scan
(closest to the PE date, as in 02_audit.py and the side project). Death
comes from INSPECT's native tte_mortality (minutes converted to days) and
its censoring flag, followed to 366 days. Patients who leave Stanford care
are censored when last seen, so censoring is real.
Frame M: the MIMIC partner. First admission per patient, three-modality
subset (nested, as in Frame B), clock at admission, death to 365 days,
event-free patients censored at 366. Follow-up is complete (state
registry).
Both frames use EHR and CTPA only (INSPECT has no ECG), on the columns the
two datasets share.

Time convention (D1). MIMIC is day-level: an event on day d sits at
d + 0.5, event-free patients are censored at tau + 1, horizon h is read as
time <= h + 0.5 and the C-index is truncated at h + 1. INSPECT is
minute-level: horizon h is read as time <= h, which reproduces INSPECT's
own 1-, 6- and 12-month labels (02_audit.py), and the C-index is truncated
just past h. frames_config.csv records this for every frame.

Also: age, condition flags and EHR missingness are compared between
persons who died within a year, were censored within a year and were known
alive at one year, and one-year INSPECT mortality is bounded by assuming
that every early-censored person survived, or died.

Checks (nothing is written if one fails): 3,300 index scans; Kaplan-Meier
risks reproduce 02_audit.py at all five horizons; times obey D1; shared
feature columns exist and none is a label or timing column; each person
is in one fold per seed.

Inputs:  data/inputs/*.csv, data/inputs/inspect/ (copied by 02_audit.py),
         data/processed/frame_A.csv, feature_whitelist.csv,
         results/audit_horizons.csv
Outputs: results/frames_im_summary.csv, results/frames_im_events.csv,
         results/frames_im_log.txt; data/processed/frame_I.csv,
         frame_M.csv, folds_IM.csv, feature_whitelist_IM.csv,
         frames_config.csv (patient level, not in the repo)
Usage:   python scripts/04b_inspect_frame.py
"""
import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from sklearn.model_selection import (
    StratifiedGroupKFold)
from sksurv.nonparametric import (
    kaplan_meier_estimator)

T0 = time.time()
pd.set_option("display.width", 200)

DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
P2 = os.path.join(DS, "SurvivalProject")
INP = os.path.join(P2, "data", "inputs")
INS = os.path.join(INP, "inspect")
PROC = os.path.join(P2, "data", "processed")
RES = os.path.join(P2, "results")

KEY = ["subject_id", "hadm_id"]
H = [10, 30, 90, 180, 365]
TAU = 365
TRUNC = TAU + 1.0      # deviation D1
GAP_D = 30
NFOLD = 5
WORDS = ("days", "death", "_first", "composite",
         "_30d", "inhosp", "mace", "label",
         "time_", "event_", "mortality",
         "readmission", "censor", "tte")
SUM, EVT, FOLDS = [], [], []


def banner(t):
    print("", flush=True)
    print("=" * 70)
    print(t)
    print("=" * 70, flush=True)


def rec(sec, fr, item, value, status="INFO"):
    SUM.append({"section": sec, "frame": fr,
                "item": item, "value": value,
                "status": status})
    print("  %-4s %-2s %-52s %s"
          % (status, fr, item, value), flush=True)


def flt(s):
    return pd.to_numeric(
        s, errors="coerce").astype("float64")


def latest(folder, prefix):
    hits = sorted(glob.glob(os.path.join(
        folder, prefix + "*.tsv")))
    return hits[-1] if hits else None


def km_curve(t, e):
    x, s = kaplan_meier_estimator(
        np.asarray(e, bool), np.asarray(t, float))
    return x, s


def km_risk(x, s, h):
    i = np.searchsorted(x, h, side="right") - 1
    return 0.0 if i < 0 else 1.0 - float(s[i])


def make_folds(ids, groups, y, fr):
    ids = np.asarray(ids)
    groups = np.asarray(groups)
    y = np.asarray(y)
    for sd in SEEDS:
        sgk = StratifiedGroupKFold(
            n_splits=NFOLD, shuffle=True,
            random_state=sd)
        f_ = np.full(len(ids), -1)
        for k, (_, te) in enumerate(sgk.split(
                np.zeros(len(ids)), y,
                groups=groups)):
            f_[te] = k
        FOLDS.append(pd.DataFrame({
            "frame": fr, "seed": sd, "id": ids,
            "group_id": groups, "fold": f_}))
        t_ = pd.DataFrame({"g": groups, "f": f_,
                           "y": y})
        ok = bool((f_ >= 0).all() and t_.groupby(
            "g")["f"].nunique().max() == 1)
        ev = (t_[t_["y"] == 1].groupby("f").size()
              .reindex(range(NFOLD), fill_value=0))
        rec("folds", fr,
            "seed %d: one fold/person; deaths/fold"
            % sd,
            "%s; %d-%d" % ("yes" if ok else "NO",
                           ev.min(), ev.max()),
            "PASS" if ok else "FAIL")


cfg = pd.read_csv(os.path.join(
    INP, "side_config.csv")).set_index("key")
SEEDS = [int(s) for s in
         str(cfg.loc["seeds", "value"]).split(",")]

# ============ FRAME I ============
banner("FRAME I: INSPECT, clock at the index scan")
p_lab = latest(INS, "labels_")
p_map = latest(INS, "study_mapping_")
p_hl = os.path.join(INS, "inspect_labels_final.csv")
p_ic = os.path.join(INS,
                    "inspect_imp_features_v2.csv")
need = [p_lab, p_map, p_hl, p_ic]
if not all(p and os.path.exists(p) for p in need):
    print("  missing INSPECT inputs:", need)
    print("  re-run 02_audit.py, which copies them")
    sys.exit(1)

il = pd.read_csv(p_lab, sep="\t")
if "person_id" in il.columns:
    il = il.drop(columns=["person_id"])
for c in ("tte_mortality", "is_censored_mortality"):
    if c not in il.columns:
        print("  column not found:", c)
        sys.exit(1)
mp = pd.read_csv(p_map, sep="\t")[
    ["impression_id", "person_id",
     "procedure_DATETIME"]]
hl = pd.read_csv(p_hl)[["person_id", "pe_date"]]
ie = pd.read_csv(os.path.join(INP, "inspect_ehr.csv"))
ic = pd.read_csv(p_ic)

g = il.merge(mp, on="impression_id", how="inner")
g["pdt"] = pd.to_datetime(g["procedure_DATETIME"],
                          errors="coerce")
hl["pdate"] = pd.to_datetime(hl["pe_date"],
                             errors="coerce")
g = g[g["person_id"].isin(set(ie["person_id"]))]
g = g.merge(hl[["person_id", "pdate"]],
            on="person_id", how="left")
g["gap"] = (g["pdt"] - g["pdate"]).abs()
g = g[g["gap"] <= pd.Timedelta(days=GAP_D)]
g = (g.sort_values(["person_id", "gap"],
                   kind="mergesort")
     .drop_duplicates("person_id")
     .reset_index(drop=True))
rec("check", "I", "persons with an index scan",
    len(g), "PASS" if len(g) == len(ie) == 3300
    else "FAIL")

t = (flt(g["tte_mortality"]) / 1440.0).to_numpy()
cz = g["is_censored_mortality"].astype(
    str).str.strip().str.lower()
cens = cz.isin(["true", "1", "1.0", "yes"]).to_numpy()
ev = ~cens
bad = np.isnan(t) | (t < 0)
rec("check", "I", "tte missing or negative",
    int(bad.sum()), "PASS" if not bad.any()
    else "FAIL")
dead_tr = ev & (t <= TRUNC)
I = pd.DataFrame({
    "person_id": g["person_id"].values,
    "time_death": np.where(dead_tr, t,
                           np.minimum(t, TRUNC)),
    "event_death": dead_tr.astype(int)})
ok = bool((I["time_death"] <= TRUNC).all()
          and (I.loc[I["event_death"] == 0,
                     "time_death"] <= TRUNC).all())
rec("check", "I", "all times <= %.0f days" % TRUNC,
    "yes" if ok else "NO",
    "PASS" if ok else "FAIL")
rec("info", "I", "deaths between day 365 and 366",
    int((ev & (t > TAU) & (t <= TRUNC)).sum()))

x, s = km_curve(I["time_death"], I["event_death"])
ref = pd.read_csv(os.path.join(
    RES, "audit_horizons.csv"))
ref = ref[(ref["dataset"] == "INSPECT")
          & (ref["method"] == "kaplan_meier")]
n = len(I)
for h in H:
    r_ = km_risk(x, s, h)
    rr = ref.loc[ref["horizon_days"] == h, "risk"]
    ok = len(rr) == 1 and abs(
        round(r_, 4) - float(rr.iloc[0])) < 1e-9
    rec("check", "I",
        "KM risk at %d d vs stage 2" % h,
        "%.4f vs %s" % (r_, rr.iloc[0]
                        if len(rr) else "none"),
        "PASS" if ok else "FAIL")
    EVT.append({"frame": "I", "cohort": "INSPECT",
                "endpoint": "death",
                "horizon_days": h, "n": n,
                "events": int((ev & (t <= h)).sum()),
                "censored_by_h":
                int((cens & (t <= h)).sum()),
                "risk": round(r_, 4),
                "method": "kaplan_meier"})

# informative censoring, the clean comparison
dead1 = ev & (t <= TAU)
grp = np.where(dead1, "died <1y",
               np.where(cens & (t <= TAU),
                        "censored <1y",
                        "alive at 1y"))
wl = pd.read_csv(os.path.join(
    PROC, "feature_whitelist.csv"))
ehr_all = wl.loc[wl["modality"] == "ehr",
                 "column"].tolist()
ehr_c = [c for c in ehr_all if c in ie.columns]
gg = pd.DataFrame({"person_id": I["person_id"],
                   "grp": grp}).merge(
    ie, on="person_id", how="left")
gg["ehr_missing"] = gg[ehr_c].isna().mean(axis=1)
order = ["died <1y", "censored <1y", "alive at 1y"]
rec("censor", "I", "n: " + " / ".join(order),
    " / ".join(str(int((grp == o).sum()))
               for o in order))
for v in ("age", "cancer", "heart_failure",
          "copd", "afib", "ehr_missing"):
    if v not in gg.columns:
        continue
    mv = [flt(gg.loc[gg["grp"] == o, v]).mean()
          for o in order]
    rec("censor", "I", "mean %s: died/cens/alive" % v,
        " / ".join("%.3f" % m_ for m_ in mv))

nd = int(dead1.sum())
nc = int((cens & (t <= TAU)).sum())
rec("bounds", "I",
    "1y mortality: all survive / KM / all die",
    "%.4f / %.4f / %.4f"
    % (nd / n, km_risk(x, s, TAU), (nd + nc) / n))

# ============ FRAME M ============
banner("FRAME M: MIMIC partner, clock at admission")
fa = pd.read_csv(os.path.join(PROC, "frame_A.csv"))
lab = pd.read_csv(os.path.join(
    INP, "mimic_labels_harmonised.csv"))
for c in ("dod_days", "inhosp_days"):
    lab[c] = flt(lab[c])
M = fa[(fa["first_adm"] == 1)
       & (fa["in_three_mod"] == 1)][KEY].merge(
    lab[KEY + ["dod_days", "inhosp_days"]],
    on=KEY, how="left")
td = np.fmin(M["dod_days"], M["inhosp_days"]
             ).to_numpy(float)
dm = (td >= 0) & (td <= TAU)
M["time_death"] = np.where(dm, td + 0.5, TAU + 1.0)
M["event_death"] = dm.astype(int)
M = M[KEY + ["time_death", "event_death"]]
ok = bool((M.loc[M["event_death"] == 1,
                 "time_death"] <= TAU + 0.5).all()
          and (M.loc[M["event_death"] == 0,
                     "time_death"]
               == TAU + 1.0).all())
rec("check", "M", "events <= %.1f and censoring"
    " at %.1f" % (TAU + 0.5, TAU + 1.0),
    "yes" if ok else "NO",
    "PASS" if ok else "FAIL")
rec("info", "M", "first admissions, 3-mod (nested)",
    len(M))
for h in H:
    v = int(((M["event_death"] == 1)
             & (M["time_death"] <= h + 0.5)).sum())
    EVT.append({"frame": "M", "cohort": "3-mod",
                "endpoint": "death",
                "horizon_days": h, "n": len(M),
                "events": v, "censored_by_h": 0,
                "risk": round(v / len(M), 4),
                "method": "empirical"})

# ============ SHARED FEATURES ============
banner("SHARED FEATURES (EHR + CTPA)")
ct = wl[wl["modality"] == "ctpa"]
ct_c = [c for c in ct["column"] if c in ic.columns]
ct_miss = ct[~ct["column"].isin(ic.columns)]
rec("check", "IM", "EHR columns shared",
    "%d of %d" % (len(ehr_c), len(ehr_all)),
    "PASS" if len(ehr_c) == len(ehr_all)
    else "FAIL")
rec("info", "IM", "CTPA columns shared",
    "%d of %d" % (len(ct_c), len(ct)))
if len(ct_miss):
    rec("info", "IM", "CTPA columns not in INSPECT,"
        " by block",
        ct_miss["block"].value_counts().to_dict())
shared = ehr_c + ct_c
hits = [c for c in shared
        if any(w in c.lower() for w in WORDS)]
rec("check", "IM",
    "label or timing names among shared",
    len(hits), "PASS" if not hits else "FAIL")
if hits:
    print("    offending:", hits[:10])

for nm, df_, key in (("EHR", ie, ehr_c),
                     ("CTPA", ic, ct_c)):
    jj = I[["person_id"]].merge(
        df_[["person_id"] + key]
        .drop_duplicates("person_id"),
        on="person_id", how="left",
        indicator=True)
    has = jj["_merge"] == "both"
    rec("coverage", "I",
        "%s: rows present / cells missing" % nm,
        "%.4f / %.4f"
        % (has.mean(),
           jj.loc[has, key].isna().mean().mean()),
        "PASS" if has.all() else "WARN")
mats = {"EHR": ("ehr28.csv", ehr_c),
        "CTPA": ("ctpa46.csv", ct_c)}
for nm, (fn, key) in mats.items():
    mm = pd.read_csv(os.path.join(INP, fn)
                     ).drop_duplicates("hadm_id")
    jj = M[["hadm_id"]].merge(
        mm[["hadm_id"] + key], on="hadm_id",
        how="left", indicator=True)
    has = jj["_merge"] == "both"
    rec("coverage", "M",
        "%s: rows present / cells missing" % nm,
        "%.4f / %.4f"
        % (has.mean(),
           jj.loc[has, key].isna().mean().mean()),
        "PASS" if has.all() else "WARN")

wim = wl[wl["column"].isin(shared)].copy()
wim["inspect_file"] = np.where(
    wim["modality"] == "ehr", "inspect_ehr.csv",
    "inspect\\inspect_imp_features_v2.csv")
wim["mimic_file"] = np.where(
    wim["modality"] == "ehr", "ehr28.csv",
    "ctpa46.csv")

# ============ FOLDS ============
banner("FOLDS")
make_folds(I["person_id"], I["person_id"],
           I["event_death"], "I")
make_folds(M["hadm_id"], M["subject_id"],
           M["event_death"], "M")

# ============ WRITE ============
banner("WRITE")
sm = pd.DataFrame(SUM)
ev_ = pd.DataFrame(EVT)
sm.to_csv(os.path.join(
    RES, "frames_im_summary.csv"), index=False)
ev_.to_csv(os.path.join(
    RES, "frames_im_events.csv"), index=False)
nf = int((sm["status"] == "FAIL").sum())
if nf:
    banner("STOPPED: %d hard check(s) failed,"
           " no frame written" % nf)
    print(sm[sm["status"] == "FAIL"][
        ["frame", "item", "value"]]
        .to_string(index=False))
    sys.exit(1)

# evaluate horizon h at h + auc_cut_offset;
# truncate the C-index strictly below
# h + c_tau_offset
conf = pd.DataFrame([
    ("A", "admission", 30, 31.0, 0.5, 1.0,
     "MIMIC day-level"),
    ("B", "discharge", 365, 366.0, 0.5, 1.0,
     "MIMIC day-level"),
    ("M", "admission", 365, 366.0, 0.5, 1.0,
     "MIMIC day-level"),
    ("I", "ctpa_scan", 365, 366.0, 0.0, 0.001,
     "INSPECT minute-level")],
    columns=["frame", "clock", "tau",
             "censor_time", "auc_cut_offset",
             "c_tau_offset", "note"])
out = {"frame_I.csv": I, "frame_M.csv": M,
       "folds_IM.csv": pd.concat(
           FOLDS, ignore_index=True),
       "feature_whitelist_IM.csv": wim,
       "frames_config.csv": conf}
for nm, df_ in out.items():
    df_.to_csv(os.path.join(PROC, nm), index=False)
    print("  wrote %-26s %7d x %d"
          % (nm, len(df_), df_.shape[1]))

print("")
print(ev_.pivot_table(
    index=["frame", "cohort", "method"],
    columns="horizon_days", values="risk")
    .to_string())
print("")
print("elapsed %.1f s" % (time.time() - T0))
print("")
print("OUTPUT FILES")
for nm in out:
    print("  " + os.path.join(PROC, nm)
          + "  (local only)")
print("  " + os.path.join(
    RES, "frames_im_summary.csv"))
print("  " + os.path.join(
    RES, "frames_im_events.csv"))
print("  " + os.path.join(
    RES, "frames_im_log.txt"),
    "(terminal log)")
print("")
print("Outputs: frames_im_summary.csv and"
      " frames_im_events.csv")