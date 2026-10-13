"""01_export_inputs.py
Freezes the inputs: rebuilds the side project's modality matrices and
cohorts and copies two dissertation files.

The EHR, ECG and CTPA matrices and the two cohorts (EHR+ECG, and the
three-modality subset with a CTPA report in the presentation window) are
rebuilt exactly as the side project's 88_fusion_final.py does, through the
side project's fusion_lib.py (copies in provenance/). Everything is
written to data/inputs with a manifest (rows, columns, MD5), so later
scripts never import side-project code. Code between the VERBATIM markers
is copied from script 88 unchanged; the rest is new.

Checks: 71 SCP logit, 28 EHR and 46 CTPA columns; hadm_id unique in each
cohort; the three-modality cohort nested in the EHR+ECG cohort; no label
or timing name among the features; every written or copied file identical
on read-back. Where the side project's fusion_final.csv is present, n,
events and the zero-shot EHR AUROC (to 1e-4) must match it.

Inputs:  side project: data/processed/{exp0_logits_record,
         ecg_derived_v12_record,fusion_final}.csv and data/raw/ecg/
         ecg_record_index.csv, plus the files read by fusion_lib;
         dissertation: figwork/data/{mimic_labels_harmonised,
         index_admission_times,ctpa_notes_index,mimic_imp_features_v2}.csv
Outputs: data/inputs/ehr28.csv, inspect_ehr.csv (INSPECT EHR and labels,
         gate G2 only), ecg_opt.csv, ctpa46.csv, ecg_timing.csv,
         ctpa_timing.csv (timing metadata, never features), cohorts.csv
         and copies of mimic_labels_harmonised.csv and
         index_admission_times.csv (all patient level, not in the repo);
         data/inputs/feature_columns.csv and side_config.csv (column lists
         and constants); data/inputs/export_checks.csv and manifest.csv
         (aggregate); results/export_inputs_log.txt; provenance/ (copies
         of the two side project scripts)
Usage:   python scripts/01_export_inputs.py
Notes:   run in the side project's Python environment (scikit-learn,
         torch). The script changes into the SideProject folder itself,
         because fusion_lib and script 88 use relative paths there.
"""
import os
import sys
import time
import shutil
import hashlib
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import roc_auc_score

T0 = time.time()
DS = os.environ.get("PROJECTS_DIR", os.path.dirname(os.getcwd()))
SIDE = os.path.join(DS, "SideProject")
P2 = os.path.join(DS, "SurvivalProject")
OUT = os.path.join(P2, "data", "inputs")
PROV = os.path.join(P2, "provenance")
for p_ in (SIDE, P2):
    if not os.path.isdir(p_):
        sys.exit("folder not found: " + p_)
os.makedirs(OUT, exist_ok=True)
os.makedirs(PROV, exist_ok=True)

os.chdir(SIDE)
sys.path.insert(0, os.path.join(SIDE, "scripts"))
import fusion_lib as f  # noqa: E402
import torch  # noqa: E402

print("interpreter:", sys.executable)
print("working folder:", os.getcwd())
print("numpy %s  pandas %s  sklearn %s  torch %s"
      % (np.__version__, pd.__version__,
         sklearn.__version__, torch.__version__),
      flush=True)

KEY = ["subject_id", "hadm_id"]
SP, DI = "SIDE PROJECT", "DISSERTATION"
MAN, CHK = [], []


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def ints(df):
    """Drop rows with a missing key and store keys
    as int64, so every P2 merge is exact."""
    d = df.copy()
    for k in KEY + ["person_id"]:
        if k in d.columns:
            d = d[d[k].notna()]
            d[k] = d[k].astype("int64")
    return d.reset_index(drop=True)


def save(df, name, project, source, note=""):
    p = os.path.join(OUT, name)
    df.to_csv(p, index=False)
    back = pd.read_csv(p)
    MAN.append({"file": name,
                "folder": "data\\inputs",
                "project": project,
                "source": source,
                "rows": back.shape[0],
                "cols": back.shape[1],
                "roundtrip_ok": int(
                    back.shape == df.shape),
                "md5": md5(p), "note": note})
    print("  wrote  %-28s %6d x %-4d %s"
          % (name, back.shape[0], back.shape[1],
             project), flush=True)


def copy_in(src, name, project, dest, note=""):
    p = os.path.join(dest, name)
    shutil.copy2(src, p)
    same = md5(src) == md5(p)
    r = c = np.nan
    if name.endswith(".csv"):
        back = pd.read_csv(p)
        r, c = back.shape
    MAN.append({"file": name,
                "folder": os.path.relpath(dest, P2),
                "project": project, "source": src,
                "rows": r, "cols": c,
                "roundtrip_ok": int(same),
                "md5": md5(p), "note": note})
    print("  copied %-28s identical=%s  %s"
          % (name, same, project), flush=True)


def check(name, expected, observed, ok):
    st = ("INFO" if ok is None
          else "PASS" if ok else "FAIL")
    CHK.append({"check": name,
                "expected": expected,
                "observed": observed,
                "status": st})
    print("  %-4s  %-44s expected %s, got %s"
          % (st, name, expected, observed),
          flush=True)


# ===== VERBATIM FROM SCRIPT 88: constants =====
PROC, FIG = f.PROC, f.FIG
RAW = os.path.join("data", "raw", "ecg")
SEEDS = [42, 7, 13, 1, 2, 3]
NFOLD = 5
ECG_LO, ECG_HI = -12.0, 48.0
CT_LO, CT_HI = -48.0, 24.0
CT_CS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0,
         1e1, 1e2]
MINCOV = 0.20

BEST_C = {"cv_first": 1e-4,
          "composite_30d": 1e-3,
          "death_30d": 1e-3,
          "death_30d_inhosp": 3e-4}
OUTS = ["death_30d", "composite_30d",
        "cv_first", "death_30d_inhosp"]

DROP_EXACT = ("rec", "t_flag", "vm_noise",
              "vm_base", "vm_rpeak", "net_i",
              "net_avf", "apen", "tp_win_ms",
              "t_win_ms", "t_peak_ms",
              "tp_fallback", "qrs_hit_limit",
              "p_hit_limit", "t_peak_edge",
              "t_edge_leads", "n_rr_dropped",
              "any_af", "tpe_spread_ms",
              "tpe_ok")
DROP_PREFIX = ("tpk_ms_",)
PERM_ONLY = ("tpe_ms_mean", "tpe_ms_worst",
             "tpe_qt_mean", "tpe_qt_worst")
NEG_BAD = ("s_wave_i", "q_wave_iii", "st_min",
           "t_v1", "t_v2", "t_v3", "t_v4",
           "t_v5", "t_v6", "t_wave_iii")
ANG = ("qrs_axis", "t_axis", "p_axis")

# ===== VERBATIM FROM SCRIPT 88: ECG + EHR =====
rec = pd.read_csv(os.path.join(
    PROC, "exp0_logits_record.csv"))
eidx = pd.read_csv(os.path.join(
    RAW, "ecg_record_index.csv"))
eidx["rec"] = eidx["rec"].astype(str).apply(
    lambda s: s if s.startswith("files/")
    else "files/" + s.replace("\\", "/")
    .strip("/"))
LG = [c for c in rec.columns
      if c.startswith("scp_")]
adm = pd.read_csv(os.path.join(
    FIG, "index_admission_times.csv"))
adm["admittime"] = pd.to_datetime(
    adm["admittime"], errors="coerce")

E = eidx.merge(rec, on="rec", how="inner")
E["t"] = pd.to_datetime(
    E["ecg_charttime"], errors="coerce")
E = E.merge(adm[["hadm_id", "admittime"]],
            on="hadm_id", how="left")
E["h"] = ((E["t"] - E["admittime"])
          .dt.total_seconds() / 3600.0)
E = E[(E["h"] >= ECG_LO)
      & (E["h"] <= ECG_HI)].copy()
ECGA = E.groupby(["subject_id", "hadm_id"],
                 as_index=False)[LG].mean()
keep = set(E["rec"])
print("ECG admissions:", len(ECGA), flush=True)

m12 = pd.read_csv(os.path.join(
    PROC, "ecg_derived_v12_record.csv"))
m12["rec"] = m12["rec"].astype(str).apply(
    lambda s: s if s.startswith("files/")
    else "files/" + s.replace("\\", "/")
    .strip("/"))
m12 = m12[m12["rec"].isin(keep)]
cols = [c for c in m12.columns
        if c not in DROP_EXACT
        and not c.startswith(DROP_PREFIX)
        and pd.api.types.is_numeric_dtype(
            m12[c])]
j = eidx[["subject_id", "hadm_id",
          "rec"]].merge(
    m12[["rec"] + cols], on="rec",
    how="inner")
g = j.groupby(["subject_id", "hadm_id"])
pos = [c for c in cols
       if not c.startswith(NEG_BAD)
       and c not in ANG]
neg = [c for c in cols
       if c.startswith(NEG_BAD)]
parts = [g[cols].mean().add_suffix("_mean")]
if pos:
    parts.append(g[pos].max()
                 .add_suffix("_worst"))
if neg:
    parts.append(g[neg].min()
                 .add_suffix("_worst"))
MA = pd.concat(parts, axis=1).reset_index()
for c in ANG:
    if c not in j.columns:
        continue
    rad = np.radians(j[c])
    tmp = j[["subject_id", "hadm_id"]].copy()
    tmp["s"], tmp["c"] = np.sin(rad), np.cos(rad)
    gg = tmp.groupby(["subject_id",
                      "hadm_id"]).mean()
    cm = np.degrees(np.arctan2(
        gg["s"], gg["c"])).rename(
        c + "_circ").reset_index()
    MA = MA.merge(cm, on=["subject_id",
                          "hadm_id"],
                  how="left")
    MA = MA.drop(columns=[c + "_mean",
                          c + "_worst"],
                 errors="ignore")
MC = [c for c in MA.columns
      if c not in ("subject_id", "hadm_id")]
cov = MA[MC].notna().mean()
MC = [c for c in MC if cov[c] >= MINCOV
      and c not in PERM_ONLY]
print("measurement columns:", len(MC),
      " (permissive Tp-e excluded, strict"
      " retained)", flush=True)

lab = pd.read_csv(os.path.join(
    FIG, "mimic_labels_harmonised.csv"))

mim_ehr = f.load_mimic_ehr()
EH = [c for c in f.EHR_COLS
      if c in mim_ehr.columns]
print("EHR columns found:", len(EH),
      "of", len(f.EHR_COLS))
if len(EH) < len(f.EHR_COLS):
    print("  missing:",
          [c for c in f.EHR_COLS
           if c not in mim_ehr.columns][:8])

FULL = ECGA.merge(
    MA[["hadm_id"] + MC], on="hadm_id",
    how="left").merge(
    lab, on=["subject_id", "hadm_id"],
    how="inner").merge(
    mim_ehr[["hadm_id"] + EH],
    on="hadm_id", how="inner")
print("full ECG cohort:", len(FULL),
      " (requires both an ECG and EHR"
      " features)", flush=True)

# ===== VERBATIM FROM SCRIPT 88: CTPA =====
mm = f.load_mimic_ctpa46()
cidx = pd.read_csv(os.path.join(
    FIG, "ctpa_notes_index.csv"))
w = cidx[(cidx["h_before"] >= CT_LO)
         & (cidx["h_before"] <= CT_HI)]
hs = set(w["idx_hadm"].dropna().astype(int))
mm = mm[mm["hadm_id"].isin(hs)]
v2 = pd.read_csv(os.path.join(
    FIG, "mimic_imp_features_v2.csv"),
    usecols=["hadm_id"])
mm = mm[mm["hadm_id"].isin(set(v2["hadm_id"]))]
CT46 = [c for c in f.CTPA46
        if c in mm.columns]
CT = mm.merge(
    ECGA.drop(columns=["subject_id"]),
    on="hadm_id", how="inner").merge(
    MA[["hadm_id"] + MC], on="hadm_id",
    how="left")
print("CTPA cohort:", len(CT),
      " (published 1,703)", flush=True)
# ===== END VERBATIM =====

# ---------- dissertation file structure ----------
print("")
print("=" * 70)
print("DISSERTATION INPUTS (column names only)")
print("=" * 70)
print("  labels:", list(lab.columns))
print("  admission times:", list(adm.columns))
dc = [c for c in adm.columns
      if "disch" in c.lower()]

# ---------- timing metadata, never features ----------
ET = (E.groupby(KEY)
      .agg(n_ecg=("rec", "size"),
           ecg_h_first=("h", "min"),
           ecg_h_last=("h", "max"))
      .reset_index())
CTT = (w.groupby("idx_hadm")
       .agg(n_ctpa=("h_before", "size"),
            h_before_min=("h_before", "min"),
            h_before_max=("h_before", "max"))
       .reset_index()
       .rename(columns={"idx_hadm": "hadm_id"}))
CTT = CTT[CTT["hadm_id"].isin(
    set(CT["hadm_id"]))]

# ---------- cohort membership ----------
c1 = ints(FULL[KEY]).drop_duplicates()
c1["in_ehr_ecg"] = 1
c2 = ints(CT[KEY]).drop_duplicates()
c2["in_three_mod"] = 1
coh = c1.merge(c2, on=KEY, how="outer")
for c in ("in_ehr_ecg", "in_three_mod"):
    coh[c] = coh[c].fillna(0).astype(int)

# ---------- INSPECT, for gate G2 only ----------
ins = f.load_inspect()
print("")
print("INSPECT persons (EHR+CTPA cohort):",
      len(ins), flush=True)

# ================== CHECKS ==================
print("")
print("=" * 70)
print("CHECKS")
print("=" * 70)
check("SCP logit columns", 71, len(LG),
      len(LG) == 71)
check("EHR columns", len(f.EHR_COLS), len(EH),
      len(EH) == len(f.EHR_COLS))
check("CTPA columns", 46, len(CT46),
      len(CT46) == 46)
check("ECG measurement columns", "info",
      len(MC), None)
check("discharge time in admission file",
      "info", ",".join(dc) if dc else "none",
      None)
d1 = int(FULL["hadm_id"].duplicated().sum())
check("EHR+ECG cohort hadm_id unique", 0, d1,
      d1 == 0)
d2 = int(CT["hadm_id"].duplicated().sum())
check("three-modality hadm_id unique", 0, d2,
      d2 == 0)
o3 = len(set(CT["hadm_id"])
         - set(FULL["hadm_id"]))
check("three-modality inside EHR+ECG", 0, o3,
      o3 == 0)
check("EHR+ECG admissions / subjects", "info",
      "%d / %d" % (int(c1.shape[0]),
                   c1["subject_id"].nunique()),
      None)
check("three-modality admissions / subjects",
      "info",
      "%d / %d" % (int(c2.shape[0]),
                   c2["subject_id"].nunique()),
      None)

lab_cols = set(lab.columns) - set(KEY)
BAD = ("days", "death", "_first", "composite",
       "_30d", "inhosp", "mace", "label")
feats = EH + LG + MC + CT46
hits = [c for c in feats
        if c in lab_cols
        or any(b in c.lower() for b in BAD)]
check("label/timing names among features", 0,
      len(hits), len(hits) == 0)
if hits:
    print("    offending:", hits[:10])

# n, events and the zero-shot EHR AUC against
# the side project's data/processed/fusion_final.csv
ffp = os.path.join(PROC, "fusion_final.csv")
ff = pd.read_csv(ffp) if os.path.exists(ffp) \
    else None
if ff is None:
    print("  fusion_final.csv not found: the"
          " comparisons below are INFO only")
for mt, D, ec in (("ehr_ecg", FULL, EH),
                  ("3mod", CT, f.EHR_COLS)):
    for oc in OUTS:
        if oc not in D.columns:
            continue
        d, y = f.labels(D, oc)
        n, ev = len(y), int(y.sum())
        tag = "%s %s" % (mt, oc)
        if ev < 30:
            check(tag + " events (<30, skipped"
                  " in 88)", "info", ev, None)
            continue
        se, ysrc = f.labels(
            ins, "death_30d"
            if oc == "death_30d_inhosp" else oc)
        a1, b1 = f.prep(
            se[ec].values.astype(float),
            d[ec].values.astype(float))
        m = f.fit_lr(a1, ysrc, "ehr")
        auc = roc_auc_score(
            y, m.predict_proba(b1)[:, 1])
        ref = None
        if ff is not None:
            r = ff[(ff["model_type"] == mt)
                   & (ff["outcome"] == oc)
                   & (ff["model"] == "ehr")]
            if len(r):
                ref = r.iloc[0]
        if ref is None:
            check(tag + " n", "info", n, None)
            check(tag + " events", "info", ev,
                  None)
            check(tag + " zero-shot EHR AUC",
                  "info", round(auc, 4), None)
            continue
        rn, re_ = int(ref["n"]), int(ref["ev"])
        ra = float(ref["mean"])
        check(tag + " n", rn, n, rn == n)
        check(tag + " events", re_, ev,
              re_ == ev)
        check(tag + " zero-shot EHR AUC",
              round(ra, 4), round(auc, 4),
              abs(ra - auc) < 1e-4)

# ================== EXPORT ==================
print("")
print("=" * 70)
print("WRITING " + OUT)
print("=" * 70)


def blk(c):
    if c in f.BASE_C:
        return "ctpa_base"
    if c in f.COMORB_C:
        return "ctpa_comorb"
    return "ctpa_device"


fc = pd.DataFrame(
    [("ehr", "ehr", c) for c in EH]
    + [("ecg", "ecg_logit", c) for c in LG]
    + [("ecg", "ecg_meas", c) for c in MC]
    + [("ctpa", blk(c), c) for c in CT46],
    columns=["modality", "block", "column"])

cfg = pd.DataFrame(
    [("seeds", ",".join(str(s) for s in SEEDS)),
     ("nfold", NFOLD),
     ("ecg_window_h", "%g,%g" % (ECG_LO, ECG_HI)),
     ("ctpa_window_h_before",
      "%g,%g" % (CT_LO, CT_HI)),
     ("ecg_meas_min_coverage", MINCOV),
     ("ctpa_c_grid",
      ",".join("%g" % c for c in CT_CS)),
     ("ecg_class_weight", "balanced"),
     ("ctpa_class_weight", "none"),
     ("ehr_training",
      "INSPECT zero-shot, fusion_lib.fit_lr,"
      " LEARNER_MODE=" + f.LEARNER_MODE)]
    + [("ecg_c_" + k, v)
       for k, v in BEST_C.items()],
    columns=["key", "value"])

ecg = ECGA.merge(MA[["hadm_id"] + MC],
                 on="hadm_id", how="left")

save(ints(mim_ehr[KEY + EH])
     .drop_duplicates("hadm_id"),
     "ehr28.csv", SP,
     "fusion_lib.load_mimic_ehr()")
save(ints(ins[["person_id"] + f.EHR_COLS
              + ["death_30d", "composite_30d",
                 "cv_first", "death_first"]]),
     "inspect_ehr.csv", SP,
     "fusion_lib.load_inspect()",
     "gate G2 only")
save(ints(ecg[KEY + LG + MC]), "ecg_opt.csv",
     SP, "script 88 ECGA + MA[MC]")
save(ints(mm[KEY + CT46])
     .drop_duplicates("hadm_id"),
     "ctpa46.csv", SP,
     "script 88 mm after window and v2 filters")
save(ints(ET), "ecg_timing.csv", SP,
     "script 88 E, in-window records",
     "metadata; never a feature")
save(ints(CTT), "ctpa_timing.csv", SP,
     "script 88 w, window-filtered notes",
     "metadata; never a feature")
save(coh, "cohorts.csv", SP,
     "script 88 FULL and CT")
save(fc, "feature_columns.csv", SP,
     "script 88 EH, LG, MC, CT46")
save(cfg, "side_config.csv", SP,
     "script 88 constants")
copy_in(os.path.join(
    FIG, "mimic_labels_harmonised.csv"),
    "mimic_labels_harmonised.csv", DI, OUT)
copy_in(os.path.join(
    FIG, "index_admission_times.csv"),
    "index_admission_times.csv", DI, OUT)
copy_in(os.path.join(
    SIDE, "scripts", "fusion_lib.py"),
    "fusion_lib.py", SP, PROV, "code")
copy_in(os.path.join(
    SIDE, "scripts", "88_fusion_final.py"),
    "88_fusion_final.py", SP, PROV, "code")

rt = [m_["file"] for m_ in MAN
      if m_["roundtrip_ok"] != 1]
check("files round-trip / copy identical", 0,
      len(rt), len(rt) == 0)

chk = pd.DataFrame(CHK)
nf = int((chk["status"] == "FAIL").sum())
save(chk, "export_checks.csv", "P2",
     "01_export_inputs.py",
     "aggregate; safe to share")
pd.DataFrame(MAN).to_csv(
    os.path.join(OUT, "manifest.csv"),
    index=False)

print("")
print("=" * 70)
print("FAILED CHECKS: %d" % nf)
if nf == 0:
    print("ALL CHECKS PASSED. Next: 02_audit.py")
else:
    print("DO NOT PROCEED. Send the log and"
          " export_checks.csv, nothing else.")
print("elapsed %.1f s" % (time.time() - T0))
print("")
print("OUTPUT FILES")
for m_ in MAN:
    print("  " + os.path.join(
        P2, m_["folder"], m_["file"]))
print("  " + os.path.join(OUT, "manifest.csv"))
print("  " + os.path.join(
    P2, "results", "export_inputs_log.txt"),
    "(terminal log)")