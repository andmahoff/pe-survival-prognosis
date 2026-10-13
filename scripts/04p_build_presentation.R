# 04p_build_presentation.R
# At-diagnosis frame (Frame P) with EHR, ECG and CTPA modalities and sPESI
# and PESI-10 scores (PREREGISTRATION_ADDENDUM_2.md A32, A36; DEVIATIONS.md
# D17, D19).
#
# Frame P takes Frame A's admissions with time zero at admission. The
# outcome is death within 365 days of admission (in-hospital deaths
# included), time = days + 0.5, censored at 366. Deaths and times must
# reproduce Frame A's to day 30 exactly, or the script stops. Version 3
# (A36, D19): chronic history also counts from the index admission, and ECG
# timing uses real acquisition times. A36 makes this the main at-diagnosis
# model.
#
# EHR: age, sex, the first vital signs and labs from 24 h before to 24 h
# after admission (stage 3e), an indicator for no presentation vital signs
# (D17), counts of earlier admissions and codes, and Charlson and Elixhauser
# categories (1% to 99% prevalence) with both scores. A chronic condition
# counts if coded on an earlier or the index admission, except categories
# the index stay can cause or that describe the PE itself (the ACUTE
# vector), which use earlier admissions only. ehrP_sens.csv applies that
# stricter rule to every category, for the sensitivity analysis.
#
# ECG (D19, replacing D18): the first ECG from 24 h before to 24 h after
# admission, timed with the real acquisition times from stage 0r (ECGs
# sharing that time are averaged). Single-ECG admissions must reproduce
# ecg_opt.csv, or the script stops. CTPA: ctpa46.csv features for reports
# from 48 h before to 24 h after admission. Scores: sPESI and PESI-10, with
# missing vital signs scored 0, on the main or the strict history.
#
# Inputs:  data/processed/frame_A.csv, folds.csv, audit/ecg_times_real.csv,
#          presentation/frame_P_outcomes.csv; the data/inputs/ tables for
#          labels, demographics, EHR, presentation vitals and labs,
#          diagnosis codes, ECG (ecg_opt.csv) and CTPA (ctpa46.csv); per-ECG
#          files of the side project (pe-fusion-transfer);
#          ctpa_notes_index.csv of the dissertation
#          (pe-multimodal-prognosis-msc)
# Outputs: results/presentation_checks.csv, presentation_summary.csv,
#          presentation_history.csv, presentation_ecg_sources.csv,
#          presentation_log.txt; data/processed/presentation/frame_P.csv,
#          frame_P_outcomes.csv (flags updated), ehrP.csv, ehrP_sens.csv,
#          ecgP.csv, ctpaP.csv, scores_P.csv (patient level, not in the
#          repo)
# Usage:   Rscript scripts/04p_build_presentation.R

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("comorbidity", quietly = TRUE)) {
  stop("package comorbidity not installed")
}
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
ROOT <- dirname(P2)
SP <- file.path(ROOT, "SideProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
OD <- file.path(PR, "presentation")
RS <- file.path(P2, "results")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
WIN_ECG <- c(-24, 24)
WIN_CTPA <- c(-48, 24)
TAU <- 365
PREV <- c(0.01, 0.99)
# categories counted from earlier admissions only (A36)
ACUTE <- c("mi", "cevd", "hp", "pud", "carit", "valv", "pcd", "para",
           "ond", "coag", "fed", "blane", "wloss")
CHK <- list()
SUMM <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}
smy <- function(block, item, value) {
  SUMM[[length(SUMM) + 1]] <<- data.table(block = block, item = item,
                                          value = as.character(value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "presentation_checks.csv"))
    cat("\nSTOP: a check failed. see: presentation_checks.csv\n")
    quit(status = 1)
  }
}
flt <- function(x) suppressWarnings(as.numeric(x))
pc <- function(x) sprintf("%.1f%%", 100 * x)
tparse <- function(x) {
  x <- sub("T", " ", as.character(x), fixed = TRUE)
  as.POSIXct(x, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%OS",
                                           "%Y-%m-%d %H:%M", "%Y-%m-%d"),
             optional = TRUE)
}
need <- c(file.path(PR, c("frame_A.csv", "folds.csv")),
          file.path(PR, "audit", "ecg_times_real.csv"),
          file.path(IN, c("mimic_labels_harmonised.csv", "demographics.csv",
                          "ehr28.csv", "ehr_extra.csv",
                          "pres_labs_first.csv", "pres_vitals_first.csv",
                          "dx_prior.csv", "dx_index_A.csv", "ecg_opt.csv",
                          "ctpa46.csv")),
          file.path(SP, "data", "processed", "exp0_logits_record.csv"),
          file.path(SP, "data", "raw", "ecg", "ecg_record_index.csv"),
          file.path(ROOT, "Dissertation", "figwork", "data",
                    "ctpa_notes_index.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()

# ---------------- frame P ----------------
cat("=== FRAME P ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
lab <- unique(fread(file.path(IN, "mimic_labels_harmonised.csv")),
              by = "hadm_id")
L <- lab[match(fa$hadm_id, hadm_id)]
dd <- flt(L$dod_days)
ih <- flt(L$inhosp_days)
d <- ifelse(is.na(dd), ih, dd)
ev <- as.integer(!is.na(d) & d >= 0 & d <= TAU)
fp <- data.table(hadm_id = fa$hadm_id, subject_id = fa$subject_id,
                 in_three_mod = fa$in_three_mod,
                 time_death = ifelse(ev == 1, d + 0.5, TAU + 1),
                 event_death = ev)
a30 <- fa$event_death == 1
p30 <- fp$event_death == 1 & fp$time_death <= 30.5
ok <- identical(a30, p30) &&
  isTRUE(all.equal(fp$time_death[a30], fa$time_death[a30]))
chk("Frame A's deaths and times to day 30 reproduced",
    sprintf("%d vs %d deaths", sum(p30), sum(a30)),
    if (ok) "PASS" else "FAIL")
stop_if_fail()
chk("deaths within 365 days of admission", sum(fp$event_death), "INFO")
smy("frame", "admissions / deaths by 30 / by 365 days",
    sprintf("%d / %d / %d", nrow(fp), sum(p30), sum(fp$event_death)))

# ---------------- EHR at presentation ----------------
cat("\n=== EHR AT PRESENTATION ===\n")
key <- fp$hadm_id
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
ex <- unique(fread(file.path(IN, "ehr_extra.csv")), by = "hadm_id")
dem <- unique(fread(file.path(IN, "demographics.csv")), by = "hadm_id")
pl <- unique(fread(file.path(IN, "pres_labs_first.csv")), by = "hadm_id")
pv <- unique(fread(file.path(IN, "pres_vitals_first.csv")),
             by = "hadm_id")
vc <- grep("_first$", names(pv), value = TRUE)
lc <- grep("_first$", names(pl), value = TRUE)
V <- as.matrix(pv[match(key, hadm_id), ..vc])
storage.mode(V) <- "double"
has_v <- rowSums(!is.na(V)) > 0
chk("admissions with presentation vital signs (D17)", pc(mean(has_v)),
    "INFO")
smy("EHR", "admissions with presentation vital signs", pc(mean(has_v)))
LB <- as.matrix(pl[match(key, hadm_id), ..lc])
storage.mode(LB) <- "double"
smy("EHR", "first labs: median coverage (registered 18 / extra)",
    sprintf("%.1f%% / %.1f%%", 100 * median(colMeans(!is.na(LB[, 1:18]))),
            100 * median(colMeans(!is.na(LB[, -(1:18), drop = FALSE])))))
dx_codes <- function(f) {
  x <- fread(f, colClasses = list(character = "icd_code"))
  x[, icd_code := gsub("[^A-Za-z0-9]", "", toupper(icd_code))]
  x[hadm_id %in% key]
}
cmb <- function(dx, map9, map10) {
  out <- NULL
  for (v in c(9, 10)) {
    d0 <- as.data.frame(dx[icd_version == v, .(hadm_id, icd_code)])
    if (!nrow(d0)) next
    r <- comorbidity::comorbidity(d0, id = "hadm_id", code = "icd_code",
                                  map = if (v == 9) map9 else map10,
                                  assign0 = FALSE)
    out <- rbind(out, as.data.frame(r))
  }
  out <- as.data.table(out)[, lapply(.SD, max), by = hadm_id]
  out <- out[match(key, hadm_id)]
  out[, hadm_id := key]
  for (k in setdiff(names(out), "hadm_id")) {
    set(out, which(is.na(out[[k]])), k, 0)
  }
  out
}
scr <- function(D, map, w) {
  as.numeric(comorbidity::score(
    structure(as.data.frame(D), map = map,
              class = c("comorbidity", "data.frame")),
    weights = w, assign0 = TRUE))
}
dxp <- dx_codes(file.path(IN, "dx_prior.csv"))
dxi <- dx_codes(file.path(IN, "dx_index_A.csv"))
# acute right heart failure on the index admission (I50.81x except
# I50.812, chronic right heart failure) is usually the PE's own effect
rhf <- dxi$icd_version == 10 & grepl("^I5081", dxi$icd_code) &
  dxi$icd_code != "I50812"
chk("index codes for acute right heart failure not counted (admissions)",
    length(unique(dxi$hadm_id[rhf])), "INFO")
dxi <- dxi[!rhf]
CHp <- cmb(dxp, "charlson_icd9_quan", "charlson_icd10_quan")
ELp <- cmb(dxp, "elixhauser_icd9_quan", "elixhauser_icd10_quan")
CHi <- cmb(dxi, "charlson_icd9_quan", "charlson_icd10_quan")
ELi <- cmb(dxi, "elixhauser_icd9_quan", "elixhauser_icd10_quan")
chk("admissions with earlier-admission codes",
    pc(mean(key %in% dxp$hadm_id)), "INFO")
chk("admissions with index-admission codes",
    pc(mean(key %in% dxi$hadm_id)), if (mean(key %in% dxi$hadm_id) >= 0.95)
      "PASS" else "WARN")
# A36: chronic categories from earlier OR index codes; ACUTE earlier only
merge36 <- function(Pp, Ii) {
  out <- copy(Pp)
  for (k in setdiff(names(Pp), "hadm_id")) {
    if (k %in% ACUTE || !(k %in% names(Ii))) next
    set(out, j = k, value = pmax(Pp[[k]], Ii[[k]]))
  }
  out
}
CHa <- merge36(CHp, CHi)
ELa <- merge36(ELp, ELi)
mk_ind <- function(CH, EL) {
  cbind(setNames(as.data.frame(CH)[, -1], paste0("hx_ch_", names(CH)[-1])),
        setNames(as.data.frame(EL)[, -1], paste0("hx_el_", names(EL)[-1])))
}
indA <- mk_ind(CHa, ELa)
indS <- mk_ind(CHp, ELp)
pvA <- colMeans(indA)
pvS <- colMeans(indS)
kp <- pvA >= PREV[1] & pvA <= PREV[2]
HXA <- cbind(as.matrix(indA[, kp, drop = FALSE]),
             hx_charlson_score = scr(CHa, "charlson_icd10_quan",
                                     "charlson"),
             hx_elixhauser_vw = scr(ELa, "elixhauser_icd10_quan", "vw"))
HXS <- cbind(as.matrix(indS[, kp, drop = FALSE]),
             hx_charlson_score = scr(CHp, "charlson_icd10_quan",
                                     "charlson"),
             hx_elixhauser_vw = scr(ELp, "elixhauser_icd10_quan", "vw"))
chk("history categories kept (1%-99%, main history)", sum(kp), "INFO")
cat_of <- sub("^hx_(ch|el)_", "", names(indA))
HIST <- data.table(
  column = names(indA),
  index = ifelse(startsWith(names(indA), "hx_ch_"), "Charlson",
                 "Elixhauser"),
  category = cat_of,
  rule = ifelse(cat_of %in% ACUTE, "earlier admissions only",
                "earlier or index admission"),
  prev_strict = round(pvS, 4), prev_main = round(pvA, 4),
  kept = kp)
fwrite(HIST, file.path(RS, "presentation_history.csv"))
age <- flt(e28$age[match(key, e28$hadm_id)])
sexf <- as.numeric(toupper(dem$gender[match(key, dem$hadm_id)]) == "F")
npr <- cbind(n_prior = flt(ex$n_prior[match(key, ex$hadm_id)]),
             n_codes = flt(ex$n_codes[match(key, ex$hadm_id)]))
EHRP <- data.table(hadm_id = key, age = age, sex_f = sexf, V,
                   no_pres_vitals = as.numeric(!has_v), LB, HXA, npr)
fwrite(EHRP, file.path(OD, "ehrP.csv"))
EHRS <- data.table(hadm_id = key, age = age, sex_f = sexf, V,
                   no_pres_vitals = as.numeric(!has_v), LB, HXS, npr)
fwrite(EHRS, file.path(OD, "ehrP_sens.csv"))
chk("EHR at presentation: features (main = strict columns)",
    sprintf("%d (%s)", ncol(EHRP) - 1,
            if (identical(names(EHRP), names(EHRS))) "same" else "DIFFER"),
    if (identical(names(EHRP), names(EHRS))) "PASS" else "FAIL")
stop_if_fail()
smy("EHR", "features", ncol(EHRP) - 1)
canc_a <- pmax(CHa$canc, CHa$metacanc)
canc_s <- pmax(CHp$canc, CHp$metacanc)
smy("EHR", "cancer / heart failure / chronic lung disease (main)",
    sprintf("%s / %s / %s", pc(mean(canc_a)), pc(mean(CHa$chf)),
            pc(mean(CHa$cpd))))
smy("EHR", "the same, strict (earlier admissions only)",
    sprintf("%s / %s / %s", pc(mean(canc_s)), pc(mean(CHp$chf)),
            pc(mean(CHp$cpd))))
smy("EHR", "Charlson score median, main / strict",
    sprintf("%.0f / %.0f", median(HXA[, "hx_charlson_score"]),
            median(HXS[, "hx_charlson_score"])))
chk("cancer / heart failure / chronic lung disease, main vs strict",
    sprintf("%s vs %s / %s vs %s / %s vs %s", pc(mean(canc_a)),
            pc(mean(canc_s)), pc(mean(CHa$chf)), pc(mean(CHp$chf)),
            pc(mean(CHa$cpd)), pc(mean(CHp$cpd))), "INFO")

# ---------------- scores ----------------
vv <- function(nm) V[, paste0(nm, "_first")]
flag <- function(x, f) as.numeric(!is.na(x) & f(x))
hr110 <- flag(vv("hr"), function(x) x >= 110)
sbp100 <- flag(vv("sbp"), function(x) x < 100)
spo90 <- flag(vv("spo2"), function(x) x < 90)
rr30 <- flag(vv("rr"), function(x) x >= 30)
t36 <- flag(vv("temp"), function(x) x < 36)
mk <- function(canc, chf, cpd) {
  spesi <- (age > 80) + canc + pmax(chf, cpd) + hr110 + sbp100 + spo90
  pesi <- age + 10 * (1 - sexf) + 30 * canc + 10 * chf + 10 * cpd +
    20 * hr110 + 30 * sbp100 + 20 * rr30 + 20 * t36 + 20 * spo90
  list(spesi = spesi, pesi = pesi)
}
s1 <- mk(canc_a, CHa$chf, CHa$cpd)
s2 <- mk(canc_s, CHp$chf, CHp$cpd)
SC <- data.table(hadm_id = key, spesi = s1$spesi, pesi10 = s1$pesi,
                 spesi_sens = s2$spesi, pesi10_sens = s2$pesi,
                 has_pres_vitals = as.numeric(has_v))
fwrite(SC, file.path(OD, "scores_P.csv"))
smy("scores", "sPESI 0 / 1 / 2+ (main)",
    sprintf("%s / %s / %s", pc(mean(SC$spesi == 0)),
            pc(mean(SC$spesi == 1)), pc(mean(SC$spesi >= 2))))
smy("scores", "sPESI 0 / 1 / 2+ (strict)",
    sprintf("%s / %s / %s", pc(mean(SC$spesi_sens == 0)),
            pc(mean(SC$spesi_sens == 1)), pc(mean(SC$spesi_sens >= 2))))
smy("scores", "PESI-10 median (IQR), main", sprintf(
  "%.0f (%.0f-%.0f)", median(SC$pesi10), quantile(SC$pesi10, 0.25),
  quantile(SC$pesi10, 0.75)))
chk("scores computed (missing)", sum(!is.finite(as.matrix(SC[, 2:5]))),
    if (all(is.finite(as.matrix(SC[, 2:5])))) "PASS" else "FAIL")

# ---------------- admission times ----------------
cand <- c(file.path(IN, "index_admission_times.csv"),
          file.path(ROOT, "Dissertation", "figwork", "data",
                    "index_admission_times.csv"))
cand <- cand[file.exists(cand)]
adm <- NULL
used <- NA_character_
for (f in cand) {
  h <- names(fread(f, nrows = 1))
  if (all(c("hadm_id", "admittime") %in% h)) {
    adm <- unique(fread(f, select = c("hadm_id", "admittime")),
                  by = "hadm_id")
    used <- sub(paste0(ROOT, "/"), "", f, fixed = TRUE)
    break
  }
}
cov_a <- if (is.null(adm)) 0 else mean(key %in% adm$hadm_id)
chk("admission dates found",
    if (is.null(adm)) {
      "no index_admission_times.csv with hadm_id and admittime"
    } else
      sprintf("%s, %.1f%% of admissions", used, 100 * cov_a),
    if (cov_a >= 0.99) "PASS" else "FAIL")
stop_if_fail()
adm[, adm_t := tparse(admittime)]

# ---------------- ECG at presentation ----------------
cat("\n=== ECG AT PRESENTATION ===\n")
eo <- unique(fread(file.path(IN, "ecg_opt.csv")), by = "hadm_id")
oc <- setdiff(names(eo), c("subject_id", "hadm_id"))
ix <- unique(fread(file.path(SP, "data", "raw", "ecg",
                             "ecg_record_index.csv")),
             by = c("rec", "hadm_id"))
ix <- ix[hadm_id %in% key]
# dates, for the comparison with D18 only
ix[, adm_d := as.Date(adm$adm_t[match(hadm_id, adm$hadm_id)])]
ix[, d_rel := as.numeric(as.Date(tparse(ecg_charttime)) - adm_d)]
# real acquisition times (stage 0r): hours from the real admission time
nrec <- function(r) {
  r <- gsub("\\\\", "/", as.character(r))
  r <- gsub("^/+|/+$", "", r)
  ifelse(startsWith(r, "files/"), r, paste0("files/", r))
}
rt <- unique(fread(file.path(PR, "audit", "ecg_times_real.csv")),
             by = c("rec", "hadm_id"))
rk <- paste(rt$hadm_id, nrec(rt$rec))
ix[, h := flt(rt$h[match(paste(hadm_id, nrec(rec)), rk)])]
tp <- mean(is.finite(ix$h))
chk("ECG records with a real acquisition time (stage 0r)",
    sprintf("%.1f%% of %d records", 100 * tp, nrow(ix)),
    if (tp >= 0.99) "PASS" else "FAIL")
stop_if_fail()
nper <- ix[, .N, by = hadm_id]
one <- nper[N == 1, hadm_id]
r1 <- ix[hadm_id %in% one]
# every per-ECG source, read once: the model outputs and each version of
# the derived measurements (newest first)
lg <- fread(file.path(SP, "data", "processed", "exp0_logits_record.csv"))
rf <- Sys.glob(file.path(SP, "data", "processed",
                         "ecg_derived_*record*.csv"))
ver <- suppressWarnings(as.numeric(sub(".*_v(\\d+)_record.*", "\\1",
                                       basename(rf))))
ver[is.na(ver)] <- 0
rf <- rf[order(-ver)]
SRC <- c(list(logits = lg), setNames(lapply(rf, fread), basename(rf)))
# each ECG column comes from the source that reproduces it on the
# single-ECG admissions (ties: model outputs, then the newest version)
same <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (!sum(ok)) return(NA_real_)
  mean(abs(a[ok] - b[ok]) <= 1e-6 * pmax(1, abs(b[ok])))
}
src <- list()
for (c0 in oc) {
  if (grepl("_worst$", c0)) next
  b <- sub("_(mean|max|min|sd|first|last|circ)$", "", c0)
  tgt <- flt(eo[[c0]][match(r1$hadm_id, eo$hadm_id)])
  best <- NULL
  for (nm in names(SRC)) {
    D <- SRC[[nm]]
    cc <- if (c0 %in% names(D)) c0 else if (b %in% names(D)) b else NA
    if (is.na(cc)) next
    sc <- same(flt(D[[cc]][match(r1$rec, D$rec)]), tgt)
    if (is.null(best) || (is.finite(sc) && sc > best$score)) {
      best <- list(file = nm, col = cc, score = ifelse(is.finite(sc), sc,
                                                       -1))
    }
  }
  if (!is.null(best)) src[[c0]] <- best
}
nwo <- sum(grepl("_worst$", oc))
um <- setdiff(oc[!grepl("_worst$", oc)], names(src))
chk("ECG columns found in a per-ECG file (_worst duplicates left out)",
    sprintf("%d of %d, %d left out%s", length(src), length(oc) - nwo, nwo,
            if (length(um)) paste0(" (missing: ",
                                   paste(head(um, 5), collapse = ", "),
                                   ")") else ""),
    if (length(um) == 0) "PASS" else "WARN")
fwrite(data.table(column = names(src),
                  source = sapply(src, function(z) z$file),
                  source_column = sapply(src, function(z) z$col),
                  share_reproduced = round(sapply(src, function(z)
                    z$score), 4)),
       file.path(RS, "presentation_ecg_sources.csv"))
cols <- names(src)
REC <- data.table(rec = lg$rec)
for (c0 in cols) {
  D <- SRC[[src[[c0]]$file]]
  set(REC, j = c0, value = flt(D[[src[[c0]]$col]][match(REC$rec, D$rec)]))
}
A <- as.matrix(REC[match(r1$rec, REC$rec), ..cols])
B <- as.matrix(eo[match(r1$hadm_id, hadm_id), ..cols])
storage.mode(B) <- "double"
okv <- !is.na(A) & !is.na(B)
dif <- abs(A - B) > 1e-6 * pmax(1, abs(B))
shd <- if (sum(okv)) mean(dif[okv]) else NA
chk("single-ECG admissions: per-ECG features reproduce ecg_opt",
    sprintf("%d admissions, %.2f%% of values differ", length(one),
            100 * shd), if (is.finite(shd) && shd <= 0.01) "PASS" else
              "FAIL")
stop_if_fail()
# D19: the first ECG by the clock from 24 h before to 24 h after admission
el <- ix[is.finite(h) & h >= WIN_ECG[1] & h <= WIN_ECG[2]]
el[, first_h := min(h), by = hadm_id]
w <- el[h == first_h]
EM <- as.data.table(REC[match(w$rec, REC$rec), ..cols])
EM[, hadm_id := w$hadm_id]
circ <- cols[grepl("_circ$", cols)]
lin <- setdiff(cols, circ)
cmean <- function(a) {
  a <- a[is.finite(a)]
  if (!length(a)) return(NA_real_)
  atan2(mean(sin(a * pi / 180)), mean(cos(a * pi / 180))) * 180 / pi
}
AG <- EM[, lapply(.SD, function(v) mean(v, na.rm = TRUE)), by = hadm_id,
         .SDcols = lin]
if (length(circ)) {
  AC <- EM[, lapply(.SD, cmean), by = hadm_id, .SDcols = circ]
  AG <- AG[AC, on = "hadm_id"]
}
for (cc in cols) set(AG, which(is.nan(AG[[cc]])), cc, NA_real_)
ECGP <- data.table(hadm_id = key)
ECGP <- cbind(ECGP, AG[match(key, AG$hadm_id), ..cols])
ECGP[, no_pres_ecg := as.numeric(!(hadm_id %in% AG$hadm_id))]
fwrite(ECGP, file.path(OD, "ecgP.csv"))
d18 <- unique(ix[d_rel >= -1 & d_rel <= 0, hadm_id])
nd <- w[, .N, by = hadm_id]
fh <- w[!duplicated(hadm_id), first_h]
chk("admissions with an ECG within 24 h of admission (D18 by dates)",
    sprintf("%s (%s)", pc(mean(key %in% AG$hadm_id)),
            pc(mean(key %in% d18))), "INFO")
smy("ECG", "admissions with an ECG within 24 h (D18 rule by dates)",
    sprintf("%s (%s)", pc(mean(key %in% AG$hadm_id)),
            pc(mean(key %in% d18))))
smy("ECG", "first ECG, hours from admission: median (IQR)",
    sprintf("%.1f (%.1f to %.1f)", median(fh), quantile(fh, 0.25),
            quantile(fh, 0.75)))
smy("ECG", "first ECG before admission / after admission",
    sprintf("%s / %s", pc(mean(fh < 0)), pc(mean(fh >= 0))))
smy("ECG", "admissions averaging ECGs with the same time",
    pc(mean(nd$N > 1)))

# ---------------- CTPA at presentation ----------------
cat("\n=== CTPA AT PRESENTATION ===\n")
ct <- unique(fread(file.path(IN, "ctpa46.csv")), by = "hadm_id")
cn <- file.path(ROOT, "Dissertation", "figwork", "data",
                "ctpa_notes_index.csv")
hn <- names(fread(cn, nrows = 1, sep = ","))
idk <- if ("idx_hadm" %in% hn) "idx_hadm" else "hadm_id"
own <- "admittime" %in% hn
nt <- fread(cn, select = c(idk, "charttime", if (own) "admittime"))
setnames(nt, idk, "ih")
if (own) {
  nt[, adm_t := tparse(admittime)]
} else {
  nt[, adm_t := adm$adm_t[match(ih, adm$hadm_id)]]
}
chk("CTPA report hours measured from",
    if (own) "ctpa_notes_index.csv's own admission times" else
      "index_admission_times.csv", "INFO")
nt[, h_rel := as.numeric(difftime(tparse(charttime), adm_t,
                                  units = "hours"))]
inw <- unique(nt[h_rel >= WIN_CTPA[1] & h_rel <= WIN_CTPA[2], ih])
cp <- ct[hadm_id %in% inw & hadm_id %in% key]
sh <- mean(ct$hadm_id[ct$hadm_id %in% key] %in% inw)
chk("ctpa46 admissions with a report in the presentation window",
    sprintf("%.1f%% (%d of %d)", 100 * sh, nrow(cp),
            sum(ct$hadm_id %in% key)), if (sh >= 0.90) "PASS" else "WARN")
fwrite(cp, file.path(OD, "ctpaP.csv"))
fp[, in_three_mod_p := as.integer(hadm_id %in% cp$hadm_id)]
fo <- fread(file.path(PR, "folds.csv"))[frame == "A" & cohort == "3-mod"]
nf <- sum(!(cp$hadm_id %in% fo$hadm_id))
chk("presentation-CTPA admissions without a three-modality fold", nf,
    if (nf == 0) "PASS" else "FAIL")
stop_if_fail()
smy("CTPA", "admissions with a presentation CTPA (deaths by 365 days)",
    sprintf("%d (%d)", sum(fp$in_three_mod_p),
            sum(fp$event_death[fp$in_three_mod_p == 1])))
fp[, has_pres_vitals := as.integer(has_v)]
fp[, has_pres_ecg := as.integer(hadm_id %in% AG$hadm_id)]
fwrite(fp, file.path(OD, "frame_P.csv"))

# ---------------- frame_P_outcomes (stage 4r) ----------------
fof <- file.path(OD, "frame_P_outcomes.csv")
if (file.exists(fof)) {
  fo2 <- fread(fof)
  m <- match(fp$hadm_id, fo2$hadm_id)
  okr <- !anyNA(m) && nrow(fo2) == nrow(fp)
  chk("frame_P_outcomes.csv: rows match frame P", sprintf(
    "%d of %d", sum(!is.na(m)), nrow(fp)), if (okr) "PASS" else "FAIL")
  stop_if_fail()
  fo2 <- fo2[m]
  okd <- identical(as.integer(fo2$event_death), fp$event_death) &&
    isTRUE(all.equal(as.numeric(fo2$time_death), fp$time_death))
  chk("frame_P_outcomes.csv: deaths and times unchanged",
      if (okd) "yes" else "no", if (okd) "PASS" else "FAIL")
  stop_if_fail()
  for (cc in c("in_three_mod_p", "has_pres_vitals", "has_pres_ecg")) {
    set(fo2, j = cc, value = fp[[cc]])
  }
  fwrite(fo2, fof)
  chk("frame_P_outcomes.csv: presentation flags updated",
      "has_pres_ecg, has_pres_vitals, in_three_mod_p", "PASS")
} else {
  chk("frame_P_outcomes.csv", "not found: run stage 4r next", "WARN")
}

ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "presentation_checks.csv"))
fwrite(rbindlist(SUMM), file.path(RS, "presentation_summary.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL")))
cat("\nOUTPUT FILES\n")
for (f in c("frame_P.csv", "frame_P_outcomes.csv", "ehrP.csv",
            "ehrP_sens.csv", "ecgP.csv", "ctpaP.csv", "scores_P.csv")) {
  cat("  ", file.path(OD, f), " (local only)\n")
}
for (f in c("presentation_checks.csv", "presentation_summary.csv",
            "presentation_history.csv", "presentation_ecg_sources.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "presentation_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: presentation_checks.csv, presentation_summary.csv,",
    "presentation_history.csv\n")