# 04r_optionB_outcomes.R
# Composite and CV-readmission outcomes for the at-diagnosis frame, and the
# discharge pipeline's three feature matrices (PREREGISTRATION_ADDENDUM_2.md
# A33, Option B).
#
# At diagnosis (Frame P, clock at admission, 365 days), the outcome is the
# first of a CV readmission (cv_days_all in followup_365.csv, days from
# admission) or death (dod_days, else inhosp_days, as in stage 4p).
# time_first = days + 0.5, censored at 366. event_cr is 1 for a CV
# readmission first (a same-day CV readmission and death counts as CV, as in
# stage 4), 2 for death first and 0 for neither, and event_comp is 1 for
# either. time_first and event_cr must reproduce Frame A's to day 30
# exactly, or the script stops.
#
# At discharge (Frame B), time_first, event_cr and event_comp are taken
# unchanged from frame_B.csv (stage 4). The EHR+ matrix is rebuilt as in
# stage 6j and must match its 139 features in name and order, or the script
# stops. The ECG and CTPA matrices use the feature whitelist (ecg_opt.csv,
# ctpa46.csv), as in stages 6a to 6i. The events file counts events per
# frame, cohort, outcome and horizon (10, 30, 90, 180 and 365 days).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, feature_whitelist.csv,
#          presentation/frame_P.csv; data/inputs/followup_365.csv,
#          mimic_labels_harmonised.csv, ehr28.csv, ehr_extra.csv,
#          dx_index.csv, ecg_opt.csv, ctpa46.csv;
#          results/ehrplus_features.csv
# Outputs: results/optionB_outcome_checks.csv, optionB_outcome_events.csv,
#          optionB_outcome_log.txt;
#          data/processed/presentation/frame_P_outcomes.csv,
#          data/processed/discharge/ehrplusB.csv, ecgB.csv, ctpaB.csv
#          (patient level, not in the repo)
# Usage:   Rscript scripts/04r_optionB_outcomes.R

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("comorbidity", quietly = TRUE)) {
  stop("package comorbidity not installed")
}
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
DD <- file.path(PR, "discharge")
RS <- file.path(P2, "results")
dir.create(DD, showWarnings = FALSE, recursive = TRUE)
TAU <- 365
HZ <- c(10, 30, 90, 180, 365)
EXCL <- c("vancomycin", "_50825_")
EXTRA <- c("spo2", "nlr", "hr_first", "sbp_first", "spo2_first", "sex_f",
           "n_prior", "n_codes")
DEST <- c("HOME HEALTH CARE" = "dest_home_health",
          "SKILLED NURSING FACILITY" = "dest_skilled_nursing",
          "REHAB" = "dest_rehab",
          "CHRONIC/LONG TERM ACUTE CARE" = "dest_long_term_care",
          "HOSPICE" = "dest_hospice")
PREV <- c(0.01, 0.99)
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "optionB_outcome_checks.csv"))
    cat("\nSTOP: a check failed. see: optionB_outcome_checks.csv\n")
    quit(status = 1)
  }
}
flt <- function(x) suppressWarnings(as.numeric(x))
need <- c(file.path(PR, c("frame_A.csv", "frame_B.csv",
                          "feature_whitelist.csv")),
          file.path(PP, "frame_P.csv"),
          file.path(IN, c("followup_365.csv", "mimic_labels_harmonised.csv",
                          "ehr28.csv", "ehr_extra.csv", "dx_index.csv",
                          "ecg_opt.csv", "ctpa46.csv")),
          file.path(RS, "ehrplus_features.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()

# ---------------- at-diagnosis outcomes ----------------
cat("=== AT-DIAGNOSIS OUTCOMES (Frame P) ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fp <- fread(file.path(PP, "frame_P.csv"))
fu <- unique(fread(file.path(IN, "followup_365.csv")), by = "hadm_id")
lab <- unique(fread(file.path(IN, "mimic_labels_harmonised.csv")),
              by = "hadm_id")
L <- lab[match(fp$hadm_id, hadm_id)]
td <- ifelse(is.na(flt(L$dod_days)), flt(L$inhosp_days), flt(L$dod_days))
tc <- flt(fu$cv_days_all[match(fp$hadm_id, fu$hadm_id)])
chk("CV readmission times found (cv_days_all)",
    sprintf("%d admissions with a later CV admission", sum(!is.na(tc))),
    if ("cv_days_all" %in% names(fu)) "PASS" else "FAIL")
stop_if_fail()
cv1 <- !is.na(tc) & tc >= 0 & tc <= TAU & (is.na(td) | tc <= td)
de1 <- !is.na(td) & td >= 0 & td <= TAU & !cv1
first <- ifelse(cv1, tc, ifelse(de1, td, NA_real_))
fp[, time_first := ifelse(cv1 | de1, first + 0.5, TAU + 1)]
fp[, event_cr := ifelse(cv1, 1L, ifelse(de1, 2L, 0L))]
fp[, event_comp := as.integer(event_cr > 0)]
A <- fa[match(fp$hadm_id, hadm_id)]
a_ev <- A$event_cr
p_ev <- ifelse(fp$time_first <= 30.5, fp$event_cr, 0L)
ok <- identical(as.integer(a_ev), as.integer(p_ev)) &&
  isTRUE(all.equal(fp$time_first[a_ev > 0], A$time_first[a_ev > 0]))
chk("Frame A's first events and times to day 30 reproduced",
    sprintf("CV %d vs %d, death %d vs %d", sum(p_ev == 1), sum(a_ev == 1),
            sum(p_ev == 2), sum(a_ev == 2)), if (ok) "PASS" else "FAIL")
stop_if_fail()
fwrite(fp, file.path(PP, "frame_P_outcomes.csv"))

# ---------------- discharge matrices ----------------
cat("\n=== DISCHARGE MATRICES (Frame B) ===\n")
fb <- fread(file.path(PR, "frame_B.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
ex <- unique(fread(file.path(IN, "ehr_extra.csv")), by = "hadm_id")
labs <- grep("^lab_", names(ex), value = TRUE)
drop <- labs[grepl(paste(EXCL, collapse = "|"), labs)]
labs <- setdiff(labs, drop)
c28 <- intersect(wl[modality == "ehr", column], names(e28))
xtra <- intersect(EXTRA, names(ex))
chk("EHR 79: registered / extra labs (excluded)",
    sprintf("%d / %d (%s)", length(c28), length(labs),
            paste(drop, collapse = ", ")),
    if (length(labs) >= 40) "PASS" else "FAIL")
chk("block 2, unused measurements found",
    sprintf("%d of %d (%s)", length(xtra), length(EXTRA),
            paste(xtra, collapse = ", ")),
    if (length(xtra) >= 6) "PASS" else "FAIL")
stop_if_fail()
# block 3: discharge destination
fu <- unique(fread(file.path(IN, "followup_365.csv")), by = "hadm_id")
loc <- toupper(trimws(as.character(
  fu$discharge_location[match(fb$hadm_id, fu$hadm_id)])))
loc[is.na(loc) | loc == "NA"] <- ""
DS <- sapply(names(DEST), function(k) as.numeric(loc == k))
colnames(DS) <- DEST
DS <- cbind(DS, dest_other = as.numeric(!(loc %in% c("HOME", names(DEST)))))
chk("block 3, discharge destination (home / other share)",
    sprintf("%.1f%% / %.1f%%", 100 * mean(loc == "HOME"),
            100 * mean(DS[, "dest_other"])), "INFO")
# block 4: comorbidity from the index admission's codes
dx <- fread(file.path(IN, "dx_index.csv"), colClasses = list(
  character = "icd_code"))
dx <- dx[hadm_id %in% fb$hadm_id]
dx[, icd_code := gsub("[^A-Za-z0-9]", "", toupper(icd_code))]
cov_dx <- mean(fb$hadm_id %in% dx$hadm_id)
chk("block 4, admissions with diagnosis codes", sprintf("%.1f%%",
                                                       100 * cov_dx),
    if (cov_dx >= 0.99) "PASS" else "FAIL")
stop_if_fail()
cmb <- function(map9, map10) {
  out <- NULL
  for (v in c(9, 10)) {
    d <- as.data.frame(dx[icd_version == v, .(hadm_id, icd_code)])
    if (!nrow(d)) next
    r <- comorbidity::comorbidity(d, id = "hadm_id", code = "icd_code",
                                  map = if (v == 9) map9 else map10,
                                  assign0 = FALSE)
    out <- rbind(out, as.data.frame(r))
  }
  out <- as.data.table(out)[, lapply(.SD, max), by = hadm_id]
  out <- out[match(fb$hadm_id, hadm_id)]
  out[is.na(hadm_id), hadm_id := fb$hadm_id[is.na(out$hadm_id)]]
  for (k in setdiff(names(out), "hadm_id")) {
    set(out, which(is.na(out[[k]])), k, 0)
  }
  out
}
CH <- cmb("charlson_icd9_quan", "charlson_icd10_quan")
EL <- cmb("elixhauser_icd9_quan", "elixhauser_icd10_quan")
sc_ch <- as.numeric(comorbidity::score(
  structure(as.data.frame(CH), map = "charlson_icd10_quan",
            class = c("comorbidity", "data.frame")),
  weights = "charlson", assign0 = TRUE))
sc_el <- as.numeric(comorbidity::score(
  structure(as.data.frame(EL), map = "elixhauser_icd10_quan",
            class = c("comorbidity", "data.frame")),
  weights = "vw", assign0 = TRUE))
ind <- cbind(setNames(as.data.frame(CH)[, -1],
                      paste0("ch_", names(CH)[-1])),
             setNames(as.data.frame(EL)[, -1],
                      paste0("el_", names(EL)[-1])))
pv <- colMeans(ind)
keep <- pv >= PREV[1] & pv <= PREV[2]
CM <- cbind(as.matrix(ind[, keep, drop = FALSE]), charlson_score = sc_ch,
            elixhauser_vw_score = sc_el)
chk("block 4, Charlson / Elixhauser categories kept (1%-99%)",
    sprintf("%d of %d / %d of %d", sum(keep[startsWith(names(keep), "ch_")]),
            ncol(CH) - 1, sum(keep[startsWith(names(keep), "el_")]),
            ncol(EL) - 1), "INFO")
chk("block 4, median Charlson / Elixhauser score",
    sprintf("%g / %g", median(sc_ch), median(sc_el)), "INFO")
FT <- rbind(
  data.table(block = "1 EHR 79", feature = c(c28, labs), prevalence = NA),
  data.table(block = "2 unused measurements", feature = xtra,
             prevalence = NA),
  data.table(block = "3 discharge destination", feature = colnames(DS),
             prevalence = colMeans(DS)),
  data.table(block = "4 comorbidity", feature = names(pv),
             prevalence = pv)[, kept := keep][],
  data.table(block = "4 comorbidity", feature = c("charlson_score",
                                                   "elixhauser_vw_score"),
             prevalence = NA), fill = TRUE)
FT[is.na(kept), kept := TRUE]
XALL <- cbind(as.matrix(e28[match(fb$hadm_id, hadm_id), ..c28]),
              as.matrix(ex[match(fb$hadm_id, hadm_id), ..labs]),
              as.matrix(ex[match(fb$hadm_id, hadm_id), ..xtra]),
              DS, CM)
storage.mode(XALL) <- "double"
ref <- fread(file.path(RS, "ehrplus_features.csv"))
ref <- ref[kept == TRUE, feature]
okf <- identical(colnames(XALL), ref)
chk("EHR+ rebuilt exactly as stage 6j (features, order)",
    sprintf("%d vs %d", ncol(XALL), length(ref)), if (okf) "PASS" else
      "FAIL")
stop_if_fail()
fwrite(data.table(hadm_id = fb$hadm_id, XALL),
       file.path(DD, "ehrplusB.csv"))
for (m in c("ecg", "ctpa")) {
  f <- if (m == "ecg") "ecg_opt.csv" else "ctpa46.csv"
  D <- unique(fread(file.path(IN, f)), by = "hadm_id")
  cc <- intersect(wl[modality == m, column], names(D))
  sub <- if (m == "ctpa") fb[in_three_mod == 1] else fb
  X <- D[match(sub$hadm_id, hadm_id), ..cc]
  fwrite(data.table(hadm_id = sub$hadm_id, X),
         file.path(DD, sprintf("%sB.csv", m)))
  chk(sprintf("%s for Frame B: features (admissions)", toupper(m)),
      sprintf("%d (%d)", length(cc), nrow(sub)),
      if (length(cc) > 0 && !anyNA(match(sub$hadm_id, D$hadm_id))) "PASS"
      else "FAIL")
}
stop_if_fail()

# ---------------- event counts ----------------
EV <- list()
cnt <- function(D, fr, co) {
  for (h in HZ) {
    EV[[length(EV) + 1]] <<- data.table(
      frame = fr, cohort = co, horizon = h, patients = nrow(D),
      death = sum(D$event_death == 1 & D$time_death <= h + 0.5),
      composite = sum(D$event_comp == 1 & D$time_first <= h + 0.5),
      cv_first = sum(D$event_cr == 1 & D$time_first <= h + 0.5),
      death_first = sum(D$event_cr == 2 & D$time_first <= h + 0.5))
  }
}
cnt(fp, "P (at diagnosis)", "EHR+ECG")
cnt(fp[in_three_mod_p == 1], "P (at diagnosis)", "3-mod")
cnt(fb, "B (at discharge)", "EHR+ECG")
cnt(fb[in_three_mod == 1], "B (at discharge)", "3-mod")
ev <- rbindlist(EV)
fwrite(ev, file.path(RS, "optionB_outcome_events.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "optionB_outcome_checks.csv"))
cat("\n=== EVENTS BY HORIZON ===\n")
print(ev, row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL")))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(PP, "frame_P_outcomes.csv"), " (local only)\n")
for (f in c("ehrplusB.csv", "ecgB.csv", "ctpaB.csv")) {
  cat("  ", file.path(DD, f), " (local only)\n")
}
cat("  ", file.path(RS, "optionB_outcome_checks.csv"), "\n")
cat("  ", file.path(RS, "optionB_outcome_events.csv"), "\n")
cat("  ", file.path(RS, "optionB_outcome_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: optionB_outcome_checks.csv, optionB_outcome_events.csv\n")