# 07b_primary_hospice.R
# Sensitivity analysis for the primary contrast, excluding hospice
# discharges (PREREGISTRATION.md, section 8).
#
# Refits the 07_primary.R model (Frame B, EHR+ECG cohort, death within 365
# days of discharge) after dropping patients discharged to hospice. The
# seed-averaged ridge-Cox scores are standardised with the full-cohort mean
# and SD, so the units match the primary. Everything else is as in
# 07_primary.R: the split at 30.5 days and the 2,000-refit patient-level
# bootstrap (seed 42). The ECG ratio is printed beside the primary result
# from results/primary.csv.
#
# Checks: 2,909 patients before exclusion, no missing score, and 177 deaths
# in days 0-30 and 443 in days 31-365 in the full cohort. The script stops
# before any fit if one fails.
#
# Inputs:  data/processed/level1/oof_linear_B_ehrecg.csv
#          data/processed/frame_B.csv, results/primary.csv
# Outputs: results/primary_hospice.csv, results/primary_hospice_boot.csv,
#          results/primary_hospice_log.txt
# Usage:   Rscript scripts/07b_primary_hospice.R [--probe-only]
# Notes:   A timing probe of 20 refits projects the bootstrap run time and
#          pauses 5 s; --probe-only stops after it. The bootstrap runs in
#          parallel (future.apply, physical cores minus 1).

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
NBOOT <- 2000L
SEED <- 42L
CUT <- 30.5
t_start <- Sys.time()

cat("=== LOAD ===\n")
oof <- fread(file.path(P2, "data", "processed", "level1",
                       "oof_linear_B_ehrecg.csv"),
             select = c("seed", "hadm_id", "ehr.cox.death.lp",
                        "ecg.cox.death.lp"))
sc <- oof[, .(ehr = mean(`ehr.cox.death.lp`),
              ecg = mean(`ecg.cox.death.lp`)), by = hadm_id]
fb <- fread(file.path(P2, "data", "processed", "frame_B.csv"))
d <- merge(fb[, .(hadm_id, hospice, time_death, event_death)],
           sc, by = "hadm_id")
full_early <- sum(d$event_death == 1 & d$time_death <= CUT)
full_late <- sum(d$event_death == 1 & d$time_death > CUT)
d[, ehr_z := (ehr - mean(ehr)) / sd(ehr)]
d[, ecg_z := (ecg - mean(ecg)) / sd(ecg)]
n_hosp <- sum(d$hospice == 1)
hd <- d[hospice == 1]
cat(sprintf("  hospice discharges %d: deaths 0-30 d %d, 31-365 d %d\n",
            n_hosp,
            sum(hd$event_death == 1 & hd$time_death <= CUT),
            sum(hd$event_death == 1 & hd$time_death > CUT)))
cat(sprintf("  mean ECG score (SD units): hospice %.2f, others %.2f\n",
            mean(hd$ecg_z), mean(d[hospice == 0]$ecg_z)))
dd <- as.data.frame(d[hospice == 0,
                      .(time_death, event_death, ehr_z, ecg_z)])
n_early <- sum(dd$event_death == 1 & dd$time_death <= CUT)
n_late <- sum(dd$event_death == 1 & dd$time_death > CUT)
ok <- nrow(d) == 2909 && nrow(dd) == 2909 - n_hosp &&
  full_early == 177 && full_late == 443 &&
  !anyNA(dd$ehr_z) && !anyNA(dd$ecg_z)
cat(sprintf("  analysed %d patients: deaths 0-30 d %d, 31-365 d %d\n",
            nrow(dd), n_early, n_late))
if (!ok) {
  cat("HARD CHECK FAILED: stopping before any fit.\n")
  quit(status = 1)
}
cat("  HARD CHECKS PASS\n")

fit_split <- function(x) {
  s <- survSplit(Surv(time_death, event_death) ~ ehr_z + ecg_z,
                 data = x, cut = CUT, episode = "period")
  s$ehr_e <- s$ehr_z * (s$period == 1)
  s$ehr_l <- s$ehr_z * (s$period == 2)
  s$ecg_e <- s$ecg_z * (s$period == 1)
  s$ecg_l <- s$ecg_z * (s$period == 2)
  coxph(Surv(tstart, time_death, event_death) ~
          ehr_e + ehr_l + ecg_e + ecg_l,
        data = s, ties = "efron")
}
boot_one <- function(i, x) {
  idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
  f <- tryCatch(fit_split(x[idx, ]), error = function(e) NULL)
  if (is.null(f)) return(rep(NA_real_, 4))
  coef(f)
}

cat("\n=== TIMING PROBE (20 real refits) ===\n")
set.seed(1)
tp <- system.time(for (i in 1:20) boot_one(i, dd))[["elapsed"]] / 20
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat(sprintf("  %.3f s per refit; parallel (%d workers) about %.1f min\n",
            tp, workers, (tp * NBOOT / workers * 1.3 + 10) / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

f0 <- fit_split(dd)
b <- coef(f0)
V <- vcov(f0)
plan(multisession, workers = workers)
B <- do.call(rbind, future_lapply(
  seq_len(NBOOT), boot_one, x = dd, future.seed = SEED,
  future.packages = "survival"))
plan(sequential)
colnames(B) <- names(b)
n_fail <- sum(!complete.cases(B))
Bc <- B[complete.cases(B), , drop = FALSE]
q <- function(v) unname(quantile(v, c(0.025, 0.975)))
r_ecg <- exp(Bc[, "ecg_e"] - Bc[, "ecg_l"])
r_ehr <- exp(Bc[, "ehr_e"] - Bc[, "ehr_l"])
ci_ecg <- q(r_ecg)
ci_ehr <- q(r_ehr)
est_ecg <- unname(exp(b["ecg_e"] - b["ecg_l"]))
est_ehr <- unname(exp(b["ehr_e"] - b["ehr_l"]))
verdict <- if (ci_ecg[1] > 1) {
  "consistent with H1: lower bound above 1"
} else if (ci_ecg[2] < 1) {
  "reversed: upper bound below 1"
} else {
  "interval includes 1"
}

rows <- list()
add <- function(quantity, est, lo, hi, method, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    quantity = quantity, estimate = est, lo = lo, hi = hi,
    method = method, note = note)
}
add("hospice discharges excluded", n_hosp, NA, NA, "count")
add("patients analysed", nrow(dd), NA, NA, "count")
add("deaths days 0-30", n_early, NA, NA, "count")
add("deaths days 31-365", n_late, NA, NA, "count")
for (v in c("ehr_e", "ehr_l", "ecg_e", "ecg_l")) {
  bq <- q(exp(Bc[, v]))
  add(sprintf("HR per SD, %s, days %s", toupper(substr(v, 1, 3)),
              ifelse(grepl("_e$", v), "0-30", "31-365")),
      unname(exp(b[v])), bq[1], bq[2], "bootstrap percentile")
}
add("SENSITIVITY: ECG HR ratio, early / late, no hospice",
    est_ecg, ci_ecg[1], ci_ecg[2],
    "bootstrap percentile, 2,000 refits", verdict)
add("EHR HR ratio, early / late, no hospice", est_ehr,
    ci_ehr[1], ci_ehr[2], "bootstrap percentile, 2,000 refits")
add("bootstrap refits that failed", n_fail, NA, NA, "count")
out <- do.call(rbind, rows)
fwrite(out, file.path(P2, "results", "primary_hospice.csv"))
fwrite(data.table(replicate = seq_len(nrow(B)), B,
                  ecg_ratio = exp(B[, "ecg_e"] - B[, "ecg_l"]),
                  ehr_ratio = exp(B[, "ehr_e"] - B[, "ehr_l"])),
       file.path(P2, "results", "primary_hospice_boot.csv"))

cat("\n=== RESULT ===\n")
print(out, row.names = FALSE, digits = 4)
pr <- fread(file.path(P2, "results", "primary.csv"))
p1 <- pr[grepl("^PRIMARY", quantity)]
cat(sprintf("\n  primary (all 2,909):       %.3f [%.3f, %.3f]\n",
            p1$estimate, p1$lo, p1$hi))
cat(sprintf("  hospice excluded (%d):   %.3f [%.3f, %.3f]  %s\n",
            nrow(dd), est_ecg, ci_ecg[1], ci_ecg[2], verdict))
cat(sprintf("  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start,
                                units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(P2, "results", "primary_hospice.csv"), "\n")
cat("  ", file.path(P2, "results", "primary_hospice_boot.csv"), "\n")
cat("  ", file.path(P2, "results", "primary_hospice_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: primary_hospice.csv\n")