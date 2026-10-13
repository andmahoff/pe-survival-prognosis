# 07c_s5_three_modality.R
# S5: the primary contrast repeated in the three-modality cohort
# (PREREGISTRATION.md, section 6). Secondary and exploratory, so intervals
# only and no decision rule.
#
# Frame B, three-modality cohort: patients discharged alive whose first
# admission had a CTPA, with death within 365 days of discharge. Two Cox
# models are fitted on the standardised, seed-averaged, out-of-fold
# ridge-Cox scores, with every coefficient allowed to differ between days
# 0-30 and 31-365 (split at 30.5): (a) EHR and ECG, the primary model in
# this smaller cohort, and (b) EHR, ECG and CTPA. Comparing (a) with the
# primary isolates the cohort, and comparing (b) with (a) isolates the
# adjustment for CTPA. Both models are refitted on the same 2,000
# patient-level bootstrap samples (seed 42), so their intervals are
# directly comparable.
#
# Checks: 1,405 patients, six seeds per patient, no missing score, and 79
# deaths in days 0-30 and 213 in days 31-365 (stage 4). The script stops
# before any fit if one fails.
#
# Inputs:  data/processed/level1/oof_linear_B_3mod.csv
#          data/processed/frame_B.csv, results/primary.csv
# Outputs: results/s5_primary.csv, results/s5_primary_boot.csv,
#          results/s5_primary_log.txt
# Usage:   Rscript scripts/07c_s5_three_modality.R [--probe-only]
# Notes:   A timing probe of 20 refits of both models projects the
#          bootstrap run time and pauses 5 s; --probe-only stops after it.
#          The bootstrap runs in parallel (future.apply, physical cores
#          minus 1).

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
cols <- c("ehr.cox.death.lp", "ecg.cox.death.lp", "ctpa.cox.death.lp")
oof <- fread(file.path(P2, "data", "processed", "level1",
                       "oof_linear_B_3mod.csv"),
             select = c("seed", "hadm_id", cols))
sc <- oof[, .(ehr = mean(`ehr.cox.death.lp`),
              ecg = mean(`ecg.cox.death.lp`),
              ctpa = mean(`ctpa.cox.death.lp`),
              n_seeds = .N), by = hadm_id]
fb <- fread(file.path(P2, "data", "processed", "frame_B.csv"))
d <- merge(fb[in_three_mod == 1,
              .(hadm_id, time_death, event_death)],
           sc, by = "hadm_id")
for (v in c("ehr", "ecg", "ctpa")) {
  d[[paste0(v, "_z")]] <- (d[[v]] - mean(d[[v]])) / sd(d[[v]])
}
dd <- as.data.frame(d[, .(time_death, event_death,
                          ehr_z, ecg_z, ctpa_z)])
n_early <- sum(dd$event_death == 1 & dd$time_death <= CUT)
n_late <- sum(dd$event_death == 1 & dd$time_death > CUT)
cat(sprintf("  patients %d, deaths 0-30 d %d, 31-365 d %d\n",
            nrow(dd), n_early, n_late))
cat("  score correlations:\n")
print(round(cor(dd[, c("ehr_z", "ecg_z", "ctpa_z")]), 3))
ok <- nrow(dd) == 1405 && all(sc$n_seeds == 6) &&
  n_early == 79 && n_late == 213 &&
  !anyNA(dd[, c("ehr_z", "ecg_z", "ctpa_z")])
if (!ok) {
  cat("HARD CHECK FAILED: stopping before any fit.\n")
  quit(status = 1)
}
cat("  HARD CHECKS PASS\n")

VA <- c("ehr_z", "ecg_z")
VB <- c("ehr_z", "ecg_z", "ctpa_z")
fit_split <- function(x, vars) {
  f0 <- as.formula(paste("Surv(time_death, event_death) ~",
                         paste(vars, collapse = " + ")))
  s <- survSplit(f0, data = x, cut = CUT, episode = "period")
  tt <- character(0)
  for (v in vars) {
    s[[paste0(v, "_e")]] <- s[[v]] * (s$period == 1)
    s[[paste0(v, "_l")]] <- s[[v]] * (s$period == 2)
    tt <- c(tt, paste0(v, "_e"), paste0(v, "_l"))
  }
  coxph(as.formula(paste("Surv(tstart, time_death, event_death) ~",
                         paste(tt, collapse = " + "))),
        data = s, ties = "efron")
}
boot_one <- function(i, x) {
  idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
  xa <- x[idx, ]
  a <- tryCatch(coef(fit_split(xa, VA)),
                error = function(e) rep(NA_real_, 4))
  b <- tryCatch(coef(fit_split(xa, VB)),
                error = function(e) rep(NA_real_, 6))
  c(setNames(a, paste0("a.", c("ehr_z_e", "ehr_z_l",
                               "ecg_z_e", "ecg_z_l"))),
    setNames(b, paste0("b.", c("ehr_z_e", "ehr_z_l", "ecg_z_e",
                               "ecg_z_l", "ctpa_z_e",
                               "ctpa_z_l"))))
}

cat("\n=== TIMING PROBE (20 real refits of both models) ===\n")
set.seed(1)
tp <- system.time(for (i in 1:20) boot_one(i, dd))[["elapsed"]] / 20
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat(sprintf("  %.3f s per replicate; parallel (%d workers) about",
            tp, workers),
    sprintf("%.1f min\n", (tp * NBOOT / workers * 1.3 + 10) / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

fa <- fit_split(dd, VA)
fbm <- fit_split(dd, VB)
plan(multisession, workers = workers)
B <- do.call(rbind, future_lapply(
  seq_len(NBOOT), boot_one, x = dd, future.seed = SEED,
  future.packages = "survival"))
plan(sequential)
n_fail <- sum(!complete.cases(B))
Bc <- B[complete.cases(B), , drop = FALSE]
q <- function(v) unname(quantile(v, c(0.025, 0.975)))

rows <- list()
add <- function(model, quantity, est, lo, hi, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    model = model, quantity = quantity, estimate = est,
    lo = lo, hi = hi, note = note)
}
add("-", "patients", nrow(dd), NA, NA)
add("-", "deaths days 0-30", n_early, NA, NA)
add("-", "deaths days 31-365", n_late, NA, NA)
for (m in list(list("a", "(a) EHR+ECG", fa, VA),
               list("b", "(b) EHR+ECG+CTPA", fbm, VB))) {
  pre <- m[[1]]
  cf <- coef(m[[3]])
  for (v in m[[4]]) {
    lab <- toupper(sub("_z", "", v))
    for (per in c("e", "l")) {
      nm <- paste0(v, "_", per)
      bq <- q(exp(Bc[, paste0(pre, ".", nm)]))
      add(m[[2]], sprintf("HR per SD, %s, days %s", lab,
                          ifelse(per == "e", "0-30", "31-365")),
          unname(exp(cf[nm])), bq[1], bq[2])
    }
    r_b <- exp(Bc[, paste0(pre, ".", v, "_e")] -
                 Bc[, paste0(pre, ".", v, "_l")])
    rq <- q(r_b)
    note <- if (rq[1] > 1) "lower bound above 1" else
      if (rq[2] < 1) "upper bound below 1" else "interval includes 1"
    add(m[[2]], sprintf("%s HR ratio, early / late", lab),
        unname(exp(cf[paste0(v, "_e")] - cf[paste0(v, "_l")])),
        rq[1], rq[2], note)
  }
}
add("-", "bootstrap replicates that failed", n_fail, NA, NA)
out <- do.call(rbind, rows)
fwrite(out, file.path(P2, "results", "s5_primary.csv"))
fwrite(data.table(replicate = seq_len(nrow(B)), B),
       file.path(P2, "results", "s5_primary_boot.csv"))

cat("\n=== RESULT ===\n")
print(out, row.names = FALSE, digits = 4)
pr <- fread(file.path(P2, "results", "primary.csv"))
p1 <- pr[grepl("^PRIMARY", quantity)]
ea <- out[out$model == "(a) EHR+ECG" &
            out$quantity == "ECG HR ratio, early / late", ]
eb <- out[out$model == "(b) EHR+ECG+CTPA" &
            out$quantity == "ECG HR ratio, early / late", ]
cat(sprintf("\n  ECG ratio, primary (EHR+ECG, n=2,909):   %.3f [%.3f, %.3f]\n",
            p1$estimate, p1$lo, p1$hi))
cat(sprintf("  ECG ratio, (a) EHR+ECG, 3-mod cohort:     %.3f [%.3f, %.3f]\n",
            ea$estimate, ea$lo, ea$hi))
cat(sprintf("  ECG ratio, (b) EHR+ECG+CTPA, 3-mod:       %.3f [%.3f, %.3f]\n",
            eb$estimate, eb$lo, eb$hi))
cat(sprintf("  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start,
                                units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(P2, "results", "s5_primary.csv"), "\n")
cat("  ", file.path(P2, "results", "s5_primary_boot.csv"), "\n")
cat("  ", file.path(P2, "results", "s5_primary_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: s5_primary.csv\n")