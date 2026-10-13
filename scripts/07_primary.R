# 07_primary.R
# Registered primary analysis (PREREGISTRATION.md, section 5).
#
# Fits a Cox model for death within 365 days of discharge (Frame B,
# EHR+ECG cohort) on the standardised out-of-fold ridge-Cox scores of the
# EHR and the ECG. Each score's hazard ratio may differ between days 0-30
# and 31-365. The estimand is the ratio of the ECG's early to late hazard
# ratio, with a 95% percentile interval from 2,000 patient-level
# bootstrap refits (seed 42). H1 is supported if the lower bound exceeds
# 1, reversed if the upper bound is below 1, and inconclusive otherwise.
#
# The scores are the seed-averaged linear predictors from stage 6a,
# standardised to mean 0 and SD 1 and then held fixed. Follow-up is split
# at 30.5 days. The EHR ratio is reported with no decision rule, and the
# period-specific hazard ratios and the cox.zph test of the unsplit model
# are reported descriptively.
#
# Checks: 2,909 patients, six seeds per patient, no missing score, and 177
# deaths in days 0-30 and 443 in days 31-365 (stage 4). The script stops
# before any fit if one fails.
#
# Inputs:  data/processed/level1/oof_linear_B_ehrecg.csv
#          data/processed/frame_B.csv
# Outputs: results/primary.csv, results/primary_boot.csv,
#          results/primary_log.txt
# Usage:   Rscript scripts/07_primary.R [--probe-only]
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
              ecg = mean(`ecg.cox.death.lp`),
              n_seeds = .N), by = hadm_id]
fb <- fread(file.path(P2, "data", "processed", "frame_B.csv"))
d <- merge(fb[, .(subject_id, hadm_id, time_death, event_death)],
           sc, by = "hadm_id")
d[, ehr_z := (ehr - mean(ehr)) / sd(ehr)]
d[, ecg_z := (ecg - mean(ecg)) / sd(ecg)]
dd <- as.data.frame(d[, .(time_death, event_death, ehr_z, ecg_z)])

n_early <- sum(dd$event_death == 1 & dd$time_death <= CUT)
n_late <- sum(dd$event_death == 1 & dd$time_death > CUT)
ok <- nrow(dd) == 2909 && all(sc$n_seeds == 6) &&
  !anyNA(dd$ehr_z) && !anyNA(dd$ecg_z) &&
  n_early == 177 && n_late == 443
cat(sprintf("  patients %d, seeds per patient %s, deaths 0-30 d %d,",
            nrow(dd), paste(unique(sc$n_seeds), collapse = "/"),
            n_early),
    sprintf("31-365 d %d\n", n_late))
cat(sprintf("  score correlation EHR vs ECG: %.3f\n",
            cor(dd$ehr_z, dd$ecg_z)))
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
serial <- tp * NBOOT
par_t <- serial / workers * 1.3 + 10
cat(sprintf("  %.3f s per refit; %d refits: serial %.1f min,",
            tp, NBOOT, serial / 60),
    sprintf("parallel (%d workers) about %.1f min\n",
            workers, par_t / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

cat("\n=== POINT ESTIMATES ===\n")
f0 <- fit_split(dd)
b <- coef(f0)
V <- vcov(f0)
wald <- function(est, se) c(exp(est), exp(est - 1.96 * se),
                            exp(est + 1.96 * se))
lr_se <- function(a, c) sqrt(V[a, a] + V[c, c] - 2 * V[a, c])

cat(sprintf("=== BOOTSTRAP: %d refits on %d workers ===\n",
            NBOOT, workers))
plan(multisession, workers = workers)
B <- do.call(rbind, future_lapply(
  seq_len(NBOOT), boot_one, x = dd, future.seed = SEED,
  future.packages = "survival"))
plan(sequential)
colnames(B) <- names(b)
n_fail <- sum(!complete.cases(B))
Bc <- B[complete.cases(B), , drop = FALSE]
r_ecg <- exp(Bc[, "ecg_e"] - Bc[, "ecg_l"])
r_ehr <- exp(Bc[, "ehr_e"] - Bc[, "ehr_l"])
q <- function(v) unname(quantile(v, c(0.025, 0.975)))
ci_ecg <- q(r_ecg)
ci_ehr <- q(r_ehr)
verdict <- if (ci_ecg[1] > 1) {
  "H1 SUPPORTED: ECG prognostic signal concentrated in days 0-30"
} else if (ci_ecg[2] < 1) {
  "REVERSED: ECG prognostic signal stronger in days 31-365"
} else {
  "INCONCLUSIVE: interval includes 1"
}

zph <- cox.zph(coxph(Surv(time_death, event_death) ~ ehr_z + ecg_z,
                     data = dd, ties = "efron"))

rows <- list()
add <- function(quantity, est, lo, hi, method, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    quantity = quantity, estimate = est, lo = lo, hi = hi,
    method = method, note = note)
}
add("patients", nrow(dd), NA, NA, "count")
add("deaths days 0-30", n_early, NA, NA, "count")
add("deaths days 31-365", n_late, NA, NA, "count")
for (v in c("ehr_e", "ehr_l", "ecg_e", "ecg_l")) {
  w <- wald(b[v], sqrt(V[v, v]))
  bq <- q(exp(Bc[, v]))
  lab <- sprintf("HR per SD, %s, days %s",
                 toupper(substr(v, 1, 3)),
                 ifelse(grepl("_e$", v), "0-30", "31-365"))
  add(lab, w[1], w[2], w[3], "Wald (model-based)")
  add(lab, w[1], bq[1], bq[2], "bootstrap percentile")
}
est_ecg <- exp(b["ecg_e"] - b["ecg_l"])
est_ehr <- exp(b["ehr_e"] - b["ehr_l"])
add("PRIMARY: ECG HR ratio, early / late", est_ecg,
    ci_ecg[1], ci_ecg[2], "bootstrap percentile, 2,000 refits",
    verdict)
w <- wald(b["ecg_e"] - b["ecg_l"], lr_se("ecg_e", "ecg_l"))
add("ECG HR ratio, early / late", w[1], w[2], w[3],
    "Wald (information only)")
add("EHR HR ratio, early / late (no decision rule)", est_ehr,
    ci_ehr[1], ci_ehr[2], "bootstrap percentile, 2,000 refits")
add("bootstrap refits that failed", n_fail, NA, NA, "count")
zt <- zph$table
for (r in rownames(zt)) {
  add(sprintf("cox.zph %s (unsplit model)", r), zt[r, "chisq"],
      NA, NA, "chi-square", sprintf("p = %.4f", zt[r, "p"]))
}

out <- do.call(rbind, rows)
fwrite(out, file.path(P2, "results", "primary.csv"))
fwrite(data.table(replicate = seq_len(nrow(B)), B,
                  ecg_ratio = exp(B[, "ecg_e"] - B[, "ecg_l"]),
                  ehr_ratio = exp(B[, "ehr_e"] - B[, "ehr_l"])),
       file.path(P2, "results", "primary_boot.csv"))

cat("\n=== RESULT ===\n")
print(out, row.names = FALSE, digits = 4)
cat(sprintf("\n  PRIMARY: ECG HR ratio early/late %.3f",
            est_ecg),
    sprintf("[%.3f, %.3f]\n  DECISION: %s\n", ci_ecg[1],
            ci_ecg[2], verdict))
cat(sprintf("  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start,
                                units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(P2, "results", "primary.csv"), "\n")
cat("  ", file.path(P2, "results", "primary_boot.csv"), "\n")
cat("  ", file.path(P2, "results", "primary_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: primary.csv\n")