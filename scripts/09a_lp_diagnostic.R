# 09a_lp_diagnostic.R
# Diagnosis of extreme level-1 scores (deviation D3), and a sensitivity
# check of the primary contrast to them.
#
# Part 1: for every level-1 ridge-Cox linear predictor (all frames,
# cohorts, modalities, targets and seeds), the range, the 0.5th, 50th and
# 99.5th percentiles, how many values sit at the +-50 cap used in stage 6a,
# and how many seed-averaged scores lie beyond |z| > 5 and |z| > 10.
# Part 2: for patients with an extreme (|z| > 10) seed-averaged score for a
# modality, which of that modality's input features are themselves extreme
# (|z| > 10 across the cohort). Counts only.
# Part 3: the primary (07_primary.R) refitted three ways, each with 2,000
# patient-level bootstrap refits (seed 42): the registered scores (as in
# stage 7), scores winsorised at their 0.5th and 99.5th percentiles before
# standardising, and rank-based normal scores. The registered primary
# stays the primary.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, feature_whitelist.csv
#          data/processed/level1/oof_linear_*.csv
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: results/lp_diag.csv, results/lp_extreme_features.csv,
#          results/primary_robust.csv, results/lp_diag_log.txt
# Usage:   Rscript scripts/09a_lp_diagnostic.R [--probe-only]
# Notes:   A timing probe of 20 refits projects the three bootstraps and
#          pauses 5 s; --probe-only stops after it. The bootstraps run in
#          parallel (future.apply, physical cores minus 1).

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
IN <- file.path(P2, "data", "inputs")
RS <- file.path(P2, "results")
NBOOT <- 2000L
SEED <- 42L
CUT <- 30.5
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

FILES <- data.table(frame = c("A", "A", "B", "B"),
                    cohort = c("EHR+ECG", "3-mod", "EHR+ECG", "3-mod"),
                    code = c("ehrecg", "3mod", "ehrecg", "3mod"))
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- list(ehr = fread(file.path(IN, "ehr28.csv")),
             ecg = fread(file.path(IN, "ecg_opt.csv")),
             ctpa = fread(file.path(IN, "ctpa46.csv")))

cat("=== PART 1 AND 2: level-1 linear predictors ===\n")
P1 <- list()
P2L <- list()
for (i in seq_len(nrow(FILES))) {
  fr <- FILES$frame[i]
  co <- FILES$cohort[i]
  path <- file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv", fr,
                                          FILES$code[i]))
  hdr <- names(fread(path, nrows = 0))
  lpc <- grep("\\.lp$", hdr, value = TRUE)
  x <- fread(path, select = c("seed", "hadm_id", lpc))
  sub <- if (co == "EHR+ECG") frames[[fr]] else
    frames[[fr]][in_three_mod == 1]
  avg <- x[, lapply(.SD, mean), by = hadm_id, .SDcols = lpc]
  avg <- avg[match(sub$hadm_id, hadm_id)]
  for (cc in lpc) {
    v <- x[[cc]]
    a <- avg[[cc]]
    z <- (a - mean(a)) / sd(a)
    q <- quantile(v, c(0.005, 0.5, 0.995), names = FALSE)
    P1[[length(P1) + 1]] <- data.table(
      frame = fr, cohort = co, column = cc, n_values = length(v),
      min = min(v), p005 = q[1], median = q[2], p995 = q[3],
      max = max(v), at_cap = sum(abs(v) >= 49.99),
      max_abs_z = max(abs(z)), n_abs_z_gt5 = sum(abs(z) > 5),
      n_abs_z_gt10 = sum(abs(z) > 10))
  }
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  for (m in mods) {
    mc <- grep(sprintf("^%s\\.cox\\..*\\.lp$", m), lpc, value = TRUE)
    zz <- sapply(mc, function(cc) {
      a <- avg[[cc]]
      (a - mean(a)) / sd(a)
    })
    ext <- which(apply(abs(zz) > 10, 1, any))
    fcols <- wl[modality == m, column]
    M <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..fcols])
    FZ <- suppressWarnings(scale(M))
    if (length(ext) == 0) {
      P2L[[length(P2L) + 1]] <- data.table(
        frame = fr, cohort = co, modality = m, n_extreme_patients = 0L,
        feature = "(none)", n_patients = 0L, max_abs_feature_z = NA_real_)
      next
    }
    E <- abs(FZ[ext, , drop = FALSE])
    cnt <- colSums(E > 10, na.rm = TRUE)
    mx <- suppressWarnings(apply(E, 2, max, na.rm = TRUE))
    keep <- order(-cnt)[seq_len(min(15, sum(cnt > 0)))]
    if (length(keep) == 0 || sum(cnt > 0) == 0) {
      P2L[[length(P2L) + 1]] <- data.table(
        frame = fr, cohort = co, modality = m,
        n_extreme_patients = length(ext),
        feature = "(no feature beyond |z| 10)", n_patients = 0L,
        max_abs_feature_z = max(mx[is.finite(mx)]))
      next
    }
    for (j in keep) {
      P2L[[length(P2L) + 1]] <- data.table(
        frame = fr, cohort = co, modality = m,
        n_extreme_patients = length(ext), feature = colnames(FZ)[j],
        n_patients = unname(cnt[j]), max_abs_feature_z = unname(mx[j]))
    }
  }
}
p1 <- rbindlist(P1)
p2 <- rbindlist(P2L)
fwrite(p1, file.path(RS, "lp_diag.csv"))
fwrite(p2, file.path(RS, "lp_extreme_features.csv"))
cat("  columns with any value at the +-50 cap or |z| > 10:\n")
print(p1[at_cap > 0 | n_abs_z_gt10 > 0,
         .(frame, cohort, column, max = round(max, 1),
           at_cap, max_abs_z = round(max_abs_z, 1), n_abs_z_gt10)],
      row.names = FALSE)

cat("\n=== PART 3: the primary with robust scores ===\n")
oof <- fread(file.path(PR, "level1", "oof_linear_B_ehrecg.csv"),
             select = c("seed", "hadm_id", "ehr.cox.death.lp",
                        "ecg.cox.death.lp"))
sc <- oof[, .(ehr = mean(`ehr.cox.death.lp`),
              ecg = mean(`ecg.cox.death.lp`)), by = hadm_id]
d <- merge(frames$B[, .(hadm_id, time_death, event_death)], sc,
           by = "hadm_id")
std <- function(x) (x - mean(x)) / sd(x)
wins <- function(x) {
  q <- quantile(x, c(0.005, 0.995), names = FALSE)
  std(pmin(pmax(x, q[1]), q[2]))
}
rank_norm <- function(x) qnorm((rank(x) - 0.5) / length(x))
VARS <- list(registered = std, winsorised = wins,
             rank_normal = rank_norm)

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

dd0 <- as.data.frame(d[, .(time_death, event_death,
                           ehr_z = std(ehr), ecg_z = std(ecg))])
set.seed(1)
tp <- system.time(for (i in 1:20) boot_one(i, dd0))[["elapsed"]] / 20
cat(sprintf("  %.3f s per refit; 3 x %d refits on %d workers: about",
            tp, NBOOT, workers),
    sprintf("%.1f min\n", (3 * tp * NBOOT / workers * 1.3 + 30) / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

plan(multisession, workers = workers)
R3 <- list()
q <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (vn in names(VARS)) {
  fn <- VARS[[vn]]
  dd <- as.data.frame(d[, .(time_death, event_death)])
  dd$ehr_z <- fn(d$ehr)
  dd$ecg_z <- fn(d$ecg)
  f0 <- fit_split(dd)
  b <- coef(f0)
  B <- do.call(rbind, future_lapply(
    seq_len(NBOOT), boot_one, x = dd, future.seed = SEED,
    future.packages = "survival"))
  colnames(B) <- names(b)
  Bc <- B[complete.cases(B), , drop = FALSE]
  add <- function(qty, est, v) {
    ci <- q(v)
    R3[[length(R3) + 1]] <<- data.table(
      variant = vn, quantity = qty, estimate = est,
      lo = ci[1], hi = ci[2],
      max_abs_z_ehr = max(abs(dd$ehr_z)),
      max_abs_z_ecg = max(abs(dd$ecg_z)))
  }
  add("ECG HR ratio, early / late", unname(exp(b["ecg_e"] - b["ecg_l"])),
      exp(Bc[, "ecg_e"] - Bc[, "ecg_l"]))
  add("EHR HR ratio, early / late", unname(exp(b["ehr_e"] - b["ehr_l"])),
      exp(Bc[, "ehr_e"] - Bc[, "ehr_l"]))
  for (v in c("ehr_e", "ehr_l", "ecg_e", "ecg_l")) {
    add(sprintf("HR per SD, %s, days %s", toupper(substr(v, 1, 3)),
                ifelse(grepl("_e$", v), "0-30", "31-365")),
        unname(exp(b[v])), exp(Bc[, v]))
  }
  cat(sprintf("  %-12s ECG ratio %.3f [%.3f, %.3f]   max|z| EHR %.1f ECG %.1f\n",
              vn, R3[[length(R3) - 5]]$estimate, R3[[length(R3) - 5]]$lo,
              R3[[length(R3) - 5]]$hi, max(abs(dd$ehr_z)),
              max(abs(dd$ecg_z))))
}
plan(sequential)
r3 <- rbindlist(R3)
fwrite(r3, file.path(RS, "primary_robust.csv"))
cat("\n")
print(r3[, .(variant, quantity, estimate = round(estimate, 3),
             lo = round(lo, 3), hi = round(hi, 3))], row.names = FALSE)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("lp_diag.csv", "lp_extreme_features.csv",
            "primary_robust.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "lp_diag_log.txt"), " (terminal log)\n")
cat("\nOutputs: lp_diag.csv, lp_extreme_features.csv, primary_robust.csv\n")