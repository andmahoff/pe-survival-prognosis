# 15a_dca.R
# Decision curve analysis of clinical utility (PREREGISTRATION_ADDENDUM.md,
# A2). Exploratory; intervals only.
#
# Net benefit at risk threshold pt, treating everyone whose predicted risk
# is at least pt, is NB = TP/n - FP/n * pt / (1 - pt). Follow-up is
# complete, so the outcome is whether the event happened by the horizon;
# for cardiovascular (CV) readmission a death first counts as no event
# (van Geloven et al., BMJ 2022). The references are treat all,
# NB = prev - (1 - prev) * pt / (1 - pt), and treat none, NB = 0.
#
# Death: Frame A at 30 days and Frame B at 30 and 365 days, in both
# cohorts, for the out-of-fold seed-averaged predictions of EHR alone
# (level-1 ridge-Cox absolute risk), the Cox stack, the binary stack and
# recalibrated late fusion (script 11b). CV readmission: the same frames,
# horizons and cohorts for the four S2 formulations (F1 conditional binary,
# the definition used in the dissertation (pe-multimodal-prognosis-msc);
# F2 unconditional binary; F3 cause-specific Cox; F4 Fine-Gray).
#
# Curves run over thresholds of 1% to 50%, with the net reduction in
# interventions per 100 patients against treat all. At 5%, 10% and 20% the
# script reports every model's net benefit and the pre-specified
# differences with 2,000-sample bootstrap intervals (patient-clustered in
# Frame A, seed 42). For death these are Cox stack minus EHR alone, binary
# stack minus EHR alone and Cox stack minus binary stack; for CV
# readmission, F2, F3 and F4 each minus F1.
#
# Checks: every prediction is present and within [0, 1], otherwise the
# script stops. A timing probe runs a parallel batch of bootstrap samples.
#
# Inputs:  data/processed/frame_A.csv, data/processed/frame_B.csv,
#          data/processed/calibration/cal_preds_*.csv,
#          data/processed/s2/s2_preds_*.csv,
#          data/processed/level1/oof_linear_*.csv
# Outputs: results/dca_curves.csv, results/dca_contrasts.csv,
#          results/dca_checks.csv, results/dca_log.txt
# Usage:   Rscript scripts/15a_dca.R [--probe-only]

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
NBOOT <- 2000L
SEED <- 42L
TH_CURVE <- seq(0.01, 0.50, by = 0.01)
TH_KEY <- c(0.05, 0.10, 0.20)
HZ <- list(A = 30, B = c(30, 365))
DM <- c(ehr = "EHR alone", coxstack = "Cox stack",
        binstack = "binary stack",
        latefusion = "late fusion, recalibrated")
FM <- c(F1 = "F1 conditional binary", F2 = "F2 unconditional binary",
        F3 = "F3 cause-specific Cox", F4 = "F4 Fine-Gray")
PAIRS <- list(
  death = list(c(DM[["coxstack"]], DM[["ehr"]]),
               c(DM[["binstack"]], DM[["ehr"]]),
               c(DM[["coxstack"]], DM[["binstack"]])),
  cv = list(c(FM[["F2"]], FM[["F1"]]), c(FM[["F3"]], FM[["F1"]]),
            c(FM[["F4"]], FM[["F1"]])))
PKG <- c("data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

nb_at <- function(y, p, th) {
  n <- length(y)
  vapply(th, function(t) {
    tr <- p >= t
    (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / n
  }, numeric(1))
}
nb_vec <- function(J, idx, th) {
  nbf <- function(y, p, th) {
    n <- length(y)
    vapply(th, function(t) {
      tr <- p >= t
      (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / n
    }, numeric(1))
  }
  v <- numeric(0)
  for (k in names(J$P)) {
    pz <- strsplit(k, "|", fixed = TRUE)[[1]]
    y <- J$Y[[paste(pz[1], pz[3], sep = "|")]][idx]
    v[paste(k, th, sep = "|")] <- nbf(y, J$P[[k]][idx], th)
  }
  for (yk in names(J$Y)) {
    pz <- strsplit(yk, "|", fixed = TRUE)[[1]]
    y <- J$Y[[yk]][idx]
    prev <- mean(y)
    v[paste(pz[1], "treat all", pz[2], th, sep = "|")] <-
      prev - (1 - prev) * th / (1 - th)
  }
  v
}
boot_nb <- function(i, J, th) {
  if (J$fr == "A") {
    idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                       length(J$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(J$n, J$n, replace = TRUE)
  }
  nb_vec(J, idx, th)
}
avg <- function(path, cols) {
  x <- fread(path, select = c("hadm_id", cols))
  x[, lapply(.SD, mean), by = hadm_id, .SDcols = cols]
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
JOBS <- list()
CHK <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  Fd <- if (fr == "A") fa else fb
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  hs <- HZ[[fr]]
  pc <- file.path(PR, "calibration",
                  sprintf("cal_preds_%s_%s.csv", fr, code))
  p2 <- file.path(PR, "s2", sprintf("s2_preds_%s_%s.csv", fr, code))
  p1 <- file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv", fr, code))
  cc <- intersect(
    as.vector(outer(c("coxstack", "binstack", "latefusion"), hs,
                    function(m, h) sprintf("death|%s|%d", m, h))),
    names(fread(pc, nrows = 0)))
  c2 <- as.vector(outer(names(FM), hs,
                        function(f, h) sprintf("%s_h%d", f, h)))
  c1 <- sprintf("ehr.cox.death.r%d", hs)
  C <- avg(pc, cc)[match(sub$hadm_id, hadm_id)]
  S2 <- avg(p2, c2)[match(sub$hadm_id, hadm_id)]
  L1 <- avg(p1, c1)[match(sub$hadm_id, hadm_id)]
  Y <- list()
  P <- list()
  for (h in hs) {
    Y[[sprintf("death|%d", h)]] <-
      as.integer(sub$event_death == 1 & sub$time_death <= h + 0.5)
    Y[[sprintf("cv|%d", h)]] <-
      as.integer(sub$event_cr == 1 & sub$time_first <= h + 0.5)
    P[[sprintf("death|%s|%d", DM[["ehr"]], h)]] <-
      L1[[sprintf("ehr.cox.death.r%d", h)]]
    for (m in c("coxstack", "binstack", "latefusion")) {
      k <- sprintf("death|%s|%d", m, h)
      if (k %in% names(C)) P[[sprintf("death|%s|%d", DM[[m]], h)]] <- C[[k]]
    }
    for (f in names(FM)) {
      P[[sprintf("cv|%s|%d", FM[[f]], h)]] <- S2[[sprintf("%s_h%d", f, h)]]
    }
  }
  pm <- do.call(cbind, P)
  nna <- sum(is.na(pm))
  rng <- range(pm, na.rm = TRUE)
  ok <- nna == 0 && rng[1] >= -1e-9 && rng[2] <= 1 + 1e-9
  CHK[[length(CHK) + 1]] <- data.table(
    frame = fr, cohort = co, prediction_series = length(P),
    missing = nna, min_risk = rng[1], max_risk = rng[2],
    status = if (ok) "PASS" else "FAIL")
  cat(sprintf("  %s %-7s n=%d  series %d  missing %d  range %.5f-%.5f  %s\n",
              fr, co, nrow(sub), length(P), nna, rng[1], rng[2],
              if (ok) "PASS" else "FAIL"))
  JOBS[[paste(fr, co)]] <- list(
    fr = fr, co = co, n = nrow(sub), hs = hs, Y = Y, P = P,
    rows_by = split(seq_len(nrow(sub)),
                    match(sub$subject_id, unique(sub$subject_id))))
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "dca_checks.csv"))
if (any(ck$status == "FAIL")) {
  cat("HARD CHECK FAILED: stopping. Outputs: dca_checks.csv\n")
  quit(status = 1)
}

# ---------------- timing probe ----------------
cat(sprintf("\n=== TIMING PROBE: bootstrap (%d workers) ===\n", workers))
plan(multisession, workers = workers)
nb <- workers * 2L
proj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), boot_nb, J = J, th = TH_KEY,
                          future.seed = 1L, future.packages = PKG))
  proj <- proj + as.numeric(difftime(Sys.time(), t0, units = "secs")) *
    NBOOT / nb
}
cat(sprintf("  projected bootstrap: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- curves, bootstrap, contrasts ----------------
CUR <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cat(sprintf("=== Frame %s %s ===\n", J$fr, J$co))
  for (k in names(J$P)) {
    pz <- strsplit(k, "|", fixed = TRUE)[[1]]
    y <- J$Y[[paste(pz[1], pz[3], sep = "|")]]
    prev <- mean(y)
    nbm <- nb_at(y, J$P[[k]], TH_CURVE)
    nba <- prev - (1 - prev) * TH_CURVE / (1 - TH_CURVE)
    CUR[[length(CUR) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, endpoint = pz[1],
      horizon = as.integer(pz[3]), model = pz[2], threshold = TH_CURVE,
      prevalence = prev, net_benefit = nbm, nb_treat_all = nba,
      net_reduction_per_100 = (nbm - nba) /
        (TH_CURVE / (1 - TH_CURVE)) * 100)
  }
  B <- do.call(rbind, future_lapply(
    seq_len(NBOOT), boot_nb, J = J, th = TH_KEY, future.seed = SEED,
    future.packages = PKG))
  pt <- nb_vec(J, seq_len(J$n), TH_KEY)
  for (k in names(J$P)) {
    pz <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (t in TH_KEY) {
      kk <- paste(k, t, sep = "|")
      q <- q2(B[, kk])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, endpoint = pz[1],
        horizon = as.integer(pz[3]), threshold = t, type = "net benefit",
        name = pz[2], estimate = pt[[kk]], lo = q[1], hi = q[2],
        excludes_zero = NA)
    }
  }
  for (yk in names(J$Y)) {
    pz <- strsplit(yk, "|", fixed = TRUE)[[1]]
    for (t in TH_KEY) {
      kk <- paste(pz[1], "treat all", pz[2], t, sep = "|")
      q <- q2(B[, kk])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, endpoint = pz[1],
        horizon = as.integer(pz[2]), threshold = t, type = "net benefit",
        name = "treat all", estimate = pt[[kk]], lo = q[1], hi = q[2],
        excludes_zero = NA)
    }
    for (pr in PAIRS[[pz[1]]]) for (t in TH_KEY) {
      a <- paste(pz[1], pr[1], pz[2], t, sep = "|")
      b <- paste(pz[1], pr[2], pz[2], t, sep = "|")
      if (!all(c(a, b) %in% names(pt))) next
      q <- q2(B[, a] - B[, b])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, endpoint = pz[1],
        horizon = as.integer(pz[2]), threshold = t, type = "difference",
        name = sprintf("%s minus %s", pr[1], pr[2]),
        estimate = pt[[a]] - pt[[b]], lo = q[1], hi = q[2],
        excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
cu <- rbindlist(CUR)
cn <- rbindlist(CON)
fwrite(cu, file.path(RS, "dca_curves.csv"))
fwrite(cn, file.path(RS, "dca_contrasts.csv"))

cat("\n=== DIFFERENCES IN NET BENEFIT (per patient) ===\n")
print(cn[type == "difference",
         .(frame, cohort, endpoint, horizon, threshold, name,
           estimate = round(estimate, 4), lo = round(lo, 4),
           hi = round(hi, 4), excludes_zero)], row.names = FALSE)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("dca_curves.csv", "dca_contrasts.csv", "dca_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "dca_log.txt"), " (terminal log)\n")
cat("\nOutputs: dca_curves.csv, dca_contrasts.csv, dca_checks.csv\n")