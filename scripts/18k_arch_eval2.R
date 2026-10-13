# 18k_arch_eval2.R
# Evaluation of the forest and intermediate fusion architectures beside
# the linear ones (PREREGISTRATION_ADDENDUM_2.md A26, A27).
#
# For each frame, cohort and modality set, scores the out-of-fold risks of
# death from EHR alone (PH and TV), Early fusion by plain concatenation
# and nested Late fusion (each PH and TV, from 18a), Early RSF and Block
# forest (18i), Intermediate PH and TV (18j) and Block forest, optimised
# (18l). Measures at every horizon: AUROC, IPA, log-likelihood skill score
# and observed/expected, plus the calibration slope (point estimate).
# Horizons are days 10 and 30 in Frame A and days 10, 30, 90, 180 and 365
# in Frame B.
#
# Contrasts (A10): each new model minus Late PH, Late TV and EHR alone PH
# (the fusion gain); Block forest minus Early RSF; Early RSF minus Early
# PH; Intermediate TV minus Intermediate PH; Block forest, optimised minus
# Block forest (A27). Intervals are 95% percentile intervals from a
# 2,000-sample bootstrap (seed 42), clustered by patient in Frame A.
#
# Version 2 (A27) adds the optimised block forest; earlier results are
# unchanged because the bootstrap samples depend only on the patients.
# The script stops if any model lacks a risk for any patient or horizon.
# A timing probe projects the bootstrap time and pauses 5 s before it.
#
# Inputs:  data/processed/frame_A.csv, data/processed/frame_B.csv,
#          data/processed/fusion_arch/arch_preds_<linear|forest|neural|
#          blockopt>_*.csv
# Outputs: results/arch2_metrics.csv, results/arch2_contrasts.csv,
#          results/arch2_checks.csv, results/arch2_log.txt
# Usage:   Rscript scripts/18k_arch_eval2.R [--probe-only] [--workers=N]
#          (--probe-only stops after the timing projection)

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
FA <- file.path(PR, "fusion_arch")
RS <- file.path(P2, "results")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
NBOOT <- 2000L
SEED <- 42L
PKG <- c("data.table")
LIN <- c("EHR alone PH", "EHR alone TV", "Early PH plain concatenation",
         "Early TV plain concatenation", "Late PH", "Late TV")
NEW <- c("Early RSF", "Block forest", "Block forest, optimised",
         "Intermediate PH", "Intermediate TV")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 220)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(J$hz)) {
    y <- J$Y[idx, hi]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    for (k in names(J$P)) {
      p <- pmin(pmax(J$P[[k]][idx, hi], 1e-6), 1 - 1e-6)
      key <- paste(k, J$hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, p)
      v[paste(key, "ipa", sep = "|")] <- 1 - mean((y - p)^2) /
        (pr * (1 - pr))
      v[paste(key, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      v[paste(key, "oe", sep = "|")] <- pr / mean(p)
    }
  }
  v
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                     length(J$rows_by), replace = TRUE)],
                use.names = FALSE)
  stat_vec(J, idx)
}

cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
JOBS <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  key <- paste(fr, co)
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") frames[[fr]] else
    frames[[fr]][in_three_mod == 1]
  hz <- HZ[[fr]]
  ff <- file.path(FA, sprintf("arch_preds_%s_%s_%s.csv",
                              c("linear", "forest", "neural", "blockopt"),
                              fr, code))
  ok <- all(file.exists(ff))
  chk(sprintf("%s: stage 18a, 18i, 18j and 18l files found", key),
      paste(basename(ff), collapse = ", "), if (ok) "PASS" else "FAIL")
  if (!ok) next
  D <- rbindlist(lapply(ff, fread))[model %in% c(LIN, NEW)]
  P <- list()
  nbad <- 0L
  for (k in unique(paste(D$model, D$set, sep = "|"))) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    x <- D[model == z[1] & set == z[2]]
    M <- sapply(hz, function(h) {
      y <- x[horizon == h]
      y$risk[match(sub$hadm_id, y$hadm_id)]
    })
    M <- matrix(M, ncol = length(hz))
    nbad <- nbad + sum(!is.finite(M))
    P[[k]] <- M
  }
  chk(sprintf("%s: models / missing risks", key),
      sprintf("%d / %d", length(P), nbad), if (nbad == 0) "PASS" else
        "FAIL")
  Y <- sapply(hz, function(h) as.integer(sub$event_death == 1 &
                                           sub$time_death <= h + 0.5))
  JOBS[[key]] <- list(fr = fr, co = co, hz = hz, P = P,
                      Y = matrix(Y, ncol = length(hz)),
                      sets = unique(D$set[D$set != "EHR"]),
                      rows_by = split(seq_len(nrow(sub)), match(
                        sub$subject_id, unique(sub$subject_id))))
}
ck <- rbindlist(CHK)
if (any(ck$status == "FAIL")) {
  fwrite(ck, file.path(RS, "arch2_checks.csv"))
  cat("\nSTOP: a check failed. see: arch2_checks.csv\n")
  quit(status = 1)
}
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
plan(multisession, workers = workers)
proj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(workers * 2), boot_rep, J = J,
                          future.seed = 1L, future.packages = PKG))
  proj <- proj + as.numeric(difftime(Sys.time(), t0, units = "secs")) *
    NBOOT / (workers * 2)
}
cat(sprintf("  projected bootstrap: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
MET <- list()
CON <- list()
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  Frame %s %s done\n", J$fr, J$co))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in names(J$P)) for (hi in seq_along(J$hz)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
    sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                   family = binomial))[2]),
                   error = function(e) NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = z[1], set = z[2],
      horizon = J$hz[hi], measure = "slope", estimate = sl, lo = NA_real_,
      hi = NA_real_)
  }
  for (sn in J$sets) {
    K <- function(m) paste(m, if (grepl("^EHR alone", m)) "EHR" else sn,
                           sep = "|")
    cds <- list()
    for (a in NEW) for (b in c("Late PH", "Late TV", "EHR alone PH")) {
      cds[[length(cds) + 1]] <- c(a, b)
    }
    cds <- c(cds, list(c("Block forest", "Early RSF"),
                       c("Block forest, optimised", "Block forest"),
                       c("Early RSF", "Early PH plain concatenation"),
                       c("Intermediate TV", "Intermediate PH")))
    for (cd in cds) for (h in J$hz) for (ms in c("auc", "ipa", "ll")) {
      ka <- paste(K(cd[1]), h, ms, sep = "|")
      kb <- paste(K(cd[2]), h, ms, sep = "|")
      if (!(ka %in% names(pt)) || !(kb %in% names(pt))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, set = sn,
        contrast = sprintf("%s - %s", cd[1], cd[2]), horizon = h,
        measure = ms, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "arch2_metrics.csv"))
fwrite(cn, file.path(RS, "arch2_contrasts.csv"))
fwrite(ck, file.path(RS, "arch2_checks.csv"))
cat("\n=== AUROC, every arm ===\n")
x <- mt[measure == "auc"]
print(dcast(x, frame + cohort + set + model ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE, nrows = 500)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("arch2_metrics.csv", "arch2_contrasts.csv",
            "arch2_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "arch2_log.txt"), " (terminal log)\n")
cat("\nOutputs: arch2_metrics.csv, arch2_contrasts.csv, arch2_checks.csv\n")