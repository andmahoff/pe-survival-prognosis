# 18c_arch_eval.R
# Evaluates every fusion architecture found in data/processed/fusion_arch
# (PREREGISTRATION_ADDENDUM_2.md, A10; also D10, D11, A11, A13).
#
# Picks up whatever prediction files are present (linear arms from 18a, 18b
# and 18d, late time-varying-coefficient stacks from 18e, and the forest and
# neural arms when they exist), so it is re-run after each new stage. Files
# in subfolders, such as archived runs, are ignored. It stops if any risk is
# missing or if a model appears in more than one prediction file.
#
# For each model, set and horizon: AUROC, IPA, log-likelihood skill score
# (1 - LL / LL of the event-rate model), observed/expected, and the
# calibration slope (point estimate). Contrasts are computed where both
# models exist:
#   architecture vs late     each fused model minus the matching late model
#                            (PH vs Late PH, TV vs Late TV, horizon-tuned
#                            vs Late PH horizon-tuned)
#   TV vs PH                 for EHR alone, early (both versions), late and
#                            intermediate
#   horizon tuning (A11)     horizon-tuned minus whole-follow-up PH, and TV
#                            minus horizon-tuned, for EHR alone, early, late
#   block-aware early        block penalties, EHR unpenalised and plain
#                            concatenation vs Early PH; plain vs block-
#                            standardised TV; horizon-tuned block penalties
#                            vs horizon-tuned Early PH
#   time-varying coefs (A13) Late TVC minus Late PH, Late TV and Late PH
#                            horizon-tuned; Late TVC log-time minus Late TVC
#   fusion gain              each fused model minus the matching EHR alone
#   dense vs sparse          three-modality cohort: (Early - Late) with
#                            EHR+ECG minus (Early - Late) with EHR+CTPA, for
#                            PH, TV, plain and horizon-tuned
# Intervals are 95% percentile intervals from 2,000 bootstrap resamples
# (seed 42), clustered by patient in Frame A. The bootstrap runs in
# parallel with future.apply after a timing probe that pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, data/processed/frame_B.csv,
#          data/processed/fusion_arch/arch_preds_*.csv
# Outputs: results/arch_metrics.csv, results/arch_contrasts.csv,
#          results/arch_checks.csv, results/arch_log.txt
# Usage:   Rscript scripts/18c_arch_eval.R [--probe-only]

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OD <- file.path(PR, "fusion_arch")
NBOOT <- 2000L
SEED <- 42L
PKG <- c("data.table")
options(width = 220)
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-60s %s\n", status, item, value))
}
TVM <- c("EHR alone TV", "Early TV", "Early TV plain concatenation",
         "Late TV", "Intermediate TV")

stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (h in J$hz) {
    y <- J$y[[as.character(h)]][idx]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    for (k in J$keys[[as.character(h)]]) {
      p <- pmin(pmax(J$P[[k]][idx], 1e-6), 1 - 1e-6)
      v[paste(k, "auc", sep = "|")] <- auc(y, p)
      v[paste(k, "ipa", sep = "|")] <- 1 - mean((y - p)^2) / (pr * (1 - pr))
      v[paste(k, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      v[paste(k, "oe", sep = "|")] <- pr / mean(p)
    }
  }
  v
}
boot_rep <- function(i, J) {
  if (J$fr == "A") {
    idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                       length(J$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(J$n, J$n, replace = TRUE)
  }
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
files <- list.files(OD, pattern = "^arch_preds_.*\\.csv$", full.names = TRUE)
chk("prediction files found", length(files),
    if (length(files) > 0) "PASS" else "FAIL")
if (!length(files)) {
  fwrite(rbindlist(CHK), file.path(RS, "arch_checks.csv"))
  quit(status = 1)
}
JOBS <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  ff <- files[grepl(sprintf("_%s_%s\\.csv$", fr, code), files)]
  if (!length(ff)) next
  L <- rbindlist(lapply(ff, function(f) fread(f)[, src := basename(f)]),
                 use.names = TRUE)
  dup <- L[, .(files = uniqueN(src)), by = .(model, set)][files > 1]
  chk(sprintf("%s %s: models found in more than one file", fr, co),
      if (nrow(dup)) paste(unique(dup$model), collapse = "; ") else "none",
      if (nrow(dup) == 0) "PASS" else "FAIL")
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  L[, key := paste(model, set, horizon, sep = "|")]
  keys <- unique(L$key)
  P <- list()
  miss <- 0L
  for (k in keys) {
    x <- L[key == k]
    P[[k]] <- x$risk[match(sub$hadm_id, x$hadm_id)]
    miss <- miss + sum(is.na(P[[k]]))
  }
  chk(sprintf("%s %s: models x sets x horizons / missing risks", fr, co),
      sprintf("%d / %d", length(keys), miss),
      if (miss == 0) "PASS" else "FAIL")
  hz <- sort(unique(L$horizon))
  J <- list(fr = fr, co = co, n = nrow(sub), hz = hz, P = P,
            y = setNames(lapply(hz, function(h) {
              as.integer(sub$event_death == 1 & sub$time_death <= h + 0.5)
            }), as.character(hz)),
            keys = setNames(lapply(hz, function(h) {
              keys[L[match(keys, key), horizon] == h]
            }), as.character(hz)),
            rows_by = split(seq_len(nrow(sub)),
                            match(sub$subject_id, unique(sub$subject_id))),
            models = unique(L[, .(model, set)]))
  JOBS[[paste(fr, co)]] <- J
  cat(sprintf("  Frame %s %-7s n=%d, %d models x sets, horizons %s\n", fr,
              co, nrow(sub), nrow(J$models), paste(hz, collapse = ", ")))
}
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "arch_checks.csv"))
  cat("\nSTOP: a check failed. see: arch_checks.csv\n")
  quit(status = 1)
}

# contrasts: list of (a, b) keys without the horizon, or a difference of
# differences (a - b) - (c - d)
HTM <- function(m) grepl("horizon-tuned", m, fixed = TRUE)
TVCM <- function(m) grepl("TVC", m, fixed = TRUE)
LATES <- c("Late PH", "Late TV", "Late PH horizon-tuned", "Late TVC",
           "Late TVC log-time")
fusion_base <- function(m) {
  if (m %in% TVM) "EHR alone TV"
  else if (HTM(m)) "EHR alone PH horizon-tuned"
  else if (TVCM(m)) "EHR alone TVC" else "EHR alone PH"
}
late_base <- function(m) {
  if (m %in% TVM) "Late TV"
  else if (HTM(m)) "Late PH horizon-tuned" else "Late PH"
}
contrast_defs <- function(M) {
  have <- function(m, s) nrow(M[model == m & set == s]) > 0
  out <- list()
  add <- function(label, type, a, b, c = NULL, d = NULL) {
    out[[length(out) + 1]] <<- list(label = label, type = type, a = a,
                                    b = b, c = c, d = d)
  }
  pair <- function(p, s, type) {
    if (have(p[1], s) && have(p[2], s)) {
      add(sprintf("%s - %s [%s]", p[1], p[2], s), type,
          paste(p[1], s, sep = "|"), paste(p[2], s, sep = "|"))
    }
  }
  # EHR alone: time-varying and horizon tuning
  pair(c("EHR alone TV", "EHR alone PH"), "EHR", "time-varying vs PH")
  pair(c("EHR alone PH horizon-tuned", "EHR alone PH"), "EHR",
       "horizon tuning")
  pair(c("EHR alone TV", "EHR alone PH horizon-tuned"), "EHR",
       "TV vs horizon-tuned")
  sets <- setdiff(unique(M$set), "EHR")
  for (s in sets) {
    for (m in unique(M[set == s, model])) {
      base <- fusion_base(m)
      if (have(base, "EHR")) {
        add(sprintf("%s [%s] - %s", m, s, base), "fusion gain",
            paste(m, s, sep = "|"), paste(base, "EHR", sep = "|"))
      }
      if (!(m %in% LATES)) {
        pair(c(m, late_base(m)), s, "architecture vs late")
      }
    }
    for (p in list(c("Late TV", "Late PH"), c("Early TV", "Early PH"),
                   c("Early TV plain concatenation",
                     "Early PH plain concatenation"),
                   c("Intermediate TV", "Intermediate PH"))) {
      pair(p, s, "time-varying vs PH")
    }
    for (p in list(c("Early PH horizon-tuned", "Early PH"),
                   c("Early PH block penalties horizon-tuned",
                     "Early PH block penalties"),
                   c("Late PH horizon-tuned", "Late PH"))) {
      pair(p, s, "horizon tuning")
    }
    for (p in list(c("Early TV", "Early PH horizon-tuned"),
                   c("Late TV", "Late PH horizon-tuned"))) {
      pair(p, s, "TV vs horizon-tuned")
    }
    for (p in list(c("Late TVC", "Late PH"), c("Late TVC", "Late TV"),
                   c("Late TVC", "Late PH horizon-tuned"),
                   c("Late TVC log-time", "Late TVC"))) {
      pair(p, s, "time-varying coefficients")
    }
    for (p in list(c("Early PH block penalties", "Early PH"),
                   c("Early PH EHR unpenalised", "Early PH"),
                   c("Early PH plain concatenation", "Early PH"),
                   c("Early TV plain concatenation", "Early TV"),
                   c("Early PH block penalties horizon-tuned",
                     "Early PH horizon-tuned"))) {
      pair(p, s, "block-aware early")
    }
  }
  for (p in list(c("Early PH", "Late PH"), c("Early TV", "Late TV"),
                 c("Early PH plain concatenation", "Late PH"),
                 c("Early TV plain concatenation", "Late TV"),
                 c("Early PH horizon-tuned", "Late PH horizon-tuned"))) {
    e <- p[1]
    l <- p[2]
    if (have(e, "EHR+ECG") && have(e, "EHR+CTPA") &&
        have(l, "EHR+ECG") && have(l, "EHR+CTPA")) {
      add(sprintf("(%s - %s) with ECG minus with CTPA", e, l),
          "dense vs sparse", paste(e, "EHR+ECG", sep = "|"),
          paste(l, "EHR+ECG", sep = "|"), paste(e, "EHR+CTPA", sep = "|"),
          paste(l, "EHR+CTPA", sep = "|"))
    }
  }
  out
}

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
plan(multisession, workers = workers)
nb <- workers * 2L
proj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), boot_rep, J = J, future.seed = 1L,
                          future.packages = PKG))
  proj <- proj + as.numeric(difftime(Sys.time(), t0, units = "secs")) *
    NBOOT / nb
}
cat(sprintf("  projected: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cat(sprintf("=== Frame %s %s: %d bootstrap samples ===\n", J$fr, J$co,
              NBOOT))
  pt <- stat_vec(J, seq_len(J$n))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in unlist(J$keys)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    h <- z[3]
    p <- pmin(pmax(J$P[[k]], 1e-6), 1 - 1e-6)
    y <- J$y[[h]]
    sl <- tryCatch(unname(coef(glm(y ~ qlogis(p), family = binomial))[2]),
                   error = function(e) NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = z[1], set = z[2],
      horizon = as.numeric(h), measure = "slope", estimate = sl,
      lo = NA_real_, hi = NA_real_)
  }
  for (cd in contrast_defs(J$models)) {
    for (h in J$hz) for (ms in c("auc", "ipa", "ll")) {
      ka <- paste(cd$a, h, ms, sep = "|")
      kb <- paste(cd$b, h, ms, sep = "|")
      if (!(ka %in% names(pt)) || !(kb %in% names(pt))) next
      if (is.null(cd$c)) {
        est <- pt[[ka]] - pt[[kb]]
        bv <- B[, ka] - B[, kb]
      } else {
        kc <- paste(cd$c, h, ms, sep = "|")
        kd <- paste(cd$d, h, ms, sep = "|")
        est <- (pt[[ka]] - pt[[kb]]) - (pt[[kc]] - pt[[kd]])
        bv <- (B[, ka] - B[, kb]) - (B[, kc] - B[, kd])
      }
      q <- q2(bv)
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, type = cd$type, contrast = cd$label,
        horizon = h, measure = ms, estimate = est, lo = q[1], hi = q[2],
        excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "arch_metrics.csv"))
fwrite(cn, file.path(RS, "arch_contrasts.csv"))
fwrite(ck, file.path(RS, "arch_checks.csv"))

cat("\n=== AUROC, Frame B, EHR+ECG cohort ===\n")
print(dcast(mt[frame == "B" & cohort == "EHR+ECG" & measure == "auc"],
            model + set ~ horizon, value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 4)
                       else x)], row.names = FALSE)
cat("\n=== CALIBRATION SLOPE, Frame B, EHR+ECG cohort ===\n")
print(dcast(mt[frame == "B" & cohort == "EHR+ECG" & measure == "slope"],
            model + set ~ horizon, value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                       else x)], row.names = FALSE)
cat("\n=== ARCHITECTURE CONTRASTS, AUROC ===\n")
print(cn[measure == "auc" & type != "fusion gain",
         .(frame, cohort, contrast, horizon, est = round(estimate, 4),
           lo = round(lo, 4), hi = round(hi, 4), excludes_zero)],
      row.names = FALSE, nrows = 2000)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("arch_metrics.csv", "arch_contrasts.csv", "arch_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "arch_log.txt"), " (terminal log)\n")
cat("\nOutputs: arch_metrics.csv, arch_contrasts.csv, arch_checks.csv\n")