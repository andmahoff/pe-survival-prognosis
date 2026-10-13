# 11a_eval_s1_s3.R
# Formal intervals for S1 and S3 (secondary analyses, PREREGISTRATION.md,
# section 6).
#
# S1: Frame A (clock at admission, 30 days), both cohorts; death,
# composite and CV readmission (competing risks, cause-specific). Every
# fused score (binary or Cox level 1, weights by AUROC or by Uno's C) and
# every single modality is scored by AUROC at 30.5 days and Uno's C
# truncated at 31. Paired contrasts are binary minus Cox (fused, and each
# modality), weights by Uno's C minus weights by AUROC, and fusion gain
# over EHR alone. The bootstrap is patient-clustered: a resampled patient
# brings all of their admissions.
#
# S3: Frame B (clock at discharge, death), both cohorts, at 10, 30, 90,
# 180 and 365 days. Reported are the AUROC of every score, Cox minus
# binary, and the fusion gain at each horizon. The horizon curve itself is
# the set of weights that best discriminate at each horizon in the whole
# cohort (on the Cox scores), with intervals, and the contrasts of the ECG
# weight at 10 days minus that at 90, 180 and 365 days (CTPA too in the
# three-modality cohort). The bootstrap resamples patients (one row each).
#
# All scores are seed-averaged out-of-fold scores, held fixed; only the
# evaluation is resampled, with 2,000 bootstrap samples (seed 42).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv
#          data/processed/level1/oof_linear_*.csv
#          data/processed/fusion/fused_*.csv
# Outputs: results/eval_s1.csv, results/eval_s3.csv,
#          results/eval_s3_weights.csv, results/eval_s1_s3_log.txt
# Usage:   Rscript scripts/11a_eval_s1_s3.R [--probe-only]
# Notes:   A timing probe of real replicates, run in parallel on physical
#          cores minus 1, projects the run time and pauses 5 s;
#          --probe-only stops after it.

suppressPackageStartupMessages({
  library(survival)
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
HB <- c(10, 30, 90, 180, 365)
PKG <- c("survival", "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

auc_fast <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cidx <- function(tt, ee, s, ymax) {
  out <- tryCatch(
    concordance(Surv(tt, ee) ~ s, reverse = TRUE,
                timewt = "n/G2", ymax = ymax)$concordance,
    error = function(err) NA_real_)
  as.numeric(out)
}
avg_scores <- function(path, cols) {
  x <- fread(path, select = c("hadm_id", cols))
  x[, lapply(.SD, mean), by = hadm_id, .SDcols = cols]
}
attach_cols <- function(X, L, cols) {
  m <- match(X$hadm_id, L$hadm_id)
  for (cc in cols) set(X, j = cc, value = L[[cc]][m])
  X
}
grid_w <- function(k, step = 0.05) {
  n <- round(1 / step)
  if (k == 2) return(cbind((n:0) / n, (0:n) / n))
  g <- NULL
  for (i in 0:n) for (j in 0:(n - i)) {
    g <- rbind(g, c(i / n, j / n, (n - i - j) / n))
  }
  g
}

fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))

make_A <- function(co) {
  # local copies travel with the closures to the parallel workers
  # (future does not look inside a list for the helpers it uses)
  auc_fast <- auc_fast
  cidx <- cidx
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fa else fa[in_three_mod == 1]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  E <- list(
    death = list(cox = "cox.death.lp", lr = "lr.death_h30",
                 t = sub$time_death, e = sub$event_death == 1),
    composite = list(cox = "cox.composite.lp",
                     lr = "lr.composite_h30",
                     t = sub$time_first, e = sub$event_comp == 1),
    cv = list(cox = "cox.cs_cv.lp", lr = "lr.cvu_h30",
              t = sub$time_first, e = sub$event_cr == 1))
  SC <- list()
  ucols <- character(0)
  fcols <- character(0)
  for (en in names(E)) {
    s <- list()
    for (l1 in c("lr", "cox")) for (cr in c("auc", "uno")) {
      cn <- sprintf("fuse.%s.%s.%s", l1, cr, en)
      s[[sprintf("fused %s %s", l1, cr)]] <- cn
      fcols <- c(fcols, cn)
    }
    for (m in mods) for (l1 in c("cox", "lr")) {
      cn <- paste0(m, ".", E[[en]][[l1]])
      s[[sprintf("%s %s", toupper(m), l1)]] <- cn
      ucols <- c(ucols, cn)
    }
    SC[[en]] <- s
  }
  ucols <- unique(ucols)
  X <- sub[, .(hadm_id, subject_id)]
  X <- attach_cols(X, avg_scores(file.path(
    PR, "level1", sprintf("oof_linear_A_%s.csv", code)), ucols), ucols)
  X <- attach_cols(X, avg_scores(file.path(
    PR, "fusion", sprintf("fused_A_%s.csv", code)), fcols), fcols)
  if (anyNA(X)) stop("missing scores in Frame A ", co)
  lab <- lapply(E, function(z) z$e & z$t <= 30.5)
  Xm <- as.list(X)
  met <- function(idx) {
    v <- numeric(0)
    for (en in names(E)) {
      y <- lab[[en]][idx]
      tt <- E[[en]]$t[idx]
      ee <- E[[en]]$e[idx]
      for (nm in names(SC[[en]])) {
        s <- Xm[[SC[[en]][[nm]]]][idx]
        v[paste(en, "auc30", nm, sep = "|")] <- auc_fast(y, s)
        v[paste(en, "unoC", nm, sep = "|")] <- cidx(tt, ee, s, 31)
      }
    }
    v
  }
  us <- unique(X$subject_id)
  rows_by <- split(seq_len(nrow(X)), match(X$subject_id, us))
  draw <- function() {
    unlist(rows_by[sample.int(length(us), length(us), replace = TRUE)],
           use.names = FALSE)
  }
  list(frame = "A", cohort = co, n = nrow(X), mods = mods, SC = SC,
       met = met, draw = draw, events = sapply(lab, sum))
}

make_B <- function(co) {
  # local copies travel with the closures to the parallel workers
  auc_fast <- auc_fast
  HB <- HB
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  fcols <- c(sprintf("fuse.cox.auc%d.death", HB),
             sprintf("fuse.lr.auc%d.death", HB))
  ucols <- c(paste0(mods, ".cox.death.lp"),
             as.vector(outer(mods, HB, function(m, h)
               sprintf("%s.lr.death_h%d", m, h))))
  X <- sub[, .(hadm_id, subject_id)]
  X <- attach_cols(X, avg_scores(file.path(
    PR, "level1", sprintf("oof_linear_B_%s.csv", code)), ucols), ucols)
  X <- attach_cols(X, avg_scores(file.path(
    PR, "fusion", sprintf("fused_B_%s.csv", code)), fcols), fcols)
  if (anyNA(X)) stop("missing scores in Frame B ", co)
  td <- sub$time_death
  ed <- sub$event_death == 1
  lab <- lapply(HB, function(h) ed & td <= h + 0.5)
  names(lab) <- HB
  R <- sapply(mods, function(m) frank(X[[paste0(m, ".cox.death.lp")]]) /
                nrow(X))
  G <- grid_w(length(mods))
  Xm <- as.list(X)
  met <- function(idx) {
    v <- numeric(0)
    Ri <- R[idx, , drop = FALSE]
    Fw <- Ri %*% t(G)
    for (h in HB) {
      y <- lab[[as.character(h)]][idx]
      sc <- c("fused cox" = sprintf("fuse.cox.auc%d.death", h),
              "fused lr" = sprintf("fuse.lr.auc%d.death", h))
      for (m in mods) {
        sc[paste(toupper(m), "cox")] <- paste0(m, ".cox.death.lp")
        sc[paste(toupper(m), "lr")] <- sprintf("%s.lr.death_h%d", m, h)
      }
      for (nm in names(sc)) {
        v[paste(h, "auc", nm, sep = "|")] <-
          auc_fast(y, Xm[[sc[[nm]]]][idx])
      }
      a <- apply(Fw, 2, function(s) auc_fast(y, s))
      w <- G[which.max(a), ]
      for (j in seq_along(mods)) {
        v[paste(h, "weight", toupper(mods[j]), sep = "|")] <- w[j]
      }
    }
    v
  }
  draw <- function() sample.int(nrow(X), nrow(X), replace = TRUE)
  list(frame = "B", cohort = co, n = nrow(X), mods = mods,
       met = met, draw = draw, events = sapply(lab, sum))
}

cat("=== BUILD ===\n")
jobs <- list(make_A("EHR+ECG"), make_A("3-mod"),
             make_B("EHR+ECG"), make_B("3-mod"))
for (J in jobs) {
  cat(sprintf("  Frame %s %-7s n=%d events %s\n", J$frame, J$cohort,
              J$n, paste(J$events, collapse = "/")))
}

cat(sprintf("\n=== TIMING PROBE (%d workers, replicates in parallel) ===\n",
            workers))
plan(multisession, workers = workers)
nb <- workers * 2L
proj <- 0
for (J in jobs) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), function(i) J$met(J$draw()),
                          future.seed = 1L, future.packages = PKG))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  proj <- proj + dt * NBOOT / nb
  cat(sprintf("  Frame %s %-7s %d replicates in %.1f s\n", J$frame,
              J$cohort, nb, dt))
}
cat(sprintf("  projected total: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

row_of <- function(point, B, a, b = NULL) {
  est <- if (is.null(b)) point[[a]] else point[[a]] - point[[b]]
  v <- if (is.null(b)) B[, a] else B[, a] - B[, b]
  q <- quantile(v, c(0.025, 0.975), na.rm = TRUE)
  data.table(estimate = est, lo = q[[1]], hi = q[[2]])
}

S1 <- list()
S3 <- list()
SW <- list()
for (J in jobs) {
  cat(sprintf("\n=== Frame %s %s: %d bootstrap samples ===\n",
              J$frame, J$cohort, NBOOT))
  point <- J$met(seq_len(J$n))
  B <- do.call(rbind, future_lapply(
    seq_len(NBOOT), function(i) J$met(J$draw()),
    future.seed = SEED, future.packages = PKG))
  if (J$frame == "A") {
    for (en in names(J$SC)) {
      nms <- names(J$SC[[en]])
      cons <- list(
        "fused: binary - Cox (AUROC weights)" =
          c("fused lr auc", "fused cox auc"),
        "fused: binary - Cox (Uno C weights)" =
          c("fused lr uno", "fused cox uno"),
        "fused Cox: Uno C - AUROC weights" =
          c("fused cox uno", "fused cox auc"),
        "fused binary: Uno C - AUROC weights" =
          c("fused lr uno", "fused lr auc"),
        "fusion gain, Cox (fused - EHR)" =
          c("fused cox auc", "EHR cox"),
        "fusion gain, binary (fused - EHR)" =
          c("fused lr auc", "EHR lr"))
      for (m in J$mods) {
        cons[[sprintf("%s: binary - Cox", toupper(m))]] <-
          c(sprintf("%s lr", toupper(m)), sprintf("%s cox", toupper(m)))
      }
      for (ms in c("auc30", "unoC")) {
        for (nm in nms) {
          k <- paste(en, ms, nm, sep = "|")
          S1[[length(S1) + 1]] <- cbind(data.table(
            cohort = J$cohort, endpoint = en, events = J$events[[en]],
            measure = ms, type = "score", name = nm),
            row_of(point, B, k))
        }
        for (cn in names(cons)) {
          a <- paste(en, ms, cons[[cn]][1], sep = "|")
          b <- paste(en, ms, cons[[cn]][2], sep = "|")
          S1[[length(S1) + 1]] <- cbind(data.table(
            cohort = J$cohort, endpoint = en, events = J$events[[en]],
            measure = ms, type = "contrast", name = cn),
            row_of(point, B, a, b))
        }
      }
    }
  } else {
    M <- toupper(J$mods)
    for (h in HB) {
      scn <- c("fused cox", "fused lr", paste(M, "cox"), paste(M, "lr"))
      for (nm in scn) {
        S3[[length(S3) + 1]] <- cbind(data.table(
          cohort = J$cohort, horizon = h,
          events = J$events[[as.character(h)]], type = "score",
          name = nm), row_of(point, B, paste(h, "auc", nm, sep = "|")))
      }
      cons <- list("fused: Cox - binary" = c("fused cox", "fused lr"),
                   "fusion gain, Cox (fused - EHR)" =
                     c("fused cox", "EHR cox"),
                   "fusion gain, binary (fused - EHR)" =
                     c("fused lr", "EHR lr"))
      for (m in M) {
        cons[[sprintf("%s: Cox - binary", m)]] <-
          c(paste(m, "cox"), paste(m, "lr"))
      }
      for (cn in names(cons)) {
        S3[[length(S3) + 1]] <- cbind(data.table(
          cohort = J$cohort, horizon = h,
          events = J$events[[as.character(h)]], type = "contrast",
          name = cn),
          row_of(point, B, paste(h, "auc", cons[[cn]][1], sep = "|"),
                 paste(h, "auc", cons[[cn]][2], sep = "|")))
      }
      for (m in M) {
        SW[[length(SW) + 1]] <- cbind(data.table(
          cohort = J$cohort, type = "weight", horizon = h,
          name = sprintf("%s weight", m)),
          row_of(point, B, paste(h, "weight", m, sep = "|")))
      }
    }
    wm <- c("ECG", if ("CTPA" %in% M) "CTPA")
    for (m in wm) for (h in c(90, 180, 365)) {
      SW[[length(SW) + 1]] <- cbind(data.table(
        cohort = J$cohort, type = "contrast", horizon = h,
        name = sprintf("%s weight: 10 d - %d d", m, h)),
        row_of(point, B, paste(10, "weight", m, sep = "|"),
               paste(h, "weight", m, sep = "|")))
    }
  }
}
plan(sequential)

flag <- function(d) {
  d[, excludes_zero := ifelse(type == "contrast",
                              lo > 0 | hi < 0, NA)]
  d
}
s1 <- flag(rbindlist(S1))
s3 <- flag(rbindlist(S3))
sw <- flag(rbindlist(SW))
fwrite(s1, file.path(RS, "eval_s1.csv"))
fwrite(s3, file.path(RS, "eval_s3.csv"))
fwrite(sw, file.path(RS, "eval_s3_weights.csv"))

cat("\n=== KEY CONTRASTS ===\n")
print(s1[type == "contrast" & measure == "auc30" &
           grepl("^fused", name),
         .(cohort, endpoint, name, estimate = round(estimate, 4),
           lo = round(lo, 4), hi = round(hi, 4), excludes_zero)],
      row.names = FALSE)
print(sw[type == "contrast",
         .(cohort, name, estimate = round(estimate, 3),
           lo = round(lo, 3), hi = round(hi, 3), excludes_zero)],
      row.names = FALSE)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("eval_s1.csv", "eval_s3.csv", "eval_s3_weights.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "eval_s1_s3_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: eval_s1.csv, eval_s3.csv and eval_s3_weights.csv\n")