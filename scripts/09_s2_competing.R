# 09_s2_competing.R
# S2: CV readmission as a competing risk (PREREGISTRATION.md, section 6;
# deviation D3).
#
# Four ways to predict the risk of a CV readmission as the first event by
# horizon h. Each is a level-2 model on the out-of-fold level-1 scores of
# every modality (stage 6a), fitted inside each outer training fold and
# applied to the held-out fold (same folds and seeds):
#   F1  conditional binary (the dissertation's definition,
#       pe-multimodal-prognosis-msc): logistic model on the conditional
#       twins (cvc_h), trained without patients who died first by h.
#   F2  unconditional binary: logistic model on the unconditional twins
#       (cvu_h), trained on everyone; a death first counts as no event.
#   F3  cause-specific Cox (riskRegression::CSC): one Cox model for CV
#       readmission on the cs_cv scores and one for death first on the
#       cs_death scores, combined into the cumulative incidence by h.
#   F4  Fine-Gray (riskRegression::FGR) on the same scores.
# Under D3 every level-2 input is winsorised at its training-fold 0.5th
# and 99.5th percentiles, then standardised, identically for all four.
# Horizons are 10 and 30 days from admission in Frame A, and 10, 30, 90,
# 180 and 365 days from discharge in Frame B, in both cohorts.
#
# Deployment evaluation uses the full cohort (a patient who dies first
# counts as no CV readmission): AUROC, Brier score, observed/expected,
# calibration slope and ICI. Conditional evaluation (the dissertation's
# question) is the AUROC among patients not dying first by h. Follow-up is
# complete, so no censoring weights are needed. Contrasts F2-F1, F3-F1,
# F4-F1 and F4-F3 have 95% percentile intervals from 2,000 bootstrap
# samples (seed 42; resampled by patient in Frame A). The per-fold AUROC
# (D2) is a sensitivity. F1 is expected to be worst on the deployment
# estimand (registered).
#
# Checks: no failed fits, every patient predicted, risks within [0, 1],
# and CSC and Fine-Gray risks never falling with the horizon. A failure
# stops the script before evaluation.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv
#          data/processed/level1/oof_linear_*.csv
# Outputs: results/s2_metrics.csv, results/s2_contrasts.csv,
#          results/s2_checks.csv, results/s2_log.txt
#          data/processed/s2/s2_preds_*.csv (patient level, not in repo)
# Usage:   Rscript scripts/09_s2_competing.R [--probe-only]
# Notes:   Parallel (future.apply). --probe-only stops after the fit probe.

suppressPackageStartupMessages({
  library(survival)
  library(prodlim)
  library(riskRegression)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
S2D <- file.path(PR, "s2")
dir.create(S2D, showWarnings = FALSE, recursive = TRUE)
NBOOT <- 2000L
SEED <- 42L
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
FL <- c(F1 = "conditional binary (DISSERTATION)",
        F2 = "unconditional binary",
        F3 = "cause-specific Cox",
        F4 = "Fine-Gray")
PKG <- c("survival", "prodlim", "riskRegression", "cmprsk",
         "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

# ---------------- level-2 fitting (one task = frame, cohort, seed) ----
fit_stacks <- function(tk, folds_to_run = 0:4) {
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  wz <- function(v, rows) {
    qv <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, qv[1]), qv[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
  n <- length(tk$fold)
  hs <- tk$hs
  mods <- tk$mods
  X <- tk$X
  tf <- tk$tf
  er <- tk$er
  out <- list()
  for (h in hs) for (f in c("F1", "F2", "F3", "F4")) {
    out[[sprintf("%s_h%d", f, h)]] <- rep(NA_real_, n)
  }
  nfail <- 0L
  maxcov <- 0
  for (k in folds_to_run) {
    tr <- tk$fold != k
    te <- tk$fold == k
    cvL <- sapply(mods, function(m)
      wz(X[[paste0(m, ".cox.cs_cv.lp")]], tr))
    dL <- sapply(mods, function(m)
      wz(X[[paste0(m, ".cox.cs_death.lp")]], tr))
    colnames(cvL) <- paste0("cv_", mods)
    colnames(dL) <- paste0("d_", mods)
    maxcov <- max(maxcov, abs(cvL), abs(dL))
    D <- data.frame(time = tf, event = er, cvL, dL)
    fcv <- paste(colnames(cvL), collapse = " + ")
    fd <- paste(colnames(dL), collapse = " + ")
    f_cv <- as.formula(paste("Hist(time, event) ~", fcv))
    f_d <- as.formula(paste("Hist(time, event) ~", fd))
    f_fg <- as.formula(paste("Hist(time, event) ~", fcv, "+", fd))
    r3 <- tryCatch({
      m3 <- CSC(list(f_cv, f_d), data = D[tr, ])
      predictRisk(m3, newdata = D[te, ], times = hs + 0.5, cause = 1)
    }, error = function(e) NULL)
    r4 <- tryCatch({
      m4 <- tryCatch(
        FGR(f_fg, data = D[tr, ], cause = 1, variance = FALSE),
        error = function(e) FGR(f_fg, data = D[tr, ], cause = 1))
      predictRisk(m4, newdata = D[te, ], times = hs + 0.5)
    }, error = function(e) NULL)
    nfail <- nfail + is.null(r3) + is.null(r4)
    for (j in seq_along(hs)) {
      h <- hs[j]
      if (!is.null(r3)) out[[sprintf("F3_h%d", h)]][te] <- r3[, j]
      if (!is.null(r4)) out[[sprintf("F4_h%d", h)]][te] <- r4[, j]
      y <- as.integer(er == 1 & tf <= h + 0.5)
      dfirst <- er == 2 & tf <= h + 0.5
      trc <- tr & !dfirst
      L1 <- as.data.frame(sapply(mods, function(m)
        wz(qlogis(clip(X[[sprintf("%s.lr.cvc_h%d", m, h)]])), trc)))
      L2 <- as.data.frame(sapply(mods, function(m)
        wz(qlogis(clip(X[[sprintf("%s.lr.cvu_h%d", m, h)]])), tr)))
      maxcov <- max(maxcov, abs(as.matrix(L1)), abs(as.matrix(L2)))
      g1 <- tryCatch(glm(y ~ ., family = binomial,
                         data = cbind(y = y[trc],
                                      L1[trc, , drop = FALSE])),
                     error = function(e) NULL)
      g2 <- tryCatch(glm(y ~ ., family = binomial,
                         data = cbind(y = y[tr],
                                      L2[tr, , drop = FALSE])),
                     error = function(e) NULL)
      nfail <- nfail + is.null(g1) + is.null(g2)
      if (!is.null(g1)) {
        out[[sprintf("F1_h%d", h)]][te] <-
          predict(g1, newdata = L1[te, , drop = FALSE], type = "response")
      }
      if (!is.null(g2)) {
        out[[sprintf("F2_h%d", h)]][te] <-
          predict(g2, newdata = L2[te, , drop = FALSE], type = "response")
      }
    }
  }
  list(fr = tk$fr, co = tk$co, seed = tk$seed,
       pred = as.data.table(out), nfail = nfail, maxcov = maxcov)
}

# ---------------- evaluation helpers ----------------
met_job <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (h in J$hs) {
    hk <- as.character(h)
    y <- J$y[[hk]][idx]
    keep <- !J$df[[hk]][idx]
    for (f in J$forms) {
      p <- J$P[[sprintf("%s_h%d", f, h)]][idx]
      v[paste(h, f, "auc_full", sep = "|")] <- auc(y, p)
      v[paste(h, f, "auc_cond", sep = "|")] <- auc(y[keep], p[keep])
      v[paste(h, f, "brier", sep = "|")] <- mean((y - p)^2)
      v[paste(h, f, "oe", sep = "|")] <- mean(y) / mean(p)
    }
  }
  v
}
draw_idx <- function(J) {
  if (J$fr == "A") {
    unlist(J$rows_by[sample.int(length(J$rows_by),
                                length(J$rows_by), replace = TRUE)],
           use.names = FALSE)
  } else {
    sample.int(J$n, J$n, replace = TRUE)
  }
}
boot_rep <- function(i, J) met_job(J, draw_idx(J))
auc1 <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ici_fun <- function(y, p) {
  f <- tryCatch(loess(y ~ p, span = 0.75, degree = 2),
                error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  mean(abs(fitted(f) - p))
}

# ---------------- build tasks ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  hs <- HZ[[fr]]
  cols <- c(paste0(mods, ".cox.cs_cv.lp"),
            paste0(mods, ".cox.cs_death.lp"),
            as.vector(outer(mods, hs, function(m, h)
              sprintf("%s.lr.cvc_h%d", m, h))),
            as.vector(outer(mods, hs, function(m, h)
              sprintf("%s.lr.cvu_h%d", m, h))))
  oof <- fread(file.path(PR, "level1",
                         sprintf("oof_linear_%s_%s.csv", fr, code)),
               select = c("seed", "hadm_id", "fold", cols))
  seeds <- sort(unique(oof$seed))
  key <- paste(fr, co)
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub,
                      mods = mods, hs = hs, seeds = seeds,
                      folds = list())
  for (s in seeds) {
    g <- oof[oof$seed == s][match(sub$hadm_id, hadm_id)]
    if (anyNA(g$fold) || anyNA(g[, ..cols])) {
      stop("missing level-1 score: ", key, " seed ", s)
    }
    META[[key]]$folds[[as.character(s)]] <- g$fold
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = g$fold,
      X = as.list(g[, ..cols]), tf = sub$time_first,
      er = sub$event_cr, hs = hs, mods = mods)
  }
  cat(sprintf("  Frame %s %-7s n=%d, CV readmissions first by %s: %s\n",
              fr, co, nrow(sub), paste(hs, collapse = "/"),
              paste(sapply(hs, function(h)
                sum(sub$event_cr == 1 & sub$time_first <= h + 0.5)),
                collapse = "/")))
}

# ---------------- timing probe: fitting ----------------
cat(sprintf("\n=== TIMING PROBE: fold 0 of %d tasks in parallel ===\n",
            workers))
plan(multisession, workers = workers)
ord <- order(-sapply(tasks, function(t) length(t$fold)))
heavy <- tasks[ord[seq_len(min(workers, length(tasks)))]]
t0 <- Sys.time()
invisible(future_lapply(heavy, fit_stacks, folds_to_run = 0L,
                        future.packages = PKG, future.seed = TRUE))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
fit_proj <- dt * 5 * ceiling(length(tasks) / workers)
cat(sprintf("  %d fold-0 fits in parallel: %.1f s; fitting all %d",
            length(heavy), dt, length(tasks)),
    sprintf("tasks about %.1f min\n", fit_proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- fit everything ----------------
cat("\n=== FITTING ===\n")
res <- future_lapply(tasks[ord], fit_stacks, future.packages = PKG,
                     future.seed = TRUE)
CHK <- list()
chk <- function(fr, co, item, value, ok) {
  st <- if (is.na(ok)) "INFO" else if (ok) "PASS" else "FAIL"
  CHK[[length(CHK) + 1]] <<- data.table(frame = fr, cohort = co,
                                        item = item,
                                        value = as.character(value),
                                        status = st)
  cat(sprintf("  %-4s %s %-7s %-46s %s\n", st, fr, co, item, value))
}
JOBS <- list()
for (key in names(META)) {
  M <- META[[key]]
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  rr <- rr[order(sapply(rr, function(r) r$seed))]
  pcols <- names(rr[[1]]$pred)
  nf <- sum(sapply(rr, function(r) r$nfail))
  chk(M$fr, M$co, "failed fits", nf, nf == 0)
  chk(M$fr, M$co, "largest level-2 input after winsorising",
      sprintf("%.2f", max(sapply(rr, function(r) r$maxcov))), NA)
  long <- rbindlist(lapply(rr, function(r) {
    cbind(data.table(seed = r$seed, hadm_id = M$sub$hadm_id,
                     fold = M$folds[[as.character(r$seed)]]), r$pred)
  }))
  fwrite(long, file.path(S2D, sprintf("s2_preds_%s_%s.csv",
                                      M$fr, M$code)))
  avg <- long[, lapply(.SD, mean), by = hadm_id, .SDcols = pcols]
  avg <- avg[match(M$sub$hadm_id, hadm_id)]
  nna <- sum(is.na(as.matrix(avg[, ..pcols])))
  chk(M$fr, M$co, "missing predictions", nna, nna == 0)
  rng <- range(as.matrix(long[, ..pcols]), na.rm = TRUE)
  chk(M$fr, M$co, "risk range (all seeds)",
      sprintf("%.6f to %.6f", rng[1], rng[2]),
      rng[1] >= -1e-9 && rng[2] <= 1 + 1e-9)
  for (f in c("F3", "F4")) {
    if (length(M$hs) < 2) next
    mm <- as.matrix(long[, sprintf("%s_h%d", f, M$hs), with = FALSE])
    viol <- sum(apply(mm, 1, function(r) any(diff(r) < -1e-9)))
    chk(M$fr, M$co, sprintf("%s risk falls as horizon lengthens", f),
        viol, viol == 0)
  }
  sub <- M$sub
  JOBS[[key]] <- list(
    fr = M$fr, co = M$co, n = nrow(sub), hs = M$hs,
    forms = names(FL),
    y = setNames(lapply(M$hs, function(h)
      sub$event_cr == 1 & sub$time_first <= h + 0.5), M$hs),
    df = setNames(lapply(M$hs, function(h)
      sub$event_cr == 2 & sub$time_first <= h + 0.5), M$hs),
    P = as.list(avg[, ..pcols]),
    rows_by = split(seq_len(nrow(sub)),
                    match(sub$subject_id, unique(sub$subject_id))),
    long = long)
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "s2_checks.csv"))
if (any(ck$status == "FAIL")) {
  plan(sequential)
  cat("\nHARD CHECK FAILED: stopping before evaluation.\n")
  cat("Outputs: s2_checks.csv\n")
  quit(status = 1)
}

# ---------------- timing probe: bootstrap ----------------
cat(sprintf("\n=== TIMING PROBE: bootstrap (%d workers) ===\n", workers))
nb <- workers * 2L
bproj <- 0
for (J in JOBS) {
  Jl <- J[setdiff(names(J), "long")]
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), boot_rep, J = Jl,
                          future.seed = 1L, future.packages = PKG))
  d_ <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  bproj <- bproj + d_ * NBOOT / nb
}
cat(sprintf("  projected bootstrap: about %.1f min\n", bproj / 60))
Sys.sleep(5)

# ---------------- bootstrap and summaries ----------------
MET <- list()
CON <- list()
pairs <- list(c("F2", "F1"), c("F3", "F1"), c("F4", "F1"),
              c("F4", "F3"))
for (J in JOBS) {
  cat(sprintf("\n=== Frame %s %s: %d bootstrap samples ===\n",
              J$fr, J$co, NBOOT))
  Jl <- J[setdiff(names(J), "long")]
  point <- met_job(Jl, seq_len(J$n))
  B <- do.call(rbind, future_lapply(
    seq_len(NBOOT), boot_rep, J = Jl, future.seed = SEED,
    future.packages = PKG))
  qq <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
  for (h in J$hs) {
    hk <- as.character(h)
    y <- J$y[[hk]]
    for (f in J$forms) {
      pc <- sprintf("%s_h%d", f, h)
      p <- J$P[[pc]]
      pf <- c()
      for (s in unique(J$long$seed)) {
        ls <- J$long[seed == s]
        for (k in 0:4) {
          ii <- which(ls$fold == k)
          pf <- c(pf, auc1(y[ii], ls[[pc]][ii]))
        }
      }
      slope <- tryCatch(unname(coef(glm(
        y ~ qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)),
        family = binomial))[2]), error = function(e) NA_real_)
      row <- data.table(
        frame = J$fr, cohort = J$co, horizon = h, formulation = f,
        label = FL[[f]], n = J$n, cv_events = sum(y),
        deaths_first = sum(J$df[[hk]]), observed_risk = mean(y),
        mean_predicted = mean(p), cal_slope = slope,
        ici = ici_fun(as.numeric(y), p),
        auc_full_perfold = mean(pf, na.rm = TRUE))
      for (ms in c("auc_full", "auc_cond", "brier", "oe")) {
        k_ <- paste(h, f, ms, sep = "|")
        q <- qq(B[, k_])
        row[[ms]] <- point[[k_]]
        row[[paste0(ms, "_lo")]] <- q[[1]]
        row[[paste0(ms, "_hi")]] <- q[[2]]
      }
      MET[[length(MET) + 1]] <- row
    }
    for (pr in pairs) for (ms in c("auc_full", "auc_cond", "brier")) {
      a <- paste(h, pr[1], ms, sep = "|")
      b <- paste(h, pr[2], ms, sep = "|")
      q <- qq(B[, a] - B[, b])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, horizon = h, measure = ms,
        contrast = sprintf("%s - %s", pr[1], pr[2]),
        label = sprintf("%s minus %s", FL[[pr[1]]], FL[[pr[2]]]),
        estimate = point[[a]] - point[[b]], lo = q[[1]], hi = q[[2]],
        excludes_zero = q[[1]] > 0 | q[[2]] < 0,
        better_if = ifelse(ms == "brier", "negative", "positive"))
    }
  }
}
plan(sequential)

mt <- rbindlist(MET, fill = TRUE)
ct <- rbindlist(CON)
fwrite(mt, file.path(RS, "s2_metrics.csv"))
fwrite(ct, file.path(RS, "s2_contrasts.csv"))

cat("\n=== KEY RESULTS (Frame A 30 d, Frame B 365 d) ===\n")
print(mt[(frame == "A" & horizon == 30) | (frame == "B" & horizon == 365),
         .(frame, cohort, horizon, formulation,
           auc_full = round(auc_full, 4), auc_cond = round(auc_cond, 4),
           brier = round(brier, 4), oe = round(oe, 3),
           slope = round(cal_slope, 3), ici = round(ici, 4))],
      row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d FAIL\n", sum(ck$status == "PASS"),
            sum(ck$status == "FAIL")))
cat(sprintf("  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("s2_metrics.csv", "s2_contrasts.csv", "s2_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", S2D, " (patient-level predictions, local only)\n")
cat("  ", file.path(RS, "s2_log.txt"), " (terminal log)\n")
cat("\nOutputs: s2_metrics.csv, s2_contrasts.csv and s2_checks.csv\n")