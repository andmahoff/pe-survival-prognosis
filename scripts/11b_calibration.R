# 11b_calibration.R
# Calibration and overall accuracy of the death and composite predictions
# (PREREGISTRATION.md section 7; DEVIATIONS.md D3, D4).
#
# Horizons are 30 days in Frame A (admission clock) and 30 and 365 days in
# Frame B (discharge clock), in both cohorts. CV readmission was calibrated
# in S2 (stage 9). All models are out of fold and averaged over seeds. The
# Cox stack is a level-2 Cox model on the modalities' winsorised ridge-Cox
# linear predictors (D3), fitted in each outer training fold, with absolute
# risk 1 - exp(-H0(h + 0.5) * exp(lp)) from the Breslow baseline hazard. The
# binary stack is a level-2 logistic model on the winsorised logits of the
# modalities' horizon-specific logistic twins (D3). Late fusion is the stage
# 8 fused rank score (Cox level 1, weights by AUROC), recalibrated in fold
# by logistic regression on its logit, since ranks are not risks
# (pre-registered). Single modalities use their level-1 ridge-Cox absolute
# risk and logistic twin as they come (point estimates only).
#
# Follow-up is complete, so the outcome at h is whether the event happened
# by h. Measures are observed/expected, calibration intercept and slope,
# ICI, E50 and E90 from a loess curve (span 0.75, degree 2), Brier score,
# IPA (reference: the observed event proportion, D4) and AUROC, plus decile
# calibration tables for the figures. Intervals for the fused models come
# from 2,000 bootstrap samples (seed 42), clustered by patient in Frame A
# (D4), including paired contrasts (Cox stack minus binary stack, Cox stack
# minus late fusion) for Brier, IPA and ICI.
#
# Checks: every patient predicted, every risk in [0, 1], and the Cox-stack
# risk never falls from 30 to 365 days. The script stops before the
# bootstrap if one fails.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv,
#          data/processed/level1/oof_linear_*.csv,
#          data/processed/fusion/fused_*.csv
# Outputs: results/cal_metrics.csv, results/cal_contrasts.csv,
#          results/cal_curves.csv, results/cal_checks.csv,
#          results/cal_log.txt (terminal log),
#          data/processed/calibration/cal_preds_*.csv (patient level, not
#          in the repo)
# Usage:   Rscript scripts/11b_calibration.R [--probe-only]
# Notes:   --probe-only stops after the timing probe. Fitting and the
#          bootstrap run in parallel (future.apply, physical cores - 1).

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
CD <- file.path(PR, "calibration")
dir.create(CD, showWarnings = FALSE, recursive = TRUE)
NBOOT <- 2000L
SEED <- 42L
HZ <- list(A = 30, B = c(30, 365))
ENDS <- c("death", "composite")
FUSED <- c("coxstack", "binstack", "latefusion")
FLAB <- c(coxstack = "Cox stack", binstack = "binary stack",
          latefusion = "late fusion, recalibrated")
PKG <- c("survival", "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

# ---------------- level-2 fitting (one task = frame, cohort, seed) ----
fit_task <- function(tk, folds_to_run = 0:4) {
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
  n <- length(tk$fold)
  out <- list()
  for (en in names(tk$E)) for (h in tk$hs) {
    out[[sprintf("%s|coxstack|%d", en, h)]] <- rep(NA_real_, n)
    out[[sprintf("%s|binstack|%d", en, h)]] <- rep(NA_real_, n)
    if (!is.null(tk$F[[paste(en, h)]])) {
      out[[sprintf("%s|latefusion|%d", en, h)]] <- rep(NA_real_, n)
    }
  }
  nfail <- 0L
  for (k in folds_to_run) {
    tr <- tk$fold != k
    te <- tk$fold == k
    for (en in names(tk$E)) {
      tt <- tk$E[[en]]$t
      ee <- as.integer(tk$E[[en]]$e)
      Z <- sapply(tk$mods, function(m)
        wz(tk$X[[sprintf("%s.cox.%s.lp", m, en)]], tr))
      D <- data.frame(tt = tt, ee = ee, Z)
      fit <- tryCatch(coxph(Surv(tt, ee) ~ ., data = D[tr, ],
                            ties = "efron"),
                      error = function(e) NULL)
      if (is.null(fit)) {
        nfail <- nfail + 1L
      } else {
        bh <- basehaz(fit, centered = FALSE)
        lp <- as.vector(as.matrix(D[te, tk$mods, drop = FALSE]) %*%
                          coef(fit)[tk$mods])
        for (h in tk$hs) {
          i <- findInterval(h + 0.5, bh$time)
          H0 <- if (i == 0) 0 else bh$hazard[i]
          out[[sprintf("%s|coxstack|%d", en, h)]][te] <-
            1 - exp(-H0 * exp(lp))
        }
      }
      for (h in tk$hs) {
        y <- as.integer(ee == 1 & tt <= h + 0.5)
        L <- as.data.frame(sapply(tk$mods, function(m)
          wz(qlogis(clip(tk$X[[sprintf("%s.lr.%s_h%d", m, en, h)]])),
             tr)))
        g <- tryCatch(glm(y ~ ., family = binomial,
                          data = cbind(y = y[tr],
                                       L[tr, , drop = FALSE])),
                      error = function(e) NULL)
        if (is.null(g)) {
          nfail <- nfail + 1L
        } else {
          out[[sprintf("%s|binstack|%d", en, h)]][te] <-
            predict(g, newdata = L[te, , drop = FALSE],
                    type = "response")
        }
        fc <- tk$F[[paste(en, h)]]
        if (!is.null(fc)) {
          dd <- data.frame(y = y, s = qlogis(clip(tk$X[[fc]])))
          g2 <- tryCatch(glm(y ~ s, family = binomial, data = dd[tr, ]),
                         error = function(e) NULL)
          if (is.null(g2)) {
            nfail <- nfail + 1L
          } else {
            out[[sprintf("%s|latefusion|%d", en, h)]][te] <-
              predict(g2, newdata = dd[te, ], type = "response")
          }
        }
      }
    }
  }
  list(fr = tk$fr, co = tk$co, seed = tk$seed,
       pred = as.data.table(out), nfail = nfail)
}

# ---------------- measures ----------------
auc1 <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cal_point <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  lg <- qlogis(p)
  prev <- mean(y)
  sl <- tryCatch(suppressWarnings(unname(glm.fit(
    cbind(1, lg), y, family = binomial())$coefficients[2])),
    error = function(e) NA_real_)
  it <- tryCatch(suppressWarnings(unname(glm.fit(
    matrix(1, length(y), 1), y, family = binomial(),
    offset = lg)$coefficients[1])), error = function(e) NA_real_)
  f <- tryCatch(fitted(loess(y ~ p, span = 0.75, degree = 2)),
                error = function(e) rep(NA_real_, length(p)))
  d <- abs(f - p)
  br <- mean((y - p)^2)
  c(oe = prev / mean(p), cal_intercept = it, cal_slope = sl,
    ici = mean(d), e50 = median(d),
    e90 = unname(quantile(d, 0.9, na.rm = TRUE)), brier = br,
    ipa = 1 - br / (prev * (1 - prev)), auc = auc1(y == 1, p))
}
boot_rep <- function(i, J) {
  if (J$fr == "A") {
    idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                       length(J$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(J$n, J$n, replace = TRUE)
  }
  v <- numeric(0)
  for (s in J$series) {
    y <- J$Y[[s$y]][idx]
    p <- pmin(pmax(J$P[[s$p]][idx], 1e-6), 1 - 1e-6)
    prev <- mean(y)
    sl <- tryCatch(suppressWarnings(unname(glm.fit(
      cbind(1, qlogis(p)), y, family = binomial())$coefficients[2])),
      error = function(e) NA_real_)
    f <- tryCatch(fitted(loess(y ~ p, span = 0.75, degree = 2)),
                  error = function(e) rep(NA_real_, length(p)))
    br <- mean((y - p)^2)
    v[paste(s$key, c("oe", "cal_slope", "ici", "brier", "ipa"),
            sep = "|")] <-
      c(prev / mean(p), sl, mean(abs(f - p)), br,
        1 - br / (prev * (1 - prev)))
  }
  v
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
  E <- list(death = list(t = sub$time_death, e = sub$event_death == 1),
            composite = list(t = sub$time_first,
                             e = sub$event_comp == 1))
  Fmap <- list()
  for (en in ENDS) for (h in hs) {
    fc <- if (fr == "A") sprintf("fuse.cox.auc.%s", en) else
      if (en == "death") sprintf("fuse.cox.auc%d.death", h) else NULL
    if (!is.null(fc)) Fmap[[paste(en, h)]] <- fc
  }
  fcols <- unique(unlist(Fmap))
  lcols <- c(
    as.vector(outer(mods, ENDS, function(m, en)
      sprintf("%s.cox.%s.lp", m, en))),
    unlist(lapply(hs, function(h) as.vector(outer(mods, ENDS,
      function(m, en) sprintf("%s.lr.%s_h%d", m, en, h))))),
    unlist(lapply(hs, function(h) as.vector(outer(mods, ENDS,
      function(m, en) sprintf("%s.cox.%s.r%d", m, en, h))))))
  p1 <- file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv", fr, code))
  p2 <- file.path(PR, "fusion", sprintf("fused_%s_%s.csv", fr, code))
  h1 <- names(fread(p1, nrows = 0))
  h2 <- names(fread(p2, nrows = 0))
  miss <- c(setdiff(c("seed", "hadm_id", "fold", lcols), h1),
            setdiff(c("seed", "hadm_id", fcols), h2))
  if (length(miss)) stop("columns not found: ", paste(miss, collapse = ", "))
  oof <- fread(p1, select = c("seed", "hadm_id", "fold", lcols))
  fu <- fread(p2, select = c("seed", "hadm_id", fcols))
  seeds <- sort(unique(oof$seed))
  key <- paste(fr, co)
  folds <- list()
  for (s in seeds) {
    g <- oof[oof$seed == s][match(sub$hadm_id, hadm_id)]
    gf <- fu[fu$seed == s][match(sub$hadm_id, hadm_id)]
    if (anyNA(g$fold) || anyNA(g[, ..lcols]) || anyNA(gf[, ..fcols])) {
      stop("missing score: ", key, " seed ", s)
    }
    folds[[as.character(s)]] <- g$fold
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = g$fold,
      X = c(as.list(g[, ..lcols]), as.list(gf[, ..fcols])),
      E = E, hs = hs, mods = mods, F = Fmap)
  }
  avg1 <- oof[, lapply(.SD, mean), by = hadm_id, .SDcols = lcols]
  avg1 <- avg1[match(sub$hadm_id, hadm_id)]
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub,
                      mods = mods, hs = hs, E = E, Fmap = Fmap,
                      folds = folds, avg1 = avg1)
  cat(sprintf("  Frame %s %-7s n=%d; deaths / composite by %s: %s\n",
              fr, co, nrow(sub), paste(hs, collapse = "/"),
              paste(sapply(hs, function(h) sprintf("%d/%d",
                sum(E$death$e & E$death$t <= h + 0.5),
                sum(E$composite$e & E$composite$t <= h + 0.5))),
                collapse = ", ")))
}

# ---------------- timing probe: fitting ----------------
cat(sprintf("\n=== TIMING PROBE: fold 0 of %d tasks in parallel ===\n",
            workers))
plan(multisession, workers = workers)
ord <- order(-sapply(tasks, function(t) length(t$fold)))
heavy <- tasks[ord[seq_len(min(workers, length(tasks)))]]
t0 <- Sys.time()
invisible(future_lapply(heavy, fit_task, folds_to_run = 0L,
                        future.packages = PKG, future.seed = TRUE))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  %d fold-0 fits in parallel: %.1f s; fitting all %d",
            length(heavy), dt, length(tasks)),
    sprintf("tasks about %.1f min\n",
            dt * 5 * ceiling(length(tasks) / workers) / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- fit everything, hard checks ----------------
cat("\n=== FITTING ===\n")
res <- future_lapply(tasks[ord], fit_task, future.packages = PKG,
                     future.seed = TRUE)
CHK <- list()
chk <- function(fr, co, item, value, ok) {
  st <- if (is.na(ok)) "INFO" else if (ok) "PASS" else "FAIL"
  CHK[[length(CHK) + 1]] <<- data.table(frame = fr, cohort = co,
                                        item = item,
                                        value = as.character(value),
                                        status = st)
  cat(sprintf("  %-4s %s %-7s %-44s %s\n", st, fr, co, item, value))
}
JOBS <- list()
for (key in names(META)) {
  M <- META[[key]]
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  rr <- rr[order(sapply(rr, function(r) r$seed))]
  pcols <- names(rr[[1]]$pred)
  nf <- sum(sapply(rr, function(r) r$nfail))
  chk(M$fr, M$co, "failed fits", nf, nf == 0)
  long <- rbindlist(lapply(rr, function(r) {
    cbind(data.table(seed = r$seed, hadm_id = M$sub$hadm_id,
                     fold = M$folds[[as.character(r$seed)]]), r$pred)
  }))
  fwrite(long, file.path(CD, sprintf("cal_preds_%s_%s.csv",
                                     M$fr, M$code)))
  avg <- long[, lapply(.SD, mean), by = hadm_id, .SDcols = pcols]
  avg <- avg[match(M$sub$hadm_id, hadm_id)]
  nna <- sum(is.na(as.matrix(avg[, ..pcols])))
  chk(M$fr, M$co, "missing predictions", nna, nna == 0)
  rng <- range(as.matrix(long[, ..pcols]), na.rm = TRUE)
  chk(M$fr, M$co, "risk range (all seeds)",
      sprintf("%.6f to %.6f", rng[1], rng[2]),
      rng[1] >= -1e-9 && rng[2] <= 1 + 1e-9)
  if (length(M$hs) > 1) {
    for (en in ENDS) {
      a <- long[[sprintf("%s|coxstack|%d", en, min(M$hs))]]
      b <- long[[sprintf("%s|coxstack|%d", en, max(M$hs))]]
      viol <- sum(b < a - 1e-9, na.rm = TRUE)
      chk(M$fr, M$co, sprintf("%s Cox-stack risk falls 30 -> 365 d", en),
          viol, viol == 0)
    }
  }
  P <- as.list(avg[, ..pcols])
  Y <- list()
  series <- list()
  for (en in ENDS) for (h in M$hs) {
    Y[[paste(en, h)]] <- as.numeric(M$E[[en]]$e &
                                      M$E[[en]]$t <= h + 0.5)
    for (m in M$mods) {
      P[[sprintf("%s|%s Cox|%d", en, toupper(m), h)]] <-
        M$avg1[[sprintf("%s.cox.%s.r%d", m, en, h)]]
      P[[sprintf("%s|%s binary|%d", en, toupper(m), h)]] <-
        M$avg1[[sprintf("%s.lr.%s_h%d", m, en, h)]]
    }
    for (f in FUSED) {
      pk <- sprintf("%s|%s|%d", en, f, h)
      if (pk %in% names(P)) {
        series[[length(series) + 1]] <- list(
          key = sprintf("%s|%d|%s", en, h, f), y = paste(en, h), p = pk)
      }
    }
  }
  sub <- M$sub
  JOBS[[key]] <- list(
    fr = M$fr, co = M$co, n = nrow(sub), hs = M$hs, mods = M$mods,
    Y = Y, P = P, series = series,
    rows_by = split(seq_len(nrow(sub)),
                    match(sub$subject_id, unique(sub$subject_id))))
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "cal_checks.csv"))
if (any(ck$status == "FAIL")) {
  plan(sequential)
  cat("\nHARD CHECK FAILED: stopping before evaluation.\n")
  cat("Outputs: cal_checks.csv\n")
  quit(status = 1)
}

# ---------------- timing probe: bootstrap ----------------
cat(sprintf("\n=== TIMING PROBE: bootstrap (%d workers) ===\n", workers))
nb <- workers * 2L
bproj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), boot_rep, J = J,
                          future.seed = 1L, future.packages = PKG))
  d_ <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  bproj <- bproj + d_ * NBOOT / nb
}
cat(sprintf("  projected bootstrap: about %.1f min\n", bproj / 60))
Sys.sleep(5)

# ---------------- bootstrap and summaries ----------------
MET <- list()
CON <- list()
CUR <- list()
qq <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cat(sprintf("\n=== Frame %s %s: %d bootstrap samples ===\n",
              J$fr, J$co, NBOOT))
  B <- do.call(rbind, future_lapply(
    seq_len(NBOOT), boot_rep, J = J, future.seed = SEED,
    future.packages = PKG))
  for (en in ENDS) for (h in J$hs) {
    y <- J$Y[[paste(en, h)]]
    mk <- c(FUSED[sprintf("%s|%s|%d", en, FUSED, h) %in% names(J$P)],
            as.vector(outer(toupper(J$mods), c("Cox", "binary"), paste)))
    pts <- list()
    for (m in mk) {
      pk <- sprintf("%s|%s|%d", en, m, h)
      p <- J$P[[pk]]
      cp <- cal_point(y, p)
      pts[[m]] <- cp
      fused <- m %in% FUSED
      row <- cbind(data.table(
        frame = J$fr, cohort = J$co, endpoint = en, horizon = h,
        model = if (fused) FLAB[[m]] else paste(m, "(level 1)"),
        type = if (fused) "fused" else "single modality",
        n = J$n, events = sum(y), observed_risk = mean(y),
        mean_predicted = mean(p)), as.data.table(as.list(cp)))
      if (fused) {
        for (ms in c("oe", "cal_slope", "ici", "brier", "ipa")) {
          q <- qq(B[, sprintf("%s|%d|%s|%s", en, h, m, ms)])
          row[[paste0(ms, "_lo")]] <- q[1]
          row[[paste0(ms, "_hi")]] <- q[2]
        }
      }
      MET[[length(MET) + 1]] <- row
      g <- ceiling(10 * frank(p, ties.method = "first") / length(p))
      cu <- data.table(pp = p, yy = y, bin = g)[
        , .(mean_predicted = mean(pp), observed = mean(yy), n = .N),
        by = bin][order(bin)]
      CUR[[length(CUR) + 1]] <- cbind(data.table(
        frame = J$fr, cohort = J$co, endpoint = en, horizon = h,
        model = row$model), cu)
    }
    pairs <- list(c("coxstack", "binstack"), c("coxstack", "latefusion"))
    for (pr in pairs) {
      if (!all(pr %in% names(pts))) next
      for (ms in c("brier", "ipa", "ici")) {
        a <- B[, sprintf("%s|%d|%s|%s", en, h, pr[1], ms)]
        b <- B[, sprintf("%s|%d|%s|%s", en, h, pr[2], ms)]
        q <- qq(a - b)
        CON[[length(CON) + 1]] <- data.table(
          frame = J$fr, cohort = J$co, endpoint = en, horizon = h,
          measure = ms,
          contrast = sprintf("%s minus %s", FLAB[[pr[1]]],
                             FLAB[[pr[2]]]),
          estimate = pts[[pr[1]]][[ms]] - pts[[pr[2]]][[ms]],
          lo = q[1], hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0,
          better_if = ifelse(ms == "ipa", "positive", "negative"))
      }
    }
  }
}
plan(sequential)

mt <- rbindlist(MET, fill = TRUE)
ct <- rbindlist(CON)
cu <- rbindlist(CUR)
fwrite(mt, file.path(RS, "cal_metrics.csv"))
fwrite(ct, file.path(RS, "cal_contrasts.csv"))
fwrite(cu, file.path(RS, "cal_curves.csv"))

cat("\n=== FUSED MODELS ===\n")
print(mt[type == "fused",
         .(frame, cohort, endpoint, horizon, model,
           oe = round(oe, 3), slope = round(cal_slope, 3),
           ici = round(ici, 4), brier = round(brier, 4),
           ipa = round(ipa, 3), auc = round(auc, 4))],
      row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d FAIL\n", sum(ck$status == "PASS"),
            sum(ck$status == "FAIL")))
cat(sprintf("  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("cal_metrics.csv", "cal_contrasts.csv", "cal_curves.csv",
            "cal_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", CD, " (patient-level predictions, local only)\n")
cat("  ", file.path(RS, "cal_log.txt"), " (terminal log)\n")
cat("\nOutputs: cal_metrics.csv, cal_contrasts.csv, cal_curves.csv,",
    "cal_checks.csv\n")