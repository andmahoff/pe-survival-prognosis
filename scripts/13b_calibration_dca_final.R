# 13b_calibration_dca_final.R
# Calibration and decision curves of the final models, for Table 2
# (stage 13c) and the results figures (stage 13d; A38).
#
# Nothing is refitted: the script reads the final models' out-of-fold
# predictions, averaged over seeds, and writes aggregates only. At
# diagnosis (Frame P, time from admission) the models are Super learner P
# for death (stage 18q; A32, A36) and the stage 18r super learners for the
# composite and CV readmission (A33, A34). At discharge (Frame B, time from
# discharge) death uses Super learner v7 (stage 18o; A30) and the other
# outcomes are as at diagnosis. All are recalibrated. Each cohort uses its
# own modality set (EHR+ECG; EHR+ECG+CTPA) at 10, 30, 90, 180 and 365 days.
#
# The outcome at h is the event by h + 0.5 days (death; death or CV
# readmission; CV readmission first, a death first counting as no event).
# Calibration is observed/expected, intercept (logit of the risk as
# offset), slope, and ICI, E50 and E90 from a loess curve (span 0.75,
# degree 2, as in stage 11b) evaluated on a grid of 60 risks from the 1st
# to the 99th percentile. Risk groups for the plots are deciles of
# predicted risk, merged from the lowest up until each holds at least 10
# events, so that no group implies a count of 1 to 9.
#
# Decision curves give net benefit NB = TP/n - FP/n x pt / (1 - pt) at
# thresholds of 0.5% to 50% for the model, treat all, treat none and, at
# diagnosis for death, the sPESI rule (1 or more = high risk) and the
# PESI-10 rule (86 or more; PESI-10 lacks altered mental status). They are
# exploratory. At 5%, 10% and 20%, each strategy's NB and the model minus
# every other strategy are tabulated.
#
# Intervals use 2,000 bootstrap samples (seed 42), clustered by patient.
# AUROC, IPA and O/E keep the intervals of stages 18q and 18r. Checks
# (FAIL): inputs found, every risk in [0, 1], scores present, and AUROC,
# IPA, O/E and slope reproducing stages 18q and 18r to within 1e-6.
#
# Inputs:  data/processed/presentation/frame_P.csv, frame_P_outcomes.csv,
#          scores_P.csv; data/processed/frame_B.csv;
#          data/processed/super_learner/slP_preds_*.csv, sl7_preds_B_*.csv,
#          optB2_preds_*.csv; results/slP_metrics.csv, optB2_sl_summary.csv
# Outputs: results/fin_metrics.csv, fin_cal_bins.csv, fin_cal_curves.csv,
#          fin_dca_curves.csv, fin_dca_key.csv, fin_checks.csv,
#          fin_log.txt (terminal log)
# Usage:   Rscript scripts/13b_calibration_dca_final.R [--probe-only]
#          [--workers=N]

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
OD <- file.path(PR, "super_learner")
RS <- file.path(P2, "results")
HZ <- c(10, 30, 90, 180, 365)
NBOOT <- 2000L
SEED <- 42L
NGRID <- 60L
TH_CURVE <- seq(0.005, 0.50, by = 0.005)
TH_KEY <- c(0.05, 0.10, 0.20)
PESI_CUT <- 86
TOL <- 1e-6
PKG <- c("data.table")
options(width = 200)
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-66s %s\n", status, item, value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "fin_checks.csv"))
    cat("\nSTOP: a check failed. see: fin_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- the final models ----------------
COH <- c("EHR+ECG" = "ehrecg", "3-mod" = "3mod")
SETN <- c("EHR+ECG" = "EHR+ECG", "3-mod" = "EHR+ECG+CTPA")
SPEC <- list(
  list(fr = "P", oc = "death", file = "slP_preds_%s.csv",
       model = "Super learner P, recalibrated"),
  list(fr = "P", oc = "comp", file = "optB2_preds_P_comp_%s.csv",
       model = "Super learner, recalibrated"),
  list(fr = "P", oc = "cv", file = "optB2_preds_P_cv_%s.csv",
       model = "Super learner, recalibrated"),
  list(fr = "B", oc = "death", file = "sl7_preds_B_%s.csv",
       model = "Super learner v7, recalibrated"),
  list(fr = "B", oc = "comp", file = "optB2_preds_B_comp_%s.csv",
       model = "Super learner, recalibrated"),
  list(fr = "B", oc = "cv", file = "optB2_preds_B_cv_%s.csv",
       model = "Super learner, recalibrated"))

# ---------------- statistics (sent to workers) ----------------
clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
auc1 <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- data.table::frank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
lo_fit <- function(y, p) {
  tryCatch(suppressWarnings(loess(
    y ~ p, span = 0.75, degree = 2,
    control = loess.control(trace.hat = "approximate"))),
    error = function(e) NULL)
}
nb <- function(y, treat, th) {
  n <- length(y)
  vapply(th, function(t) {
    (sum(treat & y == 1) - sum(treat & y == 0) * t / (1 - t)) / n
  }, numeric(1))
}
nb_model <- function(y, p, th) {
  n <- length(y)
  vapply(th, function(t) {
    tr <- p >= t
    (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / n
  }, numeric(1))
}
# one resample (or the full data): calibration, the curve on the grid,
# net benefit at the key thresholds
stat_vec <- function(J, idx) {
  v <- numeric(0)
  for (ih in seq_along(J$hz)) {
    h <- J$hz[ih]
    y <- J$Y[idx, ih]
    p <- clip(J$P[idx, ih])
    lg <- qlogis(p)
    k <- function(m) paste(h, m, sep = "|")
    ok <- sum(y) > 0 && sum(y) < length(y)
    sl <- if (ok) tryCatch(suppressWarnings(unname(glm.fit(
      cbind(1, lg), y, family = binomial())$coefficients[2])),
      error = function(e) NA_real_) else NA_real_
    it <- if (ok) tryCatch(suppressWarnings(unname(glm.fit(
      matrix(1, length(y), 1), y, family = binomial(),
      offset = lg)$coefficients[1])), error = function(e) NA_real_) else
        NA_real_
    f <- lo_fit(y, p)
    if (is.null(f)) {
      d <- rep(NA_real_, length(p))
      gr <- rep(NA_real_, NGRID)
    } else {
      d <- abs(fitted(f) - p)
      gr <- as.vector(predict(f, newdata = data.frame(p = J$grid[[ih]])))
    }
    v[k("citl")] <- it
    v[k("slope")] <- sl
    v[k("ici")] <- mean(d)
    v[k("e50")] <- unname(quantile(d, 0.5, na.rm = TRUE))
    v[k("e90")] <- unname(quantile(d, 0.9, na.rm = TRUE))
    v[paste(h, "curve", seq_len(NGRID), sep = "|")] <- gr
    nm <- nb_model(y, p, J$th_key)
    na <- nb(y, rep(TRUE, length(y)), J$th_key)
    v[paste(h, "nb", "model", J$th_key, sep = "|")] <- nm
    v[paste(h, "nb", "all", J$th_key, sep = "|")] <- na
    v[paste(h, "nbd", "all", J$th_key, sep = "|")] <- nm - na
    v[paste(h, "nbd", "none", J$th_key, sep = "|")] <- nm
    if (!is.null(J$R)) {
      for (rn in names(J$R)) {
        nr <- nb(y, J$R[[rn]][idx], J$th_key)
        v[paste(h, "nb", rn, J$th_key, sep = "|")] <- nr
        v[paste(h, "nbd", rn, J$th_key, sep = "|")] <- nm - nr
      }
    }
  }
  v
}
# decile groups of predicted risk; from the lowest, neighbouring deciles
# are merged until each group holds at least 10 events (or none), so no
# group reveals a count from 1 to 9 (MIMIC small-cell practice)
risk_groups <- function(y, p) {
  g <- ceiling(10 * frank(p, ties.method = "first") / length(p))
  ev <- tapply(y, g, sum)
  grp <- integer(10)
  cur <- 1L
  acc <- 0
  for (d in 1:10) {
    grp[d] <- cur
    acc <- acc + ev[[d]]
    if (acc >= 10 || d == 10) {
      cur <- cur + 1L
      acc <- 0
    }
  }
  # a last group with 1 to 9 events joins the one below it
  last <- max(grp)
  if (last > 1 && sum(ev[grp == last]) %in% 1:9) grp[grp == last] <- last - 1L
  D <- data.table(pp = p, yy = y, grp = grp[g], dec = g)
  out <- D[, .(deciles = if (min(dec) == max(dec)) as.character(min(dec))
                 else sprintf("%d to %d", min(dec), max(dec)),
               mean_predicted = mean(pp), observed = mean(yy), n = .N,
               ev = sum(yy)), by = grp][order(grp)]
  # only possible when the whole outcome has 1 to 9 events
  out[ev %in% 1:9, observed := NA_real_]
  setnames(out, "grp", "group")
  out[, ev := NULL]
  out[]
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by), length(J$rows_by),
                                     replace = TRUE)], use.names = FALSE)
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
FR <- list(P_death = file.path(PP, "frame_P.csv"),
           P_out = file.path(PP, "frame_P_outcomes.csv"),
           B = file.path(PR, "frame_B.csv"))
need <- c(unlist(FR), file.path(PP, "scores_P.csv"),
          file.path(RS, c("slP_metrics.csv", "optB2_sl_summary.csv")))
for (s in SPEC) need <- c(need, file.path(OD, sprintf(s$file, COH)))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else
  sprintf("all %d", length(need)), if (length(miss) == 0) "PASS" else
    "FAIL")
stop_if_fail()
FD <- lapply(FR, fread)
if (!("in_three_mod_p" %in% names(FD$B))) {
  FD$B[, in_three_mod_p := in_three_mod]
}
SCO <- unique(fread(file.path(PP, "scores_P.csv")), by = "hadm_id")
REFP <- fread(file.path(RS, "slP_metrics.csv"))[evaluation == "all"]
REFB <- fread(file.path(RS, "optB2_sl_summary.csv"))[evaluation == "all"]
JOBS <- list()
for (s in SPEC) for (co in names(COH)) {
  F0 <- if (s$fr == "B") FD$B else if (s$oc == "death") FD$P_death else
    FD$P_out
  sub <- if (co == "EHR+ECG") F0 else F0[in_three_mod_p == 1]
  n <- nrow(sub)
  tt <- if (s$oc == "death") sub$time_death else sub$time_first
  ee <- switch(s$oc, death = sub$event_death == 1,
               comp = sub$event_comp == 1, cv = sub$event_cr == 1)
  Y <- sapply(HZ, function(h) as.integer(ee & tt <= h + 0.5))
  Y <- matrix(Y, ncol = length(HZ))
  L <- fread(file.path(OD, sprintf(s$file, COH[[co]])))
  L <- L[model == s$model & set == SETN[[co]]]
  P <- sapply(HZ, function(h) {
    x <- L[horizon == h]
    x$risk[match(sub$hadm_id, x$hadm_id)]
  })
  P <- matrix(P, ncol = length(HZ))
  job <- sprintf("%s %s %s", s$fr, s$oc, co)
  bad <- sum(!is.finite(P))
  rg <- suppressWarnings(range(P, na.rm = TRUE))
  chk(sprintf("%s: risks present and in [0, 1] (missing; range)", job),
      sprintf("%d; %.5f to %.5f", bad, rg[1], rg[2]),
      if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
  stop_if_fail()
  R <- NULL
  if (s$fr == "P" && s$oc == "death") {
    sc <- SCO[match(sub$hadm_id, hadm_id)]
    ns <- sum(!is.finite(sc$spesi)) + sum(!is.finite(sc$pesi10))
    chk(sprintf("%s: sPESI and PESI-10 present (missing)", job), ns,
        if (ns == 0) "PASS" else "FAIL")
    stop_if_fail()
    R <- list(spesi = sc$spesi >= 1, pesi = sc$pesi10 >= PESI_CUT)
    chk(sprintf("%s: high risk by sPESI / by PESI-10 class III to V", job),
        sprintf("%.1f%% / %.1f%%", 100 * mean(R$spesi),
                100 * mean(R$pesi)), "INFO")
  }
  grid <- lapply(seq_along(HZ), function(ih) {
    q <- quantile(P[, ih], c(0.01, 0.99), names = FALSE)
    seq(q[1], q[2], length.out = NGRID)
  })
  JOBS[[job]] <- list(
    job = job, fr = s$fr, oc = s$oc, co = co, set = SETN[[co]],
    model = s$model, hz = HZ, n = n, Y = Y, P = P, R = R, grid = grid,
    th_key = TH_KEY,
    rows_by = split(seq_len(n), match(sub$subject_id,
                                      unique(sub$subject_id))))
  cat(sprintf("  %-18s n=%d patients=%d events by 30/365 days %d/%d\n",
              job, n, length(JOBS[[job]]$rows_by), sum(Y[, 2]),
              sum(Y[, 5])))
}

# ---------------- point estimates and the reproduction check ----------
cat("\n=== REPRODUCTION OF STAGES 18q AND 18r (point estimates) ===\n")
PT <- list()
REP <- list()
for (J in JOBS) {
  PT[[J$job]] <- stat_vec(J, seq_len(J$n))
  for (ih in seq_along(HZ)) {
    y <- J$Y[, ih]
    raw <- J$P[, ih]
    p <- clip(raw)
    pr <- mean(y)
    mine <- c(auc = auc1(y, raw), ipa = 1 - mean((y - p)^2) /
                (pr * (1 - pr)), oe = pr / mean(p),
              slope = unname(PT[[J$job]][paste(HZ[ih], "slope",
                                               sep = "|")]))
    ref <- if (J$fr == "P" && J$oc == "death") {
      REFP[cohort == J$co & model == J$model & set == J$set &
             horizon == HZ[ih]]
    } else {
      REFB[frame == J$fr & outcome == J$oc & cohort == J$co &
             model == J$model & set == J$set & horizon == HZ[ih]]
    }
    for (m in names(mine)) {
      r <- ref[measure == m]
      REP[[length(REP) + 1]] <- data.table(
        job = J$job, horizon = HZ[ih], measure = m, mine = mine[[m]],
        reference = if (nrow(r)) r$estimate[1] else NA_real_)
    }
  }
}
REP <- rbindlist(REP)
REP[, diff := abs(mine - reference)]
for (J in JOBS) {
  x <- REP[job == J$job]
  nm <- sum(is.na(x$reference))
  md <- max(x$diff, na.rm = TRUE)
  chk(sprintf("%s: AUROC, IPA, O/E, slope reproduced (missing; max diff)",
              J$job), sprintf("%d; %.1e", nm, md),
      if (nm == 0 && md <= TOL) "PASS" else "FAIL")
}
stop_if_fail()

# ---------------- timing probe ----------------
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(JOBS, function(J) J$n))
set.seed(1)
t0 <- Sys.time()
for (i in 1:20) invisible(boot_rep(i, JOBS[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 20
tot <- sum(sapply(JOBS, function(J) J$n)) / JOBS[[big]]$n * tp * NBOOT
cat(sprintf(paste("  %.3f s per resample of the largest job; %d jobs x",
                  "%d resamples: about %.1f min on %d workers\n"),
            tp, length(JOBS), NBOOT, tot * 1.25 / workers / 60, workers))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples, clustered by patient) ===\n",
            NBOOT))
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
MET <- list()
CUR <- list()
KEY <- list()
BIN <- list()
DCA <- list()
plan(multisession, workers = workers)
for (J in JOBS) {
  pt <- PT[[J$job]]
  BL <- future_lapply(seq_len(NBOOT), boot_rep, J = J,
                      future.seed = SEED, future.packages = PKG)
  B <- do.call(rbind, lapply(BL, function(v) v[names(pt)]))
  colnames(B) <- names(pt)
  nfail <- sum(!is.finite(B[, grep("\\|slope$", colnames(B))]))
  cat(sprintf("  %-18s done (%d non-finite slopes of %d)\n", J$job, nfail,
              NBOOT * length(HZ)))
  lab <- data.table(frame = J$fr, outcome = J$oc, cohort = J$co,
                    set = J$set, model = J$model)
  for (ih in seq_along(HZ)) {
    h <- HZ[ih]
    y <- J$Y[, ih]
    p <- J$P[, ih]
    ev <- sum(y)
    for (m in c("citl", "slope", "ici", "e50", "e90")) {
      k <- paste(h, m, sep = "|")
      q <- q2(B[, k])
      MET[[length(MET) + 1]] <- cbind(lab, data.table(
        horizon = h, n = J$n, events = ev, measure = m,
        estimate = pt[[k]], lo = q[1], hi = q[2],
        source = "stage 13b bootstrap"))
    }
    x <- REP[job == J$job & horizon == h]
    for (m in c("auc", "ipa", "oe")) {
      MET[[length(MET) + 1]] <- cbind(lab, data.table(
        horizon = h, n = J$n, events = ev, measure = m,
        estimate = x[measure == m]$mine, lo = NA_real_, hi = NA_real_,
        source = "point estimate (intervals: stages 18q, 18r)"))
    }
    ks <- paste(h, "curve", seq_len(NGRID), sep = "|")
    qq <- apply(B[, ks, drop = FALSE], 2, q2)
    CUR[[length(CUR) + 1]] <- cbind(lab, data.table(
      horizon = h, point = seq_len(NGRID), risk = J$grid[[ih]],
      observed = unname(pt[ks]), lo = qq[1, ], hi = qq[2, ],
      boot_ok = colSums(is.finite(B[, ks, drop = FALSE]))))
    BIN[[length(BIN) + 1]] <- cbind(lab, risk_groups(y, p)[
      , horizon := h])
    strat <- c(list(model = NULL, all = rep(TRUE, length(y))), J$R)
    for (sn in names(strat)) {
      v <- if (sn == "model") nb_model(y, clip(p), TH_CURVE) else
        nb(y, strat[[sn]], TH_CURVE)
      DCA[[length(DCA) + 1]] <- cbind(lab, data.table(
        horizon = h, prevalence = mean(y), strategy = sn,
        threshold = TH_CURVE, net_benefit = v))
    }
    for (tk in TH_KEY) {
      for (sn in names(strat)) {
        k <- paste(h, "nb", sn, tk, sep = "|")
        q <- q2(B[, k])
        KEY[[length(KEY) + 1]] <- cbind(lab, data.table(
          horizon = h, threshold = tk, quantity = "net benefit",
          strategy = sn, estimate = pt[[k]], lo = q[1], hi = q[2],
          excludes_zero = NA))
      }
      for (sn in c("none", setdiff(names(strat), "model"))) {
        k <- paste(h, "nbd", sn, tk, sep = "|")
        q <- q2(B[, k])
        KEY[[length(KEY) + 1]] <- cbind(lab, data.table(
          horizon = h, threshold = tk,
          quantity = "model minus strategy", strategy = sn,
          estimate = pt[[k]], lo = q[1], hi = q[2],
          excludes_zero = q[1] > 0 | q[2] < 0))
      }
    }
  }
}
plan(sequential)
SNAME <- c(model = "Model", all = "Treat all", none = "Treat none",
           spesi = "sPESI 1 or more",
           pesi = "PESI-10 class III to V")
mt <- rbindlist(MET)
mt[, events := ifelse(events > 0 & events < 10, "<10",
                      as.character(events))]
cu <- rbindlist(CUR)
bn <- rbindlist(BIN)
setcolorder(bn, c("frame", "outcome", "cohort", "set", "model", "horizon"))
dc <- rbindlist(DCA)
dc[, strategy := SNAME[strategy]]
ky <- rbindlist(KEY)
ky[, strategy := SNAME[strategy]]
fwrite(mt, file.path(RS, "fin_metrics.csv"))
fwrite(bn, file.path(RS, "fin_cal_bins.csv"))
fwrite(cu, file.path(RS, "fin_cal_curves.csv"))
fwrite(dc, file.path(RS, "fin_dca_curves.csv"))
fwrite(ky, file.path(RS, "fin_dca_key.csv"))
nb_curve <- cu[, sum(boot_ok < NBOOT * 0.9)]
chk("curve grid points with under 90% usable resamples", nb_curve,
    if (nb_curve == 0) "PASS" else "INFO")
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "fin_checks.csv"))

# ---------------- summary ----------------
cat("\n=== CALIBRATION, EHR+ECG COHORT (slope; ICI) ===\n")
x <- mt[cohort == "EHR+ECG" & measure %in% c("slope", "ici")]
x[, cell := sprintf("%.3f (%.3f, %.3f)", estimate, lo, hi)]
print(dcast(x, frame + outcome + measure ~ horizon, value.var = "cell"),
      row.names = FALSE)
cat("\n=== NET BENEFIT, MODEL MINUS TREAT ALL / sPESI (EHR+ECG cohort) ===\n")
x <- ky[cohort == "EHR+ECG" & quantity == "model minus strategy" &
          strategy %in% c("Treat all", "sPESI 1 or more") &
          horizon %in% c(30, 365)]
print(x[, .(frame, outcome, horizon, threshold, strategy,
            est = round(estimate, 4), lo = round(lo, 4), hi = round(hi, 4),
            excludes_zero)], row.names = FALSE, nrows = 200)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("fin_metrics.csv", "fin_cal_bins.csv", "fin_cal_curves.csv",
            "fin_dca_curves.csv", "fin_dca_key.csv", "fin_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "fin_log.txt"), " (terminal log)\n")
cat("\nOutputs: fin_metrics.csv, fin_cal_bins.csv, fin_cal_curves.csv,",
    "fin_dca_curves.csv, fin_dca_key.csv, fin_checks.csv\n")