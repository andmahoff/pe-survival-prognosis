# 14a_rp_curve.R
# ECG hazard ratio as a smooth curve over the year after discharge
# (PREREGISTRATION_ADDENDUM.md, A1.1 and A1.2). Exploratory, specified after
# the primary analysis; intervals only.
#
# Cohort: Frame B, death within 365 days of discharge. The scores are the
# seed-averaged out-of-fold ridge-Cox death scores from script 06a,
# winsorised at the 0.5th and 99.5th percentiles and standardised, so a
# hazard ratio is per 1 SD of a score.
#
# A1.1: rstpm2::stpm2 with a proportional-hazards link, a natural spline in
# log time for the baseline (4 df) and a time-varying effect (3 df) for
# every score. The curve is the hazard ratio for a 1 SD higher score at each
# day from 1 to 365 for a patient with average scores, with a 95% pointwise
# delta-method band. For the main model, 1,000 patient-level bootstrap
# refits (seed 42) give intervals for the hazard ratio at days 7, 14, 30,
# 90, 198 and 365 and for HR(day 14) / HR(day 198), the midpoints of the
# primary analysis's two periods. A likelihood-ratio test (3 df) assesses
# each score's time-varying effect. Sensitivity: baseline df 3 to 5 by
# time-varying df 1 to 4 (AIC and key hazard ratios), registered scaling,
# hospice discharges excluded, and the three-modality cohort with CTPA.
#
# A1.2: Cox models on data split at every event time, with constant effects,
# a step at day 30 (the primary model), effects linear in log time and a
# natural spline in log time (3 df), compared by AIC. A change-point profile
# steps the ECG effect at day 7, 14, 21, 30, 45, 60, 90, 120, 180 or 270
# (EHR constant) and reports the partial log-likelihood and the hazard
# ratios either side of each step.
#
# Notes:   the log-likelihood and parameter count are read from the slots of
#          the bbmle 'mle2' objects (min, coef), not through method
#          dispatch. A timing probe stops the script if no bootstrap refit
#          returns values.
#
# Inputs:  data/processed/frame_B.csv,
#          data/processed/level1/oof_linear_B_*.csv
# Outputs: results/rp_curve.csv, rp_summary.csv, rp_models.csv,
#          tvc_shape.csv, rp_log.txt
# Usage:   Rscript scripts/14a_rp_curve.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(splines)
  library(data.table)
  library(future.apply)
  library(bbmle)
  library(rstpm2)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
NBOOT <- 1000L
SEED <- 42L
DF_MAIN <- 4L
TVC_MAIN <- 3L
KEYDAYS <- c(7, 14, 30, 90, 198, 365)
DAYS <- 1:365
TAUS <- c(7, 14, 21, 30, 45, 60, 90, 120, 180, 270)
PKG <- c("survival", "splines", "bbmle", "rstpm2")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
cat(sprintf("rstpm2 %s, bbmle %s; %d workers\n",
            as.character(packageVersion("rstpm2")),
            as.character(packageVersion("bbmle")), workers))

# ---------------- data ----------------
std <- function(v) (v - mean(v)) / sd(v)
wins <- function(v) {
  q <- quantile(v, c(0.005, 0.995), names = FALSE)
  std(pmin(pmax(v, q[1]), q[2]))
}
fb <- fread(file.path(PR, "frame_B.csv"))
make_data <- function(cohort, scaling = "winsorised",
                      drop_hospice = FALSE) {
  code <- if (cohort == "EHR+ECG") "ehrecg" else "3mod"
  mods <- c("ehr", "ecg", if (cohort == "3-mod") "ctpa")
  sub <- if (cohort == "EHR+ECG") fb else fb[in_three_mod == 1]
  if (drop_hospice) sub <- sub[hospice == 0]
  lp <- sprintf("%s.cox.death.lp", mods)
  x <- fread(file.path(PR, "level1",
                       sprintf("oof_linear_B_%s.csv", code)),
             select = c("hadm_id", lp))
  x <- x[, lapply(.SD, mean), by = hadm_id, .SDcols = lp]
  d <- merge(sub[, .(hadm_id, subject_id, time_death, event_death)], x,
             by = "hadm_id")
  f <- if (scaling == "winsorised") wins else std
  for (j in seq_along(mods)) set(d, j = mods[j], value = f(d[[lp[j]]]))
  list(d = as.data.frame(d[, c("subject_id", "time_death", "event_death",
                               mods), with = FALSE]),
       mods = mods)
}

# ---------------- model functions ----------------
rp_one <- function(d, mods, df, tvc_df, tvc_drop = character(0)) {
  fml <- as.formula(paste("Surv(time_death, event_death) ~",
                          paste(mods, collapse = " + ")))
  tv <- setdiff(mods, tvc_drop)
  tvl <- NULL
  if (tvc_df > 0 && length(tv) > 0) {
    tvl <- setNames(as.list(rep(as.integer(tvc_df), length(tv))), tv)
  }
  rstpm2::stpm2(fml, data = d, df = df, tvc = tvl)
}
rp_ll <- function(f) -as.numeric(f@min)
rp_np <- function(f) length(f@coef)
rp_conv <- function(f) {
  tryCatch(isTRUE(f@details$convergence == 0), error = function(e) NA)
}
rp_hr <- function(f, mods, var, days, se = FALSE) {
  nd <- data.frame(time_death = days)
  for (m in mods) nd[[m]] <- 0
  predict(f, newdata = nd, type = "hr", var = var, se.fit = se)
}
boot_rp <- function(i, d, mods, df, tvc_df, days) {
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
  f <- tryCatch(rp_one(d[idx, , drop = FALSE], mods, df, tvc_df),
                error = function(e) NULL)
  out <- rep(NA_real_, length(mods) * length(days))
  if (is.null(f)) return(out)
  for (j in seq_along(mods)) {
    out[(j - 1) * length(days) + seq_along(days)] <- tryCatch(
      as.numeric(rp_hr(f, mods, mods[j], days)),
      error = function(e) rep(NA_real_, length(days)))
  }
  out
}
grid_rp <- function(g, d, mods, days) {
  row <- data.frame(baseline_df = g$df, tvc_df = g$tvc,
                    loglik = NA_real_, npar = NA_integer_, aic = NA_real_,
                    converged = NA, error = "")
  f <- tryCatch(rp_one(d, mods, g$df, g$tvc),
                error = function(e) conditionMessage(e))
  if (is.character(f)) {
    row$error <- substr(f, 1, 150)
    return(row)
  }
  row$loglik <- rp_ll(f)
  row$npar <- rp_np(f)
  row$aic <- -2 * row$loglik + 2 * row$npar
  row$converged <- rp_conv(f)
  for (m in mods) {
    h <- tryCatch(as.numeric(rp_hr(f, mods, m, days)),
                  error = function(e) rep(NA_real_, length(days)))
    for (i in seq_along(days)) {
      row[[sprintf("%s_hr_day%d", toupper(m), days[i])]] <- h[i]
    }
  }
  row
}
variant_rp <- function(v, df, tvc_df, days) {
  d <- v$d
  mods <- v$mods
  errs <- character(0)
  f <- tryCatch(rp_one(d, mods, df, tvc_df),
                error = function(e) conditionMessage(e))
  if (is.character(f)) {
    return(list(curve = NULL, tests = NULL,
                errors = paste(v$name, "main fit:", f)))
  }
  conv <- rp_conv(f)
  ll <- rp_ll(f)
  cur <- list()
  tst <- list()
  for (m in mods) {
    pr <- tryCatch(as.data.frame(rp_hr(f, mods, m, days, se = TRUE)),
                   error = function(e) conditionMessage(e))
    if (is.character(pr)) {
      errs <- c(errs, paste(v$name, m, "curve:", pr))
    } else {
      cur[[m]] <- data.frame(variant = v$name, modality = toupper(m),
                             day = days, hr = pr[[1]], lo = pr[[2]],
                             hi = pr[[3]])
    }
    f0 <- tryCatch(rp_one(d, mods, df, tvc_df, tvc_drop = m),
                   error = function(e) conditionMessage(e))
    if (is.character(f0)) {
      errs <- c(errs, paste(v$name, m, "reduced fit:", f0))
      lr <- NA_real_
    } else {
      lr <- 2 * (ll - rp_ll(f0))
    }
    tst[[m]] <- data.frame(variant = v$name, modality = toupper(m),
                           lr = lr, df = tvc_df,
                           p = pchisq(lr, tvc_df, lower.tail = FALSE),
                           n = nrow(d), deaths = sum(d$event_death),
                           converged = conv)
  }
  list(curve = if (length(cur)) do.call(rbind, cur) else NULL,
       tests = do.call(rbind, tst), errors = errs)
}

# ---------------- shape: step or gradual (Cox) ----------------
shape_split <- function(d) {
  ev <- sort(unique(d$time_death[d$event_death == 1]))
  s <- survSplit(Surv(time_death, event_death) ~ ehr + ecg, data = d,
                 cut = ev, episode = "iv")
  s$lt <- log(s$time_death)
  kn <- quantile(log(d$time_death[d$event_death == 1]), c(1 / 3, 2 / 3),
                 names = FALSE)
  B <- ns(s$lt, knots = kn, Boundary.knots = range(log(ev)))
  for (m in c("ehr", "ecg")) {
    s[[paste0(m, "_e")]] <- s[[m]] * (s$time_death <= 30.5)
    s[[paste0(m, "_l")]] <- s[[m]] * (s$time_death > 30.5)
    s[[paste0(m, "_lt")]] <- s[[m]] * s$lt
    for (j in 1:3) s[[sprintf("%s_b%d", m, j)]] <- s[[m]] * B[, j]
  }
  s
}
shape_fit <- function(s, rhs) {
  f <- coxph(as.formula(paste(
    "Surv(tstart, time_death, event_death) ~", rhs)),
    data = s, ties = "efron")
  list(ll = f$loglik[2], k = length(coef(f)), coef = coef(f))
}

cat("=== DATA ===\n")
V <- list(
  c(list(name = "main"), make_data("EHR+ECG")),
  c(list(name = "registered scaling"),
    make_data("EHR+ECG", "standardised")),
  c(list(name = "hospice excluded"),
    make_data("EHR+ECG", drop_hospice = TRUE)),
  c(list(name = "three-modality cohort"), make_data("3-mod")))
for (v in V) {
  cat(sprintf("  %-22s n=%d  deaths=%d  scores %s\n", v$name, nrow(v$d),
              sum(v$d$event_death), paste(v$mods, collapse = "/")))
}
main <- V[[1]]

# ---------------- timing probe ----------------
cat("\n=== TIMING PROBE (real data) ===\n")
t0 <- Sys.time()
f_probe <- rp_one(main$d, main$mods, DF_MAIN, TVC_MAIN)
t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
h_probe <- as.numeric(rp_hr(f_probe, main$mods, "ecg", c(14, 198)))
cat(sprintf("  check on the real fit: log-likelihood %.1f, %d parameters,",
            rp_ll(f_probe), rp_np(f_probe)),
    sprintf("converged %s; ECG HR day 14 %.3f, day 198 %.3f\n",
            rp_conv(f_probe), h_probe[1], h_probe[2]))
plan(multisession, workers = workers)
nb <- workers * 2L
t0 <- Sys.time()
bp <- do.call(rbind, future_lapply(
  seq_len(nb), boot_rp, d = main$d, mods = main$mods, df = DF_MAIN,
  tvc_df = TVC_MAIN, days = KEYDAYS, future.packages = PKG,
  future.seed = 1L))
t_boot <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / nb
n_bp <- sum(complete.cases(bp))
cat(sprintf("  probe bootstrap refits returning values: %d of %d\n",
            n_bp, nb))
if (n_bp == 0) {
  plan(sequential)
  cat("STOP: no bootstrap refit returned values in the workers.\n")
  quit(status = 1)
}
t0 <- Sys.time()
S <- shape_split(main$d)
invisible(shape_fit(S, "ehr + ecg"))
t_shape <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
proj <- t_boot * NBOOT + t_fit * (12 + 11) / workers * 1.5 +
  t_shape * 14
cat(sprintf("  one fit %.1f s; bootstrap %.2f s per refit in parallel;",
            t_fit, t_boot),
    sprintf("one shape fit %.1f s (%d split rows)\n", t_shape, nrow(S)))
cat(sprintf("  projected total: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- A1.1 variants: curves and LR tests ----------------
cat("\n=== A1.1 CURVES (4 variants) ===\n")
VR <- future_lapply(V, variant_rp, df = DF_MAIN, tvc_df = TVC_MAIN,
                    days = DAYS, future.packages = PKG,
                    future.seed = TRUE)
errs <- unlist(lapply(VR, function(r) r$errors))
if (length(errs)) {
  cat("  errors inside the variant fits:\n")
  for (e in errs) cat("   ", e, "\n")
}
curve <- rbindlist(lapply(VR, function(r) r$curve), fill = TRUE)
tests <- rbindlist(lapply(VR, function(r) r$tests), fill = TRUE)
if (!nrow(curve) || !"main" %in% curve$variant) {
  plan(sequential)
  cat("STOP: the main curve could not be computed (errors above).\n")
  quit(status = 1)
}
fwrite(curve, file.path(RS, "rp_curve.csv"))
cat(sprintf("  curve rows %d (expected %d)\n", nrow(curve),
            length(DAYS) * sum(sapply(V, function(v) length(v$mods)))))

# ---------------- A1.1 model grid ----------------
cat("=== A1.1 MODEL GRID (baseline df x time-varying df) ===\n")
G <- list()
for (a in 3:5) for (b in 1:4) G[[length(G) + 1]] <- list(df = a, tvc = b)
GR <- future_lapply(G, grid_rp, d = main$d, mods = main$mods,
                    days = KEYDAYS, future.packages = PKG,
                    future.seed = TRUE)
grid <- rbindlist(GR, fill = TRUE)
grid[, delta_aic := aic - min(aic, na.rm = TRUE)]
fwrite(grid, file.path(RS, "rp_models.csv"))
if (any(grid$error != "")) {
  cat("  grid fits with errors:", sum(grid$error != ""), "\n")
}

# ---------------- A1.1 bootstrap of the main model ----------------
cat(sprintf("=== A1.1 BOOTSTRAP (%d refits) ===\n", NBOOT))
B <- do.call(rbind, future_lapply(
  seq_len(NBOOT), boot_rp, d = main$d, mods = main$mods, df = DF_MAIN,
  tvc_df = TVC_MAIN, days = KEYDAYS, future.packages = PKG,
  future.seed = SEED))
plan(sequential)
n_ok <- sum(complete.cases(B))
mc <- curve[variant == "main"]
SUM <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (j in seq_along(main$mods)) {
  M <- toupper(main$mods[j])
  cb <- (j - 1) * length(KEYDAYS) + seq_along(KEYDAYS)
  for (i in seq_along(KEYDAYS)) {
    r <- mc[modality == M & day == KEYDAYS[i]]
    qb <- q2(B[, cb[i]])
    SUM[[length(SUM) + 1]] <- data.table(
      variant = "main", modality = M,
      quantity = sprintf("HR per SD at day %d", KEYDAYS[i]),
      estimate = r$hr, lo_boot = qb[1], hi_boot = qb[2],
      lo_delta = r$lo, hi_delta = r$hi, statistic = NA_real_,
      p_value = NA_real_, n_boot_ok = n_ok)
  }
  i14 <- which(KEYDAYS == 14)
  i198 <- which(KEYDAYS == 198)
  est <- mc[modality == M & day == 14, hr] /
    mc[modality == M & day == 198, hr]
  qb <- q2(B[, cb[i14]] / B[, cb[i198]])
  SUM[[length(SUM) + 1]] <- data.table(
    variant = "main", modality = M,
    quantity = "ratio HR(day 14) / HR(day 198)", estimate = est,
    lo_boot = qb[1], hi_boot = qb[2], lo_delta = NA_real_,
    hi_delta = NA_real_, statistic = NA_real_, p_value = NA_real_,
    n_boot_ok = n_ok)
}
for (i in seq_len(nrow(tests))) {
  SUM[[length(SUM) + 1]] <- data.table(
    variant = tests$variant[i], modality = tests$modality[i],
    quantity = sprintf("LR test of time-varying effect (%d df)",
                       tests$df[i]),
    estimate = NA_real_, lo_boot = NA_real_, hi_boot = NA_real_,
    lo_delta = NA_real_, hi_delta = NA_real_,
    statistic = tests$lr[i], p_value = tests$p[i], n_boot_ok = NA)
}
sm <- rbindlist(SUM)
fwrite(sm, file.path(RS, "rp_summary.csv"))

# ---------------- A1.2 shape ----------------
cat("=== A1.2 SHAPE: step or gradual ===\n")
forms <- c(
  "constant (proportional hazards)" = "ehr + ecg",
  "step at day 30 (primary model)" = "ehr_e + ehr_l + ecg_e + ecg_l",
  "linear in log time" = "ehr + ecg + ehr_lt + ecg_lt",
  "natural spline in log time, 3 df" =
    "ehr + ecg + ehr_b1 + ehr_b2 + ehr_b3 + ecg_b1 + ecg_b2 + ecg_b3")
SH <- list()
for (nm in names(forms)) {
  r <- shape_fit(S, forms[[nm]])
  SH[[length(SH) + 1]] <- data.table(
    part = "model comparison", model = nm, step_day = NA_real_,
    loglik = r$ll, k = r$k, aic = -2 * r$ll + 2 * r$k,
    ecg_hr_before = NA_real_, ecg_hr_after = NA_real_)
}
for (tau in TAUS) {
  S$c_e <- S$ecg * (S$time_death <= tau + 0.5)
  S$c_l <- S$ecg * (S$time_death > tau + 0.5)
  r <- shape_fit(S, "ehr + c_e + c_l")
  SH[[length(SH) + 1]] <- data.table(
    part = "change-point profile (ECG step, EHR constant)",
    model = sprintf("ECG step at day %d", tau), step_day = tau,
    loglik = r$ll, k = r$k, aic = -2 * r$ll + 2 * r$k,
    ecg_hr_before = exp(r$coef[["c_e"]]),
    ecg_hr_after = exp(r$coef[["c_l"]]))
}
sh <- rbindlist(SH)
sh[, delta_aic := aic - min(aic), by = part]
sh[, ratio_before_after := ecg_hr_before / ecg_hr_after]
fwrite(sh, file.path(RS, "tvc_shape.csv"))

# ---------------- print ----------------
cat("\n=== MAIN CURVE (hazard ratio per SD, delta-method 95% band) ===\n")
print(mc[day %in% KEYDAYS, .(modality, day, hr = round(hr, 3),
                             lo = round(lo, 3), hi = round(hi, 3))],
      row.names = FALSE)
cat("\n=== SUMMARY ===\n")
print(sm[, .(variant, modality, quantity, estimate = round(estimate, 3),
             lo_boot = round(lo_boot, 3), hi_boot = round(hi_boot, 3),
             statistic = round(statistic, 2),
             p_value = signif(p_value, 3))], row.names = FALSE)
cat("\n=== MODEL GRID ===\n")
print(grid[, .(baseline_df, tvc_df, aic = round(aic, 1),
               delta_aic = round(delta_aic, 1), converged,
               ECG_d14 = round(ECG_hr_day14, 3),
               ECG_d198 = round(ECG_hr_day198, 3))], row.names = FALSE)
cat("\n=== SHAPE ===\n")
print(sh[, .(part, model, delta_aic = round(delta_aic, 1),
             ecg_before = round(ecg_hr_before, 3),
             ecg_after = round(ecg_hr_after, 3))], row.names = FALSE)
cat(sprintf("\n  bootstrap refits OK: %d of %d; elapsed %.1f min\n", n_ok,
            NBOOT,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("rp_curve.csv", "rp_summary.csv", "rp_models.csv",
            "tvc_shape.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "rp_log.txt"), " (terminal log)\n")
cat("\nOutputs: rp_curve.csv, rp_summary.csv, rp_models.csv,",
    "tvc_shape.csv\n")