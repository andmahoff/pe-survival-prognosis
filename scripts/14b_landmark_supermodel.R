# 14b_landmark_supermodel.R
# Landmark supermodel for the time-varying hazard ratio of each score
# (PREREGISTRATION_ADDENDUM.md, A1.3). Exploratory, specified after the
# primary analysis; intervals only.
#
# Landmarks s are placed every 15 days from 0 to 330 after discharge. At
# each landmark the patients still alive are followed for deaths on days s
# to s + 30 and censored at s + 31 (D1 convention). The landmark data sets
# are stacked and one Cox model, stratified by landmark with a sandwich
# variance clustered by patient, gives each score's log hazard ratio as
# theta0 + theta1 (s / smax) + theta2 (s / smax)^2. Reported are the hazard
# ratio per 1 SD for death in the next 30 days as a function of time since
# discharge (95% bands), landmark-by-landmark estimates, the Wald test of
# theta1 = theta2 = 0 (2 df) and the ratio HR(s = 0) / HR(s = 180).
#
# Cohort: Frame B. Scores are the seed-averaged out-of-fold ridge-Cox death
# scores, winsorised at the 0.5th and 99.5th percentiles and standardised
# (main). Sensitivity: a 60-day window (landmarks 0 to 300), landmarks every
# 30 days, registered scaling, hospice discharges excluded, and the
# three-modality cohort with CTPA. The six variants run in parallel, after a
# timing probe of the main variant.
#
# Inputs:  data/processed/frame_B.csv,
#          data/processed/level1/oof_linear_B_*.csv
# Outputs: results/lm_curve.csv, lm_points.csv, lm_tests.csv, lm_log.txt
# Usage:   Rscript scripts/14b_landmark_supermodel.R [--probe-only]

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
PKG <- c("survival")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

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

run_lm <- function(v) {
  d <- v$d
  mods <- v$mods
  lms <- seq(0, 365 - v$w, by = v$step)
  smax <- max(lms)
  st <- list()
  for (s in lms) {
    a <- d[d$time_death > s, , drop = FALSE]
    ev <- a$event_death == 1 & a$time_death <= s + v$w + 0.5
    st[[length(st) + 1]] <- data.frame(
      subject_id = a$subject_id, lm = s, tstart = s,
      tstop = ifelse(ev, a$time_death, s + v$w + 1),
      status = as.integer(ev), a[, mods, drop = FALSE])
  }
  S <- do.call(rbind, st)
  S$g1 <- S$lm / smax
  S$g2 <- S$g1^2
  terms <- c()
  for (m in mods) {
    S[[paste0(m, "_g1")]] <- S[[m]] * S$g1
    S[[paste0(m, "_g2")]] <- S[[m]] * S$g2
    terms <- c(terms, m, paste0(m, "_g1"), paste0(m, "_g2"))
  }
  fml <- as.formula(paste("Surv(tstart, tstop, status) ~",
                          paste(terms, collapse = " + "),
                          "+ strata(lm) + cluster(subject_id)"))
  fit <- survival::coxph(fml, data = S, ties = "efron")
  b <- coef(fit)
  Vr <- fit$var
  sg <- seq(0, smax, by = 1)
  G <- cbind(1, sg / smax, (sg / smax)^2)
  cur <- list()
  tst <- list()
  for (m in mods) {
    ix <- match(c(m, paste0(m, "_g1"), paste0(m, "_g2")), names(b))
    th <- b[ix]
    Vm <- Vr[ix, ix]
    eta <- as.vector(G %*% th)
    se <- sqrt(rowSums((G %*% Vm) * G))
    cur[[m]] <- data.frame(
      variant = v$name, modality = toupper(m), window_days = v$w,
      landmark_day = sg, hr = exp(eta), lo = exp(eta - 1.96 * se),
      hi = exp(eta + 1.96 * se))
    w2 <- as.numeric(t(th[2:3]) %*% solve(Vm[2:3, 2:3]) %*% th[2:3])
    dg <- c(1, 0, 0) - c(1, 180 / smax, (180 / smax)^2)
    lr <- sum(dg * th)
    lse <- sqrt(as.numeric(t(dg) %*% Vm %*% dg))
    tst[[m]] <- data.frame(
      variant = v$name, modality = toupper(m), window_days = v$w,
      landmark_step = v$step, landmarks = length(lms),
      stacked_rows = nrow(S), n = nrow(d), deaths = sum(d$event_death),
      theta0 = th[[1]], theta1 = th[[2]], theta2 = th[[3]],
      wald_2df = w2, p_value = pchisq(w2, 2, lower.tail = FALSE),
      ratio_hr_s0_over_s180 = exp(lr),
      ratio_lo = exp(lr - 1.96 * lse), ratio_hi = exp(lr + 1.96 * lse))
  }
  pts <- list()
  f1 <- as.formula(paste("Surv(tstart, tstop, status) ~",
                         paste(mods, collapse = " + ")))
  for (s in lms) {
    Ss <- S[S$lm == s, , drop = FALSE]
    g <- tryCatch(survival::coxph(f1, data = Ss, ties = "efron"),
                  error = function(e) NULL)
    for (m in mods) {
      est <- if (is.null(g)) NA_real_ else coef(g)[[m]]
      sem <- if (is.null(g)) NA_real_ else sqrt(diag(vcov(g)))[[m]]
      pts[[length(pts) + 1]] <- data.frame(
        variant = v$name, modality = toupper(m), window_days = v$w,
        landmark_day = s, at_risk = nrow(Ss), deaths = sum(Ss$status),
        hr = exp(est), lo = exp(est - 1.96 * sem),
        hi = exp(est + 1.96 * sem))
    }
  }
  list(curve = do.call(rbind, cur), points = do.call(rbind, pts),
       tests = do.call(rbind, tst))
}

cat("=== DATA ===\n")
spec <- list(
  list(name = "main", cohort = "EHR+ECG", scaling = "winsorised",
       hosp = FALSE, step = 15, w = 30),
  list(name = "window 60 days", cohort = "EHR+ECG",
       scaling = "winsorised", hosp = FALSE, step = 15, w = 60),
  list(name = "landmarks every 30 days", cohort = "EHR+ECG",
       scaling = "winsorised", hosp = FALSE, step = 30, w = 30),
  list(name = "registered scaling", cohort = "EHR+ECG",
       scaling = "standardised", hosp = FALSE, step = 15, w = 30),
  list(name = "hospice excluded", cohort = "EHR+ECG",
       scaling = "winsorised", hosp = TRUE, step = 15, w = 30),
  list(name = "three-modality cohort", cohort = "3-mod",
       scaling = "winsorised", hosp = FALSE, step = 15, w = 30))
VARS <- lapply(spec, function(sp) {
  dd <- make_data(sp$cohort, sp$scaling, sp$hosp)
  c(sp, dd)
})
for (v in VARS) {
  cat(sprintf("  %-24s n=%d deaths=%d  window %d d, landmarks every %d d\n",
              v$name, nrow(v$d), sum(v$d$event_death), v$w, v$step))
}

cat("\n=== TIMING PROBE (main variant, real data) ===\n")
t0 <- Sys.time()
invisible(run_lm(VARS[[1]]))
t1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  main variant %.1f s; all six on %d workers about %.1f min\n",
            t1, workers, t1 * 6 / min(6, workers) * 1.5 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

plan(multisession, workers = min(workers, length(VARS)))
R <- future_lapply(VARS, run_lm, future.packages = PKG,
                   future.seed = TRUE)
plan(sequential)
curve <- rbindlist(lapply(R, function(r) r$curve))
points <- rbindlist(lapply(R, function(r) r$points))
tests <- rbindlist(lapply(R, function(r) r$tests))
fwrite(curve, file.path(RS, "lm_curve.csv"))
fwrite(points, file.path(RS, "lm_points.csv"))
fwrite(tests, file.path(RS, "lm_tests.csv"))

cat("\n=== MAIN: HR per SD for death in the next 30 days ===\n")
print(curve[variant == "main" & landmark_day %in% c(0, 30, 90, 180, 330),
            .(modality, landmark_day, hr = round(hr, 3), lo = round(lo, 3),
              hi = round(hi, 3))], row.names = FALSE)
cat("\n=== TESTS (all variants) ===\n")
print(tests[, .(variant, modality, wald_2df = round(wald_2df, 2),
                p_value = signif(p_value, 3),
                ratio_s0_s180 = round(ratio_hr_s0_over_s180, 3),
                lo = round(ratio_lo, 3), hi = round(ratio_hi, 3))],
      row.names = FALSE)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("lm_curve.csv", "lm_points.csv", "lm_tests.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "lm_log.txt"), " (terminal log)\n")
cat("\nOutputs: lm_curve.csv, lm_points.csv, lm_tests.csv\n")