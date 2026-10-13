# 16c_spesi.R
# Simplified PESI (sPESI) as a benchmark (PREREGISTRATION.md sensitivity
# list; PREREGISTRATION_ADDENDUM_2.md A6; DEVIATIONS.md D7).
#
# sPESI (0 to 6) scores age > 80, cancer, heart failure or chronic lung
# disease (COPD; one point), heart rate >= 110, systolic BP < 100 and SpO2
# < 90%. Age and comorbidities come from ehr28.csv and vital signs from
# ehr_extra.csv (script 03b). A missing item counts as absent. Three
# versions (D7): presentation vitals (main), using the first value at
# presentation and otherwise the value averaged as in ehr28 (heart rate,
# systolic BP) or script 03b (SpO2); averaged vitals, using the averaged
# values for everyone; and presentation vitals, complete cases, restricted
# to admissions with all three presentation values, with every model
# re-evaluated on them.
#
# Cells: Frame A death at 30 days and Frame B death at 30 and 365 days, in
# both cohorts. Models: sPESI as a single-covariate Cox model fitted in the
# script 04 training folds (out of fold, seed-averaged); EHR alone (level-1
# ridge-Cox risk); and the Cox stack and binary stack (script 11b
# predictions).
#
# Measures: AUROC, Brier score, IPA, observed/expected, calibration slope
# and net benefit at 5%, 10% and 20%. For sPESI >= 1 as a rule:
# sensitivity, specificity, PPV, NPV, share and death rate of the low-risk
# group, and net benefit. Matched comparisons: each model's death rate in
# its lowest-risk group of sPESI's low-risk size, and its specificity at
# sPESI's sensitivity. Intervals are from a 2,000-sample bootstrap
# (patient-clustered in Frame A, seed 42) for every measure and every
# difference of a model minus sPESI. A timing probe runs a few bootstrap
# samples per cell first.
#
# Checks:  stops if an input column is missing, an ehr28 admission has no
#          ehr_extra row, median SpO2 is not on a percent scale, or an
#          sPESI value or prediction is missing.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          level1/oof_linear_*.csv, calibration/cal_preds_*.csv;
#          data/inputs/ehr28.csv, ehr_extra.csv
# Outputs: results/spesi_distribution.csv, results/spesi_items.csv,
#          results/spesi_metrics.csv, results/spesi_contrasts.csv,
#          results/spesi_checks.csv, results/spesi_log.txt
# Usage:   Rscript scripts/16c_spesi.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
NBOOT <- 2000L
SEED <- 42L
TH <- c(0.05, 0.10, 0.20)
CELLS <- list(
  list(fr = "A", co = "EHR+ECG", h = 30),
  list(fr = "A", co = "3-mod", h = 30),
  list(fr = "B", co = "EHR+ECG", h = 30),
  list(fr = "B", co = "EHR+ECG", h = 365),
  list(fr = "B", co = "3-mod", h = 30),
  list(fr = "B", co = "3-mod", h = 365))
VERS <- c(main = "presentation vitals (main)", avg = "averaged vitals",
          cc = "presentation vitals, complete cases")
PKG <- c("data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-60s %s\n", status, item, value))
}

# ---------------- sPESI ----------------
cat("=== sPESI INPUTS ===\n")
ehr <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
ext <- unique(fread(file.path(IN, "ehr_extra.csv")), by = "hadm_id")
n28 <- c("age", "cancer", "heart_failure", "copd", "hr", "sbp")
nx <- c("hr_first", "sbp_first", "spo2_first", "spo2")
miss <- c(setdiff(n28, names(ehr)), setdiff(nx, names(ext)))
if (length(miss)) {
  chk("columns missing", paste(miss, collapse = ", "), "FAIL")
  fwrite(rbindlist(CHK), file.path(RS, "spesi_checks.csv"))
  quit(status = 1)
}
X <- merge(ehr[, c("hadm_id", n28), with = FALSE],
           ext[, c("hadm_id", nx), with = FALSE],
           by = "hadm_id", all.x = TRUE)
chk("ehr28 admissions without an ehr_extra row",
    sum(is.na(match(ehr$hadm_id, ext$hadm_id))),
    if (all(ehr$hadm_id %in% ext$hadm_id)) "PASS" else "FAIL")
chk("median SpO2 at presentation / averaged (percent scale)",
    sprintf("%.1f / %.1f", median(X$spo2_first, na.rm = TRUE),
            median(X$spo2, na.rm = TRUE)),
    if (median(X$spo2, na.rm = TRUE) > 50) "PASS" else "FAIL")
ind <- function(v) ifelse(is.na(v), 0L, as.integer(v))
X[, `:=`(
  i_age80 = ind(age > 80), i_cancer = ind(cancer > 0),
  i_cardiopulm = pmax(ind(heart_failure > 0), ind(copd > 0)),
  hr_p = fifelse(is.na(hr_first), hr, hr_first),
  sbp_p = fifelse(is.na(sbp_first), sbp, sbp_first),
  spo2_p = fifelse(is.na(spo2_first), spo2, spo2_first))]
X[, base := i_age80 + i_cancer + i_cardiopulm]
X[, `:=`(
  i_hr_main = ind(hr_p >= 110), i_sbp_main = ind(sbp_p < 100),
  i_spo2_main = ind(spo2_p < 90),
  i_hr_avg = ind(hr >= 110), i_sbp_avg = ind(sbp < 100),
  i_spo2_avg = ind(spo2 < 90),
  i_hr_cc = ind(hr_first >= 110), i_sbp_cc = ind(sbp_first < 100),
  i_spo2_cc = ind(spo2_first < 90),
  complete_pres = !is.na(hr_first) & !is.na(sbp_first) &
    !is.na(spo2_first))]
X[, `:=`(spesi_main = base + i_hr_main + i_sbp_main + i_spo2_main,
         spesi_avg = base + i_hr_avg + i_sbp_avg + i_spo2_avg,
         spesi_cc = base + i_hr_cc + i_sbp_cc + i_spo2_cc)]
chk("admissions with all three presentation vitals",
    sprintf("%.3f", mean(X$complete_pres)), "INFO")

# ---------------- statistics (bootstrap-safe) ----------------
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  nbf <- function(y, tr, t) {
    (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / length(y)
  }
  y <- J$y[idx]
  prev <- mean(y)
  v <- numeric(0)
  for (m in names(J$P)) {
    p <- J$P[[m]][idx]
    br <- mean((y - p)^2)
    v[paste(m, "auc", sep = "|")] <- auc(y, p)
    v[paste(m, "brier", sep = "|")] <- br
    v[paste(m, "ipa", sep = "|")] <- 1 - br / (prev * (1 - prev))
    v[paste(m, "oe", sep = "|")] <- prev / mean(p)
    for (t in J$th) {
      v[paste(m, sprintf("nb%02d", round(100 * t)), sep = "|")] <-
        nbf(y, p >= t, t)
    }
  }
  hi <- J$sp[idx] >= 1
  rl <- "sPESI >= 1 rule"
  for (t in J$th) {
    v[paste(rl, sprintf("nb%02d", round(100 * t)), sep = "|")] <-
      nbf(y, hi, t)
  }
  sens <- sum(hi & y == 1) / sum(y == 1)
  v[paste(rl, "sensitivity", sep = "|")] <- sens
  v[paste(rl, "specificity", sep = "|")] <- sum(!hi & y == 0) / sum(y == 0)
  v[paste(rl, "ppv", sep = "|")] <- sum(hi & y == 1) / sum(hi)
  v[paste(rl, "npv", sep = "|")] <- sum(!hi & y == 0) / sum(!hi)
  v[paste(rl, "share_low_risk", sep = "|")] <- mean(!hi)
  v[paste(rl, "death_rate_low_risk", sep = "|")] <- mean(y[!hi])
  n0 <- sum(!hi)
  for (m in J$cmp) {
    p <- J$P[[m]][idx]
    v[paste(m, "death_rate_low_risk", sep = "|")] <-
      if (n0 > 0) mean(y[order(p)[seq_len(n0)]]) else NA_real_
    k_ <- paste(m, "specificity_at_spesi_sensitivity", sep = "|")
    if (sum(y == 1) == 0 || is.na(sens)) {
      v[k_] <- NA_real_
      next
    }
    o <- order(-p)
    cs <- cumsum(y[o] == 1) / sum(y == 1)
    kk <- which(cs >= sens - 1e-12)[1]
    flag <- rep(FALSE, length(p))
    flag[o[seq_len(kk)]] <- TRUE
    v[k_] <- sum(!flag & y == 0) / sum(y == 0)
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
avg <- function(path, cols) {
  x <- fread(path, select = c("hadm_id", cols))
  x[, lapply(.SD, mean), by = hadm_id, .SDcols = cols]
}

# ---------------- build cells ----------------
cat("\n=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
JOBS <- list()
DIST <- list()
ITEMS <- list()
for (cl in CELLS) {
  code <- if (cl$co == "EHR+ECG") "ehrecg" else "3mod"
  Fd <- frames[[cl$fr]]
  sub0 <- if (cl$co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  cut <- cl$h + 0.5
  fo <- folds[frame == cl$fr & cohort == cl$co]
  seeds <- sort(unique(fo$seed))
  l1all <- avg(file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv",
                                               cl$fr, code)),
               sprintf("ehr.cox.death.r%d", cl$h))
  cc <- sprintf("death|%s|%d", c("coxstack", "binstack"), cl$h)
  cpall <- avg(file.path(PR, "calibration",
                         sprintf("cal_preds_%s_%s.csv", cl$fr, code)), cc)
  for (vk in names(VERS)) {
    xi <- match(sub0$hadm_id, X$hadm_id)
    sub <- sub0
    if (vk == "cc") {
      keep <- !is.na(xi) & X$complete_pres[xi]
      sub <- sub0[keep]
      xi <- xi[keep]
    }
    sp <- X[[paste0("spesi_", vk)]][xi]
    nmiss <- sum(is.na(sp))
    sp[is.na(sp)] <- 0L
    y <- as.integer(sub$event_death == 1 & sub$time_death <= cut)
    R <- matrix(NA_real_, nrow(sub), length(seeds))
    for (si in seq_along(seeds)) {
      fs <- fo[fo$seed == seeds[si]][match(sub$hadm_id, hadm_id), fold]
      for (k in 0:4) {
        tr <- fs != k
        te <- fs == k
        d <- data.frame(tt = sub$time_death[tr], ee = sub$event_death[tr],
                        sp = sp[tr])
        f <- coxph(Surv(tt, ee) ~ sp, data = d, ties = "efron")
        bh <- basehaz(f, centered = FALSE)
        i <- findInterval(cut, bh$time)
        H0 <- if (i == 0) 0 else bh$hazard[i]
        R[te, si] <- 1 - exp(-H0 * exp(coef(f)[[1]] * sp[te]))
      }
    }
    l1 <- l1all[match(sub$hadm_id, hadm_id)]
    cp_ <- cpall[match(sub$hadm_id, hadm_id)]
    P <- list("sPESI (Cox)" = rowMeans(R),
              "EHR alone" = l1[[2]],
              "Cox stack (fused)" = cp_[[cc[1]]],
              "binary stack (fused)" = cp_[[cc[2]]])
    key <- sprintf("%s %s %dd | %s", cl$fr, cl$co, cl$h, VERS[[vk]])
    nna <- sum(sapply(P, function(p) sum(is.na(p))))
    chk(sprintf("%s: missing sPESI / predictions", key),
        sprintf("%d / %d", nmiss, nna),
        if (nmiss == 0 && nna == 0) "PASS" else "FAIL")
    for (sc in 0:6) {
      ii <- sp == sc
      DIST[[length(DIST) + 1]] <- data.table(
        frame = cl$fr, cohort = cl$co, horizon = cl$h,
        version = VERS[[vk]], spesi = sc, n = sum(ii),
        deaths = sum(y[ii]),
        death_rate = if (any(ii)) mean(y[ii]) else NA_real_)
    }
    Xs <- X[xi]
    src <- if (vk == "main") {
      c(hr = mean(!is.na(Xs$hr_first)), sbp = mean(!is.na(Xs$sbp_first)),
        spo2 = mean(!is.na(Xs$spo2_first)))
    } else if (vk == "cc") c(hr = 1, sbp = 1, spo2 = 1) else
      c(hr = 0, sbp = 0, spo2 = 0)
    ITEMS[[length(ITEMS) + 1]] <- data.table(
      frame = cl$fr, cohort = cl$co, horizon = cl$h, version = VERS[[vk]],
      n = nrow(sub), deaths = sum(y),
      age_over_80 = mean(Xs$i_age80), cancer = mean(Xs$i_cancer),
      heart_failure_or_copd = mean(Xs$i_cardiopulm),
      hr_110_or_more = mean(Xs[[paste0("i_hr_", vk)]]),
      sbp_under_100 = mean(Xs[[paste0("i_sbp_", vk)]]),
      spo2_under_90 = mean(Xs[[paste0("i_spo2_", vk)]]),
      share_hr_from_presentation = src[["hr"]],
      share_sbp_from_presentation = src[["sbp"]],
      share_spo2_from_presentation = src[["spo2"]],
      share_spesi_0 = mean(sp == 0))
    JOBS[[key]] <- list(
      key = key, fr = cl$fr, co = cl$co, h = cl$h, version = VERS[[vk]],
      n = nrow(sub), y = y, sp = sp, P = P, th = TH,
      cmp = c("EHR alone", "Cox stack (fused)", "binary stack (fused)"),
      rows_by = split(seq_len(nrow(sub)),
                      match(sub$subject_id, unique(sub$subject_id))))
    cat(sprintf("  %-52s n=%d deaths=%d  sPESI 0: %.1f%%\n", key,
                nrow(sub), sum(y), 100 * mean(sp == 0)))
  }
}
fwrite(rbindlist(DIST), file.path(RS, "spesi_distribution.csv"))
fwrite(rbindlist(ITEMS), file.path(RS, "spesi_items.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "spesi_checks.csv"))
if (any(ck$status == "FAIL")) {
  cat("\nSTOP: a check failed. see: spesi_checks.csv\n")
  quit(status = 1)
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
cat(sprintf("  %d cells x %d samples: about %.1f min\n", length(JOBS),
            NBOOT, proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- bootstrap and summaries ----------------
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cat(sprintf("=== %s: %d bootstrap samples ===\n", J$key, NBOOT))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  pt <- stat_vec(J, seq_len(J$n))
  lab <- data.table(frame = J$fr, cohort = J$co, horizon = J$h,
                    version = J$version)
  for (k in names(pt)) {
    pz <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- cbind(lab, data.table(
      model = pz[1], measure = pz[2], estimate = pt[[k]], lo = q[1],
      hi = q[2]))
  }
  for (m in names(J$P)) {
    p <- pmin(pmax(J$P[[m]], 1e-6), 1 - 1e-6)
    sl <- tryCatch(unname(coef(glm(J$y ~ qlogis(p),
                                   family = binomial))[2]),
                   error = function(e) NA_real_)
    MET[[length(MET) + 1]] <- cbind(lab, data.table(
      model = m, measure = "calibration slope", estimate = sl,
      lo = NA_real_, hi = NA_real_))
  }
  pairs <- list()
  for (m in J$cmp) {
    for (ms in c("auc", "ipa", sprintf("nb%02d", round(100 * TH)))) {
      pairs[[length(pairs) + 1]] <- c(m, "sPESI (Cox)", ms)
    }
    for (ms in sprintf("nb%02d", round(100 * TH))) {
      pairs[[length(pairs) + 1]] <- c(m, "sPESI >= 1 rule", ms)
    }
    pairs[[length(pairs) + 1]] <- c(m, "sPESI >= 1 rule",
                                    "death_rate_low_risk")
  }
  for (pr in pairs) {
    a <- paste(pr[1], pr[3], sep = "|")
    b <- paste(pr[2], pr[3], sep = "|")
    q <- q2(B[, a] - B[, b])
    CON[[length(CON) + 1]] <- cbind(lab, data.table(
      measure = pr[3], contrast = sprintf("%s minus %s", pr[1], pr[2]),
      estimate = pt[[a]] - pt[[b]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0))
  }
  for (m in J$cmp) {
    a <- paste(m, "specificity_at_spesi_sensitivity", sep = "|")
    b <- "sPESI >= 1 rule|specificity"
    q <- q2(B[, a] - B[, b])
    CON[[length(CON) + 1]] <- cbind(lab, data.table(
      measure = "specificity at sPESI's sensitivity",
      contrast = sprintf("%s minus sPESI >= 1 rule", m),
      estimate = pt[[a]] - pt[[b]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0))
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
fwrite(mt, file.path(RS, "spesi_metrics.csv"))
fwrite(cn, file.path(RS, "spesi_contrasts.csv"))

cat("\n=== AUROC, MAIN VERSION ===\n")
print(dcast(mt[measure == "auc" & version == VERS[["main"]]],
            frame + cohort + horizon ~ model, value.var = "estimate"),
      row.names = FALSE)
cat("\n=== KEY CONTRASTS, MAIN VERSION ===\n")
print(cn[version == VERS[["main"]] &
           measure %in% c("auc", "death_rate_low_risk",
                          "specificity at sPESI's sensitivity") &
           grepl("^Cox stack", contrast),
         .(frame, cohort, horizon, measure, contrast,
           estimate = round(estimate, 4), lo = round(lo, 4),
           hi = round(hi, 4), excludes_zero)], row.names = FALSE)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("spesi_distribution.csv", "spesi_items.csv",
            "spesi_metrics.csv", "spesi_contrasts.csv",
            "spesi_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "spesi_log.txt"), " (terminal log)\n")
cat("\nOutputs: spesi_distribution.csv, spesi_items.csv,",
    "spesi_metrics.csv, spesi_contrasts.csv, spesi_checks.csv\n")