# 22a_subgroups.R
# Performance of the final models by age, sex and recorded race (A29).
#
# Frame B (death to 365 days after discharge), EHR+ECG cohort: the final
# discharge model, read from the sl7_decision.csv of 18o (A30's rule;
# predictions from sl7_preds_B_ehrecg.csv). Frame A (death by day 30 of the
# admission), EHR+ECG cohort: Super learner v4, recalibrated
# (sl5_preds_A_ehrecg.csv). Comparators are Late PH (the registered late
# fusion, arch_preds_linear of 18a) and EHR alone (the spline Cox stack,
# from the same super learner file). Version 2 (A30): Super learner v7
# replaces v6 as the final Frame B model.
#
# Subgroups are age at the index admission (the EHR feature age; under 65,
# 65-79, 80 and over), sex (MIMIC patients.gender) and recorded race in six
# groups (03c: White, Black, Hispanic or Latino, Asian, Other, Unknown).
# Race is not a model input.
#
# Per subgroup and horizon: patients, deaths, AUROC, IPA, observed/expected
# and calibration slope of the final model; the final model minus Late PH
# and minus EHR alone (AUROC); and the subgroup minus the rest of the cohort
# (AUROC and observed/expected of the final model). Estimates need at least
# 10 deaths in the subgroup at that horizon, otherwise they are reported as
# too few events. Intervals come from a 2,000-sample bootstrap of the whole
# cohort (seed 42; patient-clustered in Frame A) with every subgroup
# recomputed in each resample, so between-group differences have intervals.
#
# The script stops if inputs are missing, if the final model or a comparator
# is missing for any patient and horizon, or if demographics are found for
# fewer than 99% of patients.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv,
#          super_learner/sl7_preds_B_ehrecg.csv, sl5_preds_A_ehrecg.csv,
#          fusion_arch/arch_preds_linear_*_ehrecg.csv;
#          results/sl7_decision.csv; data/inputs/ehr28.csv,
#          demographics.csv
# Outputs: results/subgroup_{metrics,contrasts,checks}.csv,
#          results/subgroup_log.txt
# Usage:   Rscript scripts/22a_subgroups.R [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
OD <- file.path(PR, "super_learner")
RS <- file.path(P2, "results")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MINEV <- 10L
NBOOT <- 2000L
SEED <- 42L
PKG <- c("data.table")
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
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "subgroup_checks.csv"))
    cat("\nSTOP: a check failed. see: subgroup_checks.csv\n")
    quit(status = 1)
  }
}
need <- c(file.path(PR, c("frame_A.csv", "frame_B.csv")),
          file.path(OD, c("sl7_preds_B_ehrecg.csv",
                          "sl5_preds_A_ehrecg.csv")),
          file.path(PR, "fusion_arch", c("arch_preds_linear_A_ehrecg.csv",
                                         "arch_preds_linear_B_ehrecg.csv")),
          file.path(RS, "sl7_decision.csv"),
          file.path(IN, c("ehr28.csv", "demographics.csv")))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
dec <- fread(file.path(RS, "sl7_decision.csv"))
FINAL_B <- dec$final_model_frame_B[1]
chk("final Frame B model (from stage 18o)", FINAL_B, "INFO")
dem <- unique(fread(file.path(IN, "demographics.csv")), by = "hadm_id")
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
if (!("age" %in% names(e28))) {
  chk("ehr28.csv has an age column", "no", "FAIL")
  stop_if_fail()
}

# ---------------- statistics for one resample ----------------
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (gv in names(J$G)) {
    gr <- J$G[[gv]][idx]
    for (lv in J$levels[[gv]]) {
      inn <- which(gr == lv)
      out <- which(gr != lv)
      for (hi in seq_along(J$hz)) {
        y <- J$Y[idx, hi]
        key <- paste(gv, lv, J$hz[hi], sep = "|")
        ok <- sum(y[inn]) >= J$minev && sum(1 - y[inn]) >= J$minev
        for (m in names(J$P)) {
          p <- pmin(pmax(J$P[[m]][idx, hi], 1e-6), 1 - 1e-6)
          v[paste(key, m, "auc", sep = "|")] <- if (ok) {
            auc(y[inn], p[inn])
          } else NA_real_
        }
        p <- pmin(pmax(J$P[["final"]][idx, hi], 1e-6), 1 - 1e-6)
        pr <- mean(y[inn])
        v[paste(key, "final", "oe", sep = "|")] <- if (ok) {
          pr / mean(p[inn])
        } else NA_real_
        v[paste(key, "final", "ipa", sep = "|")] <- if (ok) {
          1 - mean((y[inn] - p[inn])^2) / (pr * (1 - pr))
        } else NA_real_
        okr <- sum(y[out]) >= J$minev && sum(1 - y[out]) >= J$minev
        v[paste(key, "rest", "auc", sep = "|")] <- if (okr) {
          auc(y[out], p[out])
        } else NA_real_
        v[paste(key, "rest", "oe", sep = "|")] <- if (okr) {
          mean(y[out]) / mean(p[out])
        } else NA_real_
      }
    }
  }
  v
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
getP <- function(D, mdl, st, sub, hz) {
  x <- D[D$model == mdl & D$set == st]
  if (!nrow(x)) return(NULL)
  matrix(sapply(hz, function(h) {
    y <- x[horizon == h]
    y$risk[match(sub$hadm_id, y$hadm_id)]
  }), ncol = length(hz))
}
JOBS <- list()
for (fr in c("B", "A")) {
  sub <- fread(file.path(PR, sprintf("frame_%s.csv", fr)))
  hz <- HZ[[fr]]
  f <- if (fr == "B") file.path(OD, "sl7_preds_B_ehrecg.csv") else
    file.path(OD, "sl5_preds_A_ehrecg.csv")
  fm <- if (fr == "B") FINAL_B else "Super learner v4, recalibrated"
  D <- fread(f)
  D18 <- fread(file.path(PR, "fusion_arch",
                         sprintf("arch_preds_linear_%s_ehrecg.csv", fr)))
  P <- list(final = getP(D, fm, "EHR+ECG", sub, hz),
            late = getP(D18, "Late PH", "EHR+ECG", sub, hz),
            ehr = getP(D, "SL candidate: EHR alone, spline Cox", "EHR",
                       sub, hz))
  rm(D, D18)
  nb <- sum(sapply(P, function(m) if (is.null(m)) Inf else
    sum(!is.finite(m))))
  chk(sprintf("Frame %s: %s, Late PH and EHR alone present (missing)", fr,
              fm), nb, if (nb == 0) "PASS" else "FAIL")
  stop_if_fail()
  d <- dem[match(sub$hadm_id, hadm_id)]
  fd <- mean(!is.na(d$hadm_id))
  chk(sprintf("Frame %s: demographics found", fr), sprintf("%.1f%%",
                                                          100 * fd),
      if (fd >= 0.99) "PASS" else "FAIL")
  stop_if_fail()
  age <- e28$age[match(sub$hadm_id, e28$hadm_id)]
  G <- list(
    age = ifelse(is.na(age), "missing", ifelse(age < 65, "under 65",
                                               ifelse(age < 80, "65-79",
                                                      "80 and over"))),
    sex = ifelse(toupper(d$gender) == "F", "female",
                 ifelse(toupper(d$gender) == "M", "male", "missing")),
    race = ifelse(is.na(d$race_group), "Unknown", d$race_group))
  lev <- list(age = c("under 65", "65-79", "80 and over"),
              sex = c("female", "male"),
              race = c("White", "Black", "Hispanic or Latino", "Asian",
                       "Other", "Unknown"))
  for (gv in names(lev)) lev[[gv]] <- intersect(lev[[gv]], unique(G[[gv]]))
  Y <- sapply(hz, function(h) as.integer(sub$event_death == 1 &
                                           sub$time_death <= h + 0.5))
  Y <- matrix(Y, ncol = length(hz))
  for (gv in names(lev)) {
    sz <- sapply(lev[[gv]], function(l) sum(G[[gv]] == l))
    ev <- sapply(lev[[gv]], function(l) sum(Y[G[[gv]] == l, length(hz)]))
    chk(sprintf("Frame %s %s: patients (deaths by %d days)", fr, gv,
                max(hz)), paste(sprintf("%s %d (%d)", names(sz), sz, ev),
                                collapse = "; "), "INFO")
  }
  JOBS[[fr]] <- list(fr = fr, hz = hz, P = P, Y = Y, G = G, levels = lev,
                     minev = MINEV, model = fm,
                     rows_by = split(seq_len(nrow(sub)),
                                     match(sub$subject_id,
                                           unique(sub$subject_id))))
}
boot_rep <- function(i, J) {
  idx <- if (J$fr == "A") {
    unlist(J$rows_by[sample.int(length(J$rows_by), length(J$rows_by),
                                replace = TRUE)], use.names = FALSE)
  } else {
    sample.int(nrow(J$Y), nrow(J$Y), replace = TRUE)
  }
  stat_vec(J, idx)
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
cat(sprintf("  projected: about %.1f min\n", proj / 60))
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
  cat(sprintf("  Frame %s done\n", J$fr))
  for (gv in names(J$G)) for (lv in J$levels[[gv]]) {
    inn <- J$G[[gv]] == lv
    for (hi in seq_along(J$hz)) {
      key <- paste(gv, lv, J$hz[hi], sep = "|")
      y <- J$Y[inn, hi]
      p <- pmin(pmax(J$P$final[inn, hi], 1e-6), 1 - 1e-6)
      sl <- if (sum(y) >= MINEV) {
        tryCatch(unname(coef(glm(y ~ qlogis(p), family = binomial))[2]),
                 error = function(e) NA_real_)
      } else NA_real_
      for (ms in c("auc", "ipa", "oe")) {
        k <- paste(key, "final", ms, sep = "|")
        q <- q2(B[, k])
        MET[[length(MET) + 1]] <- data.table(
          frame = J$fr, model = J$model, variable = gv, group = lv,
          horizon = J$hz[hi], patients = sum(inn), deaths = sum(y),
          measure = ms, estimate = pt[[k]], lo = q[1], hi = q[2])
      }
      MET[[length(MET) + 1]] <- data.table(
        frame = J$fr, model = J$model, variable = gv, group = lv,
        horizon = J$hz[hi], patients = sum(inn), deaths = sum(y),
        measure = "slope", estimate = sl, lo = NA_real_, hi = NA_real_)
      for (cp in list(c("final", "late", "final model - Late PH"),
                      c("final", "ehr", "final model - EHR alone"),
                      c("final", "rest", "group - rest of cohort"))) {
        ka <- paste(key, cp[1], "auc", sep = "|")
        kb <- paste(key, cp[2], "auc", sep = "|")
        q <- q2(B[, ka] - B[, kb])
        CON[[length(CON) + 1]] <- data.table(
          frame = J$fr, model = J$model, variable = gv, group = lv,
          horizon = J$hz[hi], contrast = cp[3], measure = "auc",
          estimate = pt[[ka]] - pt[[kb]], lo = q[1], hi = q[2],
          excludes_zero = q[1] > 0 | q[2] < 0)
      }
      ka <- paste(key, "final", "oe", sep = "|")
      kb <- paste(key, "rest", "oe", sep = "|")
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, model = J$model, variable = gv, group = lv,
        horizon = J$hz[hi], contrast = "group - rest of cohort",
        measure = "oe", estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
cn[!is.finite(estimate), excludes_zero := NA]
fwrite(mt, file.path(RS, "subgroup_metrics.csv"))
fwrite(cn, file.path(RS, "subgroup_contrasts.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "subgroup_checks.csv"))
cat("\n=== FINAL MODEL AUROC BY SUBGROUP ===\n")
x <- mt[measure == "auc"]
print(dcast(x, frame + variable + group ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE, nrows = 200)
cat("\n=== GROUP MINUS REST OF COHORT, significant cells ===\n")
print(cn[contrast == "group - rest of cohort" & excludes_zero == TRUE,
         .(frame, variable, group, horizon, measure,
           est = round(estimate, 3), lo = round(lo, 3), hi = round(hi, 3))],
      row.names = FALSE, nrows = 200)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("subgroup_metrics.csv", "subgroup_contrasts.csv",
            "subgroup_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "subgroup_log.txt"), " (terminal log)\n")
cat("\nOutputs: subgroup_metrics.csv, subgroup_contrasts.csv,",
    "subgroup_checks.csv\n")