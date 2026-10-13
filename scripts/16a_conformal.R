# 16a_conformal.R
# Conformal lower predictive bounds on survival time
# (PREREGISTRATION_ADDENDUM.md, A5; DEVIATIONS.md D8). Exploratory;
# intervals only.
#
# Version 2 (D8): a survivor's CDR score is alpha - 1, because
# Y = min(T, 366) has a point mass at 366 where its CDF is 1 (version 1 used
# alpha - F(366 | x), valid but far too conservative). Only CDR and the two
# Mondrian CDR methods change; CQR, the uncalibrated Cox quantile and the
# covariate-free bound use the same random splits and must reproduce version
# 1 exactly (results/conformal_v1/).
#
# Setting: Frame B, death within one year of discharge, one row per patient
# (Candes, Lei and Ren, JRSSB 2023). Every patient is followed to day 366
# (Type I administrative censoring), so Algorithm 1 applies with c0 = 366
# and all weights 1. The target is Y = min(T, 366), and a lower bound on Y
# is a lower bound on T. Base models: the Cox stack (level-2 Cox on the
# winsorised ridge-Cox scores of every modality) and an EHR-only Cox model,
# fitted on the four training folds of each seed, with a Breslow baseline.
# Inside each held-out fold, 10 random half splits each calibrate the other
# half (600 calibration/test splits per cohort and level).
#
# Scores: CDR (main), V = alpha - F(Y | x), bound = earliest time the
# predicted risk reaches alpha - eta; CQR, V = q_alpha(x) - Y, bound =
# q_alpha(x) - eta. Here eta is the ceil((1 - alpha)(n + 1))-th smallest
# score, and bounds are clipped to [0, 366]. Comparators: the Cox model's
# own alpha-quantile (uncalibrated), a covariate-free conformal bound, and
# Mondrian CDR within tertiles of predicted one-year risk (cut-points from
# the training folds) and within age bands. Levels: alpha = 0.10 (main),
# 0.05 and 0.20.
#
# Checks:  FAIL if a held-out patient is not predicted, a bound leaves
#          [0, 366] or an unchanged method differs from version 1; WARN if
#          mean coverage is below 1 - alpha - 0.005.
#
# Inputs:  data/processed/frame_B.csv,
#          data/processed/level1/oof_linear_B_*.csv, data/inputs/ehr28.csv,
#          results/conformal_v1/conformal_summary.csv
# Outputs: results/conformal_summary.csv, conformal_groups.csv,
#          conformal_contrast.csv, conformal_checks.csv, conformal_log.txt
# Usage:   Rscript scripts/16a_conformal.R [--probe-only]

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
CDIR <- file.path(PR, "conformal")
dir.create(CDIR, showWarnings = FALSE, recursive = TRUE)
ALPHAS <- c(0.10, 0.05, 0.20)
AMAIN <- 0.10
NREP <- 10L
C0 <- 366
NBOOT <- 2000L
PKG <- c("survival")
CONFORMAL <- c("CDR", "CQR", "CDR, Mondrian by risk tertile",
               "CDR, Mondrian by age band", "conformal (no covariates)")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

# ---------------- one task = cohort, seed, held-out fold ----------------
conf_task <- function(tk) {
  c0 <- tk$c0
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
  qconf <- function(V, a) {
    n <- length(V)
    k <- ceiling((1 - a) * (n + 1))
    if (n == 0 || k > n) Inf else sort(V)[k]
  }
  first_reach <- function(Fm, g, thr) {
    idx <- rowSums(Fm < thr) + 1
    out <- ifelse(idx > length(g), c0, g[pmin(idx, length(g))])
    out[rep_len(thr, nrow(Fm)) <= 0] <- 0
    out
  }
  tr <- which(tk$fold != tk$k)
  te <- which(tk$fold == tk$k)
  n <- length(tk$y)
  nte <- length(te)
  models <- list()
  n_na <- 0L
  for (mn in names(tk$sets)) {
    mods <- tk$sets[[mn]]
    Z <- sapply(mods, function(m) wz(tk$LP[[m]], tr))
    if (is.null(dim(Z))) {
      Z <- matrix(Z, ncol = 1, dimnames = list(NULL, mods))
    }
    D <- data.frame(tt = tk$t, ee = tk$e, Z)
    f <- coxph(Surv(tt, ee) ~ ., data = D[tr, ], ties = "efron")
    bh <- basehaz(f, centered = FALSE)
    keep <- bh$time <= c0 - 0.5
    g <- bh$time[keep]
    H <- bh$hazard[keep]
    lp <- as.vector(as.matrix(D[, mods, drop = FALSE]) %*% coef(f)[mods])
    Fm <- 1 - exp(-outer(exp(lp[te]), H))
    n_na <- n_na + sum(is.na(Fm))
    iy <- findInterval(tk$y[te], g)
    Fy <- ifelse(iy == 0, 0, Fm[cbind(seq_len(nte), pmax(iy, 1))])
    # D8: Y = min(T, c0) has a point mass at c0, so its CDF there is 1
    Fy[tk$y[te] >= c0] <- 1
    r1_tr <- 1 - exp(-H[length(H)] * exp(lp[tr]))
    cuts <- quantile(r1_tr, c(1 / 3, 2 / 3), names = FALSE)
    r1_te <- Fm[, ncol(Fm)]
    models[[mn]] <- list(g = g, Fm = Fm, Fy = Fy, r1 = r1_te,
                         tert = findInterval(r1_te, cuts) + 1)
  }
  age_te <- tk$ageband[te]
  qa_all <- lapply(models, function(M) {
    setNames(lapply(tk$alphas, function(a) first_reach(M$Fm, M$g, a)),
             as.character(tk$alphas))
  })
  rows <- list()
  acc <- list()
  cnt <- numeric(n)
  lpb_min <- Inf
  lpb_max <- -Inf
  for (r in seq_len(tk$nrep)) {
    set.seed(tk$taskseed * 1000 + r)
    perm <- sample.int(nte)
    h1 <- perm[seq_len(floor(nte / 2))]
    h2 <- perm[-seq_len(floor(nte / 2))]
    for (dir in 1:2) {
      cal <- if (dir == 1) h1 else h2
      tst <- if (dir == 1) h2 else h1
      yc <- tk$y[te][cal]
      yt <- tk$y[te][tst]
      gidx <- te[tst]
      cnt[gidx] <- cnt[gidx] + 1
      for (a in tk$alphas) {
        L <- list(list(model = "covariate-free",
                       method = "conformal (no covariates)",
                       lpb = rep(-qconf(-yc, a), length(tst))))
        for (mn in names(models)) {
          M <- models[[mn]]
          eta <- qconf(a - M$Fy[cal], a)
          L[[length(L) + 1]] <- list(
            model = mn, method = "CDR",
            lpb = first_reach(M$Fm[tst, , drop = FALSE], M$g, a - eta))
          qa <- qa_all[[mn]][[as.character(a)]]
          etaq <- qconf(qa[cal] - yc, a)
          L[[length(L) + 1]] <- list(model = mn, method = "CQR",
                                     lpb = qa[tst] - etaq)
          L[[length(L) + 1]] <- list(model = mn,
                                     method = "Cox quantile, uncalibrated",
                                     lpb = qa[tst])
          for (gv in c("risk tertile", "age band")) {
            grp <- if (gv == "risk tertile") M$tert else age_te
            lpg <- numeric(length(tst))
            for (lv in unique(grp[tst])) {
              ic <- cal[grp[cal] == lv]
              it <- which(grp[tst] == lv)
              eg <- qconf(a - M$Fy[ic], a)
              lpg[it] <- first_reach(M$Fm[tst[it], , drop = FALSE], M$g,
                                     a - eg)
            }
            L[[length(L) + 1]] <- list(
              model = mn, method = paste("CDR, Mondrian by", gv),
              lpb = lpg)
          }
        }
        for (el in L) {
          lpb <- pmin(c0, pmax(0, el$lpb))
          lpb_min <- min(lpb_min, lpb)
          lpb_max <- max(lpb_max, lpb)
          cov <- yt >= lpb
          rows[[length(rows) + 1]] <- data.frame(
            model = el$model, method = el$method, alpha = a,
            seed = tk$seed, fold = tk$k, rep = r, dir = dir,
            n_cal = length(cal), n_test = length(tst),
            coverage = mean(cov), mean_lpb = mean(lpb),
            median_lpb = median(lpb), share_366 = mean(lpb >= c0),
            share_lt30 = mean(lpb < 30))
          if (a == tk$amain) {
            key <- paste(el$model, el$method, sep = "|")
            if (is.null(acc[[key]])) {
              acc[[key]] <- list(lpb = numeric(n), cov = numeric(n))
            }
            v <- acc[[key]]
            v$lpb[gidx] <- v$lpb[gidx] + lpb
            v$cov[gidx] <- v$cov[gidx] + cov
            acc[[key]] <- v
          }
        }
      }
    }
  }
  risk_sum <- numeric(n)
  risk_sum[te] <- models[[tk$fused]]$r1
  list(co = tk$co, rows = do.call(rbind, rows), acc = acc, cnt = cnt,
       risk_sum = risk_sum, risk_cnt = as.numeric(seq_len(n) %in% te),
       n_na = n_na, lpb_min = lpb_min, lpb_max = lpb_max)
}

# ---------------- build tasks ----------------
cat("=== BUILD ===\n")
fb <- fread(file.path(PR, "frame_B.csv"))
ehrm <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
nm <- names(ehrm)
agec <- nm[tolower(nm) %in% c("age", "age_at_admit", "anchor_age")][1]
cat(sprintf("  age column: %s\n", ifelse(is.na(agec), "NOT FOUND", agec)))
tasks <- list()
META <- list()
for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  fused <- paste(toupper(mods), collapse = " + ")
  sets <- setNames(list("ehr", mods), c("EHR only", fused))
  lpc <- sprintf("%s.cox.death.lp", mods)
  oof <- fread(file.path(PR, "level1", sprintf("oof_linear_B_%s.csv", code)),
               select = c("seed", "hadm_id", "fold", lpc))
  seeds <- sort(unique(oof$seed))
  a <- if (is.na(agec)) rep(NA_real_, nrow(sub)) else
    ehrm[match(sub$hadm_id, hadm_id)][[agec]]
  ageband <- ifelse(is.na(a), "age missing",
                    ifelse(a < 65, "<65", ifelse(a < 80, "65-79", "80+")))
  META[[co]] <- list(code = code, sub = sub, fused = fused,
                     ageband = ageband)
  for (s in seeds) {
    g <- oof[oof$seed == s][match(sub$hadm_id, hadm_id)]
    if (anyNA(g$fold) || anyNA(g[, ..lpc])) stop("missing score: ", co, s)
    for (k in 0:4) {
      tasks[[length(tasks) + 1]] <- list(
        co = co, seed = s, k = k, fold = g$fold,
        LP = setNames(lapply(lpc, function(cc) g[[cc]]), mods),
        t = sub$time_death, e = as.integer(sub$event_death == 1),
        y = pmin(sub$time_death, C0), ageband = ageband, sets = sets,
        fused = fused, alphas = ALPHAS, amain = AMAIN, nrep = NREP,
        c0 = C0, taskseed = s * 10 + k)
    }
  }
  cat(sprintf("  %-7s n=%d  one-year deaths=%d  seeds %d  models: %s\n",
              co, nrow(sub), sum(sub$event_death), length(seeds),
              paste(names(sets), collapse = " / ")))
}

cat("\n=== TIMING PROBE ===\n")
t0 <- Sys.time()
invisible(conf_task(tasks[[1]]))
t1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  one task %.1f s; %d tasks on %d workers: about %.1f min\n",
            t1, length(tasks), workers,
            t1 * length(tasks) / workers * 1.3 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, conf_task, future.packages = PKG,
                     future.seed = TRUE)
plan(sequential)

# ---------------- assemble ----------------
SUM <- list()
GRP <- list()
CON <- list()
CHK <- list()
chk <- function(co, item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(cohort = co, item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-4s %-7s %-60s %s\n", status, co, item, value))
}
cat("\n=== CHECKS ===\n")
for (co in names(META)) {
  M <- META[[co]]
  rr <- Filter(function(r) r$co == co, res)
  rows <- rbindlist(lapply(rr, function(r) r$rows))
  nna <- sum(sapply(rr, function(r) r$n_na))
  lo <- min(sapply(rr, function(r) r$lpb_min))
  hi <- max(sapply(rr, function(r) r$lpb_max))
  chk(co, "missing predicted risks", nna, if (nna == 0) "PASS" else "FAIL")
  chk(co, "bounds within [0, 366]", sprintf("%.1f to %.1f", lo, hi),
      if (lo >= 0 && hi <= C0) "PASS" else "FAIL")
  s <- rows[, .(n_splits = .N, coverage_mean = mean(coverage),
                coverage_p05 = quantile(coverage, 0.05, names = FALSE),
                coverage_p95 = quantile(coverage, 0.95, names = FALSE),
                mean_lpb = mean(mean_lpb), median_lpb = mean(median_lpb),
                share_366 = mean(share_366), share_lt30 = mean(share_lt30)),
            by = .(alpha, model, method)]
  s[, cohort := co]
  SUM[[co]] <- s
  for (i in which(s$method %in% CONFORMAL)) {
    ok <- s$coverage_mean[i] >= 1 - s$alpha[i] - 0.005
    chk(co, sprintf("coverage, alpha %.2f, %s, %s", s$alpha[i],
                    s$model[i], s$method[i]),
        sprintf("%.4f", s$coverage_mean[i]), if (ok) "PASS" else "WARN")
  }
  cnt <- Reduce(`+`, lapply(rr, function(r) r$cnt))
  risk <- Reduce(`+`, lapply(rr, function(r) r$risk_sum)) /
    Reduce(`+`, lapply(rr, function(r) r$risk_cnt))
  keys <- names(rr[[1]]$acc)
  P <- data.table(hadm_id = M$sub$hadm_id,
                  died_1y = as.integer(M$sub$event_death == 1),
                  risk_1y_fused = risk, age_band = M$ageband,
                  appearances = cnt)
  for (kk in keys) {
    P[[paste0("lpb|", kk)]] <-
      Reduce(`+`, lapply(rr, function(r) r$acc[[kk]]$lpb)) / cnt
    P[[paste0("cov|", kk)]] <-
      Reduce(`+`, lapply(rr, function(r) r$acc[[kk]]$cov)) / cnt
  }
  fwrite(P, file.path(CDIR, sprintf("lpb_patients_%s.csv", M$code)))
  dec <- ceiling(10 * frank(risk, ties.method = "first") / length(risk))
  ter <- ceiling(3 * frank(risk, ties.method = "first") / length(risk))
  groupings <- list(
    "fused risk decile" = paste("decile", dec),
    "fused risk tertile" = paste("tertile", ter),
    "age band" = M$ageband)
  for (kk in keys) {
    pz <- strsplit(kk, "|", fixed = TRUE)[[1]]
    lp <- P[[paste0("lpb|", kk)]]
    cv <- P[[paste0("cov|", kk)]]
    for (gn in names(groupings)) {
      gg <- groupings[[gn]]
      for (lv in sort(unique(gg))) {
        ii <- which(gg == lv)
        GRP[[length(GRP) + 1]] <- data.table(
          cohort = co, alpha = AMAIN, model = pz[1], method = pz[2],
          grouping = gn, level = lv, n = length(ii),
          deaths_1y = sum(P$died_1y[ii]),
          mean_risk_1y = mean(risk[ii]), coverage = mean(cv[ii]),
          mean_lpb = mean(lp[ii]), median_lpb = median(lp[ii]))
      }
    }
  }
  set.seed(42)
  for (mth in c("CDR", "CQR")) {
    a1 <- P[[paste0("lpb|", M$fused, "|", mth)]]
    a0 <- P[[paste0("lpb|EHR only|", mth)]]
    d <- a1 - a0
    bs <- replicate(NBOOT, mean(d[sample.int(length(d), replace = TRUE)]))
    q <- quantile(bs, c(0.025, 0.975), names = FALSE)
    CON[[length(CON) + 1]] <- data.table(
      cohort = co, alpha = AMAIN, method = mth,
      contrast = sprintf("%s minus EHR only", M$fused),
      mean_lpb_fused = mean(a1), mean_lpb_ehr = mean(a0),
      mean_difference_days = mean(d), lo = q[1], hi = q[2],
      share_patients_higher = mean(d > 0),
      share_patients_lower = mean(d < 0))
  }
}
sm <- rbindlist(SUM)
setcolorder(sm, c("cohort", "alpha", "model", "method"))
gp <- rbindlist(GRP)
cn <- rbindlist(CON)
OLD <- file.path(RS, "conformal_v1", "conformal_summary.csv")
SAME <- c("CQR", "Cox quantile, uncalibrated", "conformal (no covariates)")
if (file.exists(OLD)) {
  old <- fread(OLD)
  a <- merge(sm[method %in% SAME], old[method %in% SAME],
             by = c("cohort", "alpha", "model", "method"),
             suffixes = c("", ".v1"))
  d <- max(abs(a$coverage_mean - a$coverage_mean.v1),
           abs(a$mean_lpb - a$mean_lpb.v1))
  nexp <- nrow(sm[method %in% SAME])
  chk("both", "unchanged methods reproduce the first run (max diff)",
      sprintf("%.1e over %d of %d rows", d, nrow(a), nexp),
      if (nrow(a) == nexp && d < 1e-6) "PASS" else "FAIL")
} else {
  chk("both", "first run not found in results/conformal_v1", "skipped",
      "WARN")
}
ck <- rbindlist(CHK)
fwrite(sm, file.path(RS, "conformal_summary.csv"))
fwrite(gp, file.path(RS, "conformal_groups.csv"))
fwrite(cn, file.path(RS, "conformal_contrast.csv"))
fwrite(ck, file.path(RS, "conformal_checks.csv"))

cat("\n=== SUMMARY, alpha 0.10 ===\n")
print(sm[alpha == AMAIN, .(cohort, model, method,
                           coverage = round(coverage_mean, 3),
                           p05 = round(coverage_p05, 3),
                           mean_lpb = round(mean_lpb, 1),
                           share_366 = round(share_366, 3),
                           share_lt30 = round(share_lt30, 3))],
      row.names = FALSE)
cat("\n=== FUSED MINUS EHR ONLY (days) ===\n")
print(cn[, .(cohort, method, diff = round(mean_difference_days, 1),
             lo = round(lo, 1), hi = round(hi, 1))], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("conformal_summary.csv", "conformal_groups.csv",
            "conformal_contrast.csv", "conformal_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", CDIR, " (patient-level bounds, local only)\n")
cat("  ", file.path(RS, "conformal_log.txt"), " (terminal log)\n")
cat("\nOutputs: conformal_summary.csv, conformal_groups.csv,",
    "conformal_contrast.csv, conformal_checks.csv\n")