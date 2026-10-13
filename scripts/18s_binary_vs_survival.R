# 18s_binary_vs_survival.R
# Matched binary super learners, and their comparison with the survival
# super learners (A33).
#
# Targets are death, the composite, CV readmission and a dissertation-style
# CV target (cv_excl), at diagnosis and at discharge, in the primary cohort
# (EHR+ECG). Level 1 is the out-of-fold binary predictions of 06r, per
# target, horizon, modality and seed: ridge logistic and probability forest.
# Per target and horizon the candidates are the prevalence, EHR alone for
# each family (logistic regression on the normal-scored prediction) and
# EHR+ECG for each family (logistic stack on both modalities' normal-scored
# predictions). Non-negative weights summing to one come from grouped inner
# 3-fold out-of-fold predictions within each training fold, and the result
# is recalibrated by logistic regression on the logit. SL-EHR uses the
# EHR-alone candidates only, and seeds are averaged. In cv_excl, patients
# who died first by the horizon are left out of that horizon's fitting and
# evaluation.
#
# Survival counterparts (same patients, EHR+ECG, recalibrated): at
# diagnosis, death Super learner P (18q) and the composite and CV super
# learners of 18r; at discharge, death v7-core (the matched four families;
# v7 is also reported) and the same 18r super learners. For cv_excl the
# survival CV super learner is evaluated among patients who did not die
# first, on AUROC only, since its probabilities refer to readmission before
# death.
#
# AUROC, IPA and calibration slope, with a 2,000-sample patient-clustered
# bootstrap (seed 42). Contrasts: survival minus binary (AUROC and IPA;
# AUROC only for cv_excl) and binary SL minus binary SL-EHR (AUROC). A33
# expected the two within 0.01 AUROC at most horizons. The script stops on
# missing inputs, misaligned seeds, rows or folds, predictions outside [0,
# 1], or missing survival predictions.
#
# Inputs:  data/processed/frame_B.csv, presentation/frame_P_outcomes.csv,
#          level1/oof_bin_*.csv, super_learner/slP_preds_ehrecg.csv,
#          optB2_preds_*_ehrecg.csv, sl7_preds_B_ehrecg.csv (all under
#          data/processed/)
# Outputs: results/binsurv_{summary,contrasts,metrics,weights,checks}.csv,
#          results/binsurv_log.txt
# Usage:   Rscript scripts/18s_binary_vs_survival.R [--probe-only]
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
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
OD <- file.path(PR, "super_learner")
HZ <- c(10, 30, 90, 180, 365)
TARGETS <- c("death", "comp", "cv", "cv_excl")
FAMS <- c(ridge = "lp", rf = "prob")
NBOOT <- 2000L
SEED <- 42L
PKG <- c("data.table")
options(width = 220)
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
  cat(sprintf("  %-5s %-64s %s\n", status, item, value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "binsurv_checks.csv"))
    cat("\nSTOP: a check failed. see: binsurv_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- helpers (sent to workers) ----------------
nscore <- function(x, ref) {
  qnorm((findInterval(x, sort(ref)) + 0.5) / (length(ref) + 1))
}
proj_simplex <- function(v) {
  u <- sort(v, decreasing = TRUE)
  cs <- cumsum(u)
  r <- max(which(u - (cs - 1) / seq_along(u) > 0))
  pmax(v - (cs[r] - 1) / r, 0)
}
sl_weights <- function(X, y) {
  J <- ncol(X)
  a <- rep(1 / J, J)
  XtX <- crossprod(X)
  Xty <- as.vector(crossprod(X, y))
  L <- max(eigen(XtX, symmetric = TRUE, only.values = TRUE)$values)
  for (it in seq_len(20000)) {
    an <- proj_simplex(a - (as.vector(XtX %*% a) - Xty) / L)
    if (max(abs(an - a)) < 1e-12) {
      a <- an
      break
    }
    a <- an
  }
  a
}
logit_fit <- function(Za, ya, Zb) {
  g <- suppressWarnings(glm.fit(cbind(1, Za), ya, family = binomial()))
  cf <- g$coefficients
  cf[is.na(cf)] <- 0
  plogis(as.vector(cbind(1, Zb) %*% cf))
}
# candidates for one horizon, fitted on rows a and predicting rows b
bin_candidates <- function(S, a, b, y) {
  Z <- function(fam, ms, rows) {
    matrix(sapply(ms, function(m) nscore(S[[fam]][[m]][rows],
                                         S[[fam]][[m]][a])),
           ncol = length(ms))
  }
  out <- list()
  out[["EHR|prevalence"]] <- rep(mean(y[a]), length(b))
  for (fam in names(S)) {
    out[[paste0("EHR|EHR alone, ", fam)]] <-
      logit_fit(Z(fam, "ehr", a), y[a], Z(fam, "ehr", b))
    out[[paste0("EHR+ECG|EHR+ECG stack, ", fam)]] <-
      logit_fit(Z(fam, c("ehr", "ecg"), a), y[a],
                Z(fam, c("ehr", "ecg"), b))
  }
  out
}
# one task = frame x target x seed; every horizon, every outer fold
bin_task <- function(tk) {
  n <- length(tk$subj)
  out <- list()
  wts <- list()
  put <- function(k, idx, hi, v) {
    if (is.null(out[[k]])) out[[k]] <<- matrix(NA_real_, n, length(tk$hz))
    out[[k]][idx, hi] <<- v
  }
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  for (hi in seq_along(tk$hz)) {
    y <- tk$Y[, hi]
    S <- tk$S[[hi]]
    ok <- !is.na(y)
    for (k in sort(unique(tk$fold))) {
      tr <- which(tk$fold != k & ok)
      te <- which(tk$fold == k & ok)
      inner <- grouped(tr, tk$seed * 100 + k, 3)
      PI <- NULL
      for (j in 1:3) {
        a <- tr[inner != j]
        b <- tr[inner == j]
        cj <- bin_candidates(S, a, b, y)
        if (is.null(PI)) {
          PI <- lapply(cj, function(v) rep(NA_real_, length(tr)))
        }
        for (nm in names(cj)) PI[[nm]][inner == j] <- cj[[nm]]
      }
      CT <- bin_candidates(S, tr, te, y)
      for (nm in names(CT)) {
        z <- strsplit(nm, "|", fixed = TRUE)[[1]]
        put(paste0("Binary candidate: ", z[2], "|", z[1]), te, hi, CT[[nm]])
      }
      libs <- list("Binary super learner-EHR|EHR" =
                     grep("^EHR\\|", names(CT), value = TRUE),
                   "Binary super learner|EHR+ECG" = names(CT))
      for (ln in names(libs)) {
        cn <- libs[[ln]]
        X <- do.call(cbind, PI[cn])
        w <- sl_weights(X, y[tr])
        ens <- as.vector(do.call(cbind, CT[cn]) %*% w)
        ein <- as.vector(X %*% w)
        g <- tryCatch(suppressWarnings(glm.fit(
          cbind(1, qlogis(clip(ein))), y[tr], family = binomial())),
          error = function(e) NULL)
        cf <- if (is.null(g)) c(NA, NA) else g$coefficients
        rc <- if (any(!is.finite(cf))) ens else
          plogis(cf[1] + cf[2] * qlogis(clip(ens)))
        z <- strsplit(ln, "|", fixed = TRUE)[[1]]
        put(paste0(z[1], "|", z[2]), te, hi, ens)
        put(paste0(z[1], ", recalibrated|", z[2]), te, hi, rc)
        wts[[length(wts) + 1]] <- data.table(
          model = z[1], horizon = tk$hz[hi], fold = k,
          candidate = sub("^[^|]*\\|", "", cn), weight = w)
      }
    }
  }
  list(job = tk$job, seed = tk$seed, out = out,
       wts = rbindlist(wts)[, `:=`(job = tk$job, seed = tk$seed)])
}
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(J$hz)) {
    ii <- idx[J$keep[idx, hi]]
    y <- J$Y[ii, hi]
    pr <- mean(y)
    for (k in names(J$P)) {
      raw <- J$P[[k]][ii, hi]
      key <- paste(k, J$hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, raw)
      p <- pmin(pmax(raw, 1e-6), 1 - 1e-6)
      v[paste(key, "ipa", sep = "|")] <- if (k %in% J$auc_only) NA_real_
      else 1 - mean((y - p)^2) / (pr * (1 - pr))
    }
  }
  v
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by), length(J$rows_by),
                                     replace = TRUE)], use.names = FALSE)
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat("=== BUILD (primary cohort, EHR+ECG) ===\n")
FRM <- list(P = file.path(PP, "frame_P_outcomes.csv"),
            B = file.path(PR, "frame_B.csv"))
SURV <- list(
  P = list(death = c("slP_preds_ehrecg.csv", "Super learner P, recalibrated"),
           comp = c("optB2_preds_P_comp_ehrecg.csv",
                    "Super learner, recalibrated"),
           cv = c("optB2_preds_P_cv_ehrecg.csv",
                  "Super learner, recalibrated")),
  B = list(death = c("optB2_preds_B_death_ehrecg.csv",
                     "Super learner, recalibrated"),
           comp = c("optB2_preds_B_comp_ehrecg.csv",
                    "Super learner, recalibrated"),
           cv = c("optB2_preds_B_cv_ehrecg.csv",
                  "Super learner, recalibrated")))
need <- c(unlist(FRM), file.path(OD, "sl7_preds_B_ehrecg.csv"))
for (fr in names(FRM)) {
  need <- c(need, file.path(L1, sprintf("oof_bin_%s_%s.csv", fr, TARGETS)),
            file.path(OD, sapply(SURV[[fr]], `[`, 1)))
}
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else
  sprintf("all %d", length(need)), if (length(miss) == 0) "PASS" else
    "FAIL")
stop_if_fail()
target_mat <- function(fp, tg) {
  sapply(HZ, function(h) {
    cu <- h + 0.5
    switch(tg,
           death = as.integer(fp$event_death == 1 & fp$time_death <= cu),
           comp = as.integer(fp$event_comp == 1 & fp$time_first <= cu),
           cv = as.integer(fp$event_cr == 1 & fp$time_first <= cu),
           cv_excl = ifelse(fp$event_cr == 2 & fp$time_first <= cu,
                            NA_integer_,
                            as.integer(fp$event_cr == 1 &
                                         fp$time_first <= cu)))
  })
}
tasks <- list()
META <- list()
for (fr in names(FRM)) {
  fp <- fread(FRM[[fr]])
  for (tg in TARGETS) {
    job <- paste(fr, tg)
    D <- fread(file.path(L1, sprintf("oof_bin_%s_%s.csv", fr, tg)))
    Y <- target_mat(fp, tg)
    seeds <- sort(unique(D$seed))
    for (s in seeds) {
      x <- D[D$seed == s][match(fp$hadm_id, hadm_id)]
      if (anyNA(x$hadm_id) || anyNA(x$fold)) {
        chk(sprintf("%s seed %s: rows and folds complete", job, s), "no",
            "FAIL")
        stop_if_fail()
      }
      SH <- lapply(HZ, function(h) {
        setNames(lapply(names(FAMS), function(fam) {
          setNames(lapply(c("ehr", "ecg"), function(m) {
            x[[sprintf("%s.%s.h%d.%s", m, fam, h, FAMS[[fam]])]]
          }), c("ehr", "ecg"))
        }), names(FAMS))
      })
      bad <- sum(sapply(seq_along(HZ), function(hi) {
        sum(sapply(SH[[hi]], function(z) sum(sapply(z, function(v)
          sum(!is.finite(v[!is.na(Y[, hi])]))))))
      }))
      if (bad > 0) {
        chk(sprintf("%s seed %s: missing level-1 predictions", job, s), bad,
            "FAIL")
        stop_if_fail()
      }
      tasks[[length(tasks) + 1]] <- list(job = job, seed = s, hz = HZ,
                                         Y = Y, S = SH, fold = x$fold,
                                         subj = fp$subject_id)
    }
    META[[job]] <- list(fr = fr, tg = tg, fp = fp, Y = Y)
    cat(sprintf("  %-10s n=%d  events by 10/30/365 days %s  seeds %d\n",
                job, nrow(fp),
                paste(colSums(Y, na.rm = TRUE)[c(1, 2, 5)], collapse = "/"),
                length(seeds)))
  }
}
cat(sprintf("  %d tasks (frame x target x seed)\n", length(tasks)))
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
t0 <- Sys.time()
invisible(bin_task(tasks[[1]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  one task %.0f s; %d tasks on %d workers: about %.1f",
                  "min, plus the bootstrap\n"), tp, length(tasks), workers,
            tp * ceiling(length(tasks) / workers) * 1.2 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, bin_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)

# ---------------- assemble ----------------
cat("\n=== ASSEMBLE ===\n")
V7 <- fread(file.path(OD, "sl7_preds_B_ehrecg.csv"))
V7 <- V7[model == "Super learner v7, recalibrated" & set == "EHR+ECG"]
JOBS <- list()
for (job in names(META)) {
  M <- META[[job]]
  n <- nrow(M$fp)
  rv <- Filter(function(r) r$job == job, res)
  P <- list()
  for (k in names(rv[[1]]$out)) {
    P[[k]] <- Reduce(`+`, lapply(rv, function(r) r$out[[k]])) / length(rv)
  }
  keep <- !is.na(M$Y)
  bad <- sum(sapply(P, function(m) sum(!is.finite(m[keep]))))
  rg <- range(unlist(lapply(P, function(m) m[keep])))
  chk(sprintf("%s: binary risks missing / range", job),
      sprintf("%d / %.5f to %.5f", bad, rg[1], rg[2]),
      if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
  LONG <- rbindlist(lapply(names(P), function(k) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = rep(M$fp$hadm_id, length(HZ)), model = z[1],
               set = z[2], horizon = rep(HZ, each = n),
               risk = as.vector(P[[k]]))
  }))
  fwrite(LONG, file.path(OD, sprintf("bin_preds_%s_%s.csv", M$fr, M$tg)))
  # survival counterpart on the same patients and horizons
  base <- if (M$tg == "cv_excl") "cv" else M$tg
  sv <- SURV[[M$fr]][[base]]
  SP <- fread(file.path(OD, sv[1]))[model == sv[2] & set == "EHR+ECG"]
  sm <- sapply(HZ, function(h) {
    y <- SP[horizon == h]
    y$risk[match(M$fp$hadm_id, y$hadm_id)]
  })
  nm <- sum(is.na(sm))
  chk(sprintf("%s: survival predictions present (missing)", job), nm,
      if (nm == 0) "PASS" else "FAIL")
  P[["Survival super learner|EHR+ECG"]] <- sm
  if (M$fr == "B" && M$tg == "death") {
    P[["Survival super learner v7|EHR+ECG"]] <- sapply(HZ, function(h) {
      y <- V7[horizon == h]
      y$risk[match(M$fp$hadm_id, y$hadm_id)]
    })
  }
  keepP <- c(grep("^(Binary super learner|Survival)", names(P),
                  value = TRUE),
             grep("^Binary candidate", names(P), value = TRUE))
  JOBS[[job]] <- list(job = job, M = M, hz = HZ, P = P[keepP],
                      Y = replace(M$Y, is.na(M$Y), 0L), keep = keep,
                      auc_only = if (M$tg == "cv_excl")
                        grep("^Survival", keepP, value = TRUE) else
                          character(0),
                      rows_by = split(seq_len(n),
                                      match(M$fp$subject_id,
                                            unique(M$fp$subject_id))))
}
stop_if_fail()
W <- rbindlist(lapply(res, function(r) r$wts))
wsum <- W[, .(mean_weight = mean(weight), share_nonzero = mean(weight > 1e-6)),
          by = .(job, model, horizon, candidate)]
setorder(wsum, job, model, horizon, -mean_weight)
fwrite(wsum, file.path(RS, "binsurv_weights.csv"))

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples, clustered by patient) ===\n",
            NBOOT))
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
plan(multisession, workers = workers)
for (J in JOBS) {
  M <- J$M
  JB <- J
  JB$P <- J$P[grep("^(Binary super learner|Survival)", names(J$P),
                   value = TRUE)]
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  bn <- names(stat_vec(JB, seq_len(nrow(J$Y))))
  BL <- future_lapply(seq_len(NBOOT), boot_rep, J = JB,
                      future.seed = SEED, future.packages = PKG)
  B <- do.call(rbind, lapply(BL, function(v) v[bn]))
  colnames(B) <- bn
  cat(sprintf("  %s done\n", J$job))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- if (k %in% colnames(B)) q2(B[, k]) else c(NA_real_, NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      frame = M$fr, target = M$tg, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in names(JB$P)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      ii <- which(J$keep[, hi])
      p <- pmin(pmax(J$P[[k]][ii, hi], 1e-6), 1 - 1e-6)
      sl <- if (k %in% J$auc_only) NA_real_ else
        tryCatch(unname(coef(glm(J$Y[ii, hi] ~ qlogis(p),
                                 family = binomial))[2]),
                 error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        frame = M$fr, target = M$tg, model = z[1], set = z[2],
        horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  pairs <- list(c("Survival super learner", "EHR+ECG",
                  "Binary super learner, recalibrated", "EHR+ECG"),
                c("Binary super learner, recalibrated", "EHR+ECG",
                  "Binary super learner-EHR, recalibrated", "EHR"))
  if (M$fr == "B" && M$tg == "death") {
    pairs[[3]] <- c("Survival super learner v7", "EHR+ECG",
                    "Binary super learner, recalibrated", "EHR+ECG")
  }
  for (cd in pairs) for (h in J$hz) for (m0 in c("auc", "ipa")) {
    ka <- paste(cd[1], cd[2], h, m0, sep = "|")
    kb <- paste(cd[3], cd[4], h, m0, sep = "|")
    if (!(ka %in% colnames(B)) || !(kb %in% colnames(B))) next
    d <- B[, ka] - B[, kb]
    if (all(is.na(d))) next
    q <- q2(d)
    CON[[length(CON) + 1]] <- data.table(
      frame = M$fr, target = M$tg,
      contrast = sprintf("%s - %s", cd[1], cd[3]), horizon = h,
      measure = m0, estimate = pt[[ka]] - pt[[kb]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0)
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
fwrite(mt, file.path(RS, "binsurv_metrics.csv"))
fwrite(cn, file.path(RS, "binsurv_contrasts.csv"))

# ---------------- summary ----------------
pv <- function(x) sprintf("%.3f", x)
s1 <- mt[measure == "auc" & model %in% c("Survival super learner",
                                         "Binary super learner, recalibrated"),
         .(frame, target, horizon, model, auc = estimate)]
s1 <- dcast(s1, frame + target + horizon ~ model, value.var = "auc")
setnames(s1, c("Binary super learner, recalibrated",
               "Survival super learner"), c("auc_binary", "auc_survival"))
s2 <- cn[measure == "auc" & grepl("^Survival super learner - ", contrast),
         .(frame, target, horizon, diff = estimate, lo, hi,
           excludes_zero)]
smy <- merge(s1, s2, by = c("frame", "target", "horizon"))
smy[, within_001 := abs(diff) <= 0.01]
setorder(smy, frame, target, horizon)
fwrite(smy, file.path(RS, "binsurv_summary.csv"))
chk("A33: survival and binary within 0.01 AUROC (cells)",
    sprintf("%d of %d", sum(smy$within_001), nrow(smy)), "INFO")
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "binsurv_checks.csv"))
cat("\n=== SURVIVAL vs BINARY SUPER LEARNER (AUROC) ===\n")
print(smy[, .(frame, target, horizon, survival = round(auc_survival, 3),
              binary = round(auc_binary, 3), diff = round(diff, 3),
              lo = round(lo, 3), hi = round(hi, 3), excludes_zero)],
      row.names = FALSE, nrows = 200)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (bin_preds_*.csv, patient-level, local only)\n")
for (f in c("binsurv_summary.csv", "binsurv_contrasts.csv",
            "binsurv_metrics.csv", "binsurv_weights.csv",
            "binsurv_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "binsurv_log.txt"), " (terminal log)\n")
cat("\nOutputs: binsurv_summary.csv, binsurv_contrasts.csv,",
    "binsurv_metrics.csv, binsurv_checks.csv\n")