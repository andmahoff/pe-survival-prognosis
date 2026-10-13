# 18b_early_blockstd.R
# Registered block-standardised early fusion arms
# (PREREGISTRATION_ADDENDUM_2.md, A10; DEVIATIONS.md D10, D11).
#
# For each set, fits Early PH (ridge-Cox) and Early TV (ridge logistic per
# horizon) on the block-standardised concatenation: each block's
# standardised columns are divided by the square root of the block's size,
# and glmnet's own standardisation is switched off so that the scaling
# survives (D10). The penalty path has 80 values down to 1e-7 of the
# largest (D11).
#
# Frames, cohorts, sets, horizons, folds (every seed), preprocessing in each
# training fold and tuning are as in 18a_early_linear.R: median imputation,
# winsorising of non-binary features at the training 0.5th and 99.5th
# percentiles, standardisation, inner 3-fold subject-grouped folds with the
# same seeds, deviance tuning and Breslow risks on centred linear
# predictors. The 18a arm "Early PH plain concatenation" is refitted here
# exactly as in 18a and must reproduce its predictions (largest absolute
# risk difference below 1e-6), which confirms that both scripts see
# identical data and folds. If it does not, nothing is written.
#
# Checks:  18a reproduced; block scaling changes the model (rank
#          correlation with the plain concatenation below 0.9999); no
#          missing risks (all FAIL). Risks of exactly 0 or 1, and over 20%
#          of penalty choices at the least-penalised end (WARN).
# Notes:   stops if 18a predictions are missing or from an earlier version,
#          or if 18b predictions from an earlier run remain. Tasks (frame,
#          cohort, seed) run in parallel with future.apply after a timing
#          probe that pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, fusion_arch/arch_preds_linear_*.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: data/processed/fusion_arch/arch_preds_blockstd_*.csv (patient
#          level, not in the repo); results/blockstd_lambda.csv,
#          results/blockstd_checks.csv, results/blockstd_warnings.csv,
#          results/blockstd_log.txt
# Usage:   Rscript scripts/18b_early_blockstd.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OD <- file.path(PR, "fusion_arch")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
NLAM <- 80L
LMR <- 1e-7
PKG <- c("survival", "glmnet", "data.table")
options(width = 220)
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}

# ---------------- one task = frame x cohort x seed (all 5 folds) --------
bs_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  out <- list()
  edges <- list()
  put <- function(model, set, h, idx, val) {
    k <- paste(model, set, h, sep = "|")
    if (is.null(out[[k]])) out[[k]] <<- rep(NA_real_, n)
    out[[k]][idx] <<- val
  }
  rec <- function(model, set, h, cv) {
    K <- length(cv$lambda)
    i <- which(cv$lambda == cv$lambda.min)[1]
    edges[[length(edges) + 1]] <<- data.table(
      model = model, set = set, horizon = h, nlam = K,
      least = i == K, most = i == 1)
  }
  breslow <- function(tt, ee, lp, cut) {
    o <- order(tt)
    tt <- tt[o]
    ee <- ee[o]
    w <- exp(lp[o] - max(lp))
    rs <- rev(cumsum(rev(w)))
    ut <- unique(tt[ee == 1 & tt <= cut])
    if (!length(ut)) return(0)
    d <- sapply(ut, function(u) sum(ee[tt == u]))
    sum(d / rs[match(ut, tt)])
  }
  prep <- function(X, tr) {
    for (j in seq_len(ncol(X))) {
      m <- median(X[tr, j], na.rm = TRUE)
      if (!is.finite(m)) m <- 0
      X[is.na(X[, j]), j] <- m
      if (length(unique(X[tr, j])) > 2) {
        q <- quantile(X[tr, j], c(0.005, 0.995), names = FALSE)
        X[, j] <- pmin(pmax(X[, j], q[1]), q[2])
      }
    }
    mu <- colMeans(X[tr, , drop = FALSE])
    s <- apply(X[tr, , drop = FALSE], 2, sd)
    keep <- is.finite(s) & s > 0
    X <- sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
    X
  }
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  warn <- character(0)
  notes <- 0L
  for (k in sort(unique(tk$fold))) {
    t_fold <- Sys.time()
    nw0 <- length(warn)
    withCallingHandlers({
      tr <- which(tk$fold != k)
      te <- which(tk$fold == k)
      inner <- grouped(tr, tk$seed * 100 + k, 3)
      B <- lapply(tk$X, prep, tr = tr)
      tt <- tk$t
      ee <- tk$e
      ph_risk <- function(lp_tr, lp_te) {
        c0 <- max(lp_tr)
        sapply(H, function(h) {
          hs <- breslow(tt[tr], ee[tr], lp_tr, h + 0.5)
          if (hs == 0) return(rep(0, length(lp_te)))
          -expm1(-hs * exp(lp_te - c0))
        })
      }
      ph_model <- function(Xm, model, set) {
        cv <- cv.glmnet(Xm[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                        family = "cox", alpha = 0, foldid = inner,
                        nlambda = NLAM, lambda.min.ratio = LMR,
                        type.measure = "deviance", standardize = FALSE)
        if (!is.null(model)) rec(model, set, NA_real_, cv)
        lp <- function(rows) {
          as.vector(predict(cv, Xm[rows, , drop = FALSE], s = "lambda.min"))
        }
        matrix(ph_risk(lp(tr), lp(te)), ncol = length(H))
      }
      tv_model <- function(Xm, model, set) {
        matrix(sapply(H, function(h) {
          y <- as.numeric(ee == 1 & tt <= h + 0.5)
          cv <- suppressWarnings(cv.glmnet(
            Xm[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0,
            foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
            type.measure = "deviance", standardize = FALSE))
          rec(model, set, h, cv)
          as.vector(predict(cv, Xm[te, , drop = FALSE], s = "lambda.min",
                            type = "response"))
        }), ncol = length(H))
      }
      for (sn in names(tk$sets)) {
        ms <- tk$sets[[sn]]
        Xc <- do.call(cbind, lapply(ms, function(m) {
          B[[m]] / sqrt(ncol(B[[m]]))
        }))
        rb <- ph_model(Xc, "Early PH", sn)
        tv <- tv_model(Xc, "Early TV", sn)
        rp <- ph_model(do.call(cbind, B[ms]), NULL, sn)
        for (i in seq_along(H)) {
          put("Early PH", sn, H[i], te, rb[, i])
          put("Early TV", sn, H[i], te, tv[, i])
          put("check 18a plain", sn, H[i], te, rp[, i])
        }
      }
    }, warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("tie-handling", msg, fixed = TRUE)) {
        notes <<- notes + 1L
      } else {
        warn <<- c(warn, msg)
      }
      invokeRestart("muffleWarning")
    })
    if (!is.null(tk$progress)) {
      cat(sprintf("%s  %s %-7s seed %s fold %d done in %.0f s, %d warnings\n",
                  format(Sys.time(), "%H:%M:%S"), tk$fr, tk$co, tk$seed, k,
                  as.numeric(difftime(Sys.time(), t_fold, units = "secs")),
                  length(warn) - nw0),
          file = tk$progress, append = TRUE)
    }
  }
  list(fr = tk$fr, co = tk$co, seed = tk$seed, out = out, warn = warn,
       notes = notes,
       edges = rbindlist(edges)[, `:=`(frame = tk$fr, cohort = tk$co)])
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- list(ehr = fread(file.path(IN, "ehr28.csv")),
             ecg = fread(file.path(IN, "ecg_opt.csv")),
             ctpa = fread(file.path(IN, "ctpa46.csv")))
mats <- lapply(mats, function(x) unique(x, by = "hadm_id"))
COH <- list("EHR+ECG" = list(code = "ehrecg",
                             sets = list("EHR+ECG" = c("ehr", "ecg"))),
            "3-mod" = list(code = "3mod",
                           sets = list("EHR+ECG" = c("ehr", "ecg"),
                                       "EHR+CTPA" = c("ehr", "ctpa"),
                                       "EHR+ECG+CTPA" =
                                         c("ehr", "ecg", "ctpa"))))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in names(COH)) {
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  lf <- file.path(OD, sprintf("arch_preds_linear_%s_%s.csv", fr,
                              COH[[co]]$code))
  ok <- file.exists(lf)
  labs <- if (ok) unique(fread(lf, select = "model")$model) else character(0)
  chk(sprintf("%s %s: stage 18a predictions found, D11 version", fr, co),
      basename(lf),
      if (ok && "Early PH plain concatenation" %in% labs &&
          !("Early PH" %in% labs)) "PASS" else "FAIL")
  X <- list()
  for (m in unique(unlist(COH[[co]]$sets))) {
    cols <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cols])
    storage.mode(Xm) <- "double"
    X[[m]] <- Xm
  }
  fo <- folds[frame == fr & cohort == co]
  META[[paste(fr, co)]] <- list(fr = fr, co = co, sub = sub, lf = lf,
                                code = COH[[co]]$code)
  for (s in sort(unique(fo$seed))) {
    fs <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(fs)) stop("fold missing: ", fr, " ", co, " seed ", s)
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = fs, subj = sub$subject_id,
      t = sub$time_death, e = as.integer(sub$event_death == 1),
      X = X, sets = COH[[co]]$sets, hz = HZ[[fr]])
  }
}
old <- list.files(OD, pattern = "^arch_preds_blockstd_.*\\.csv$")
chk("previous stage 18b prediction files (archive them to v1)",
    length(old), if (length(old) == 0) "PASS" else "FAIL")
cat(sprintf("  %d tasks (frame x cohort x seed), 5 folds each\n",
            length(tasks)))
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "blockstd_checks.csv"))
  cat("\nSTOP: a check failed. see: blockstd_checks.csv\n")
  quit(status = 1)
}

cat("\n=== TIMING PROBE: the slowest kinds of fit, per frame x cohort ===\n")
prep_fold <- function(X, tr) {
  for (j in seq_len(ncol(X))) {
    m <- median(X[tr, j], na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    X[is.na(X[, j]), j] <- m
    if (length(unique(X[tr, j])) > 2) {
      q <- quantile(X[tr, j], c(0.005, 0.995), names = FALSE)
      X[, j] <- pmin(pmax(X[, j], q[1]), q[2])
    }
  }
  mu <- colMeans(X[tr, , drop = FALSE])
  s <- apply(X[tr, , drop = FALSE], 2, sd)
  keep <- is.finite(s) & s > 0
  sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
}
probe_kind <- function(tk) {
  k <- sort(unique(tk$fold))[1]
  tr <- which(tk$fold != k)
  us <- unique(tk$subj[tr])
  set.seed(tk$seed * 100 + k)
  f <- sample(rep_len(1:3, length(us)))
  inner <- f[match(tk$subj[tr], us)]
  B <- lapply(tk$X, prep_fold, tr = tr)
  Xs <- do.call(cbind, B[tk$sets[[length(tk$sets)]]])
  y <- as.numeric(tk$e == 1 & tk$t <= min(tk$hz) + 0.5)
  tm <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    proc.time()[["elapsed"]] - t0
  }
  cx <- function(X) {
    suppressWarnings(cv.glmnet(
      X[tr, , drop = FALSE], Surv(tk$t[tr], tk$e[tr]), family = "cox",
      alpha = 0, foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
      type.measure = "deviance", standardize = FALSE))
  }
  c(big = tm(cx(Xs)), small = tm(cx(B$ehr)),
    bin = tm(suppressWarnings(cv.glmnet(
      Xs[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0,
      foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
      type.measure = "deviance", standardize = FALSE))))
}
fits_per_fold <- function(tk) {
  se <- sum(sapply(tk$sets, function(s) "ecg" %in% s))
  ns <- length(tk$sets)
  c(big = 2 * se, small = 2 * (ns - se), bin = length(tk$hz) * ns)
}
kind <- sapply(tasks, function(t) paste(t$fr, t$co))
reps <- tasks[!duplicated(kind)]
names(reps) <- unique(kind)
plan(multisession, workers = min(workers, length(reps)))
t0 <- Sys.time()
pt <- future_lapply(reps, probe_kind, future.packages = PKG,
                    future.seed = TRUE, future.chunk.size = 1)
plan(sequential)
fold_s <- c()
for (kd in names(reps)) {
  u <- pt[[kd]]
  nf <- fits_per_fold(reps[[kd]])
  fold_s[kd] <- sum(u * nf[names(u)])
  cat(sprintf(paste("  %-10s wide Cox fit %5.1f s (x%d), narrow %4.1f s",
                    "(x%d), logistic %4.1f s (x%d): a fold %4.1f min\n"),
              kd, u[["big"]], nf[["big"]], u[["small"]], nf[["small"]],
              u[["bin"]], nf[["bin"]], fold_s[kd] / 60))
}
cost <- fold_s[kind] * 5
tasks <- tasks[order(-cost)]
cost <- sort(cost, decreasing = TRUE)
cat(sprintf(paste("  probe took %.1f min; %d tasks on %d workers: about",
                  "%.1f hours (longest task %.1f hours)\n"),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            length(tasks), workers,
            max(sum(cost) / workers, max(cost)) * 1.2 / 3600,
            max(cost) / 3600))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
PROG <- file.path(OD, "progress_18b.txt")
cat(sprintf("%s  started, %d tasks x 5 folds\n",
            format(Sys.time(), "%H:%M:%S"), length(tasks)), file = PROG)
for (i in seq_along(tasks)) tasks[[i]]$progress <- PROG
cat("  progress: one line per finished fold in", PROG, "\n")
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, bs_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)

# ---------------- assemble and check ----------------
cat("\n=== CHECKS ===\n")
NEW <- list()
for (key in names(META)) {
  M <- META[[key]]
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  rr <- rr[order(sapply(rr, function(r) r$seed))]
  keys <- names(rr[[1]]$out)
  L <- rbindlist(lapply(keys, function(k) {
    v <- Reduce(`+`, lapply(rr, function(r) r$out[[k]])) / length(rr)
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = M$sub$hadm_id, model = z[1], set = z[2],
               horizon = as.numeric(z[3]), risk = v)
  }))
  old <- fread(M$lf)
  lab <- "Early PH plain concatenation"
  cmp <- merge(L[model == "check 18a plain",
                 .(hadm_id, set, horizon, new = risk)],
               old[model == lab, .(hadm_id, set, horizon, old = risk)],
               by = c("hadm_id", "set", "horizon"))
  dmax <- max(abs(cmp$new - cmp$old))
  chk(sprintf("%s: stage 18a plain concatenation reproduced (rows, max)",
              key), sprintf("%d, %.1e", nrow(cmp), dmax),
      if (nrow(cmp) == nrow(old[model == lab]) && dmax < 1e-6) "PASS"
      else "FAIL")
  h <- max(HZ[[M$fr]])
  for (sn in unique(L$set)) {
    a <- L[model == "Early PH" & set == sn & horizon == h][
      match(M$sub$hadm_id, hadm_id), risk]
    b <- cmp[set == sn & horizon == h][match(M$sub$hadm_id, hadm_id), new]
    rho <- cor(a, b, method = "spearman")
    chk(sprintf("%s %s: block scaling changes the model (rho at %g d)", key,
                sn, h), sprintf("%.4f", rho),
        if (rho < 0.9999) "PASS" else "FAIL")
  }
  L <- L[model != "check 18a plain"]
  nna <- sum(is.na(L$risk))
  nsat <- sum(L$risk <= 0 | L$risk >= 1, na.rm = TRUE)
  rg <- range(L$risk, na.rm = TRUE)
  chk(sprintf("%s: missing predictions", key), nna,
      if (nna == 0) "PASS" else "FAIL")
  chk(sprintf("%s: risks of exactly 0 or 1 (range)", key),
      sprintf("%d (%.2e to %.6f)", nsat, rg[1], rg[2]),
      if (nsat == 0) "PASS" else "WARN")
  NEW[[key]] <- L
}
EDG <- rbindlist(lapply(res, function(r) r$edges))
esum <- EDG[, .(fits = .N, median_path_length = median(nlam),
                share_least_penalised = mean(least),
                share_most_penalised = mean(most)),
            by = .(frame, cohort, model, set, horizon)]
setorder(esum, frame, cohort, model, set, horizon)
fwrite(esum, file.path(RS, "blockstd_lambda.csv"))
eall <- EDG[, .(share = mean(least)), by = model]
for (i in seq_len(nrow(eall))) {
  chk(sprintf("%s: choices at the least-penalised end", eall$model[i]),
      sprintf("%.1f%%", 100 * eall$share[i]),
      if (eall$share[i] <= 0.20) "PASS" else "WARN")
}
WRN <- rbindlist(lapply(res, function(r) {
  if (!length(r$warn)) return(NULL)
  data.table(frame = r$fr, cohort = r$co, message = substr(r$warn, 1, 150))
}))
if (nrow(WRN)) {
  wsum <- WRN[, .(count = .N), by = .(frame, cohort, message)]
  setorder(wsum, frame, cohort, -count)
} else {
  wsum <- data.table(frame = character(0), cohort = character(0),
                     message = character(0), count = integer(0))
}
fwrite(wsum, file.path(RS, "blockstd_warnings.csv"))
chk("warnings from model fitting (see blockstd_warnings.csv)",
    sum(wsum$count), "INFO")
chk("glmnet tie-handling notices (harmless; Breslow ties throughout)",
    sum(sapply(res, function(r) r$notes)), "INFO")
ck <- rbindlist(CHK)
if (any(ck$status == "FAIL")) {
  fwrite(ck, file.path(RS, "blockstd_checks.csv"))
  cat("\nSTOP: a check failed; nothing written.",
      "Outputs: blockstd_checks.csv\n")
  quit(status = 1)
}
for (key in names(META)) {
  M <- META[[key]]
  fwrite(NEW[[key]], file.path(OD, sprintf("arch_preds_blockstd_%s_%s.csv",
                                           M$fr, M$code)))
}
fwrite(ck, file.path(RS, "blockstd_checks.csv"))
cat("\n=== PENALTY CHOICES AT EITHER END OF THE PATH ===\n")
print(EDG[, .(fits = .N, least = round(mean(least), 3),
              most = round(mean(most), 3)), by = model], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (arch_preds_blockstd_*.csv, patient-level, local only)\n")
for (f in c("blockstd_lambda.csv", "blockstd_checks.csv",
            "blockstd_warnings.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "blockstd_log.txt"), " (terminal log)\n")
cat("\nOutputs: blockstd_lambda.csv, blockstd_checks.csv,",
    "blockstd_warnings.csv\n")