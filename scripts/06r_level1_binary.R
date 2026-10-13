# 06r_level1_binary.R
# Matched binary level-1 models, out of fold, for the comparison of survival
# and binary super learners (PREREGISTRATION_ADDENDUM_2.md, A33 Option B).
#
# --frame=P is at diagnosis (stage 4r's frame_P_outcomes.csv; EHR at
# presentation and first ECG; Frame A's folds) and --frame=B at discharge
# (frame_B.csv; EHR+ and ECG as stage 4r wrote them; Frame B's folds). Only
# the primary EHR+ECG cohort is fitted, as A33 specifies.
#
# For each horizon h (10, 30, 90, 180 and 365 days) there is a separate
# yes/no target by h:
#   death    death by h
#   comp     CV readmission or death by h
#   cv       CV readmission first by h, death first counting as no (A33's
#            primary)
#   cv_excl  as cv, but patients who died first by h are left out of that
#            horizon's training and evaluation and score NA (A33's
#            secondary; the CV target of the dissertation,
#            pe-multimodal-prognosis-msc)
#
# Per modality, target and horizon: ridge logistic regression (glmnet,
# alpha 0, 80 penalties down to 1e-4 of the largest, grouped inner 3-fold
# CV on deviance, path extended as in 06i) and a probability forest
# (randomForestSRC classification, 500 trees, nsplit 10) tuned by
# out-of-bag Brier score with tune_expand.R (node size 5, 20, 60; mtry
# sqrt(p) x 0.5, 1, 2; extended past the edges).
#
# Checks: a missing prediction for a row with a target is FAIL. D2's
# per-fold rule on AUROC (below 0.40) is only WARN, since early CV horizons
# have few events per fold. Forest failures and grid edges are WARN too.
#
# Inputs:  data/processed/presentation/frame_P_outcomes.csv, ehrP.csv,
#          ecgP.csv (frame P); data/processed/frame_B.csv and
#          discharge/ehrplusB.csv, ecgB.csv (frame B);
#          data/processed/folds.csv; scripts/tune_expand.R
# Outputs: results/bin_<frame>_level1_checks.csv,
#          results/bin_<frame>_tuning_summary.csv,
#          results/bin_<frame>_level1_log.txt
#          data/processed/level1/oof_bin_<frame>_<target>.csv (patient
#          level, not in the repo)
# Usage:   Rscript scripts/06r_level1_binary.R --frame=P|B [--probe-only]
#          [--workers=N]
# Notes:   --probe-only stops after the timing probe (fold 0, 30 days).

suppressPackageStartupMessages({
  library(glmnet)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
fa <- sub("^--frame=", "", grep("^--frame=", args, value = TRUE))
if (length(fa) != 1 || !(fa %in% c("P", "B"))) stop("give --frame=P or B")
FR <- fa
TAG <- paste0("bin_", FR)
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CKPT <- file.path(L1, paste0(TAG, "_ckpt"))
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
source(file.path(P2, "scripts", "tune_expand.R"))
HZ <- c(10, 30, 90, 180, 365)
TARGETS <- c("death", "comp", "cv", "cv_excl")
NLAM <- 80L
LMR <- 1e-4
NTREE <- 500L
NSPLIT <- 10L
NODES <- c(5, 20, 60)
MTRYF <- c(0.5, 1, 2)
TOL <- 0.001
MAXEXT <- 3L
NFOLD <- 5L
PKG <- c("glmnet", "data.table", "randomForestSRC")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 200)
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
    fwrite(ck, file.path(RS, paste0(TAG, "_level1_checks.csv")))
    cat(sprintf("\nSTOP: a check failed. see: %s_level1_checks.csv\n",
                TAG))
    quit(status = 1)
  }
}

# ---------------- one task = target x horizon x modality x seed ------------
run_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                     probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  n <- length(tk$y)
  ok <- !is.na(tk$y)
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
    sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
  }
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  out <- list(ridge = rep(NA_real_, n), rf = rep(NA_real_, n))
  prm <- list()
  exts <- list()
  fails <- list()
  for (k in folds_to_run) {
    tr <- which(tk$fold != k & ok)
    te <- which(tk$fold == k & ok)
    y <- tk$y[tr]
    Xs <- prep(tk$X, tr)
    # ridge logistic, path extended if the smallest penalty is chosen
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    lmr <- tk$lmr
    for (pass in 1:3) {
      cv <- cv.glmnet(Xs[tr, , drop = FALSE], y, family = "binomial",
                      alpha = 0, foldid = inner, nlambda = tk$nlam,
                      lambda.min.ratio = lmr, type.measure = "deviance",
                      standardize = FALSE)
      at_end <- which(cv$lambda == cv$lambda.min) == length(cv$lambda)
      if (!at_end || pass == 3) break
      lmr <- lmr / 100
    }
    out$ridge[te] <- as.vector(predict(cv, Xs[te, , drop = FALSE],
                                       s = "lambda.min", type = "link"))
    # probability forest, Brier-tuned out of bag
    X <- tk$X
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- med[na[, 2]]
    colnames(X) <- sprintf("f%d", seq_len(ncol(X)))
    dtr <- data.frame(y = factor(y, levels = c(0, 1)),
                      X[tr, , drop = FALSE])
    dte <- as.data.frame(X[te, , drop = FALSE])
    p <- ncol(X)
    gi <- 0L
    rf_eval <- function(par) {
      gi <<- gi + 1L
      f <- rfsrc(y ~ ., data = dtr, ntree = tk$ntree, nodesize = par$node,
                 mtry = par$mtry, nsplit = tk$nsplit, save.memory = TRUE,
                 seed = -as.integer(tk$seed * 1000 + k * 100 + gi))
      po <- f$predicted.oob[, "1"]
      br <- mean((y - po)^2, na.rm = TRUE)
      pr <- predict(f, newdata = dte)$predicted[, "1"]
      rm(f)
      invisible(gc(verbose = FALSE))
      list(score = br, pred = pr)
    }
    rr <- tune_expand(
      list(node = tk$nodes,
           mtry = unique(pmax(1, pmin(p, ceiling(tk$mtryf * sqrt(p)))))),
      rf_eval,
      list(node = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) v * 2, lo = 1,
                       hi = floor(length(tr) / 4)),
           mtry = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) ceiling(v * 1.75), lo = 1, hi = p)),
      tol = tk$tol, max_ext = tk$maxext)
    if (!is.null(rr$fit)) out$rf[te] <- rr$fit$pred
    pth <- rr$path
    if (!is.null(pth) && nrow(pth)) {
      bad <- pth[pth$error != ""]
      if (nrow(bad)) {
        fails[[length(fails) + 1]] <- data.table(
          fold = k, setting = sprintf("node %g, mtry %g", bad$node,
                                      bad$mtry), message = bad$error)
      }
    }
    prm[[length(prm) + 1]] <- data.table(
      fold = k, ridge_path_passes = pass,
      rf_node = if (is.null(rr$best)) NA else rr$best$node,
      rf_mtry = if (is.null(rr$best)) NA else rr$best$mtry,
      rf_edge_node = if (is.null(rr$edge)) NA else rr$edge[["node"]],
      rf_edge_mtry = if (is.null(rr$edge)) NA else rr$edge[["mtry"]],
      rf_fails = rr$nfail, p = p)
    if (nrow(rr$ext)) {
      exts[[length(exts) + 1]] <- cbind(data.table(fold = k), rr$ext)
    }
  }
  res <- list(tg = tk$tg, h = tk$h, mod = tk$mod, seed = tk$seed,
              out = out, ok = ok, prm = rbindlist(prm),
              ext = rbindlist(exts, fill = TRUE),
              fails = rbindlist(fails, fill = TRUE),
              sec = proc.time()[["elapsed"]] - t_task)
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, tk$ckpt)
  }
  res
}

# ---------------- build ----------------
cat(sprintf("=== BUILD (frame %s, primary cohort EHR+ECG) ===\n", FR))
DIR <- if (FR == "P") PP else file.path(PR, "discharge")
fframe <- if (FR == "P") file.path(PP, "frame_P_outcomes.csv") else
  file.path(PR, "frame_B.csv")
fm <- if (FR == "P") c(ehr = "ehrP.csv", ecg = "ecgP.csv") else
  c(ehr = "ehrplusB.csv", ecg = "ecgB.csv")
need <- c(fframe, file.path(DIR, fm), file.path(PR, "folds.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- fread(fframe)
MX <- list()
for (m in names(fm)) {
  D <- unique(fread(file.path(DIR, fm[[m]])), by = "hadm_id")
  cc <- setdiff(names(D), c("hadm_id", "subject_id"))
  X <- as.matrix(D[match(fp$hadm_id, hadm_id), ..cc])
  storage.mode(X) <- "double"
  nm <- sum(is.na(match(fp$hadm_id, D$hadm_id)))
  chk(sprintf("%s: features (admissions missing)", toupper(m)),
      sprintf("%d (%d)", ncol(X), nm), if (nm == 0) "PASS" else "FAIL")
  MX[[m]] <- X
}
stop_if_fail()
target <- function(tg, h) {
  cu <- h + 0.5
  switch(tg,
         death = as.integer(fp$event_death == 1 & fp$time_death <= cu),
         comp = as.integer(fp$event_comp == 1 & fp$time_first <= cu),
         cv = as.integer(fp$event_cr == 1 & fp$time_first <= cu),
         cv_excl = ifelse(fp$event_cr == 2 & fp$time_first <= cu,
                          NA_integer_,
                          as.integer(fp$event_cr == 1 &
                                       fp$time_first <= cu)))
}
folds <- fread(file.path(PR, "folds.csv"))[frame == (if (FR == "P") "A"
                                                     else "B") &
                                             cohort == "EHR+ECG"]
seeds <- sort(unique(folds$seed))
fmap <- list()
for (s in seeds) {
  f <- folds[folds$seed == s][match(fp$hadm_id, hadm_id), fold]
  if (anyNA(f)) {
    chk(sprintf("seed %s: folds complete", s), sum(is.na(f)), "FAIL")
    stop_if_fail()
  }
  fmap[[as.character(s)]] <- as.integer(f)
}
EV <- list()
tasks <- list()
for (tg in TARGETS) for (h in HZ) {
  y <- target(tg, h)
  EV[[length(EV) + 1]] <- data.table(target = tg, horizon = h,
                                     patients = sum(!is.na(y)),
                                     events = sum(y, na.rm = TRUE))
  for (m in names(MX)) for (s in seeds) {
    tasks[[length(tasks) + 1]] <- list(
      tg = tg, h = h, mod = m, seed = s, X = MX[[m]], y = y,
      subj = fp$subject_id, fold = fmap[[as.character(s)]], nlam = NLAM,
      lmr = LMR, ntree = NTREE, nsplit = NSPLIT, nodes = NODES,
      mtryf = MTRYF, tol = TOL, maxext = MAXEXT,
      ckpt = file.path(CKPT, sprintf("%s_h%d_%s_%d.rds", tg, h, m, s)))
  }
}
ev <- rbindlist(EV)
cat("\n  events per target and horizon (patients in brackets):\n")
print(dcast(ev[, .(target, horizon, v = sprintf("%d (%d)", events,
                                                patients))],
            target ~ horizon, value.var = "v"), row.names = FALSE)
cat(sprintf("\n  %d admissions, %d seeds, %d tasks\n", nrow(fp),
            length(seeds), length(tasks)))
kind <- sapply(tasks, function(t) paste(t$tg, t$mod))
if (file.exists(PLANF)) {
  PL <- readRDS(PLANF)
  costs <- PL$costs
  cat(sprintf("\n=== RESUMING (plan of %s); probe skipped ===\n", PL$when))
} else {
  cat("\n=== TIMING PROBE (real data, fold 0 at 30 days, in parallel) ===\n")
  i30 <- which(sapply(tasks, function(t) t$h == 30))
  reps <- tasks[i30[!duplicated(kind[i30])]]
  names(reps) <- sapply(reps, function(t) paste(t$tg, t$mod))
  plan(multisession, workers = min(workers, length(reps)))
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  for (kd in names(pres)) {
    cat(sprintf("  %-14s %6.1f min per fold\n", kd, pres[[kd]]$sec / 60))
  }
  costs <- sapply(pres, function(r) r$sec)[kind] * NFOLD
  if (probe_only) {
    cat(sprintf("  %d tasks on %d workers: about %.1f hours\n",
                length(tasks), workers,
                max(sum(costs) / workers, max(costs)) * 1.1 / 3600))
    cat("--probe-only: stopping here (nothing saved).\n")
    quit(status = 0)
  }
  saveRDS(list(costs = costs, when = format(Sys.time(), "%Y-%m-%d %H:%M")),
          PLANF)
}
ok_ck <- function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(r) && all(sapply(r$out, function(v) all(is.finite(v[r$ok]))))
}
done <- which(sapply(tasks, function(t) ok_ck(t$ckpt)))
todo <- setdiff(seq_along(tasks), done)
proj <- max(sum(costs[todo]) / workers, max(c(0, costs[todo]))) * 1.1
cat(sprintf(paste("  %d of %d tasks already saved; %d to fit, about %.1f",
                  "hours on %d workers\n"), length(done), length(tasks),
            length(todo), proj / 3600, workers))
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
if (length(todo)) {
  plan(multisession, workers = workers)
  ord <- todo[order(-costs[todo])]
  invisible(future_lapply(tasks[ord], run_task, future.packages = PKG,
                          future.seed = TRUE, future.chunk.size = 1L))
  plan(sequential)
}
res <- lapply(tasks, function(t) readRDS(t$ckpt))

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE AND CHECK ===\n")
auc <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- data.table::frank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
PRM <- list()
EXT <- list()
FL <- list()
for (tg in TARGETS) {
  parts <- list()
  nmiss <- 0L
  nbad <- 0L
  for (s in seeds) {
    f <- fmap[[as.character(s)]]
    d <- data.table(seed = s, hadm_id = fp$hadm_id, fold = f)
    for (r in Filter(function(r) r$tg == tg && r$seed == s, res)) {
      y <- target(tg, r$h)
      for (fam in c("ridge", "rf")) {
        cn <- sprintf("%s.%s.h%d.%s", r$mod, fam, r$h,
                      if (fam == "ridge") "lp" else "prob")
        v <- r$out[[fam]]
        d[[cn]] <- v
        nmiss <- nmiss + sum(!is.finite(v[r$ok]))
        a <- sapply(sort(unique(f)), function(k) {
          ii <- which(f == k & r$ok)
          auc(y[ii], v[ii])
        })
        if (mean(a, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
      }
      tagd <- data.table(target = tg, horizon = r$h, modality = r$mod,
                         seed = s)
      PRM[[length(PRM) + 1]] <- cbind(tagd, r$prm)
      if (nrow(r$ext)) EXT[[length(EXT) + 1]] <- cbind(tagd, r$ext)
      if (nrow(r$fails)) FL[[length(FL) + 1]] <- cbind(tagd, r$fails)
    }
    parts[[length(parts) + 1]] <- d
  }
  fwrite(rbindlist(parts), file.path(L1, sprintf("oof_bin_%s_%s.csv", FR,
                                                 tg)))
  chk(sprintf("%s: missing predictions (rows with a target)", tg), nmiss,
      if (nmiss == 0) "PASS" else "FAIL")
  chk(sprintf("%s: predictions failing the per-fold rule (D2)", tg), nbad,
      if (nbad == 0) "PASS" else "WARN")
}
P <- rbindlist(PRM, fill = TRUE)
XT <- rbindlist(EXT, fill = TRUE)
FLS <- rbindlist(FL, fill = TRUE)
fwrite(P, file.path(L1, paste0("params_", TAG, ".csv")))
chk("forest fits that failed (best remaining one used)", sum(P$rf_fails),
    if (sum(P$rf_fails) == 0) "PASS" else "WARN")
chk("ridge: folds whose penalty path was extended",
    sum(P$ridge_path_passes > 1), "INFO")
for (pr in c("rf_edge_node", "rf_edge_mtry")) {
  sh <- mean(P[[pr]], na.rm = TRUE)
  chk(sprintf("share of folds still at an edge after extension: %s",
              sub("rf_edge_", "forest ", pr)), sprintf("%.0f%%", 100 * sh),
      if (sh <= 0.20) "PASS" else "WARN")
}
cs <- P[, .(setting = sprintf("node %g, mtry %g", rf_node, rf_mtry)),
        by = .(target, modality, horizon)]
cs <- cs[, .N, by = .(target, modality, setting)][order(target, modality,
                                                        -N)]
cs[, share := round(N / sum(N), 3), by = .(target, modality)]
fwrite(rbind(cbind(type = "chosen setting", cs),
             if (nrow(XT)) cbind(type = "extensions",
                                 XT[, .(N = .N, share = round(mean(kept), 3)),
                                    by = .(target, modality,
                                           setting = paste(param,
                                                           direction))]),
             fill = TRUE),
       file.path(RS, paste0(TAG, "_tuning_summary.csv")))
if (nrow(FLS)) {
  cat("\n  fit failures:\n")
  print(FLS, row.names = FALSE)
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, paste0(TAG, "_level1_checks.csv")))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (tg in TARGETS) {
  cat("  ", file.path(L1, sprintf("oof_bin_%s_%s.csv", FR, tg)),
      " (local only)\n")
}
cat("  ", file.path(L1, paste0("params_", TAG, ".csv")), " (local only)\n")
cat("  ", file.path(RS, paste0(TAG, "_level1_checks.csv")), "\n")
cat("  ", file.path(RS, paste0(TAG, "_tuning_summary.csv")), "\n")
cat("  ", file.path(RS, paste0(TAG, "_level1_log.txt")),
    " (terminal log)\n")
cat(sprintf("\nOutputs: %s_level1_checks.csv, %s_tuning_summary.csv\n",
            TAG, TAG))