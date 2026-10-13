# 16d_transfer_forests.R
# Zero-shot survival transfer between INSPECT and MIMIC with random survival
# forests, along follow-up time (PREREGISTRATION_ADDENDUM_2.md, A12).
#
# Version 2: each task restores the rf.cores and mc.cores options it sets.
#
# Data: Frames I and M, the 28 EHR and 38 CTPA features shared by both,
# prepared as in 16b (D12), the folds of 04b and their six seeds, with
# seed-averaged predictions. Learners, per modality (EHR, CTPA): ridge-Cox
# as in 16b; a random survival forest (randomForestSRC, 500 trees, log-rank
# splitting, nsplit 10, node size 5, 15, 40, 80 or 160 and mtry
# ceiling(f sqrt(p)) with f = 0.5, 1 or 3, chosen on out-of-bag Harrell's C,
# the D11 grid); and the rank average of the two scores' training-data
# percentiles. Each learner feeds Cox stacks EHR, CTPA and EHR + CTPA on
# rank-based normal scores, with the transform fitted on the training data's
# out-of-sample scores (inner CV for ridge-Cox, out-of-bag for the forest)
# and centred Breslow risks at 10, 30, 90, 180 and 365 days.
#
# Designs: internal (5-fold CV within each dataset) and transfer (trained on
# all of one dataset, applied unchanged to the other), in both directions.
# Target measures at each horizon: IPCW AUC (Uno et al. 2007), Brier score
# and IPA (Kaplan-Meier null model), Kaplan-Meier observed / expected and
# calibration slope, with Kaplan-Meier censoring weights (1 in MIMIC).
# Contrasts: RSF and rank average minus ridge, in performance and in
# transfer loss (transfer minus internal); the EHR minus CTPA gap in AUC and
# its change from 10 to 365 days. Intervals: 2,000 bootstrap resamples of
# the target patients (seed 42), censoring weights re-estimated, predictions
# fixed.
#
# Checks:  stops if the 16b results are missing or a task fails. FAIL if
#          predictions are missing or outside [0, 1], or 16b's internal EHR
#          stack is not reproduced to 1e-6 (AUC at 30 and 365 days). WARN on
#          IPCW differences over 1e-4 from riskRegression::Score and on
#          penalty or node-size choices at the edge of the grid.
# Notes:   all tasks, probe tasks included, run in the parallel workers as
#          in 16b, so that its inner folds are reproduced.
#
# Inputs:  data/processed/frame_I.csv, frame_M.csv, folds_IM.csv,
#          feature_whitelist_IM.csv; data/inputs/inspect_ehr.csv,
#          inspect/inspect_imp_features_v2.csv, ehr28.csv, ctpa46.csv;
#          results/s6_performance.csv (the D12 rerun of 16b)
# Outputs: results/a12_performance.csv, a12_contrasts.csv, a12_gap.csv,
#          a12_tuning.csv, a12_checks.csv, a12_warnings.csv, a12_log.txt
# Usage:   Rscript scripts/16d_transfer_forests.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(prodlim)
  library(riskRegression)
  library(glmnet)
  library(randomForestSRC)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OD <- file.path(PR, "transfer_ext")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
NLAM <- 80L
LMR <- 1e-7
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5L, 15L, 40L, 80L, 160L)
MTRYF <- c(0.5, 1, 3)
NBOOT <- 2000L
SEED <- 42L
STACKS <- list("EHR" = "EHR", "CTPA" = "CTPA",
               "EHR + CTPA" = c("EHR", "CTPA"))
LRN <- c("ridge-Cox" = "ridge", "RSF" = "rsf", "rank average" = "rank")
PKG <- c("survival", "glmnet", "randomForestSRC", "data.table")
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

# ---------------- modelling (sent to workers) ----------------
impute_med <- function(Xs, Xt) {
  meds <- apply(Xs, 2, median, na.rm = TRUE)
  meds[!is.finite(meds)] <- 0
  fill <- function(X) {
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- meds[na[, 2]]
    X
  }
  list(s = fill(Xs), t = fill(Xt))
}
winsorise <- function(imp) {
  for (j in seq_len(ncol(imp$s))) {
    if (length(unique(imp$s[, j])) > 2) {
      q <- quantile(imp$s[, j], c(0.005, 0.995), names = FALSE)
      imp$s[, j] <- pmin(pmax(imp$s[, j], q[1]), q[2])
      imp$t[, j] <- pmin(pmax(imp$t[, j], q[1]), q[2])
    }
  }
  imp
}
pct <- function(x, ref) {
  (findInterval(x, sort(ref)) + 0.5) / (length(ref) + 1)
}
nscore <- function(x, ref) qnorm(pct(x, ref))
bres_cuts <- function(tt, ee, lp, cuts) {
  o <- order(tt)
  tt <- tt[o]
  ee <- ee[o]
  w <- exp(lp[o] - max(lp))
  rs <- rev(cumsum(rev(w)))
  ut <- unique(tt[ee == 1])
  if (!length(ut)) return(rep(0, length(cuts)))
  d <- tabulate(match(tt[ee == 1], ut), length(ut))
  ch <- cumsum(d / rs[match(ut, tt)])
  sapply(cuts, function(cu) {
    j <- sum(ut <= cu)
    if (j == 0) 0 else ch[j]
  })
}
fit_task <- function(tk) {
  # one core per task; the old settings come back afterwards, so a probe
  # task run in the main session cannot block the parallel workers
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t0 <- proc.time()[["elapsed"]]
  imp <- winsorise(impute_med(tk$Xs, tk$Xt))
  ttr <- tk$ts
  etr <- tk$es
  str_ <- list()
  ste <- list()
  tune <- list()
  # ridge-Cox, in S6's order so the random inner folds match S6
  set.seed(tk$seed)
  for (m in names(tk$groups)) {
    cols <- intersect(tk$groups[[m]], colnames(imp$s))
    A <- imp$s[, cols, drop = FALSE]
    B <- imp$t[, cols, drop = FALSE]
    oof <- numeric(nrow(A))
    for (k in sort(unique(tk$ftr))) {
      tr <- tk$ftr != k
      inner <- sample(rep_len(1:3, sum(tr)))
      cv <- cv.glmnet(A[tr, , drop = FALSE], Surv(ttr[tr], etr[tr]),
                      family = "cox", alpha = 0, foldid = inner,
                      type.measure = "C", nlambda = NLAM,
                      lambda.min.ratio = LMR)
      oof[!tr] <- as.vector(predict(cv, A[!tr, , drop = FALSE],
                                    s = "lambda.min"))
    }
    cvf <- cv.glmnet(A, Surv(ttr, etr), family = "cox", alpha = 0,
                     foldid = tk$ftr, type.measure = "C", nlambda = NLAM,
                     lambda.min.ratio = LMR)
    i <- which(cvf$lambda == cvf$lambda.min)[1]
    tune[[length(tune) + 1]] <- data.table(
      group = m, learner = "ridge-Cox", nodesize = NA_integer_,
      mtry = NA_integer_, mtry_level = NA_real_,
      ridge_least = i == length(cvf$lambda), ridge_most = i == 1)
    str_[[m]] <- list(ridge = oof)
    ste[[m]] <- list(ridge = as.vector(predict(cvf, B, s = "lambda.min")))
  }
  # S6's own EHR stack (winsorised z-scores), kept only as a check that
  # this stage reproduces S6's ridge-Cox exactly
  v <- str_[["EHR"]]$ridge
  q <- quantile(v, c(0.005, 0.995), names = FALSE)
  w <- pmin(pmax(v, q[1]), q[2])
  z6 <- matrix((w - mean(w)) / sd(w), ncol = 1)
  z6t <- matrix((pmin(pmax(ste[["EHR"]]$ridge, q[1]), q[2]) - mean(w)) /
                  sd(w), ncol = 1)
  f6 <- coxph(Surv(ttr, etr) ~ z6, ties = "efron")
  bh <- basehaz(f6, centered = FALSE)
  lp6 <- as.vector(z6t %*% coef(f6))
  check6 <- sapply(tk$cuts, function(ct) {
    i <- findInterval(ct, bh$time)
    H0 <- if (i == 0) 0 else bh$hazard[i]
    1 - exp(-H0 * exp(lp6))
  })
  # random survival forests, tuned on out-of-bag C
  for (m in names(tk$groups)) {
    cols <- intersect(tk$groups[[m]], colnames(imp$s))
    A <- imp$s[, cols, drop = FALSE]
    B <- imp$t[, cols, drop = FALSE]
    colnames(A) <- colnames(B) <- sprintf("f%d", seq_along(cols))
    p <- ncol(A)
    grid <- expand.grid(nodesize = NODES, mf = MTRYF)
    grid$mtry <- pmax(1L, pmin(p, ceiling(grid$mf * sqrt(p))))
    d <- data.frame(time = ttr, status = etr, A)
    best <- -Inf
    keep <- NULL
    for (g in seq_len(nrow(grid))) {
      f <- rfsrc(Surv(time, status) ~ ., data = d, ntree = NTREE,
                 nodesize = grid$nodesize[g], mtry = grid$mtry[g],
                 nsplit = NSPLIT, ntime = NTIME, splitrule = "logrank",
                 seed = -as.integer(tk$seed * 1000 + g))
      cc <- as.numeric(concordance(Surv(ttr, etr) ~ f$predicted.oob,
                                   reverse = TRUE)$concordance)
      if (is.finite(cc) && cc > best) {
        best <- cc
        keep <- list(g = g, oob = f$predicted.oob,
                     te = predict(f, newdata = as.data.frame(B))$predicted)
      }
      rm(f)
    }
    tune[[length(tune) + 1]] <- data.table(
      group = m, learner = "RSF", nodesize = grid$nodesize[keep$g],
      mtry = grid$mtry[keep$g], mtry_level = grid$mf[keep$g],
      ridge_least = NA, ridge_most = NA)
    str_[[m]]$rsf <- keep$oob
    ste[[m]]$rsf <- keep$te
    str_[[m]]$rank <- (pct(str_[[m]]$ridge, str_[[m]]$ridge) +
                         pct(str_[[m]]$rsf, str_[[m]]$rsf)) / 2
    ste[[m]]$rank <- (pct(ste[[m]]$ridge, str_[[m]]$ridge) +
                        pct(ste[[m]]$rsf, str_[[m]]$rsf)) / 2
  }
  # Cox stacks on rank-based normal scores
  out <- list()
  for (ln in names(LRN)) for (sn in names(STACKS)) {
    lk <- LRN[[ln]]
    ms <- STACKS[[sn]]
    Ztr <- sapply(ms, function(m) nscore(str_[[m]][[lk]], str_[[m]][[lk]]))
    Zte <- sapply(ms, function(m) nscore(ste[[m]][[lk]], str_[[m]][[lk]]))
    Ztr <- matrix(Ztr, ncol = length(ms))
    Zte <- matrix(Zte, ncol = length(ms))
    b <- coef(coxph(Surv(ttr, etr) ~ Ztr, ties = "efron"))
    b[is.na(b)] <- 0
    lt <- as.vector(Ztr %*% b)
    lv <- as.vector(Zte %*% b)
    Hc <- bres_cuts(ttr, etr, lt, tk$cuts)
    out[[paste(sn, ln, sep = "|")]] <- sapply(seq_along(tk$cuts), function(i)
      -expm1(-Hc[i] * exp(lv - max(lt))))
  }
  if (!is.null(tk$progress)) {
    cat(sprintf("%s  %-8s %-8s seed %s fold %s done in %.1f min\n",
                format(Sys.time(), "%H:%M:%S"), tk$design, tk$target,
                tk$seed, tk$fold,
                (proc.time()[["elapsed"]] - t0) / 60),
        file = tk$progress, append = TRUE)
  }
  list(design = tk$design, target = tk$target, seed = tk$seed,
       rows = tk$rows, out = out,
       check6 = matrix(check6, ncol = length(tk$cuts)),
       tune = rbindlist(tune)[, `:=`(design = tk$design,
                                     target = tk$target)])
}
run_task <- function(tk) {
  warn <- character(0)
  res <- withCallingHandlers(
    tryCatch(c(list(ok = TRUE, msg = ""), fit_task(tk)),
             error = function(e) list(ok = FALSE, msg = conditionMessage(e),
                                      design = tk$design,
                                      target = tk$target, seed = tk$seed)),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  res$warn <- warn
  res
}

# ---------------- evaluation (IPCW, bootstrap) ----------------
km_steps <- function(t, e, cens) {
  ut <- sort(unique(t))
  nr <- length(t) - findInterval(ut, sort(t), left.open = TRUE)
  ev <- if (cens) e == 0 else e == 1
  nd <- tabulate(match(t[ev], ut), length(ut))
  list(ut = ut, S = cumprod(1 - nd / nr))
}
km_at <- function(K, x, left = FALSE) {
  c(1, K$S)[findInterval(x, K$ut, left.open = left) + 1]
}
stat_vec <- function(J, idx) {
  t <- J$t[idx]
  e <- J$e[idx]
  G <- km_steps(t, e, TRUE)
  S <- km_steps(t, e, FALSE)
  v <- numeric(0)
  for (hi in seq_along(J$cuts)) {
    h <- J$cuts[hi]
    case <- e == 1 & t <= h
    ctrl <- t > h
    wc <- 1 / km_at(G, t[case], left = TRUE)
    gh <- km_at(G, h)
    obs <- 1 - km_at(S, h)
    b0 <- mean(ifelse(case, (1 - obs)^2 / km_at(G, t, left = TRUE), 0) +
                 ifelse(ctrl, obs^2 / gh, 0))
    for (k in names(J$P)) {
      r <- J$P[[k]][idx, hi]
      rc <- sort(r[ctrl])
      lo <- findInterval(r[case], rc, left.open = TRUE)
      up <- findInterval(r[case], rc)
      au <- sum(wc * (lo + 0.5 * (up - lo)) / length(rc)) / sum(wc)
      br <- mean(ifelse(case, (1 - r)^2 / km_at(G, t, left = TRUE), 0) +
                   ifelse(ctrl, r^2 / gh, 0))
      key <- paste(k, J$hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- au
      v[paste(key, "ipa", sep = "|")] <- 1 - br / b0
      v[paste(key, "brier", sep = "|")] <- br
      v[paste(key, "oe", sep = "|")] <- obs / mean(r)
    }
  }
  v
}
boot_rep <- function(i, J) {
  stat_vec(J, sample.int(length(J$t), length(J$t), replace = TRUE))
}

# ---------------- data ----------------
cat("=== DATA ===\n")
fi <- fread(file.path(PR, "frame_I.csv"))
fm <- fread(file.path(PR, "frame_M.csv"))
wl <- fread(file.path(PR, "feature_whitelist_IM.csv"))
fo <- fread(file.path(PR, "folds_IM.csv"))
ehr_c <- wl[modality == "ehr", column]
ct_c <- wl[modality == "ctpa", column]
ie <- unique(fread(file.path(IN, "inspect_ehr.csv")), by = "person_id")
ic <- unique(fread(file.path(IN, "inspect", "inspect_imp_features_v2.csv")),
             by = "person_id")
me <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
mc <- unique(fread(file.path(IN, "ctpa46.csv")), by = "hadm_id")
num <- function(dt, key, ids, cols) {
  m <- as.matrix(dt[match(ids, dt[[key]]), ..cols])
  storage.mode(m) <- "double"
  m
}
XI <- cbind(num(ie, "person_id", fi$person_id, ehr_c),
            num(ic, "person_id", fi$person_id, ct_c))
XM <- cbind(num(me, "hadm_id", fm$hadm_id, ehr_c),
            num(mc, "hadm_id", fm$hadm_id, ct_c))
colnames(XI) <- colnames(XM) <- c(ehr_c, ct_c)
DS <- list(
  I = list(name = "INSPECT", X = XI, t = fi$time_death,
           e = as.integer(fi$event_death), id = fi$person_id, off = 0),
  M = list(name = "MIMIC", X = XM, t = fm$time_death,
           e = as.integer(fm$event_death), id = fm$hadm_id, off = 0.5))
GROUPS <- list(EHR = ehr_c, CTPA = ct_c)
for (f in c("I", "M")) {
  chk(sprintf("%s: patients / deaths within 366 d", DS[[f]]$name),
      sprintf("%d / %d", nrow(DS[[f]]$X), sum(DS[[f]]$e)), "INFO")
}
seeds <- sort(unique(fo$seed))
fold_of <- function(f, s) {
  x <- fo[fo$frame == f & fo$seed == s]
  x$fold[match(DS[[f]]$id, x$id)]
}
for (f in c("I", "M")) for (s in seeds) {
  if (anyNA(fold_of(f, s))) stop("fold missing: ", f, " seed ", s)
}
s6f <- file.path(RS, "s6_performance.csv")
chk("S6 rerun (D12) results found", basename(s6f),
    if (file.exists(s6f)) "PASS" else "FAIL")
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "a12_checks.csv"))
  cat("\nSTOP: a check failed. see: a12_checks.csv\n")
  quit(status = 1)
}

# ---------------- tasks (same random draws as S6) ----------------
tasks <- list()
for (f in c("I", "M")) for (s in seeds) {
  fold <- fold_of(f, s)
  for (k in 0:4) {
    tr <- which(fold != k)
    te <- which(fold == k)
    set.seed(s * 100 + k)
    tasks[[length(tasks) + 1]] <- list(
      design = "internal", target = f, seed = s, fold = k, rows = te,
      Xs = DS[[f]]$X[tr, , drop = FALSE], ts = pmax(DS[[f]]$t[tr], 1e-3),
      es = DS[[f]]$e[tr], ftr = sample(rep_len(1:5, length(tr))),
      Xt = DS[[f]]$X[te, , drop = FALSE], cuts = HZ + DS[[f]]$off,
      groups = GROUPS)
  }
}
for (tgt in c("I", "M")) {
  src <- setdiff(c("I", "M"), tgt)
  for (s in seeds) {
    tasks[[length(tasks) + 1]] <- list(
      design = "transfer", target = tgt, seed = s, fold = "all",
      rows = seq_len(nrow(DS[[tgt]]$X)), Xs = DS[[src]]$X,
      ts = pmax(DS[[src]]$t, 1e-3), es = DS[[src]]$e,
      ftr = fold_of(src, s) + 1L, Xt = DS[[tgt]]$X,
      cuts = HZ + DS[[tgt]]$off, groups = GROUPS)
  }
}
cat(sprintf("  %d tasks (%d internal, %d transfer)\n", length(tasks),
            sum(sapply(tasks, function(t) t$design == "internal")),
            sum(sapply(tasks, function(t) t$design == "transfer"))))

cat(sprintf("\n=== TIMING PROBE (real data; %d workers) ===\n", workers))
i1 <- which(sapply(tasks, function(t) t$design == "internal" &&
                     t$target == "I"))[1]
i2 <- which(sapply(tasks, function(t) t$design == "transfer" &&
                     t$target == "M"))[1]
t0 <- Sys.time()
p1 <- run_task(tasks[[i1]])
ta <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
t0 <- Sys.time()
p2 <- run_task(tasks[[i2]])
tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (!p1$ok || !p2$ok) {
  cat("STOP: probe task failed:", p1$msg, p2$msg, "\n")
  quit(status = 1)
}
pw <- unique(c(p1$warn, p2$warn))
pw <- pw[!grepl("tie-handling", pw, fixed = TRUE)]
cat(sprintf("  probe warnings: %d distinct%s\n", length(pw),
            if (length(pw)) paste0("; first: ", substr(pw[1], 1, 100))
            else ""))
rest <- seq_along(tasks)
cost <- sapply(tasks[rest], function(t) {
  (if (t$design == "transfer") tb else ta) * nrow(t$Xs) /
    nrow(tasks[[if (t$design == "transfer") i2 else i1]]$Xs)
})
rest <- rest[order(-cost)]
cat(sprintf(paste("  internal task %.1f min, transfer task %.1f min;",
                  "all %d tasks: about %.1f hours\n"),
            ta / 60, tb / 60, length(rest),
            max(sum(cost) / workers, max(cost)) * 1.2 / 3600))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
PROG <- file.path(OD, "progress_16d.txt")
cat(sprintf("%s  started, %d tasks\n", format(Sys.time(), "%H:%M:%S"),
            length(rest)), file = PROG)
for (i in rest) tasks[[i]]$progress <- PROG
cat("  progress: one line per finished task in", PROG, "\n")
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks[rest], run_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)
nfail <- sum(!sapply(res, function(r) r$ok))
chk("tasks that failed", nfail, if (nfail == 0) "PASS" else "FAIL")
if (nfail > 0) {
  for (r in Filter(function(r) !r$ok, res)) cat("   ", r$msg, "\n")
  fwrite(rbindlist(CHK), file.path(RS, "a12_checks.csv"))
  quit(status = 1)
}

# ---------------- assemble ----------------
cat("\n=== ASSEMBLE ===\n")
keys <- names(res[[1]]$out)
PRED <- list()
for (f in c("I", "M")) {
  n <- nrow(DS[[f]]$X)
  for (dsg in c("internal", "transfer")) {
    rr <- Filter(function(r) r$design == dsg && r$target == f, res)
    for (k in keys) {
      arr <- array(NA_real_, c(n, length(HZ), length(seeds)))
      for (r in rr) arr[r$rows, , match(r$seed, seeds)] <- r$out[[k]]
      PRED[[f]][[paste(dsg, k, sep = "|")]] <- apply(arr, c(1, 2), mean)
    }
  }
}
CHK6 <- list()
for (f in c("I", "M")) {
  arr <- array(NA_real_, c(nrow(DS[[f]]$X), length(HZ), length(seeds)))
  for (r in Filter(function(r) r$design == "internal" && r$target == f,
                   res)) {
    arr[r$rows, , match(r$seed, seeds)] <- r$check6
  }
  CHK6[[f]] <- apply(arr, c(1, 2), mean)
}
nna <- sum(sapply(PRED, function(x) sum(sapply(x, function(m)
  sum(!is.finite(m))))))
rg <- range(unlist(PRED))
chk("missing predictions", nna, if (nna == 0) "PASS" else "FAIL")
chk("risk range", sprintf("%.5f to %.5f", rg[1], rg[2]),
    if (rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
WRN <- rbindlist(lapply(res, function(r) {
  if (!length(r$warn)) return(NULL)
  data.table(design = r$design, target = r$target,
             message = substr(r$warn, 1, 150))
}))
if (nrow(WRN)) {
  notes <- WRN[grepl("tie-handling", message, fixed = TRUE), .N]
  WRN <- WRN[!grepl("tie-handling", message, fixed = TRUE)]
  wsum <- WRN[, .(count = .N), by = .(design, target, message)]
  setorder(wsum, design, target, -count)
} else {
  notes <- 0L
  wsum <- data.table(design = character(0), target = character(0),
                     message = character(0), count = integer(0))
}
fwrite(wsum, file.path(RS, "a12_warnings.csv"))
chk("warnings from model fitting (see a12_warnings.csv)", sum(wsum$count),
    "INFO")
chk("glmnet tie-handling notices (harmless)", notes, "INFO")
TUN <- rbindlist(lapply(res, function(r) r$tune))
rl <- TUN[learner == "ridge-Cox", mean(ridge_least)]
chk("ridge-Cox choices at the least-penalised end",
    sprintf("%.1f%%", 100 * rl), if (rl <= 0.20) "PASS" else "WARN")
for (g in names(GROUPS)) {
  q <- TUN[learner == "RSF" & group == g]
  chk(sprintf("RSF %s: node size 160 / 5 share; smallest mtry share", g),
      sprintf("%.2f / %.2f; %.2f", mean(q$nodesize == 160),
              mean(q$nodesize == 5), mean(q$mtry_level == 0.5)),
      if (mean(q$nodesize == 160) > 0.5) "WARN" else "PASS")
}
tsum <- TUN[, .(fits = .N,
                share_node5 = mean(nodesize == 5),
                share_node160 = mean(nodesize == 160),
                share_mtry_half = mean(mtry_level == 0.5),
                share_mtry_3 = mean(mtry_level == 3),
                ridge_least = mean(ridge_least),
                ridge_most = mean(ridge_most)),
            by = .(learner, group, design, target)]
fwrite(tsum, file.path(RS, "a12_tuning.csv"))
LONG <- rbindlist(lapply(c("I", "M"), function(f) {
  rbindlist(lapply(names(PRED[[f]]), function(k) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(target = DS[[f]]$name, id = rep(DS[[f]]$id, length(HZ)),
               design = z[1], unit = z[2], learner = z[3],
               horizon = rep(HZ, each = nrow(DS[[f]]$X)),
               risk = as.vector(PRED[[f]][[k]]))
  }))
}))
fwrite(LONG, file.path(OD, "a12_preds.csv"))

# check 1: this stage's IPCW AUC and Brier against riskRegression::Score
D <- DS$I
k <- "transfer|EHR + CTPA|ridge-Cox"
J0 <- list(t = D$t, e = D$e, cuts = HZ + D$off, hz = HZ,
           P = setNames(list(PRED$I[[k]]), k))
mine <- stat_vec(J0, seq_along(D$t))
dmax <- 0
for (hi in c(2, 5)) {
  sc <- Score(list(m = PRED$I[[k]][, hi]), formula = Hist(time, status) ~ 1,
              data = data.frame(time = D$t, status = D$e),
              times = HZ[hi] + D$off, metrics = c("auc", "brier"),
              cens.model = "km", null.model = FALSE)
  a <- as.data.table(sc$AUC$score)$AUC
  b <- as.data.table(sc$Brier$score)[model == "m", Brier]
  dmax <- max(dmax, abs(a - mine[[paste(k, HZ[hi], "auc", sep = "|")]]),
              abs(b - mine[[paste(k, HZ[hi], "brier", sep = "|")]]))
}
chk("IPCW AUC and Brier match riskRegression (INSPECT, max |diff|)",
    sprintf("%.1e", dmax), if (dmax <= 1e-4) "PASS" else "WARN")
# check 2: S6's internal EHR stack reproduced from this stage's ridge-Cox
s6 <- fread(s6f)
for (f in c("I", "M")) {
  D <- DS[[f]]
  for (h in c(30, 365)) {
    k <- "S6 EHR"
    J0 <- list(t = D$t, e = D$e, cuts = h + D$off, hz = h,
               P = setNames(list(CHK6[[f]][, match(h, HZ), drop = FALSE]),
                            k))
    a1 <- stat_vec(J0, seq_along(D$t))[[paste(k, h, "auc", sep = "|")]]
    a0 <- s6[target == D$name & horizon == h & modality == "EHR" &
               model == "internal (refit in target)", auc]
    chk(sprintf("%s: S6 ridge-Cox EHR AUC at %d d reproduced", D$name, h),
        sprintf("%.6f vs %.6f", a1, a0),
        if (length(a0) == 1 && abs(a1 - a0) <= 1e-6) "PASS" else "FAIL")
  }
}

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples per target) ===\n", NBOOT))
MET <- list()
CON <- list()
GAP <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
add_con <- function(lab, est, bv, ...) {
  q <- q2(bv)
  CON[[length(CON) + 1]] <<- data.table(..., contrast = lab,
                                        estimate = est, lo = q[1],
                                        hi = q[2],
                                        excludes_zero = q[1] > 0 | q[2] < 0)
}
plan(multisession, workers = workers)
for (f in c("I", "M")) {
  D <- DS[[f]]
  src <- setdiff(c("I", "M"), f)
  dirl <- sprintf("%s -> %s", DS[[src]]$name, D$name)
  J <- list(t = D$t, e = D$e, cuts = HZ + D$off, hz = HZ, P = PRED[[f]])
  pt <- stat_vec(J, seq_along(D$t))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  target %s done\n", D$name))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      target = D$name, direction = dirl, design = z[1], unit = z[2],
      learner = z[3], horizon = as.numeric(z[4]), measure = z[5],
      estimate = pt[[k]], lo = q[1], hi = q[2])
  }
  for (k in names(PRED[[f]])) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(HZ)) {
      cut <- HZ[hi] + D$off
      r <- pmin(pmax(PRED[[f]][[k]][, hi], 1e-6), 1 - 1e-6)
      th <- pmin(D$t, cut)
      eh <- as.integer(D$e == 1 & D$t <= cut)
      sl <- tryCatch(unname(coef(coxph(Surv(th, eh) ~ log(-log(1 - r))))),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        target = D$name, direction = dirl, design = z[1], unit = z[2],
        learner = z[3], horizon = HZ[hi], measure = "slope",
        estimate = sl, lo = NA_real_, hi = NA_real_)
    }
  }
  kk <- function(dsg, u, l, h, ms) paste(dsg, u, l, h, ms, sep = "|")
  for (u in names(STACKS)) for (h in HZ) for (ms in c("auc", "ipa")) {
    for (l in c("RSF", "rank average")) {
      for (dsg in c("internal", "transfer")) {
        a <- kk(dsg, u, l, h, ms)
        b <- kk(dsg, u, "ridge-Cox", h, ms)
        add_con(sprintf("%s - ridge-Cox, %s", l, dsg), pt[[a]] - pt[[b]],
                B[, a] - B[, b], target = D$name, direction = dirl,
                type = "learner vs ridge-Cox", unit = u, horizon = h,
                measure = ms)
      }
      at <- kk("transfer", u, l, h, ms)
      ai <- kk("internal", u, l, h, ms)
      bt <- kk("transfer", u, "ridge-Cox", h, ms)
      bi <- kk("internal", u, "ridge-Cox", h, ms)
      add_con(sprintf("%s - ridge-Cox, in transfer loss", l),
              (pt[[at]] - pt[[ai]]) - (pt[[bt]] - pt[[bi]]),
              (B[, at] - B[, ai]) - (B[, bt] - B[, bi]),
              target = D$name, direction = dirl,
              type = "difference in transfer loss", unit = u,
              horizon = h, measure = ms)
    }
    for (l in names(LRN)) {
      a <- kk("transfer", u, l, h, ms)
      b <- kk("internal", u, l, h, ms)
      add_con(sprintf("transfer loss, %s", l), pt[[a]] - pt[[b]],
              B[, a] - B[, b], target = D$name, direction = dirl,
              type = "transfer loss", unit = u, horizon = h, measure = ms)
    }
  }
  for (dsg in c("internal", "transfer")) for (l in names(LRN)) {
    for (h in HZ) {
      a <- kk(dsg, "EHR", l, h, "auc")
      b <- kk(dsg, "CTPA", l, h, "auc")
      q <- q2(B[, a] - B[, b])
      GAP[[length(GAP) + 1]] <- data.table(
        target = D$name, direction = dirl, design = dsg, learner = l,
        horizon = as.character(h), gap_ehr_minus_ctpa = pt[[a]] - pt[[b]],
        lo = q[1], hi = q[2])
    }
    a1 <- kk(dsg, "EHR", l, 365, "auc")
    b1 <- kk(dsg, "CTPA", l, 365, "auc")
    a0 <- kk(dsg, "EHR", l, 10, "auc")
    b0 <- kk(dsg, "CTPA", l, 10, "auc")
    bv <- (B[, a1] - B[, b1]) - (B[, a0] - B[, b0])
    q <- q2(bv)
    GAP[[length(GAP) + 1]] <- data.table(
      target = D$name, direction = dirl, design = dsg, learner = l,
      horizon = "change, 365 minus 10",
      gap_ehr_minus_ctpa = (pt[[a1]] - pt[[b1]]) - (pt[[a0]] - pt[[b0]]),
      lo = q[1], hi = q[2])
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
gp <- rbindlist(GAP)
fwrite(mt, file.path(RS, "a12_performance.csv"))
fwrite(cn, file.path(RS, "a12_contrasts.csv"))
fwrite(gp, file.path(RS, "a12_gap.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "a12_checks.csv"))

cat("\n=== FUSED (EHR + CTPA) AUC, internal and transfer ===\n")
print(dcast(mt[unit == "EHR + CTPA" & measure == "auc"],
            direction + learner + design ~ horizon,
            value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                       else x)], row.names = FALSE)
cat("\n=== RSF - ridge-Cox IN TRANSFER LOSS (AUC) ===\n")
print(cn[type == "difference in transfer loss" & measure == "auc" &
           grepl("^RSF", contrast),
         .(direction, unit, horizon, est = round(estimate, 4),
           lo = round(lo, 4), hi = round(hi, 4), excludes_zero)],
      row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(OD, "a12_preds.csv"), " (patient-level, local only)\n")
for (f in c("a12_performance.csv", "a12_contrasts.csv", "a12_gap.csv",
            "a12_tuning.csv", "a12_checks.csv", "a12_warnings.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "a12_log.txt"), " (terminal log)\n")
cat("\nOutputs: a12_performance.csv, a12_contrasts.csv, a12_gap.csv,",
    "a12_tuning.csv, a12_checks.csv, a12_warnings.csv\n")