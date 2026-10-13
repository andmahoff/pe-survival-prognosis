# 06h_level1_enet_cboost.R
# Elastic-net Cox and longer CoxBoost death scores
# (PREREGISTRATION_ADDENDUM_2.md A22, items 2 and 4).
#
# Elastic-net Cox death scores, and CoxBoost with up to 500 steps where
# stage 6d's CoxBoost hit its 300-step cap. Every frame, cohort, modality
# and seed, out of fold over the stage 4 folds, with stage 18a's
# preprocessing and inner folds (as stage 6d).
#
# Elastic-net Cox (item 2): glmnet Cox at alpha 0 (ridge), 0.25, 0.5, 0.75
# and 1 (lasso), each over 80 penalties down to 1e-7 of the largest. The
# alpha and lambda with the lowest inner cross-validated partial-likelihood
# deviance are kept, with the linear predictor and risks at every horizon
# from the centred Breslow baseline. Where alpha 0 is chosen, the fit is
# exactly stage 6d's ridge-Cox.
#
# CoxBoost 500 (item 4): for every fold in which stage 6d's CoxBoost chose
# 300 steps (the cap), CoxBoost is refitted with steps 0 to 500 chosen by
# the same inner cross-validated partial log-likelihood. Every other fold
# keeps stage 6d's fit, because its optimum lay inside the old range. No
# performance numbers are written. The checks FAIL if a patient is not
# scored, if the score where alpha 0 was chosen or in an unchanged CoxBoost
# fold differs from stage 6d (beyond 1e-9) or if the D2 per-fold rule fails.
# They also report alpha shares, and WARN when over 20% of CoxBoost refits
# are still at 500 steps.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_extra_*.csv,
#          level1/params_extra.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/level1_enet_checks.csv, results/level1_enet_log.txt;
#          data/processed/level1/oof_enet_*.csv and params_enet.csv (not in
#          the repo)
# Usage:   Rscript scripts/06h_level1_enet_cboost.R [--probe-only]
#          [--workers=N]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. --workers=N sets the worker
#          count (default: physical cores - 1). Tasks run in the workers, as
#          in stage 6d, so the random inner folds match. Finished tasks are
#          saved in level1/enet_ckpt/ and skipped on a re-run.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(future.apply)
})
if (!requireNamespace("CoxBoost", quietly = TRUE)) stop("CoxBoost missing")
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
OUTD <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CKPT <- file.path(OUTD, "enet_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
NLAM <- 80L
LMR <- 1e-7
ALPHAS <- c(0, 0.25, 0.5, 0.75, 1)
CB_OLD <- 300L
CB_NEW <- 500L
PKG <- c("survival", "glmnet", "data.table", "CoxBoost")
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

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  tt <- tk$t
  ee <- tk$e
  m <- tk$mod
  prep <- function(X, tr) {
    for (j in seq_len(ncol(X))) {
      md <- median(X[tr, j], na.rm = TRUE)
      if (!is.finite(md)) md <- 0
      X[is.na(X[, j]), j] <- md
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
  cols <- list()
  cols[[sprintf("%s.enet.death.lp", m)]] <- rep(NA_real_, n)
  for (h in H) cols[[sprintf("%s.enet.death.r%d", m, h)]] <- rep(NA_real_, n)
  cols[[sprintf("%s.cb500.death.lp", m)]] <- tk$cb_old
  prm <- list()
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    X <- prep(tk$X, tr)
    best <- NULL
    bdev <- Inf
    ba <- NA_real_
    for (a in tk$alphas) {
      cv <- cv.glmnet(X[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                      family = "cox", alpha = a, foldid = inner,
                      nlambda = NLAM, lambda.min.ratio = LMR,
                      type.measure = "deviance",
                      penalty.factor = rep(1, ncol(X)), standardize = FALSE)
      dv <- min(cv$cvm)
      if (dv < bdev) {
        bdev <- dv
        best <- cv
        ba <- a
      }
    }
    l_tr <- as.vector(predict(best, X[tr, , drop = FALSE], s = "lambda.min"))
    l_te <- as.vector(predict(best, X[te, , drop = FALSE], s = "lambda.min"))
    cols[[sprintf("%s.enet.death.lp", m)]][te] <- l_te
    c0 <- max(l_tr)
    for (h in H) {
      hs <- breslow(tt[tr], ee[tr], l_tr, h + 0.5)
      cols[[sprintf("%s.enet.death.r%d", m, h)]][te] <-
        if (hs == 0) 0 else -expm1(-hs * exp(l_te - c0))
    }
    nz <- sum(as.vector(coef(best, s = "lambda.min")) != 0)
    st <- NA_integer_
    if (k %in% tk$cap_folds) {
      pen <- 9 * sum(ee[tr])
      fl <- lapply(1:3, function(j) which(inner == j))
      cvb <- CoxBoost::cv.CoxBoost(time = tt[tr], status = ee[tr],
                                   x = X[tr, , drop = FALSE],
                                   maxstepno = tk$cb_new, K = 3,
                                   type = "verweij", folds = fl,
                                   penalty = pen, trace = FALSE)
      st <- cvb$optimal.step
      fb <- CoxBoost::CoxBoost(time = tt[tr], status = ee[tr],
                               x = X[tr, , drop = FALSE], stepno = st,
                               penalty = pen)
      cols[[sprintf("%s.cb500.death.lp", m)]][te] <-
        as.vector(predict(fb, newdata = X[te, , drop = FALSE], type = "lp"))
    }
    prm[[length(prm) + 1]] <- data.table(
      fold = k, alpha = ba, nonzero = nz, p = ncol(X),
      cb_refit = k %in% tk$cap_folds, cb_steps = st)
  }
  out <- list(fr = tk$fr, co = tk$co, mod = m, seed = tk$seed,
              cols = as.data.table(cols), prm = rbindlist(prm))
  tmp <- paste0(tk$ckpt, ".tmp")
  saveRDS(out, tmp)
  file.rename(tmp, tk$ckpt)
  out
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
P6 <- fread(file.path(OUTD, "params_extra.csv"))
mats <- lapply(MATS, function(f) unique(fread(file.path(IN, f)),
                                        by = "hadm_id"))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  Fd <- if (fr == "A") fa else fb
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  X6 <- fread(file.path(OUTD, sprintf("oof_extra_%s_%s.csv", fr, code)))
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  fmap <- list()
  for (s in seeds) {
    fmap[[as.character(s)]] <- as.integer(
      fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold])
  }
  for (m in mods) {
    cc <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cc])
    storage.mode(Xm) <- "double"
    for (s in seeds) {
      x6 <- X6[X6$seed == s][match(sub$hadm_id, hadm_id)]
      capf <- P6[frame == fr & cohort == co & modality == m & seed == s &
                   cboost_steps >= CB_OLD, fold]
      tasks[[length(tasks) + 1]] <- list(
        fr = fr, co = co, mod = m, seed = s, X = Xm, t = sub$time_death,
        e = as.integer(sub$event_death == 1), subj = sub$subject_id,
        fold = fmap[[as.character(s)]], hz = HZ[[fr]], alphas = ALPHAS,
        cb_old = x6[[sprintf("%s.cboost.death.lp", m)]],
        cap_folds = capf, cb_new = CB_NEW,
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", fr, code, m, s)))
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                seeds = seeds, mods = mods, fmap = fmap,
                                X6 = X6)
}
ncap <- sum(sapply(tasks, function(t) length(t$cap_folds)))
cat(sprintf("  %d tasks; %d folds where stage 6d's CoxBoost hit %d steps\n",
            length(tasks), ncap, CB_OLD))
done <- which(sapply(tasks, function(t) file.exists(t$ckpt)))
todo <- setdiff(seq_along(tasks), done)

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
if (length(todo)) {
  cost <- sapply(tasks[todo], function(t) {
    nrow(t$X) * ncol(t$X) * (1 + 3 * length(t$cap_folds))
  })
  big <- todo[which.max(cost)]
  tb <- tasks[[big]]
  tb$ckpt <- tempfile(fileext = ".rds")
  t0 <- Sys.time()
  invisible(run_task(tb))
  tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf(paste("  largest task %.0f s; %d of %d tasks to fit on %d",
                    "workers: at most about %.1f min\n"), tp, length(todo),
              length(tasks), workers,
              tp * ceiling(length(todo) / workers) * 1.2 / 60))
} else {
  cat("  every task already saved\n")
}
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
if (length(todo)) {
  plan(multisession, workers = workers)
  invisible(future_lapply(tasks[todo], run_task, future.packages = PKG,
                          future.seed = TRUE, future.chunk.size = 1L))
  plan(sequential)
}
res <- lapply(tasks, function(t) readRDS(t$ckpt))

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE AND CHECK ===\n")
cidx_c <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
PRM <- list()
for (key in names(META)) {
  M <- META[[key]]
  parts <- list()
  for (s in M$seeds) {
    rr <- Filter(function(r) r$fr == M$fr && r$co == M$co &&
                   r$seed == s, res)
    rr <- rr[match(M$mods, sapply(rr, function(r) r$mod))]
    parts[[length(parts) + 1]] <- cbind(
      data.table(seed = s, subject_id = M$sub$subject_id,
                 hadm_id = M$sub$hadm_id,
                 fold = M$fmap[[as.character(s)]]),
      do.call(cbind, lapply(rr, function(r) r$cols)))
    for (r in rr) {
      PRM[[length(PRM) + 1]] <- cbind(
        data.table(frame = r$fr, cohort = r$co, modality = r$mod,
                   seed = r$seed), r$prm)
    }
  }
  out <- rbindlist(parts)
  pc <- setdiff(names(out), c("seed", "subject_id", "hadm_id", "fold"))
  fwrite(out, file.path(OUTD, sprintf("oof_enet_%s_%s.csv", M$fr, M$code)))
  nmiss <- sum(!is.finite(as.matrix(out[, ..pc])))
  chk(sprintf("%s: rows / columns / missing cells", key),
      sprintf("%d / %d / %d", nrow(out), length(pc), nmiss),
      if (nrow(out) == nrow(M$sub) * length(M$seeds) && nmiss == 0)
        "PASS" else "FAIL")
  x6 <- M$X6[match(paste(out$seed, out$hadm_id),
                   paste(M$X6$seed, M$X6$hadm_id))]
  pk <- rbindlist(PRM)[frame == M$fr & cohort == M$co]
  d0 <- 0
  n0 <- 0L
  dc <- 0
  for (m in M$mods) {
    for (s in M$seeds) {
      a0 <- pk[modality == m & seed == s & alpha == 0, fold]
      kp <- pk[modality == m & seed == s & cb_refit == FALSE, fold]
      ii <- which(out$seed == s & out$fold %in% a0)
      if (length(ii)) {
        d0 <- max(d0, max(abs(out[[sprintf("%s.enet.death.lp", m)]][ii] -
                                x6[[sprintf("%s.cox2.death.lp", m)]][ii])))
        n0 <- n0 + length(ii)
      }
      jj <- which(out$seed == s & out$fold %in% kp)
      if (length(jj)) {
        dc <- max(dc, max(abs(out[[sprintf("%s.cb500.death.lp", m)]][jj] -
                                x6[[sprintf("%s.cboost.death.lp", m)]][jj])))
      }
    }
  }
  chk(sprintf("%s: alpha-0 folds equal stage 6d ridge-Cox (rows, max)",
              key), sprintf("%d, %.1e", n0, d0),
      if (n0 == 0) "WARN" else if (d0 <= 1e-9) "PASS" else "FAIL")
  chk(sprintf("%s: unchanged CoxBoost folds equal stage 6d (max)", key),
      sprintf("%.1e", dc), if (dc <= 1e-9) "PASS" else "FAIL")
  td <- M$sub$time_death
  ed <- as.integer(M$sub$event_death == 1)
  nbad <- 0L
  for (cc in grep("\\.lp$", pc, value = TRUE)) {
    v <- c()
    for (s in M$seeds) {
      os <- out[out$seed == s]
      for (k in sort(unique(os$fold))) {
        ii <- which(os$fold == k)
        v <- c(v, cidx_c(td[ii], ed[ii], os[[cc]][ii]))
      }
    }
    if (!is.finite(mean(v, na.rm = TRUE)) || mean(v, na.rm = TRUE) < 0.40) {
      nbad <- nbad + 1L
    }
  }
  chk(sprintf("%s: scores failing the per-fold rule (D2)", key), nbad,
      if (nbad == 0) "PASS" else "FAIL")
}
P <- rbindlist(PRM)
fwrite(P, file.path(OUTD, "params_enet.csv"))
for (m in c("ehr", "ecg", "ctpa")) {
  q <- P[modality == m]
  if (!nrow(q)) next
  sh <- sapply(ALPHAS, function(a) mean(q$alpha == a))
  chk(sprintf("%s elastic net: share of alpha 0 / .25 / .5 / .75 / 1",
              toupper(m)), paste(sprintf("%.2f", sh), collapse = " / "),
      "INFO")
  chk(sprintf("%s elastic net: median non-zero coefficients of p",
              toupper(m)), sprintf("%g of %g", median(q$nonzero),
                                   median(q$p)), "INFO")
  r <- q[cb_refit == TRUE]
  if (nrow(r)) {
    s5 <- mean(r$cb_steps >= CB_NEW)
    chk(sprintf("%s CoxBoost: folds refitted / median steps / at %d",
                toupper(m), CB_NEW),
        sprintf("%d / %g / %.0f%%", nrow(r), median(r$cb_steps), 100 * s5),
        if (s5 <= 0.20) "PASS" else "WARN")
  }
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_enet_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(OUTD, sprintf("oof_enet_%s_%s.csv", META[[key]]$fr,
                                    META[[key]]$code)), " (local only)\n")
}
cat("  ", file.path(OUTD, "params_enet.csv"), " (local only)\n")
cat("  ", file.path(RS, "level1_enet_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_enet_log.txt"), " (terminal log)\n")
cat("\nOutputs: level1_enet_checks.csv\n")