# 06d_level1_extra.R
# Extra level-1 death scores: retuned ridge-Cox, oblique random survival
# forest and CoxBoost (PREREGISTRATION_ADDENDUM_2.md A16).
#
# Out-of-fold death scores over the stage 4 folds, for every frame, cohort,
# modality and seed. All three learners use stage 18a's preprocessing:
# median imputation, winsorising of non-binary features at the training
# fold's 0.5th and 99.5th percentiles, standardisation, and inner 3-fold
# subject-grouped folds with 18a's seeds.
#
# Ridge-Cox: exactly stage 18a's EHR alone (glmnet, alpha 0, 80 penalties
# down to 1e-7 of the largest (D11), lambda by inner cross-validated
# partial-likelihood deviance, standardize = FALSE), with risks at every
# horizon from the centred Breslow baseline. Oblique RSF: aorsf::orsf, 500
# trees grown in 10 batches of 50 to limit memory, tuned on out-of-bag
# Harrell's C over minimum leaf size 5, 20, 60 x mtry ceiling(sqrt(p) / 2),
# ceiling(sqrt(p)), ceiling(2 sqrt(p)); the score is the out-of-bag or
# predicted mortality, averaged over batches. CoxBoost: CoxBoost::CoxBoost
# with the default penalty (9 x events); steps (0 to 300) chosen by the
# inner cross-validated partial log-likelihood (Verweij), then refitted on
# the whole training fold; the score is the linear predictor.
#
# Version 3: each fold records the number of failed forest fits and the
# first error message, and a saved task with any failed fit or missing score
# is refitted on the next run (same seeds). No performance numbers are
# written. The checks FAIL on a missing score, a D2 per-fold mean C below
# 0.40, a failed forest fit, or EHR ridge-Cox risks that differ from stage
# 18a's EHR alone PH (beyond 1e-6). They WARN when over 20% of ridge-Cox
# choices are at the least-penalised end or CoxBoost choices at 300 steps,
# or over 50% of oblique-forest choices at a grid corner.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, fusion_arch/arch_preds_linear_*.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: results/level1_extra_checks.csv, results/level1_extra_log.txt;
#          data/processed/level1/oof_extra_*.csv and params_extra.csv (not
#          in the repo)
# Usage:   Rscript scripts/06d_level1_extra.R [--probe-only] [--workers=N]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. --workers=N sets the worker
#          count (default: physical cores - 1). Tasks run in the workers, as
#          in stage 18a, so the random inner folds match. Finished tasks are
#          saved in level1/extra_ckpt/ and skipped on a re-run.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(future.apply)
})
for (p in c("aorsf", "CoxBoost")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("package ", p, " is not installed: see the install step")
  }
}
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OUTD <- file.path(PR, "level1")
CKPT <- file.path(OUTD, "extra_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
PROG <- file.path(OUTD, "extra_progress.txt")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
NLAM <- 80L
LMR <- 1e-7
NTREE <- 500L
BATCH <- 50L
LEAF <- c(5L, 20L, 60L)
MTRYF <- c(0.5, 1, 2)
CB_MAX <- 300L
NFOLD <- 5L
PKG <- c("survival", "glmnet", "data.table", "aorsf", "CoxBoost")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 200)
t_start <- Sys.time()
cat(sprintf("aorsf %s, CoxBoost %s; %d workers\n",
            as.character(packageVersion("aorsf")),
            as.character(packageVersion("CoxBoost")), workers))

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                     probe = FALSE) {
  t_task <- proc.time()[["elapsed"]]
  n <- length(tk$t)
  H <- tk$hz
  tt <- tk$t
  ee <- tk$e
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
  cidx <- function(t, e, s) {
    as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                          concordance, error = function(err) NA_real_))
  }
  m <- tk$mod
  cols <- list()
  cols[[sprintf("%s.cox2.death.lp", m)]] <- rep(NA_real_, n)
  for (h in H) cols[[sprintf("%s.cox2.death.r%d", m, h)]] <- rep(NA_real_, n)
  cols[[sprintf("%s.orsf.death.score", m)]] <- rep(NA_real_, n)
  cols[[sprintf("%s.cboost.death.lp", m)]] <- rep(NA_real_, n)
  prm <- list()
  tm <- c(ridge = 0, orsf = 0, cboost = 0)
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    X <- prep(tk$X, tr)
    p <- ncol(X)
    # ridge-Cox, as stage 18a's EHR alone
    t0 <- proc.time()[["elapsed"]]
    cv <- cv.glmnet(X[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                    family = "cox", alpha = 0, foldid = inner,
                    nlambda = NLAM, lambda.min.ratio = LMR,
                    type.measure = "deviance",
                    penalty.factor = rep(1, p), standardize = FALSE)
    l_tr <- as.vector(predict(cv, X[tr, , drop = FALSE], s = "lambda.min"))
    l_te <- as.vector(predict(cv, X[te, , drop = FALSE], s = "lambda.min"))
    cols[[sprintf("%s.cox2.death.lp", m)]][te] <- l_te
    c0 <- max(l_tr)
    for (h in H) {
      hs <- breslow(tt[tr], ee[tr], l_tr, h + 0.5)
      cols[[sprintf("%s.cox2.death.r%d", m, h)]][te] <-
        if (hs == 0) 0 else -expm1(-hs * exp(l_te - c0))
    }
    K <- length(cv$lambda)
    il <- which(cv$lambda == cv$lambda.min)[1]
    tm[["ridge"]] <- tm[["ridge"]] + proc.time()[["elapsed"]] - t0
    # oblique random survival forest, tuned on out-of-bag C
    t0 <- proc.time()[["elapsed"]]
    colnames(X) <- sprintf("f%d", seq_len(p))
    dtr <- data.frame(time = tt[tr], status = ee[tr], X[tr, , drop = FALSE])
    dte <- as.data.frame(X[te, , drop = FALSE])
    grid <- expand.grid(leaf = LEAF,
                        mtry = unique(pmax(1L, pmin(p, ceiling(MTRYF *
                                                                 sqrt(p))))))
    best <- -Inf
    gi <- NA_integer_
    oc <- rep(NA_real_, nrow(grid))
    nfail <- 0L
    emsg <- ""
    for (g in seq_len(nrow(grid))) {
      ob <- list()
      pt <- list()
      for (bb in seq_len(max(1L, NTREE %/% BATCH))) {
        f <- tryCatch(aorsf::orsf(
          dtr, Surv(time, status) ~ ., n_tree = min(NTREE, BATCH),
          mtry = grid$mtry[g], leaf_min_obs = grid$leaf[g],
          split_min_obs = 2L * grid$leaf[g], oobag_pred_type = "mort",
          tree_seeds = as.integer(tk$seed * 100000 + k * 1000 + g * 20 +
                                    bb), n_thread = 1),
          error = function(e) conditionMessage(e))
        if (is.character(f)) {
          nfail <- nfail + 1L
          if (!nzchar(emsg)) emsg <- substr(f, 1, 120)
          break
        }
        ob[[bb]] <- as.vector(f$pred_oobag)
        pt[[bb]] <- as.vector(predict(f, new_data = dte, pred_type = "mort"))
        rm(f)
        invisible(gc(verbose = FALSE))
      }
      if (length(pt) < max(1L, NTREE %/% BATCH)) next
      oo <- rowMeans(do.call(cbind, ob), na.rm = TRUE)
      oc[g] <- cidx(tt[tr], ee[tr], oo)
      if (is.finite(oc[g]) && oc[g] > best) {
        best <- oc[g]
        gi <- g
        cols[[sprintf("%s.orsf.death.score", m)]][te] <-
          rowMeans(do.call(cbind, pt))
      }
    }
    tm[["orsf"]] <- tm[["orsf"]] + proc.time()[["elapsed"]] - t0
    # CoxBoost, steps by inner cross-validated partial likelihood
    t0 <- proc.time()[["elapsed"]]
    pen <- 9 * sum(ee[tr])
    fl <- lapply(1:3, function(j) which(inner == j))
    cvb <- CoxBoost::cv.CoxBoost(time = tt[tr], status = ee[tr],
                                 x = X[tr, , drop = FALSE],
                                 maxstepno = CB_MAX, K = 3,
                                 type = "verweij", folds = fl,
                                 penalty = pen, trace = FALSE)
    st <- cvb$optimal.step
    fb <- CoxBoost::CoxBoost(time = tt[tr], status = ee[tr],
                             x = X[tr, , drop = FALSE], stepno = st,
                             penalty = pen)
    cols[[sprintf("%s.cboost.death.lp", m)]][te] <-
      as.vector(predict(fb, newdata = X[te, , drop = FALSE], type = "lp"))
    tm[["cboost"]] <- tm[["cboost"]] + proc.time()[["elapsed"]] - t0
    prm[[length(prm) + 1]] <- data.table(
      fold = k, ridge_least = il == K, ridge_most = il == 1,
      orsf_leaf = if (is.na(gi)) NA_integer_ else grid$leaf[gi],
      orsf_mtry = if (is.na(gi)) NA_integer_ else grid$mtry[gi],
      orsf_mtry_level = if (is.na(gi)) NA_real_ else
        match(grid$mtry[gi], sort(unique(grid$mtry))),
      orsf_oob = paste(sprintf("%.4f", oc), collapse = ";"),
      orsf_failed = nfail, orsf_error = emsg, cboost_steps = st)
  }
  out <- list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
              cols = as.data.table(cols), prm = rbindlist(prm), tm = tm)
  if (!probe) {
    cat(sprintf("%s  done  %s %-7s %-4s seed %2d  %5.1f min\n",
                format(Sys.time(), "%H:%M:%S"), tk$fr, tk$co, tk$mod,
                tk$seed, (proc.time()[["elapsed"]] - t_task) / 60),
        file = tk$prog, append = TRUE)
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, tk$ckpt)
  }
  out
}

# ---------------- build tasks ----------------
cat("=== BUILD ===\n")
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
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
  fmap <- list()
  for (s in seeds) {
    fs <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(fs)) stop("fold missing: ", fr, " ", co, " seed ", s)
    fmap[[as.character(s)]] <- as.integer(fs)
  }
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  for (m in mods) {
    cc <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cc])
    storage.mode(Xm) <- "double"
    for (s in seeds) {
      tasks[[length(tasks) + 1]] <- list(
        fr = fr, co = co, mod = m, seed = s, X = Xm,
        t = sub$time_death, e = as.integer(sub$event_death == 1),
        subj = sub$subject_id, fold = fmap[[as.character(s)]],
        hz = HZ[[fr]], prog = PROG,
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", fr, code, m, s)))
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                seeds = seeds, mods = mods, fmap = fmap)
  cat(sprintf("  %s %-7s n=%d  modalities %s  seeds %d\n", fr, co,
              nrow(sub), paste(mods, collapse = "/"), length(seeds)))
}

# ---------------- timing probe, or resume ----------------
RESUME <- file.exists(PLANF)
if (RESUME) {
  PL <- readRDS(PLANF)
  costs <- PL$costs
  cat(sprintf("\n=== RESUMING (plan of %s); probe skipped ===\n", PL$when))
} else {
  cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
  kind <- sapply(tasks, function(t) paste(t$fr, t$co, t$mod))
  reps <- tasks[!duplicated(kind)]
  names(reps) <- unique(kind)
  plan(multisession, workers = min(workers, length(reps)))
  t0 <- Sys.time()
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  fold_s <- sapply(pres, function(r) sum(r$tm))
  for (kd in names(pres)) {
    r <- pres[[kd]]
    cat(sprintf(paste("  %-16s ridge %5.1f s, oblique forest %6.1f s,",
                      "CoxBoost %6.1f s per fold\n"), kd,
                r$tm[["ridge"]], r$tm[["orsf"]], r$tm[["cboost"]]))
  }
  costs <- fold_s[kind] * NFOLD
  hrs <- max(sum(costs) / workers, max(costs)) * 1.1 / 3600
  cat(sprintf("  probe took %.1f min; %d tasks on %d workers: about %.1f",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              length(tasks), workers, hrs), "hours\n")
  if (probe_only) {
    cat("--probe-only: stopping here (nothing saved).\n")
    quit(status = 0)
  }
  saveRDS(list(costs = costs, when = format(Sys.time(), "%Y-%m-%d %H:%M")),
          PLANF)
}
# a saved task counts as finished only if every score is present and no
# oblique-forest fit failed (older saves have no failure count: a missing
# out-of-bag C then marks a failed fit)
ck_ok <- function(f) {
  if (!file.exists(f)) return(NA)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(FALSE)
  nf <- if ("orsf_failed" %in% names(r$prm)) sum(r$prm$orsf_failed) else
    sum(grepl("NA", r$prm$orsf_oob))
  all(sapply(r$cols, function(v) all(is.finite(v)))) && nf == 0
}
okv <- sapply(tasks, function(t) ck_ok(t$ckpt))
bad <- which(okv %in% FALSE)
if (length(bad)) {
  cat(sprintf("  %d saved task(s) had failed forest fits or missing scores",
              length(bad)), "and will be refitted:\n")
  for (i in bad) {
    cat(sprintf("    %s %s %s seed %d\n", tasks[[i]]$fr, tasks[[i]]$co,
                tasks[[i]]$mod, tasks[[i]]$seed))
  }
}
done <- which(okv %in% TRUE)
todo <- setdiff(seq_along(tasks), done)
proj <- max(sum(costs[todo]) / workers, max(c(0, costs[todo]))) * 1.1
cat(sprintf(paste("  %d of %d tasks already saved; %d to fit, about %.1f",
                  "hours on %d workers\n"), length(done), length(tasks),
            length(todo), proj / 3600, workers))
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- full run ----------------
cat("\n=== FITTING (progress: extra_progress.txt) ===\n")
cat(sprintf("%s  %s %d tasks, %d workers, about %.1f hours\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            if (length(done)) "resumed with" else "started", length(todo),
            workers, proj / 3600), file = PROG, append = length(done) > 0)
if (length(todo)) {
  ord <- todo[order(-costs[todo])]
  plan(multisession, workers = workers)
  invisible(future_lapply(tasks[ord], run_task, future.packages = PKG,
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
  fwrite(out, file.path(OUTD, sprintf("oof_extra_%s_%s.csv", M$fr,
                                      M$code)))
  nmiss <- sum(!is.finite(as.matrix(out[, ..pc])))
  chk(sprintf("%s: rows / columns / missing cells", key),
      sprintf("%d / %d / %d", nrow(out), length(pc), nmiss),
      if (nrow(out) == nrow(M$sub) * length(M$seeds) && nmiss == 0)
        "PASS" else "FAIL")
  if (nmiss > 0) {
    mc <- sapply(pc, function(cc) sum(!is.finite(out[[cc]])))
    chk(sprintf("%s: columns with missing cells", key),
        paste(sprintf("%s %d", names(mc)[mc > 0], mc[mc > 0]),
              collapse = "; "), "FAIL")
  }
  # D2 per-fold rule on every score
  td <- M$sub$time_death
  ed <- as.integer(M$sub$event_death == 1)
  sc <- grep("\\.(lp|score)$", pc, value = TRUE)
  nbad <- 0L
  for (cc in sc) {
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
  # reproduction of stage 18a's EHR alone PH
  lf <- file.path(PR, "fusion_arch", sprintf("arch_preds_linear_%s_%s.csv",
                                             M$fr, M$code))
  if (!file.exists(lf)) {
    chk(sprintf("%s: stage 18a file for the reproduction check", key),
        "missing", "WARN")
  } else {
    a18 <- fread(lf)[model == "EHR alone PH" & set == "EHR"]
    av <- out[, lapply(.SD, mean), by = hadm_id,
              .SDcols = sprintf("ehr.cox2.death.r%d", HZ[[M$fr]])][
                match(M$sub$hadm_id, hadm_id)]
    dm <- max(sapply(HZ[[M$fr]], function(h) {
      x <- a18[horizon == h]
      max(abs(x$risk[match(M$sub$hadm_id, x$hadm_id)] -
                av[[sprintf("ehr.cox2.death.r%d", h)]]))
    }))
    chk(sprintf("%s: EHR ridge-Cox reproduces stage 18a (max |diff|)", key),
        sprintf("%.1e", dm), if (is.finite(dm) && dm < 1e-6) "PASS" else
          "FAIL")
  }
}
P <- rbindlist(PRM, fill = TRUE)
fwrite(P, file.path(OUTD, "params_extra.csv"))
nff <- sum(P$orsf_failed, na.rm = TRUE)
em <- unique(P$orsf_error[!is.na(P$orsf_error) & nzchar(P$orsf_error)])
chk("oblique forest fits that failed (first message)",
    sprintf("%d%s", nff, if (length(em)) paste0(" (", em[1], ")") else ""),
    if (nff == 0) "PASS" else "FAIL")
for (m in c("ehr", "ecg", "ctpa")) {
  q <- P[modality == m]
  if (!nrow(q)) next
  rl <- mean(q$ridge_least)
  chk(sprintf("%s ridge-Cox: choices at the least-penalised end",
              toupper(m)), sprintf("%.1f%%", 100 * rl),
      if (rl <= 0.20) "PASS" else "WARN")
  cm <- mean(q$cboost_steps >= CB_MAX)
  chk(sprintf("%s CoxBoost: optimal steps at %d / at 0 / median",
              toupper(m), CB_MAX),
      sprintf("%.1f%% / %.1f%% / %g", 100 * cm,
              100 * mean(q$cboost_steps == 0), median(q$cboost_steps)),
      if (cm <= 0.20) "PASS" else "WARN")
  e1 <- max(mean(q$orsf_leaf == min(LEAF), na.rm = TRUE),
            mean(q$orsf_leaf == max(LEAF), na.rm = TRUE))
  e2 <- max(mean(q$orsf_mtry_level == 1, na.rm = TRUE),
            mean(q$orsf_mtry_level == 3, na.rm = TRUE))
  chk(sprintf("%s oblique forest: leaf 5 / 60, smallest / largest mtry",
              toupper(m)),
      sprintf("%.2f / %.2f, %.2f / %.2f",
              mean(q$orsf_leaf == min(LEAF), na.rm = TRUE),
              mean(q$orsf_leaf == max(LEAF), na.rm = TRUE),
              mean(q$orsf_mtry_level == 1, na.rm = TRUE),
              mean(q$orsf_mtry_level == 3, na.rm = TRUE)),
      if (max(e1, e2) > 0.5) "WARN" else "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_extra_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(OUTD, sprintf("oof_extra_%s_%s.csv", META[[key]]$fr,
                                    META[[key]]$code)), " (local only)\n")
}
cat("  ", file.path(OUTD, "params_extra.csv"), " (local only)\n")
cat("  ", file.path(RS, "level1_extra_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_extra_log.txt"), " (terminal log)\n")
cat("\nOutputs: level1_extra_checks.csv\n")