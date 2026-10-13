# 06g_level1_forests_proper.R
# Death scores from forests tuned on a proper score
# (PREREGISTRATION_ADDENDUM_2.md A22, items 1 and 4).
#
# The random survival forest (RSF) and the oblique RSF are tuned on the
# out-of-bag Brier score summed over the frame's horizons as well as on
# out-of-bag Harrell's C, and the oblique forest also gets a wider mtry
# grid. Every frame, cohort, modality and seed, out of fold over the stage 4
# folds.
#
# RSF is exactly stage 6c's death forest (same learner, grid and seeds),
# with each grid point scored twice: rsfc is chosen by out-of-bag C
# (reproduces stage 6c exactly) and rsfb by out-of-bag Brier. Oblique RSF is
# exactly stage 6d's (same learner, preprocessing and seeds), with mtry
# extended to ceiling(4 sqrt(p)) because the EHR forests chose the largest
# value in 91% of folds; leaf size is 5, 20 or 60. Its choices are orsfc9
# (out-of-bag C over stage 6d's 9 grid points, reproduces stage 6d exactly),
# orsfc (out-of-bag C over all 12 points, item 4) and orsfb (out-of-bag
# Brier over all 12, items 1 and 4).
#
# The out-of-bag Brier score takes each training patient's out-of-bag risk
# of death by h + 0.5 days at every horizon h of the frame against whether
# they died by then (follow-up is complete), with squared errors summed over
# the horizons. No performance numbers are written. The checks FAIL on a
# missing score, a failed forest fit, a D2 per-fold rule failure, or rsfc or
# orsfc9 differing from stage 6c's or 6d's score (beyond 1e-9). They also
# report how often the Brier and C choices differ, and WARN when choices sit
# at a grid edge more than half the time.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_forest_wide_*.csv,
#          level1/oof_extra_*.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/level1_proper_checks.csv, results/level1_proper_log.txt;
#          data/processed/level1/oof_proper_*.csv and params_proper.csv (not
#          in the repo)
# Usage:   Rscript scripts/06g_level1_forests_proper.R [--probe-only]
#          [--workers=N]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. --workers=N sets the worker
#          count (default: physical cores - 1). Finished tasks are saved in
#          level1/proper_ckpt/ and skipped on a re-run.

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf not installed")
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
OUTD <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CKPT <- file.path(OUTD, "proper_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
PROG <- file.path(OUTD, "proper_progress.txt")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5L, 15L, 40L, 80L, 160L)
MTRYF <- c(0.5, 1, 3)
LEAF <- c(5L, 20L, 60L)
OMTRYF <- c(0.5, 1, 2, 4)
ONTREE <- 500L
BATCH <- 50L
NFOLD <- 5L
PKG <- c("survival", "data.table", "randomForestSRC", "aorsf")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 200)
t_start <- Sys.time()
cat(sprintf("randomForestSRC %s, aorsf %s; %d workers\n",
            as.character(packageVersion("randomForestSRC")),
            as.character(packageVersion("aorsf")), workers))

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                     probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  cidx <- function(t, e, s) {
    as.numeric(tryCatch(survival::concordance(
      survival::Surv(t, e) ~ s, reverse = TRUE)$concordance,
      error = function(err) NA_real_))
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
  n <- length(tk$t)
  H <- tk$hz
  tt <- tk$t
  ee <- tk$e
  m <- tk$mod
  CH <- c("rsfc", "rsfb", "orsfc9", "orsfc", "orsfb")
  cols <- list()
  for (ch in CH) {
    cols[[sprintf("%s.%s.death.score", m, ch)]] <- rep(NA_real_, n)
    for (h in H) {
      cols[[sprintf("%s.%s.death.r%d", m, ch, h)]] <- rep(NA_real_, n)
    }
  }
  put <- function(ch, te, sc, rk) {
    cols[[sprintf("%s.%s.death.score", m, ch)]][te] <<- sc
    for (i in seq_along(H)) {
      cols[[sprintf("%s.%s.death.r%d", m, ch, H[i])]][te] <<- rk[, i]
    }
  }
  prm <- list()
  nfail <- 0L
  tm <- c(rsf = 0, orsf = 0)
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    Ytr <- sapply(H, function(h) as.numeric(ee[tr] == 1 &
                                               tt[tr] <= h + 0.5))
    Ytr <- matrix(Ytr, ncol = length(H))
    # ---- RSF, exactly as stage 6c's death forests ----
    t0 <- proc.time()[["elapsed"]]
    X <- tk$X
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    Xi <- X
    na <- which(is.na(Xi), arr.ind = TRUE)
    if (nrow(na) > 0) Xi[na] <- med[na[, 2]]
    Xte <- as.data.frame(Xi[te, , drop = FALSE])
    ev <- sort(unique(tt[tr][ee[tr] == 1]))
    at_h <- sapply(H, function(h) {
      v <- ev[ev <= h + 0.5]
      if (length(v)) max(v) else NA_real_
    })
    pick <- ev[unique(round(seq(1, length(ev),
                                length.out = min(tk$ntime, length(ev)))))]
    nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
    dtr <- data.frame(time = tt[tr], status = ee[tr],
                      Xi[tr, , drop = FALSE])
    p <- ncol(X)
    grid <- expand.grid(nodesize = tk$nodes,
                        mtry = unique(pmin(p, ceiling(tk$mtryf * sqrt(p)))))
    oc <- rep(NA_real_, nrow(grid))
    ob <- rep(NA_real_, nrow(grid))
    PT <- list()
    for (g in seq_len(nrow(grid))) {
      f <- tryCatch(
        rfsrc(Surv(time, status) ~ ., data = dtr, ntree = tk$ntree,
              nodesize = grid$nodesize[g], mtry = grid$mtry[g],
              nsplit = tk$nsplit, ntime = nt, splitrule = "logrank",
              save.memory = TRUE,
              seed = -as.integer(tk$seed * 10000 + k * 100 + g)),
        error = function(e) NULL)
      if (is.null(f)) {
        nfail <- nfail + 1L
        next
      }
      oc[g] <- cidx(tt[tr], ee[tr], f$predicted.oob)
      ro <- sapply(H, function(h) {
        i <- findInterval(h + 0.5, f$time.interest)
        if (i == 0) rep(0, length(tr)) else 1 - f$survival.oob[, i]
      })
      ob[g] <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2))
      pr <- tryCatch(predict(f, newdata = Xte), error = function(e) NULL)
      if (is.null(pr)) {
        nfail <- nfail + 1L
        oc[g] <- NA_real_
        ob[g] <- NA_real_
      } else {
        rk <- sapply(H, function(h) {
          i <- findInterval(h + 0.5, pr$time.interest)
          if (i == 0) rep(0, length(te)) else 1 - pr$survival[, i]
        })
        PT[[g]] <- list(sc = pr$predicted,
                        rk = matrix(rk, ncol = length(H)))
      }
      rm(f)
      invisible(gc(verbose = FALSE))
    }
    gc_ <- if (any(is.finite(oc))) which.max(oc) else NA_integer_
    gb_ <- if (any(is.finite(ob))) which.min(ob) else NA_integer_
    if (!is.na(gc_)) put("rsfc", te, PT[[gc_]]$sc, PT[[gc_]]$rk)
    if (!is.na(gb_)) put("rsfb", te, PT[[gb_]]$sc, PT[[gb_]]$rk)
    tm[["rsf"]] <- tm[["rsf"]] + proc.time()[["elapsed"]] - t0
    # ---- oblique RSF, as stage 6d, wider mtry grid ----
    t0 <- proc.time()[["elapsed"]]
    Xs <- prep(tk$X, tr)
    ps <- ncol(Xs)
    colnames(Xs) <- sprintf("f%d", seq_len(ps))
    otr <- data.frame(time = tt[tr], status = ee[tr], Xs[tr, , drop = FALSE])
    ote <- as.data.frame(Xs[te, , drop = FALSE])
    og <- expand.grid(leaf = tk$leaf,
                      mtry = unique(pmax(1L, pmin(ps, ceiling(tk$omtryf *
                                                                sqrt(ps))))))
    n6 <- length(tk$leaf) * length(unique(pmax(1L, pmin(
      ps, ceiling(tk$omtryf[1:3] * sqrt(ps))))))
    occ <- rep(NA_real_, nrow(og))
    obb <- rep(NA_real_, nrow(og))
    OT <- list()
    nb <- max(1L, tk$ontree %/% tk$batch)
    for (g in seq_len(nrow(og))) {
      mo <- list()
      ro <- list()
      ms <- list()
      rs <- list()
      for (bb in seq_len(nb)) {
        f <- tryCatch(aorsf::orsf(
          otr, Surv(time, status) ~ ., n_tree = min(tk$ontree, tk$batch),
          mtry = og$mtry[g], leaf_min_obs = og$leaf[g],
          split_min_obs = 2L * og$leaf[g], oobag_pred_type = "mort",
          tree_seeds = as.integer(tk$seed * 100000 + k * 1000 + g * 20 +
                                    bb), n_thread = 1),
          error = function(e) NULL)
        if (is.null(f)) {
          nfail <- nfail + 1L
          break
        }
        mo[[bb]] <- as.vector(f$pred_oobag)
        ro[[bb]] <- matrix(predict(f, oobag = TRUE, pred_type = "risk",
                                   pred_horizon = H + 0.5,
                                   boundary_checks = FALSE),
                           ncol = length(H))
        ms[[bb]] <- as.vector(predict(f, new_data = ote, pred_type = "mort"))
        rs[[bb]] <- matrix(predict(f, new_data = ote, pred_type = "risk",
                                   pred_horizon = H + 0.5,
                                   boundary_checks = FALSE),
                           ncol = length(H))
        rm(f)
        invisible(gc(verbose = FALSE))
      }
      if (length(ms) < nb) next
      occ[g] <- cidx(tt[tr], ee[tr],
                     rowMeans(do.call(cbind, mo), na.rm = TRUE))
      rbar <- Reduce(`+`, ro) / nb
      obb[g] <- sum(colMeans((Ytr - rbar)^2, na.rm = TRUE))
      OT[[g]] <- list(sc = rowMeans(do.call(cbind, ms)),
                      rk = Reduce(`+`, rs) / nb)
    }
    g9 <- if (any(is.finite(occ[seq_len(n6)]))) {
      which.max(occ[seq_len(n6)])
    } else NA_integer_
    gcw <- if (any(is.finite(occ))) which.max(occ) else NA_integer_
    gbw <- if (any(is.finite(obb))) which.min(obb) else NA_integer_
    if (!is.na(g9)) put("orsfc9", te, OT[[g9]]$sc, OT[[g9]]$rk)
    if (!is.na(gcw)) put("orsfc", te, OT[[gcw]]$sc, OT[[gcw]]$rk)
    if (!is.na(gbw)) put("orsfb", te, OT[[gbw]]$sc, OT[[gbw]]$rk)
    tm[["orsf"]] <- tm[["orsf"]] + proc.time()[["elapsed"]] - t0
    lv <- function(v, x) if (is.na(x)) NA_real_ else v[x]
    mlev <- function(x) {
      if (is.na(x)) return(NA_real_)
      match(og$mtry[x], sort(unique(og$mtry)))
    }
    prm[[length(prm) + 1]] <- data.table(
      fold = k,
      rsfc_node = lv(grid$nodesize, gc_), rsfc_mtry = lv(grid$mtry, gc_),
      rsfb_node = lv(grid$nodesize, gb_), rsfb_mtry = lv(grid$mtry, gb_),
      orsfc9_leaf = lv(og$leaf, g9), orsfc_leaf = lv(og$leaf, gcw),
      orsfc_mlevel = mlev(gcw), orsfb_leaf = lv(og$leaf, gbw),
      orsfb_mlevel = mlev(gbw), rsf_differ = !identical(gc_, gb_),
      orsf_differ = !identical(gcw, gbw))
  }
  out <- list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
              cols = as.data.table(cols), prm = rbindlist(prm),
              nfail = nfail, tm = tm)
  if (!probe) {
    cat(sprintf("%s  done  %s %-7s %-4s seed %2d  %5.1f min  failed %d\n",
                format(Sys.time(), "%H:%M:%S"), tk$fr, tk$co, tk$mod,
                tk$seed, (proc.time()[["elapsed"]] - t_task) / 60, nfail),
        file = tk$prog, append = TRUE)
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, tk$ckpt)
  }
  out
}

# ---------------- build ----------------
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
    cc <- wl[modality == m, column]
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cc])
    storage.mode(Xm) <- "double"
    colnames(Xm) <- sprintf("f%d", seq_len(ncol(Xm)))
    for (s in seeds) {
      tasks[[length(tasks) + 1]] <- list(
        fr = fr, co = co, mod = m, seed = s, X = Xm, t = sub$time_death,
        e = as.integer(sub$event_death == 1),
        fold = fmap[[as.character(s)]], hz = HZ[[fr]], ntree = NTREE,
        nsplit = NSPLIT, ntime = NTIME, nodes = NODES, mtryf = MTRYF,
        leaf = LEAF, omtryf = OMTRYF, ontree = ONTREE, batch = BATCH,
        prog = PROG,
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", fr, code, m, s)))
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                seeds = seeds, mods = mods, fmap = fmap)
  cat(sprintf("  %s %-7s n=%d  modalities %s  seeds %d\n", fr, co,
              nrow(sub), paste(mods, collapse = "/"), length(seeds)))
}

# ---------------- timing probe, or resume ----------------
kind <- sapply(tasks, function(t) paste(t$fr, t$co, t$mod))
if (file.exists(PLANF)) {
  PL <- readRDS(PLANF)
  costs <- PL$costs
  cat(sprintf("\n=== RESUMING (plan of %s); probe skipped ===\n", PL$when))
} else {
  cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
  reps <- tasks[!duplicated(kind)]
  names(reps) <- unique(kind)
  plan(multisession, workers = min(workers, length(reps)))
  t0 <- Sys.time()
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  for (kd in names(pres)) {
    cat(sprintf("  %-16s RSF %6.1f s, oblique forest %6.1f s per fold\n",
                kd, pres[[kd]]$tm[["rsf"]], pres[[kd]]$tm[["orsf"]]))
  }
  costs <- sapply(pres, function(r) sum(r$tm))[kind] * NFOLD
  cat(sprintf("  probe took %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
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
  !is.null(r) && r$nfail == 0 &&
    all(sapply(r$cols, function(v) all(is.finite(v))))
}
done <- which(sapply(tasks, function(t) ok_ck(t$ckpt)))
todo <- setdiff(seq_along(tasks), done)
proj <- max(sum(costs[todo]) / workers, max(c(0, costs[todo]))) * 1.1
cat(sprintf(paste("  %d of %d tasks already saved; %d to fit, about %.1f",
                  "hours on %d workers\n"), length(done), length(tasks),
            length(todo), proj / 3600, workers))
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- full run ----------------
cat("\n=== FITTING (progress: proper_progress.txt) ===\n")
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
nf <- sum(sapply(res, function(r) r$nfail))
chk("forest fits that failed", nf, if (nf == 0) "PASS" else "FAIL")
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
  fwrite(out, file.path(OUTD, sprintf("oof_proper_%s_%s.csv", M$fr,
                                      M$code)))
  nmiss <- sum(!is.finite(as.matrix(out[, ..pc])))
  chk(sprintf("%s: rows / columns / missing cells", key),
      sprintf("%d / %d / %d", nrow(out), length(pc), nmiss),
      if (nrow(out) == nrow(M$sub) * length(M$seeds) && nmiss == 0)
        "PASS" else "FAIL")
  # reproduction of stage 6c (RSF, C) and stage 6d (oblique RSF, C)
  key2 <- paste(out$seed, out$hadm_id)
  rep_chk <- function(fn, pairs, what) {
    if (!file.exists(fn)) {
      chk(sprintf("%s: %s file for the reproduction check", key, what),
          "missing", "WARN")
      return(invisible())
    }
    R <- fread(fn)
    rk <- match(key2, paste(R$seed, R$hadm_id))
    dm <- 0
    nn <- 0L
    for (pp in pairs) {
      if (!(pp[2] %in% names(R))) next
      dm <- max(dm, max(abs(out[[pp[1]]] - R[[pp[2]]][rk])))
      nn <- nn + 1L
    }
    chk(sprintf("%s: reproduces %s (columns, max |diff|)", key, what),
        sprintf("%d, %.1e", nn, dm),
        if (nn > 0 && is.finite(dm) && dm <= 1e-9) "PASS" else "FAIL")
  }
  pr6c <- list()
  pr6d <- list()
  for (m in M$mods) {
    pr6c[[length(pr6c) + 1]] <- c(sprintf("%s.rsfc.death.score", m),
                                  sprintf("%s.rsf.death.score", m))
    for (h in HZ[[M$fr]]) {
      pr6c[[length(pr6c) + 1]] <- c(sprintf("%s.rsfc.death.r%d", m, h),
                                    sprintf("%s.rsf.death.r%d", m, h))
    }
    pr6d[[length(pr6d) + 1]] <- c(sprintf("%s.orsfc9.death.score", m),
                                  sprintf("%s.orsf.death.score", m))
  }
  rep_chk(file.path(OUTD, sprintf("oof_forest_wide_%s_%s.csv", M$fr,
                                  M$code)), pr6c, "stage 6c (RSF)")
  rep_chk(file.path(OUTD, sprintf("oof_extra_%s_%s.csv", M$fr, M$code)),
          pr6d, "stage 6d (oblique RSF)")
  td <- M$sub$time_death
  ed <- as.integer(M$sub$event_death == 1)
  nbad <- 0L
  for (cc in grep("\\.score$", pc, value = TRUE)) {
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
P <- rbindlist(PRM, fill = TRUE)
fwrite(P, file.path(OUTD, "params_proper.csv"))
for (m in c("ehr", "ecg", "ctpa")) {
  q <- P[modality == m]
  if (!nrow(q)) next
  chk(sprintf("%s: Brier and C choose different settings (RSF / oblique)",
              toupper(m)),
      sprintf("%.0f%% / %.0f%%", 100 * mean(q$rsf_differ),
              100 * mean(q$orsf_differ)), "INFO")
  e1 <- c(mean(q$rsfb_node == min(NODES)), mean(q$rsfb_node == max(NODES)))
  chk(sprintf("%s RSF (Brier): node size 5 / 160 share", toupper(m)),
      sprintf("%.2f / %.2f", e1[1], e1[2]),
      if (max(e1) > 0.5) "WARN" else "INFO")
  e2 <- c(mean(q$orsfc_mlevel == 4, na.rm = TRUE),
          mean(q$orsfb_mlevel == 4, na.rm = TRUE))
  chk(sprintf("%s oblique RSF: largest mtry share, C / Brier", toupper(m)),
      sprintf("%.2f / %.2f", e2[1], e2[2]),
      if (max(e2) > 0.5) "WARN" else "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_proper_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(OUTD, sprintf("oof_proper_%s_%s.csv", META[[key]]$fr,
                                    META[[key]]$code)), " (local only)\n")
}
cat("  ", file.path(OUTD, "params_proper.csv"), " (local only)\n")
cat("  ", file.path(RS, "level1_proper_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_proper_log.txt"), " (terminal log)\n")
cat("\nOutputs: level1_proper_checks.csv\n")