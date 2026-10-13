# 06c_level1_forests_wide.R
# Level-1 random forests on the wide grid, out of fold, for every frame,
# cohort, modality and seed (DEVIATIONS.md D11, D13).
#
# Sensitivity analysis S4. The registered stage 6b forests chose the most
# smoothed corner of their grid in most folds, so the grid is extended at
# both ends: node size 5, 15, 40, 80 or 160 x mtry ceiling(sqrt(p) / 2),
# ceiling(sqrt(p)) or ceiling(3 sqrt(p)), 15 points against 6. Everything
# else is as in stage 6b (D5): the same survival and classification forests,
# targets and labels, 500 trees, nsplit = 10, out-of-bag tuning inside each
# outer training fold, median imputation, the stage 4 folds and six seeds.
# Version 2: forests are grown with save.memory = TRUE (identical
# predictions, about a quarter of the peak memory per worker:
# 1.08 GB to 0.26 GB for a 2,800 x 264 node-size-5 forest), and memory is
# released after every forest.
#
# Scope (D13): the grid is never trimmed. If the full run is projected above
# LIMIT_H (10) hours, the scope is reduced to the largest of: (1) all four
# survival targets and every label; (2) all four survival targets and the
# death labels; (3) the death target and the death labels; (4) the death
# target only. No performance numbers are written; the checks file holds
# completeness and pass/fail checks (D2 per-fold rule), the scope used, and
# where the tuning choices fell.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/level1_forest_wide_checks.csv,
#          results/level1_forest_wide_log.txt;
#          data/processed/level1/oof_forest_wide_*.csv and
#          params_forest_wide.csv (not in the repo)
# Usage:   Rscript scripts/06c_level1_forests_wide.R [--probe-only]
#          [--workers=N] [--scope=N]
# Notes:   --probe-only prints the projection for every scope and stops;
#          otherwise the run starts after a 5 s pause. --workers=N lowers
#          the worker count from its default of physical cores - 1.
#          --scope=N (1 to 4) forces a scope. Finished tasks are saved in
#          level1/forest_wide_ckpt/ (plan.rds holds the scope and worker
#          count); re-running the same command skips the probe and fits only
#          the missing tasks. Delete the folder to start afresh.

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OUTD <- file.path(PR, "level1")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
PROG <- file.path(OUTD, "forest_wide_progress.txt")
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5L, 15L, 40L, 80L, 160L)
MTRYF <- c(0.5, 1, 3)
LIMIT_H <- 10
sc_arg <- grep("^--scope=", args, value = TRUE)
FORCE <- if (length(sc_arg)) as.integer(sub("^--scope=", "", sc_arg)) else NA
wk_arg <- grep("^--workers=", args, value = TRUE)
WARG <- if (length(wk_arg)) as.integer(sub("^--workers=", "", wk_arg)) else NA
CKPT <- file.path(OUTD, "forest_wide_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
NFOLD <- 5L
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
PKG <- c("survival", "data.table", "randomForestSRC")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
cat(sprintf("randomForestSRC %s; up to %d workers\n",
            as.character(packageVersion("randomForestSRC")), workers))

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk, folds_to_run = 0:4, probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  cidx <- function(t, e, s) {
    out <- tryCatch(survival::concordance(
      survival::Surv(t, e) ~ s, reverse = TRUE)$concordance,
      error = function(err) NA_real_)
    as.numeric(out)
  }
  auc <- function(y, s) {
    y <- as.logical(y)
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  X <- tk$X
  n <- nrow(X)
  p <- ncol(X)
  grid <- expand.grid(
    nodesize = tk$nodes,
    mtry = unique(pmin(p, ceiling(tk$mtryf * sqrt(p)))))
  sv <- if (probe) "death" else names(tk$ST)
  lb <- if (probe) sprintf("death_h%d", max(tk$hs)) else names(tk$BT)
  cols <- list()
  for (tg in sv) {
    cols[[sprintf("%s.rsf.%s.score", tk$mod, tg)]] <- rep(NA_real_, n)
    if (tg %in% c("death", "composite")) {
      for (h in tk$hs) {
        cols[[sprintf("%s.rsf.%s.r%d", tk$mod, tg, h)]] <-
          rep(NA_real_, n)
      }
    }
  }
  for (b in lb) cols[[sprintf("%s.rf.%s", tk$mod, b)]] <- rep(NA_real_, n)
  prm <- list()
  nfail <- 0L
  t_surv <- 0
  t_class <- 0
  peak <- 0
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    Xi <- X
    na <- which(is.na(Xi), arr.ind = TRUE)
    if (nrow(na) > 0) Xi[na] <- med[na[, 2]]
    Xte <- as.data.frame(Xi[te, , drop = FALSE])
    # ---- survival forests ----
    for (tg in sv) {
      t0 <- proc.time()[["elapsed"]]
      tt <- tk$ST[[tg]]$t
      ee <- as.integer(tk$ST[[tg]]$e)
      ev <- sort(unique(tt[tr][ee[tr] == 1]))
      at_h <- sapply(tk$hs, function(h) {
        v <- ev[ev <= h + 0.5]
        if (length(v)) max(v) else NA_real_
      })
      pick <- ev[unique(round(seq(1, length(ev),
                                  length.out = min(tk$ntime,
                                                   length(ev)))))]
      nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
      dtr <- data.frame(time = tt[tr], status = ee[tr],
                        Xi[tr, , drop = FALSE])
      best <- -Inf
      bi <- NA_integer_
      sc <- rep(NA_real_, nrow(grid))
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
        if (probe) peak <- max(peak, as.numeric(object.size(f)) / 2^20)
        sc[g] <- cidx(tt[tr], ee[tr], f$predicted.oob)
        if (is.finite(sc[g]) && sc[g] > best) {
          pr <- tryCatch(predict(f, newdata = Xte),
                         error = function(e) NULL)
          if (is.null(pr)) {
            nfail <- nfail + 1L
          } else {
            best <- sc[g]
            bi <- g
            cols[[sprintf("%s.rsf.%s.score", tk$mod, tg)]][te] <-
              pr$predicted
            if (tg %in% c("death", "composite")) {
              for (h in tk$hs) {
                i <- findInterval(h + 0.5, pr$time.interest)
                cols[[sprintf("%s.rsf.%s.r%d", tk$mod, tg, h)]][te] <-
                  if (i == 0) 0 else 1 - pr$survival[, i]
              }
            }
          }
        }
        rm(f)
        invisible(gc(verbose = FALSE))
      }
      prm[[length(prm) + 1]] <- data.frame(
        fold = k, learner = "rsf", target = tg,
        nodesize = if (is.na(bi)) NA_integer_ else grid$nodesize[bi],
        mtry = if (is.na(bi)) NA_integer_ else grid$mtry[bi],
        oob = paste(sprintf("%.4f", sc), collapse = ";"))
      t_surv <- t_surv + proc.time()[["elapsed"]] - t0
    }
    # ---- classification forests ----
    for (b in lb) {
      t0 <- proc.time()[["elapsed"]]
      y <- tk$BT[[b]]$y
      keep <- tk$BT[[b]]$keep
      trl <- tr[keep[tr]]
      dtr <- data.frame(
        y = factor(as.integer(y[trl]), levels = c(0, 1)),
        Xi[trl, , drop = FALSE])
      best <- -Inf
      bi <- NA_integer_
      sc <- rep(NA_real_, nrow(grid))
      for (g in seq_len(nrow(grid))) {
        f <- tryCatch(
          rfsrc(y ~ ., data = dtr, ntree = tk$ntree,
                nodesize = grid$nodesize[g], mtry = grid$mtry[g],
                nsplit = tk$nsplit, splitrule = "gini",
                save.memory = TRUE,
                seed = -as.integer(tk$seed * 10000 + k * 100 + g)),
          error = function(e) NULL)
        if (is.null(f)) {
          nfail <- nfail + 1L
          next
        }
        sc[g] <- auc(y[trl], f$predicted.oob[, "1"])
        if (is.finite(sc[g]) && sc[g] > best) {
          pr <- tryCatch(predict(f, newdata = Xte),
                         error = function(e) NULL)
          if (is.null(pr)) {
            nfail <- nfail + 1L
          } else {
            best <- sc[g]
            bi <- g
            cols[[sprintf("%s.rf.%s", tk$mod, b)]][te] <-
              pr$predicted[, "1"]
          }
        }
        rm(f)
        invisible(gc(verbose = FALSE))
      }
      prm[[length(prm) + 1]] <- data.frame(
        fold = k, learner = "rf", target = b,
        nodesize = if (is.na(bi)) NA_integer_ else grid$nodesize[bi],
        mtry = if (is.na(bi)) NA_integer_ else grid$mtry[bi],
        oob = paste(sprintf("%.4f", sc), collapse = ";"))
      t_class <- t_class + proc.time()[["elapsed"]] - t0
    }
    invisible(gc())
  }
  mins <- (proc.time()[["elapsed"]] - t_task) / 60
  if (!probe) {
    cat(sprintf("%s  done  %s %-7s %-4s seed %2d  %5.1f min  failed %d\n",
                format(Sys.time(), "%H:%M:%S"), tk$fr, tk$co, tk$mod,
                tk$seed, mins, nfail),
        file = tk$prog, append = TRUE)
  }
  out <- list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
              cols = as.data.table(cols),
              prm = if (length(prm)) rbindlist(prm) else NULL,
              nfail = nfail, t_surv = t_surv, t_class = t_class,
              peak = peak)
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, tk$ckpt)
  }
  out
}

# ---------------- helpers for the checks ----------------
cidx_c <- function(t, e, s) {
  out <- tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$concordance,
                  error = function(err) NA_real_)
  as.numeric(out)
}
auc_c <- function(y, s) {
  y <- as.logical(y)
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
surv_targets <- function(s) {
  list(death = list(t = s$time_death, e = s$event_death == 1),
       composite = list(t = s$time_first, e = s$event_comp == 1),
       cs_cv = list(t = s$time_first, e = s$event_cr == 1),
       cs_death = list(t = s$time_first, e = s$event_cr == 2))
}
bin_targets <- function(s, hs) {
  out <- list()
  allr <- rep(TRUE, nrow(s))
  for (h in hs) {
    cut <- h + 0.5
    cvb <- s$event_cr == 1 & s$time_first <= cut
    dfb <- s$event_cr == 2 & s$time_first <= cut
    out[[sprintf("death_h%d", h)]] <-
      list(y = s$event_death == 1 & s$time_death <= cut, keep = allr)
    out[[sprintf("composite_h%d", h)]] <-
      list(y = s$event_comp == 1 & s$time_first <= cut, keep = allr)
    out[[sprintf("cvu_h%d", h)]] <- list(y = cvb, keep = allr)
    out[[sprintf("cvc_h%d", h)]] <- list(y = cvb, keep = !dfb)
  }
  out
}

# ---------------- build tasks ----------------
cat("=== BUILD ===\n")
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
  ST <- surv_targets(sub)
  BT <- bin_targets(sub, HZ[[fr]])
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
        fr = fr, co = co, mod = m, seed = s, X = Xm, ST = ST, BT = BT,
        fold = fmap[[as.character(s)]], hs = HZ[[fr]], ntree = NTREE,
        nsplit = NSPLIT, ntime = NTIME, nodes = NODES, mtryf = MTRYF,
        prog = PROG, cost = 0,
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", fr, code, m, s)))
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                ST = ST, BT = BT, fmap = fmap,
                                seeds = seeds, mods = mods)
  cat(sprintf("  %s %-7s n=%d  modalities %s  binary labels %d\n", fr, co,
              nrow(sub), paste(mods, collapse = "/"), length(BT)))
}

# ---------------- timing probe, workers and scope ----------------
# scopes, largest first (D13): survival targets and binary labels kept
scope_keep <- function(sc, tk) {
  st <- names(tk$ST)
  lb <- names(tk$BT)
  dl <- grep("^death_h", lb, value = TRUE)
  switch(sc,
         list(st = st, lb = lb),
         list(st = st, lb = dl),
         list(st = "death", lb = dl),
         list(st = "death", lb = character(0)))
}
SCOPE_TXT <- c("all four survival targets and every binary label",
               "all four survival targets and the death labels",
               "the death survival target and the death labels",
               "the death survival target only")
RESUME <- file.exists(PLANF)
if (RESUME) {
  PL <- readRDS(PLANF)
  SCOPE <- PL$scope
  WORKERS <- if (is.na(WARG)) PL$workers else max(1L, min(workers, WARG))
  if (!is.na(FORCE) && FORCE != SCOPE) {
    stop("checkpoints were made with scope ", SCOPE, "; delete ", CKPT,
         " to start afresh with another scope")
  }
  cat(sprintf(paste("\n=== RESUMING: scope %d (%s), %d workers; probe",
                    "skipped ===\n"), SCOPE, SCOPE_TXT[SCOPE], WORKERS))
  for (i in seq_along(tasks)) tasks[[i]]$cost <- PL$costs[[i]]
} else {
  cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
  reps <- list()
  for (tk in tasks) {
    key <- paste(tk$fr, tk$co, tk$mod)
    if (is.null(reps[[key]])) reps[[key]] <- tk
  }
  plan(multisession, workers = if (is.na(WARG)) workers else
    max(1L, min(workers, WARG)))
  t0 <- Sys.time()
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  dtp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  plan(sequential)
  invisible(gc())
  unit <- list()
  for (key in names(pres)) {
    r <- pres[[key]]
    unit[[key]] <- c(ts = r$t_surv, tc = r$t_class)
    cat(sprintf(paste("  %-16s %6.1f s per survival target, %5.1f s per",
                      "label, largest forest %4.0f MB\n"),
                key, r$t_surv, r$t_class, r$peak))
  }
  WORKERS <- if (is.na(WARG)) workers else max(1L, min(workers, WARG))
  cat(sprintf("  WORKERS %d%s\n", WORKERS,
              if (is.na(WARG)) "" else " (set by --workers)"))
  proj_of <- function(sc) {
    cs <- sapply(tasks, function(t) {
      u <- unit[[paste(t$fr, t$co, t$mod)]]
      kp <- scope_keep(sc, t)
      NFOLD * (length(kp$st) * u[["ts"]] + length(kp$lb) * u[["tc"]])
    })
    list(costs = cs,
         hours = max(sum(cs) / WORKERS, max(cs)) * 1.1 / 3600)
  }
  PJ <- lapply(1:4, proj_of)
  for (sc in 1:4) {
    cat(sprintf("  scope %d (%s): about %.1f hours on %d workers\n", sc,
                SCOPE_TXT[sc], PJ[[sc]]$hours, WORKERS))
  }
  if (is.na(FORCE)) {
    ok <- which(sapply(PJ, function(z) z$hours) <= LIMIT_H)
    SCOPE <- if (length(ok)) min(ok) else 4L
  } else {
    SCOPE <- FORCE
  }
  cat(sprintf("  probe took %.1f min; SCOPE %d chosen (%s)%s\n", dtp / 60,
              SCOPE, SCOPE_TXT[SCOPE],
              if (is.na(FORCE)) sprintf(", limit %g hours", LIMIT_H)
              else ", forced by --scope"))
  for (i in seq_along(tasks)) tasks[[i]]$cost <- PJ[[SCOPE]]$costs[i]
  if (!probe_only) {
    saveRDS(list(scope = SCOPE, workers = WORKERS,
                 costs = sapply(tasks, function(t) t$cost)), PLANF)
  }
}
for (i in seq_along(tasks)) {
  kp <- scope_keep(SCOPE, tasks[[i]])
  tasks[[i]]$ST <- tasks[[i]]$ST[kp$st]
  tasks[[i]]$BT <- tasks[[i]]$BT[kp$lb]
}
costs <- sapply(tasks, function(t) t$cost)
if (probe_only) {
  cat("--probe-only: stopping here (nothing saved).\n")
  quit(status = 0)
}
done <- which(file.exists(sapply(tasks, function(t) t$ckpt)))
todo <- setdiff(seq_along(tasks), done)
proj <- max(sum(costs[todo]) / WORKERS, max(c(0, costs[todo]))) * 1.1
cat(sprintf(paste("  %d of %d tasks already saved; %d to fit, about %.1f",
                  "hours on %d workers\n"), length(done), length(tasks),
            length(todo), proj / 3600, WORKERS))
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- full run ----------------
cat("\n=== FITTING (progress: forest_wide_progress.txt) ===\n")
cat(sprintf("%s  %s %d tasks, scope %d, %d workers, about %.1f hours\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            if (length(done)) "resumed with" else "started", length(todo),
            SCOPE, WORKERS, proj / 3600), file = PROG,
    append = length(done) > 0)
res_new <- list()
if (length(todo)) {
  ord <- todo[order(-costs[todo])]
  plan(multisession, workers = WORKERS)
  res_new <- future_lapply(tasks[ord], run_task, future.packages = PKG,
                           future.seed = TRUE, future.chunk.size = 1L)
  plan(sequential)
}
res <- c(lapply(tasks[done], function(t) readRDS(t$ckpt)), res_new)

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE AND CHECK ===\n")
CHK <- list()
chk <- function(fr, co, item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(frame = fr, cohort = co,
                                        item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-4s %s %-7s %-46s %s\n", status, fr, co, item, value))
}
PRM <- list()
for (r in res) {
  if (!is.null(r$prm)) {
    PRM[[length(PRM) + 1]] <- cbind(
      data.table(frame = r$fr, cohort = r$co, modality = r$mod,
                 seed = r$seed), r$prm)
  }
}
nfail <- sum(sapply(res, function(r) r$nfail))
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
  }
  out <- rbindlist(parts)
  pc <- sort(setdiff(names(out), c("seed", "subject_id", "hadm_id",
                                   "fold")))
  out <- out[, c("seed", "subject_id", "hadm_id", "fold", pc),
             with = FALSE]
  fwrite(out, file.path(OUTD, sprintf("oof_forest_wide_%s_%s.csv",
                                      M$fr, M$code)))
  chk(M$fr, M$co, "rows (admissions x seeds)",
      sprintf("%d = %d x %d", nrow(out), nrow(M$sub), length(M$seeds)),
      if (nrow(out) == nrow(M$sub) * length(M$seeds)) "PASS" else "FAIL")
  nmiss <- sum(is.na(as.matrix(out[, ..pc])))
  chk(M$fr, M$co, "prediction columns / missing cells",
      sprintf("%d / %d", length(pc), nmiss),
      if (nmiss == 0) "PASS" else "FAIL")
  sc_cols <- pc[!grepl("\\.r[0-9]+$", pc)]
  nbad <- 0L
  for (cc in sc_cols) {
    pz <- strsplit(cc, ".", fixed = TRUE)[[1]]
    vals <- c()
    for (s in M$seeds) {
      os <- out[out$seed == s]
      v <- os[[cc]]
      for (k in 0:(NFOLD - 1)) {
        ii <- which(os$fold == k)
        if (pz[2] == "rsf") {
          vals <- c(vals, cidx_c(M$ST[[pz[3]]]$t[ii],
                                 M$ST[[pz[3]]]$e[ii], v[ii]))
        } else {
          vals <- c(vals, auc_c(M$BT[[pz[3]]]$y[ii], v[ii]))
        }
      }
    }
    if (!is.finite(mean(vals, na.rm = TRUE)) ||
        mean(vals, na.rm = TRUE) < 0.40) {
      nbad <- nbad + 1L
    }
  }
  chk(M$fr, M$co, "scores failing the per-fold rule (D2)", nbad,
      if (nbad == 0) "PASS" else "FAIL")
}
chk("-", "-", "fits that failed", nfail,
    if (nfail == 0) "PASS" else "FAIL")
chk("-", "-", "scope (D13)", sprintf("%d: %s", SCOPE, SCOPE_TXT[SCOPE]),
    "INFO")
chk("-", "-", "workers / tasks resumed from checkpoints",
    sprintf("%d / %d", WORKERS, length(done)), "INFO")
pr_all <- rbindlist(PRM, fill = TRUE)
fwrite(pr_all, file.path(OUTD, "params_forest_wide.csv"))
for (lr in c("rsf", "rf")) for (m in c("ehr", "ecg", "ctpa")) {
  q <- pr_all[learner == lr & modality == m]
  if (!nrow(q)) next
  lo <- mean(q$nodesize == min(NODES), na.rm = TRUE)
  hi <- mean(q$nodesize == max(NODES), na.rm = TRUE)
  ml <- q[, .(low = mtry == min(mtry, na.rm = TRUE),
              high = mtry == max(mtry, na.rm = TRUE)),
          by = .(frame, cohort)]
  chk("-", "-", sprintf("%s %s: node size 5 / 160 share", lr, m),
      sprintf("%.2f / %.2f", lo, hi),
      if (max(lo, hi) > 0.5) "WARN" else "INFO")
  chk("-", "-", sprintf("%s %s: smallest / largest mtry share", lr, m),
      sprintf("%.2f / %.2f", mean(ml$low, na.rm = TRUE),
              mean(ml$high, na.rm = TRUE)),
      if (max(mean(ml$low, na.rm = TRUE),
              mean(ml$high, na.rm = TRUE)) > 0.5) "WARN" else "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_forest_wide_checks.csv"))
cat(sprintf("\n  FAIL %d   WARN %d   elapsed %.1f hours\n",
            sum(ck$status == "FAIL"), sum(ck$status == "WARN"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(OUTD, sprintf("oof_forest_wide_%s_%s.csv",
                                    META[[key]]$fr, META[[key]]$code)),
      " (local only)\n")
}
cat("  ", file.path(OUTD, "params_forest_wide.csv"), " (local only)\n")
cat("  ", PROG, " (local only)\n")
cat("  ", file.path(RS, "level1_forest_wide_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_forest_wide_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: level1_forest_wide_checks.csv\n")