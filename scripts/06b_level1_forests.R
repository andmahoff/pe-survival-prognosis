# 06b_level1_forests.R
# Level-1 random forests, out of fold, for every frame, cohort, modality and
# seed (PREREGISTRATION.md section 4; DEVIATIONS.md D5).
#
# Two forests are fitted per modality, with the targets and labels of stage
# 6a. The survival forest (randomForestSRC, log-rank splitting) is fitted
# for the four ridge-Cox targets (death, composite, cs_cv, cs_death). It
# stores the ensemble mortality (a risk score; higher = worse) and, for
# death and composite, the risk at every horizon, 1 - S(h + 0.5), from a
# grid of 60 training event times plus, for each horizon, the last training
# event time at or before it. The classification forest (randomForestSRC,
# Gini splitting) is fitted for every label at every horizon (death_h,
# composite_h, cvu_h, cvc_h).
#
# Budgets are matched: 500 trees, nsplit = 10 random split points and the
# same 6-point grid (node size 5, 15 or 40 x mtry ceiling(sqrt(p)) or
# ceiling(3 sqrt(p))), tuned on out-of-bag Harrell's C (survival) or
# out-of-bag AUROC (classification) inside each outer training fold (D5).
# Median imputation is fitted inside each training fold; folds are the stage
# 4 folds, with six seeds. The ridge-Cox + forest rank average (section 4)
# is formed at evaluation. No performance numbers are written; the checks
# file holds completeness and pass/fail checks only (D2 per-fold rule: a
# score fails if its per-fold average C or AUROC is below 0.40).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/level1_forest_checks.csv, results/level1_forest_log.txt;
#          data/processed/level1/oof_forest_*.csv and params_forest.csv (not
#          in the repo)
# Usage:   Rscript scripts/06b_level1_forests.R [--probe-only]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. Tasks (frame, cohort, modality,
#          seed) run in parallel with future.apply, one core per forest.

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
PROG <- file.path(OUTD, "forest_progress.txt")
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5L, 15L, 40L)
MTRYF <- c(1, 3)
NFOLD <- 5L
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
PKG <- c("survival", "data.table", "randomForestSRC")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
cat(sprintf("randomForestSRC %s; %d workers\n",
            as.character(packageVersion("randomForestSRC")), workers))

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk, folds_to_run = 0:4, probe = FALSE) {
  options(rf.cores = 1L, mc.cores = 1L)
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
                seed = -as.integer(tk$seed * 10000 + k * 100 + g)),
          error = function(e) NULL)
        if (is.null(f)) {
          nfail <- nfail + 1L
          next
        }
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
  list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
       cols = as.data.table(cols),
       prm = if (length(prm)) rbindlist(prm) else NULL,
       nfail = nfail, t_surv = t_surv, t_class = t_class)
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
        prog = PROG, cost = 0)
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                ST = ST, BT = BT, fmap = fmap,
                                seeds = seeds, mods = mods)
  cat(sprintf("  %s %-7s n=%d  modalities %s  binary labels %d\n", fr, co,
              nrow(sub), paste(mods, collapse = "/"), length(BT)))
}

# ---------------- timing probe ----------------
cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
reps <- list()
for (tk in tasks) {
  key <- paste(tk$fr, tk$co, tk$mod)
  if (is.null(reps[[key]])) reps[[key]] <- tk
}
plan(multisession, workers = workers)
t0 <- Sys.time()
pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                      future.packages = PKG, future.seed = TRUE,
                      future.chunk.size = 1L)
dtp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
unit <- list()
for (key in names(pres)) {
  r <- pres[[key]]
  unit[[key]] <- c(ts = r$t_surv, tc = r$t_class)
  cat(sprintf("  %-16s %6.1f s per survival target, %5.1f s per label",
              key, r$t_surv, r$t_class),
      "(one fold, six grid points)\n")
}
for (i in seq_along(tasks)) {
  u <- unit[[paste(tasks[[i]]$fr, tasks[[i]]$co, tasks[[i]]$mod)]]
  tasks[[i]]$cost <- NFOLD * (length(tasks[[i]]$ST) * u[["ts"]] +
                                length(tasks[[i]]$BT) * u[["tc"]])
}
costs <- sapply(tasks, function(t) t$cost)
proj <- max(sum(costs) / workers, max(costs)) * 1.1
cat(sprintf("  probe took %.1f min; %d tasks; projected full run about",
            dtp / 60, length(tasks)),
    sprintf("%.1f hours on %d workers\n", proj / 3600, workers))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- full run ----------------
cat("\n=== FITTING (progress: forest_progress.txt) ===\n")
ord <- order(-costs)
cat(sprintf("%s  started %d tasks, projected %.1f hours\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(tasks),
            proj / 3600), file = PROG)
res <- future_lapply(tasks[ord], run_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1L)
plan(sequential)

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
  fwrite(out, file.path(OUTD, sprintf("oof_forest_%s_%s.csv",
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
pr_all <- rbindlist(PRM, fill = TRUE)
fwrite(pr_all, file.path(OUTD, "params_forest.csv"))
for (lr in c("rsf", "rf")) for (m in c("ehr", "ecg", "ctpa")) {
  q <- pr_all[learner == lr & modality == m]
  if (!nrow(q)) next
  lo <- mean(q$nodesize == min(NODES), na.rm = TRUE)
  hi <- mean(q$nodesize == max(NODES), na.rm = TRUE)
  mx <- mean(q$mtry == max(q$mtry, na.rm = TRUE), na.rm = TRUE)
  chk("-", "-", sprintf("%s %s: node size 5 / 40 share", lr, m),
      sprintf("%.2f / %.2f", lo, hi),
      if (max(lo, hi) > 0.5) "WARN" else "INFO")
  chk("-", "-", sprintf("%s %s: share at the larger mtry", lr, m),
      sprintf("%.2f", mx), "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_forest_checks.csv"))
cat(sprintf("\n  FAIL %d   WARN %d   elapsed %.1f hours\n",
            sum(ck$status == "FAIL"), sum(ck$status == "WARN"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(OUTD, sprintf("oof_forest_%s_%s.csv",
                                    META[[key]]$fr, META[[key]]$code)),
      " (local only)\n")
}
cat("  ", file.path(OUTD, "params_forest.csv"), " (local only)\n")
cat("  ", PROG, " (local only)\n")
cat("  ", file.path(RS, "level1_forest_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_forest_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: level1_forest_checks.csv\n")