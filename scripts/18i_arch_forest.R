# 18i_arch_forest.R
# Survival forests on the concatenated modalities, the forest arms of A10
# (PREREGISTRATION_ADDENDUM_2.md, A26).
#
# For every frame, cohort and modality set of A10 (EHR+ECG; in the
# three-modality cohort also EHR+CTPA and EHR+ECG+CTPA), and every seed and
# fold of 04_frames.py, two forests are fitted with the data handling and
# forest of 06c_level1_forests_wide.R: training-fold median imputation, no
# scaling, randomForestSRC with 500 trees, log-rank splitting, nsplit 10
# and save.memory. Early RSF is the registered random survival forest on
# the concatenation (PREREGISTRATION.md section 4), with node size 10 or 40
# and mtry ceiling(sqrt(p)) or ceiling(2 sqrt(p)) chosen on out-of-bag
# Brier score summed over the frame's horizons (the A22 rule). Block forest
# is a survival block forest (Hornung and Wright 2019, 'VarProb'): at every
# split each feature is offered with probability proportional to its
# modality's weight divided by the modality's number of features (xvar.wt),
# so a modality's chance of being offered no longer grows with its size. The
# EHR weight is 1 and the other weights are chosen on out-of-bag Brier from
# 0.25, 0.5 or 1 (two-modality sets) or (ECG, CTPA) in {0.25, 1} x {0.25, 1}
# (EHR+ECG+CTPA). Node size is 10 and mtry ceiling(sqrt(p)).
#
# Each forest's risk of death by h + 0.5 at every horizon is averaged over
# seeds and written in the layout of 18a (hadm_id, model, set, horizon,
# risk) for 18k_arch_eval2.R.
#
# Checks:  every patient predicted, risks in [0, 1] and never falling with
#          the horizon, no failed fits, D2 per-fold rule (all FAIL); where
#          the tuning choices fell (INFO).
# Notes:   each finished task is saved in fusion_arch/forest_ckpt/, and a
#          rerun resumes from them. A timing probe (fold 0 of one task per
#          frame, cohort and set) projects the run time, then the script
#          pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: data/processed/fusion_arch/arch_preds_forest_*.csv, with
#          params_arch_forest.csv and forest_ckpt/ (not in the repo);
#          results/arch_forest_checks.csv, results/arch_forest_log.txt
# Usage:   Rscript scripts/18i_arch_forest.R [--probe-only] [--workers=N]

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
FA <- file.path(PR, "fusion_arch")
RS <- file.path(P2, "results")
CKPT <- file.path(FA, "forest_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
SETS <- list("EHR+ECG" = list("EHR+ECG" = c("ehr", "ecg")),
             "3-mod" = list("EHR+ECG" = c("ehr", "ecg"),
                            "EHR+CTPA" = c("ehr", "ctpa"),
                            "EHR+ECG+CTPA" = c("ehr", "ecg", "ctpa")))
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(10L, 40L)
MTRYF <- c(1, 2)
BNODE <- 10L
W2 <- c(0.25, 0.5, 1)
W3 <- c(0.25, 1)
NFOLD <- 5L
PKG <- c("survival", "data.table", "randomForestSRC")
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

# ---------------- one task = frame, cohort, set, seed ----------------
run_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                     probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  n <- length(tk$t)
  H <- tk$hz
  tt <- tk$t
  ee <- tk$e
  p <- ncol(tk$X)
  colnames(tk$X) <- sprintf("f%d", seq_len(p))
  R <- list("Early RSF" = matrix(NA_real_, n, length(H)),
            "Block forest" = matrix(NA_real_, n, length(H)))
  prm <- list()
  nfail <- 0L
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- tk$X
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- med[na[, 2]]
    ev <- sort(unique(tt[tr][ee[tr] == 1]))
    at_h <- sapply(H, function(h) {
      v <- ev[ev <= h + 0.5]
      if (length(v)) max(v) else NA_real_
    })
    pick <- ev[unique(round(seq(1, length(ev),
                                length.out = min(tk$ntime, length(ev)))))]
    nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
    dtr <- data.frame(time = tt[tr], status = ee[tr], X[tr, , drop = FALSE])
    dte <- as.data.frame(X[te, , drop = FALSE])
    Ytr <- matrix(sapply(H, function(h) as.numeric(ee[tr] == 1 &
                                                      tt[tr] <= h + 0.5)),
                  ncol = length(H))
    grow <- function(node, mtry, wt, sd) {
      f <- tryCatch(rfsrc(Surv(time, status) ~ ., data = dtr,
                          ntree = tk$ntree, nodesize = node, mtry = mtry,
                          nsplit = tk$nsplit, ntime = nt,
                          splitrule = "logrank", xvar.wt = wt,
                          save.memory = TRUE, seed = -as.integer(sd)),
                    error = function(e) NULL)
      if (is.null(f)) return(NULL)
      ro <- sapply(H, function(h) {
        i <- findInterval(h + 0.5, f$time.interest)
        if (i == 0) rep(0, length(tr)) else 1 - f$survival.oob[, i]
      })
      br <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2,
                         na.rm = TRUE))
      pr <- predict(f, newdata = dte)
      rk <- sapply(H, function(h) {
        i <- findInterval(h + 0.5, pr$time.interest)
        if (i == 0) rep(0, length(te)) else 1 - pr$survival[, i]
      })
      list(brier = br, risk = matrix(rk, ncol = length(H)))
    }
    base <- tk$seed * 10000 + k * 100
    # early RSF
    best <- NULL
    g <- 0L
    for (node in tk$nodes) for (mf in tk$mtryf) {
      g <- g + 1L
      r <- grow(node, min(p, ceiling(mf * sqrt(p))), rep(1, p), base + g)
      if (is.null(r)) {
        nfail <- nfail + 1L
        next
      }
      if (is.finite(r$brier) && (is.null(best) || r$brier < best$brier)) {
        best <- c(r, list(node = node, mf = mf))
      }
    }
    if (!is.null(best)) R[["Early RSF"]][te, ] <- best$risk
    # block forest (VarProb); modalities in the set's order, EHR first
    lev <- unique(tk$blk)
    bs <- table(factor(tk$blk, levels = lev))
    W <- if (length(bs) == 2) {
      lapply(tk$w2, function(w) c(1, w))
    } else {
      g3 <- expand.grid(a = tk$w3, b = tk$w3)
      lapply(seq_len(nrow(g3)), function(i) c(1, g3$a[i], g3$b[i]))
    }
    bb <- NULL
    for (wi in seq_along(W)) {
      wv <- W[[wi]][match(tk$blk, lev)] / as.numeric(bs[tk$blk])
      r <- grow(tk$bnode, min(p, ceiling(sqrt(p))), wv, base + 50 + wi)
      if (is.null(r)) {
        nfail <- nfail + 1L
        next
      }
      if (is.finite(r$brier) && (is.null(bb) || r$brier < bb$brier)) {
        bb <- c(r, list(w = paste(W[[wi]][-1], collapse = "/")))
      }
    }
    if (!is.null(bb)) R[["Block forest"]][te, ] <- bb$risk
    prm[[length(prm) + 1]] <- data.table(
      fold = k, rsf_node = if (is.null(best)) NA else best$node,
      rsf_mtryf = if (is.null(best)) NA else best$mf,
      block_weights = if (is.null(bb)) NA else bb$w)
    invisible(gc(verbose = FALSE))
  }
  out <- list(fr = tk$fr, co = tk$co, set = tk$set, seed = tk$seed, R = R,
              prm = rbindlist(prm), nfail = nfail,
              sec = proc.time()[["elapsed"]] - t_task)
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, tk$ckpt)
  }
  out
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- lapply(MATS, function(f) unique(fread(file.path(IN, f)),
                                        by = "hadm_id"))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in names(SETS)) {
  Fd <- if (fr == "A") fa else fb
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  fmap <- lapply(seeds, function(s) {
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(f)) stop("fold missing: ", fr, " ", co, " seed ", s)
    as.integer(f)
  })
  names(fmap) <- seeds
  XM <- lapply(c("ehr", "ecg", "ctpa"), function(m) {
    cc <- intersect(wl[modality == m, column], names(mats[[m]]))
    X <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cc])
    storage.mode(X) <- "double"
    X
  })
  names(XM) <- c("ehr", "ecg", "ctpa")
  for (sn in names(SETS[[co]])) {
    ms <- SETS[[co]][[sn]]
    X <- do.call(cbind, XM[ms])
    blk <- rep(ms, sapply(XM[ms], ncol))
    for (s in seeds) {
      tasks[[length(tasks) + 1]] <- list(
        fr = fr, co = co, set = sn, seed = s, X = X, blk = blk,
        t = sub$time_death, e = as.integer(sub$event_death == 1),
        fold = fmap[[as.character(s)]], hz = HZ[[fr]], ntree = NTREE,
        nsplit = NSPLIT, ntime = NTIME, nodes = NODES, mtryf = MTRYF,
        bnode = BNODE, w2 = W2, w3 = W3,
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", fr, code,
                                       gsub("[^A-Za-z]+", "", sn), s)))
    }
  }
  META[[paste(fr, co)]] <- list(fr = fr, co = co, code = code, sub = sub,
                                seeds = seeds, fmap = fmap,
                                sets = names(SETS[[co]]))
  cat(sprintf("  %s %-7s n=%d  sets %s  seeds %d\n", fr, co, nrow(sub),
              paste(names(SETS[[co]]), collapse = ", "), length(seeds)))
}
kind <- sapply(tasks, function(t) paste(t$fr, t$co, t$set))
if (file.exists(PLANF)) {
  PL <- readRDS(PLANF)
  costs <- PL$costs
  cat(sprintf("\n=== RESUMING (plan of %s); probe skipped ===\n", PL$when))
} else {
  cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
  reps <- tasks[!duplicated(kind)]
  names(reps) <- unique(kind)
  plan(multisession, workers = min(workers, length(reps)))
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  for (kd in names(pres)) {
    cat(sprintf("  %-22s %6.1f s per fold\n", kd, pres[[kd]]$sec))
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
  !is.null(r) && r$nfail == 0 &&
    all(sapply(r$R, function(m) all(is.finite(m))))
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
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
nf <- sum(sapply(res, function(r) r$nfail))
chk("forest fits that failed", nf, if (nf == 0) "PASS" else "FAIL")
PRM <- list()
for (key in names(META)) {
  M <- META[[key]]
  hz <- HZ[[M$fr]]
  n <- nrow(M$sub)
  LONG <- list()
  nbad <- 0L
  mono <- 0L
  for (sn in M$sets) {
    rr <- Filter(function(r) r$fr == M$fr && r$co == M$co &&
                   r$set == sn, res)
    for (r in rr) {
      PRM[[length(PRM) + 1]] <- cbind(data.table(
        frame = M$fr, cohort = M$co, set = sn, seed = r$seed), r$prm)
    }
    for (am in c("Early RSF", "Block forest")) {
      P <- Reduce(`+`, lapply(rr, function(r) r$R[[am]])) / length(rr)
      mono <- mono + sum(apply(P, 1, function(v) any(diff(v) < -1e-9)))
      for (r in rr) {
        f <- M$fmap[[as.character(r$seed)]]
        v <- sapply(sort(unique(f)), function(k) {
          ii <- which(f == k)
          cidx(M$sub$time_death[ii], M$sub$event_death[ii] == 1,
               r$R[[am]][ii, length(hz)])
        })
        if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
      }
      LONG[[length(LONG) + 1]] <- data.table(
        hadm_id = rep(M$sub$hadm_id, length(hz)), model = am, set = sn,
        horizon = rep(hz, each = n), risk = as.vector(P))
    }
  }
  L <- rbindlist(LONG)
  fwrite(L, file.path(FA, sprintf("arch_preds_forest_%s_%s.csv", M$fr,
                                  M$code)))
  rg <- range(L$risk)
  chk(sprintf("%s: missing / risk range", key),
      sprintf("%d / %.4f to %.4f", sum(!is.finite(L$risk)), rg[1], rg[2]),
      if (all(is.finite(L$risk)) && rg[1] >= 0 && rg[2] <= 1) "PASS" else
        "FAIL")
  chk(sprintf("%s: patients whose risk falls with the horizon", key),
      mono, if (mono == 0) "PASS" else "FAIL")
  chk(sprintf("%s: seed-level scores failing the per-fold rule (D2)",
              key), nbad, if (nbad == 0) "PASS" else "FAIL")
}
P <- rbindlist(PRM)
fwrite(P, file.path(FA, "params_arch_forest.csv"))
chk("Early RSF: node size 10 / 40, mtry factor 1 / 2 (shares)",
    sprintf("%.2f / %.2f, %.2f / %.2f", mean(P$rsf_node == 10),
            mean(P$rsf_node == 40), mean(P$rsf_mtryf == 1),
            mean(P$rsf_mtryf == 2)), "INFO")
tb <- P[, .N, by = .(set, block_weights)][order(set, -N)]
for (i in seq_len(nrow(tb))) {
  chk(sprintf("Block forest %s: weights %s (other modalities)",
              tb$set[i], tb$block_weights[i]), tb$N[i], "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "arch_forest_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(FA, sprintf("arch_preds_forest_%s_%s.csv",
                                  META[[key]]$fr, META[[key]]$code)),
      " (local only)\n")
}
cat("  ", file.path(RS, "arch_forest_checks.csv"), "\n")
cat("  ", file.path(RS, "arch_forest_log.txt"), " (terminal log)\n")
cat("\nOutputs: arch_forest_checks.csv\n")