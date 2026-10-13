# 18l_block_forest_opt.R
# Survival block forest of 18i tuned on a wider grid
# (PREREGISTRATION_ADDENDUM_2.md A27).
#
# In 18i the block forest was the only alternative architecture to match
# late fusion, but its weights sat at the edge of a small grid (the lowest
# weight for the other modalities in 53 to 63% of folds) and node size and
# mtry were fixed. The model is as in 18i: randomForestSRC on the
# concatenated modalities, training-fold median imputation, 500 trees,
# log-rank splitting, nsplit 10, save.memory. Each feature is offered at
# every split with probability proportional to its modality's weight
# divided by the modality's number of features (xvar.wt); the EHR weight
# is 1.
#
# Tuning has two steps on every training fold, both on out-of-bag Brier
# score summed over the frame's horizons (the A22 rule). Step 1 sets the
# weights at node size 10 and mtry ceiling(sqrt(p)): the other modality's
# weight in {0.05, 0.1, 0.25, 0.5, 1} for two-modality sets, and (ECG,
# CTPA) in {0.05, 0.1, 0.25} x {0.1, 0.25, 0.5} for EHR+ECG+CTPA. Step 2
# keeps those weights and searches node size {5, 10, 20} x mtry
# ceiling(f sqrt(p)) with f in {1, 2, 3}. The best forest over both steps
# is kept, for every frame, cohort and modality set of A10, every stage 4
# seed and every fold.
#
# Out-of-fold risk of death by h + 0.5 at every horizon is saved per seed
# (blockopt_oof_*.csv, used by 18m) and averaged over seeds in 18a's layout
# as "Block forest, optimised" (arch_preds_blockopt_*.csv, used by 18k).
# Checks: failed fits, missing predictions, risks outside [0, 1] or falling
# with the horizon, and the D2 per-fold rule (all FAIL); share of choices
# at each edge of the grid (WARN above 50%). Finished tasks are saved in
# fusion_arch/blockopt_ckpt/ and skipped on a rerun.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/blockopt_checks.csv, results/blockopt_params_summary.csv,
#          results/blockopt_log.txt;
#          data/processed/fusion_arch/blockopt_oof_*.csv and
#          arch_preds_blockopt_*.csv (patient level, not in the repo)
# Usage:   Rscript scripts/18l_block_forest_opt.R [--probe-only]
#          [--workers=N]   (--probe-only stops after the timing probe)

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
CKPT <- file.path(FA, "blockopt_ckpt")
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
W2 <- c(0.05, 0.1, 0.25, 0.5, 1)
W3_ECG <- c(0.05, 0.1, 0.25)
W3_CTPA <- c(0.1, 0.25, 0.5)
NODE1 <- 10L
NODES <- c(5L, 10L, 20L)
MTRYF <- c(1, 2, 3)
NFOLD <- 5L
MODEL <- "Block forest, optimised"
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
  R <- matrix(NA_real_, n, length(H))
  prm <- list()
  nfail <- 0L
  lev <- unique(tk$blk)
  bs <- table(factor(tk$blk, levels = lev))
  W <- if (length(bs) == 2) {
    lapply(tk$w2, function(w) c(1, w))
  } else {
    g3 <- expand.grid(a = tk$w3e, b = tk$w3c)
    lapply(seq_len(nrow(g3)), function(i) c(1, g3$a[i], g3$b[i]))
  }
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
      rm(f, pr)
      invisible(gc(verbose = FALSE))
      list(brier = br, risk = matrix(rk, ncol = length(H)))
    }
    base <- tk$seed * 10000 + k * 100
    wvec <- function(w) w[match(tk$blk, lev)] / as.numeric(bs[tk$blk])
    # step 1: weights at node size 10, mtry ceiling(sqrt(p))
    best <- NULL
    for (wi in seq_along(W)) {
      r <- grow(tk$node1, min(p, ceiling(sqrt(p))), wvec(W[[wi]]),
                base + wi)
      if (is.null(r)) {
        nfail <- nfail + 1L
        next
      }
      if (is.finite(r$brier) && (is.null(best) || r$brier < best$brier)) {
        best <- c(r, list(wi = wi, node = tk$node1, mf = 1))
      }
    }
    # step 2: node size and mtry at step 1's weights
    if (!is.null(best)) {
      g2 <- expand.grid(node = tk$nodes, mf = tk$mtryf)
      g2 <- g2[!(g2$node == tk$node1 & g2$mf == 1), ]
      for (gi in seq_len(nrow(g2))) {
        r <- grow(g2$node[gi], min(p, ceiling(g2$mf[gi] * sqrt(p))),
                  wvec(W[[best$wi]]), base + 50 + gi)
        if (is.null(r)) {
          nfail <- nfail + 1L
          next
        }
        if (is.finite(r$brier) && r$brier < best$brier) {
          best <- c(r, list(wi = best$wi, node = g2$node[gi],
                            mf = g2$mf[gi]))
        }
      }
      R[te, ] <- best$risk
    }
    prm[[length(prm) + 1]] <- data.table(
      fold = k,
      weights = if (is.null(best)) NA_character_ else
        paste(W[[best$wi]][-1], collapse = "/"),
      node = if (is.null(best)) NA_real_ else best$node,
      mtryf = if (is.null(best)) NA_real_ else best$mf)
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
        nsplit = NSPLIT, ntime = NTIME, node1 = NODE1, nodes = NODES,
        mtryf = MTRYF, w2 = W2, w3e = W3_ECG, w3c = W3_CTPA,
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
  !is.null(r) && r$nfail == 0 && all(is.finite(r$R))
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
  WIDE <- list()
  nbad <- 0L
  mono <- 0L
  for (sn in M$sets) {
    rr <- Filter(function(r) r$fr == M$fr && r$co == M$co &&
                   r$set == sn, res)
    for (r in rr) {
      PRM[[length(PRM) + 1]] <- cbind(data.table(
        frame = M$fr, cohort = M$co, set = sn, seed = r$seed), r$prm)
      f <- M$fmap[[as.character(r$seed)]]
      w <- data.table(seed = r$seed, hadm_id = M$sub$hadm_id, fold = f,
                      set = sn)
      for (hi in seq_along(hz)) w[[sprintf("r%d", hz[hi])]] <- r$R[, hi]
      WIDE[[length(WIDE) + 1]] <- w
      v <- sapply(sort(unique(f)), function(k) {
        ii <- which(f == k)
        cidx(M$sub$time_death[ii], M$sub$event_death[ii] == 1,
             r$R[ii, length(hz)])
      })
      if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
    }
    P <- Reduce(`+`, lapply(rr, function(r) r$R)) / length(rr)
    mono <- mono + sum(apply(P, 1, function(v) any(diff(v) < -1e-9)))
    LONG[[length(LONG) + 1]] <- data.table(
      hadm_id = rep(M$sub$hadm_id, length(hz)), model = MODEL, set = sn,
      horizon = rep(hz, each = n), risk = as.vector(P))
  }
  L <- rbindlist(LONG)
  fwrite(L, file.path(FA, sprintf("arch_preds_blockopt_%s_%s.csv", M$fr,
                                  M$code)))
  fwrite(rbindlist(WIDE), file.path(FA, sprintf("blockopt_oof_%s_%s.csv",
                                                M$fr, M$code)))
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
fwrite(P, file.path(FA, "params_blockopt.csv"))
# where the choices fell
two <- P[set != "EHR+ECG+CTPA"]
thr <- P[set == "EHR+ECG+CTPA"]
thr[, `:=`(w_ecg = as.numeric(sub("/.*$", "", weights)),
           w_ctpa = as.numeric(sub("^.*/", "", weights)))]
edges <- list(
  list("two-modality sets: weight at 0.05 (lowest)",
       mean(as.numeric(two$weights) == min(W2))),
  list("two-modality sets: weight at 1 (highest)",
       mean(as.numeric(two$weights) == max(W2))),
  list("EHR+ECG+CTPA: ECG weight at 0.05 (lowest)",
       mean(thr$w_ecg == min(W3_ECG))),
  list("EHR+ECG+CTPA: ECG weight at 0.25 (highest)",
       mean(thr$w_ecg == max(W3_ECG))),
  list("EHR+ECG+CTPA: CTPA weight at 0.1 (lowest)",
       mean(thr$w_ctpa == min(W3_CTPA))),
  list("EHR+ECG+CTPA: CTPA weight at 0.5 (highest)",
       mean(thr$w_ctpa == max(W3_CTPA))),
  list("node size at 5 (lowest)", mean(P$node == min(NODES))),
  list("node size at 20 (highest)", mean(P$node == max(NODES))),
  list("mtry factor at 1 (lowest)", mean(P$mtryf == min(MTRYF))),
  list("mtry factor at 3 (highest)", mean(P$mtryf == max(MTRYF))))
for (e in edges) {
  if (!is.finite(e[[2]])) next
  chk(sprintf("share of folds with %s", e[[1]]),
      sprintf("%.0f%%", 100 * e[[2]]), if (e[[2]] <= 0.5) "PASS" else
        "WARN")
}
sm <- rbind(
  P[, .(choice = "weights (other modalities)", value = weights)],
  P[, .(choice = "node size", value = as.character(node))],
  P[, .(choice = "mtry factor", value = as.character(mtryf))])
sm <- cbind(P[, .(frame, cohort, set)], sm)
sm <- sm[, .N, by = .(set, choice, value)][order(set, choice, -N)]
sm[, share := round(N / sum(N), 3), by = .(set, choice)]
fwrite(sm, file.path(RS, "blockopt_params_summary.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "blockopt_checks.csv"))
cat("\n=== CHOICES (share of folds) ===\n")
print(sm, row.names = FALSE, nrows = 200)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  for (pf in c("arch_preds_blockopt", "blockopt_oof")) {
    cat("  ", file.path(FA, sprintf("%s_%s_%s.csv", pf, META[[key]]$fr,
                                    META[[key]]$code)), " (local only)\n")
  }
}
cat("  ", file.path(RS, "blockopt_checks.csv"), "\n")
cat("  ", file.path(RS, "blockopt_params_summary.csv"), "\n")
cat("  ", file.path(RS, "blockopt_log.txt"), " (terminal log)\n")
cat("\nOutputs: blockopt_checks.csv, blockopt_params_summary.csv\n")