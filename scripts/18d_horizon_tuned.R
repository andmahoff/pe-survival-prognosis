# 18d_horizon_tuned.R
# Horizon-tuned proportional-hazards models (PREREGISTRATION_ADDENDUM_2.md,
# A11). Writes out-of-fold predictions only; 18c_arch_eval.R evaluates them
# with every other model.
#
# Each model is a Cox model fitted on every training death, as in the PH
# models, but with its tuning chosen separately for each horizon: the
# penalty (and, for the block-penalty model, the block weights) that gives
# the best inner-cross-validated log-likelihood of the predicted risk at
# that horizon. The same scoring rule tunes the TV models, so horizon-tuned
# minus PH is the gain from targeting the horizon, and TV minus
# horizon-tuned is the effect of fitting on the horizon's label (only
# deaths before the horizon count) instead of every death.
#
# Frames, cohorts, sets, horizons, folds (every seed), preprocessing, inner
# 3-fold subject-grouped folds, ridge penalty, standardize = FALSE, the
# penalty path (80 values down to 1e-7 of the largest, D11) and Breslow
# risks on centred linear predictors are as in 18a_early_linear.R and
# 18b_early_blockstd.R. Models, each tuned separately at every horizon: EHR
# alone PH (ridge-Cox on the EHR features); Early PH (block-standardised,
# per set); Early PH with block penalties (full set, block weights from the
# 18a grid, 1 to 256); and Late PH (nested, per set: modality models tuned
# at the horizon, then a Cox stack on their inner out-of-fold scores).
#
# Checks:  18b's Early PH reproduced with cv.glmnet on the same data
#          (largest absolute risk difference below 1e-6); penalty path
#          identical to cv.glmnet's; no missing risks (all FAIL). Risks of
#          exactly 0 or 1, and over 20% of choices at the least-penalised
#          end of the path (WARN).
# Notes:   stops if 18b predictions are missing or 18d predictions from an
#          earlier run remain. Tasks (frame, cohort, seed) run in parallel
#          with future.apply after a timing probe that pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, fusion_arch/arch_preds_blockstd_*.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: data/processed/fusion_arch/arch_preds_htuned_*.csv (patient
#          level, not in the repo); results/htune_checks.csv,
#          results/htune_blocks.csv, results/htune_lambda.csv,
#          results/htune_warnings.csv, results/htune_log.txt
# Usage:   Rscript scripts/18d_horizon_tuned.R [--probe-only]

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
IPF <- c(1, 4, 16, 64, 256)
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
ht_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  nh <- length(H)
  tt <- tk$t
  ee <- tk$e
  out <- list()
  blocks <- list()
  lams <- list()
  put <- function(model, set, h, idx, val) {
    k <- paste(model, set, h, sep = "|")
    if (is.null(out[[k]])) out[[k]] <<- rep(NA_real_, n)
    out[[k]][idx] <<- val
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
  # Breslow baseline at every cut for every column of L (one column per
  # penalty value), centred at each column's maximum
  bres_mat <- function(t0, e0, L, cuts) {
    o <- order(t0)
    t0 <- t0[o]
    e0 <- e0[o]
    L <- L[o, , drop = FALSE]
    cm <- apply(L, 2, max)
    W <- exp(sweep(L, 2, cm))
    nr <- nrow(W)
    RSm <- apply(W[nr:1, , drop = FALSE], 2, cumsum)
    RSm <- matrix(RSm, nrow = nr)[nr:1, , drop = FALSE]
    ut <- unique(t0[e0 == 1])
    Hm <- matrix(0, length(cuts), ncol(L))
    if (length(ut)) {
      d <- tabulate(match(t0[e0 == 1], ut), length(ut))
      inc <- d / RSm[match(ut, t0), , drop = FALSE]
      CH <- matrix(apply(inc, 2, cumsum), nrow = length(ut))
      for (i in seq_along(cuts)) {
        j <- sum(ut <= cuts[i])
        if (j > 0) Hm[i, ] <- CH[j, ]
      }
    }
    list(H = Hm, cm = cm)
  }
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
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
  # penalty path on `rows`, and the inner-CV log-likelihood loss of the
  # risk at every horizon for every penalty value
  ht_path <- function(Xm, rows, fid, pf) {
    fit <- glmnet(Xm[rows, , drop = FALSE], Surv(tt[rows], ee[rows]),
                  family = "cox", alpha = 0, nlambda = NLAM,
                  lambda.min.ratio = LMR, standardize = FALSE,
                  penalty.factor = pf)
    lam <- fit$lambda
    loss <- matrix(0, nh, length(lam))
    for (j in sort(unique(fid))) {
      a <- rows[fid != j]
      b <- rows[fid == j]
      fj <- glmnet(Xm[a, , drop = FALSE], Surv(tt[a], ee[a]),
                   family = "cox", alpha = 0, lambda = lam,
                   standardize = FALSE, penalty.factor = pf)
      La <- predict(fj, Xm[a, , drop = FALSE], s = lam)
      Lb <- predict(fj, Xm[b, , drop = FALSE], s = lam)
      br <- bres_mat(tt[a], ee[a], La, H + 0.5)
      Eb <- exp(sweep(Lb, 2, br$cm))
      for (i in seq_len(nh)) {
        P <- -expm1(-sweep(Eb, 2, br$H[i, ], "*"))
        P <- pmin(pmax(P, 1e-12), 1 - 1e-12)
        y <- as.numeric(ee[b] == 1 & tt[b] <= H[i] + 0.5)
        loss[i, ] <- loss[i, ] - colSums(y * log(P) + (1 - y) * log(1 - P))
      }
    }
    list(fit = fit, lam = lam, loss = loss, X = Xm, rows = rows)
  }
  lp_at <- function(obj, k, rows) {
    as.vector(predict(obj$fit, obj$X[rows, , drop = FALSE],
                      s = obj$lam[k]))
  }
  risk_h <- function(lp_tr, lp_te, rows_tr, h) {
    hs <- breslow(tt[rows_tr], ee[rows_tr], lp_tr, h + 0.5)
    if (hs == 0) return(rep(0, length(lp_te)))
    -expm1(-hs * exp(lp_te - max(lp_tr)))
  }
  warn <- character(0)
  notes <- 0L
  folds_run <- if (is.null(tk$only)) sort(unique(tk$fold)) else tk$only
  for (k in folds_run) {
    t_fold <- Sys.time()
    nw0 <- length(warn)
    withCallingHandlers({
      tr <- which(tk$fold != k)
      te <- which(tk$fold == k)
      inner <- grouped(tr, tk$seed * 100 + k, 3)
      B <- lapply(tk$X, prep, tr = tr)
      bs <- function(ms) {
        do.call(cbind, lapply(ms, function(m) B[[m]] / sqrt(ncol(B[[m]]))))
      }
      # a horizon-tuned model: choose per horizon, predict at the horizon
      ht_model <- function(obj, model, set) {
        for (i in seq_len(nh)) {
          kk <- which.min(obj$loss[i, ])
          r <- risk_h(lp_at(obj, kk, tr), lp_at(obj, kk, te), tr, H[i])
          put(model, set, H[i], te, r)
          lams[[length(lams) + 1]] <<- data.table(
            model = model, set = set, horizon = H[i], fold = k,
            log10_lambda = log10(obj$lam[kk]),
            least = kk == length(obj$lam), most = kk == 1)
        }
      }
      # EHR alone
      pe <- ht_path(B$ehr, tr, inner, rep(1, ncol(B$ehr)))
      ht_model(pe, "EHR alone PH horizon-tuned", "EHR")
      # early, block-standardised, per set (+ 18b reproduction check)
      for (sn in names(tk$sets)) {
        Xc <- bs(tk$sets[[sn]])
        po <- ht_path(Xc, tr, inner, rep(1, ncol(Xc)))
        ht_model(po, "Early PH horizon-tuned", sn)
        cv <- cv.glmnet(Xc[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                        family = "cox", alpha = 0, foldid = inner,
                        nlambda = NLAM, lambda.min.ratio = LMR,
                        type.measure = "deviance", standardize = FALSE)
        same <- length(cv$lambda) == length(po$lam) &&
          isTRUE(all.equal(cv$lambda, po$lam))
        put("check path identical", sn, H[1], te, as.numeric(same))
        l18 <- function(rows) {
          as.vector(predict(cv, Xc[rows, , drop = FALSE], s = "lambda.min"))
        }
        for (i in seq_len(nh)) {
          put("check 18b Early PH", sn, H[i], te,
              risk_h(l18(tr), l18(te), tr, H[i]))
          lams[[length(lams) + 1]] <- data.table(
            model = "Early PH (whole follow-up, 18b)", set = sn,
            horizon = H[i], fold = k,
            log10_lambda = log10(cv$lambda.min),
            least = cv$lambda.min == min(cv$lambda),
            most = cv$lambda.min == max(cv$lambda))
        }
      }
      # early, block penalties, full set
      full <- tk$sets[[length(tk$sets)]]
      fsn <- names(tk$sets)[length(tk$sets)]
      Xs <- do.call(cbind, B[full])
      blk <- rep(full, sapply(B[full], ncol))
      oth <- setdiff(full, "ehr")
      grid <- as.matrix(expand.grid(rep(list(IPF), length(oth))))
      colnames(grid) <- oth
      paths <- lapply(seq_len(nrow(grid)), function(g) {
        ht_path(Xs, tr, inner,
                ifelse(blk == "ehr", 1, grid[g, match(blk, oth)]))
      })
      for (i in seq_len(nh)) {
        best <- sapply(paths, function(p) min(p$loss[i, ]))
        g <- which.min(best)
        kk <- which.min(paths[[g]]$loss[i, ])
        r <- risk_h(lp_at(paths[[g]], kk, tr), lp_at(paths[[g]], kk, te),
                    tr, H[i])
        put("Early PH block penalties horizon-tuned", fsn, H[i], te, r)
        blocks[[length(blocks) + 1]] <- data.table(
          frame = tk$fr, cohort = tk$co, horizon = H[i], fold = k,
          choice = paste(sprintf("%s x%g", toupper(oth), grid[g, ]),
                         collapse = ", "))
        lams[[length(lams) + 1]] <- data.table(
          model = "Early PH block penalties horizon-tuned", set = fsn,
          horizon = H[i], fold = k,
          log10_lambda = log10(paths[[g]]$lam[kk]),
          least = kk == length(paths[[g]]$lam), most = kk == 1)
      }
      # late, nested: modality models tuned at each horizon
      mods <- unique(unlist(tk$sets))
      oof <- list()
      lte <- list()
      for (m in mods) {
        V <- matrix(NA_real_, n, nh)
        for (j in 1:3) {
          a <- tr[inner != j]
          b <- tr[inner == j]
          pa <- ht_path(B[[m]], a,
                        grouped(a, tk$seed * 1000 + k * 10 + j, 3),
                        rep(1, ncol(B[[m]])))
          for (i in seq_len(nh)) {
            V[b, i] <- lp_at(pa, which.min(pa$loss[i, ]), b)
          }
        }
        oof[[m]] <- V
        pm <- if (m == "ehr") pe else {
          ht_path(B[[m]], tr, inner, rep(1, ncol(B[[m]])))
        }
        lte[[m]] <- sapply(seq_len(nh), function(i) {
          lp_at(pm, which.min(pm$loss[i, ]), te)
        })
        lte[[m]] <- matrix(lte[[m]], ncol = nh)
      }
      for (sn in names(tk$sets)) {
        ms <- tk$sets[[sn]]
        for (i in seq_len(nh)) {
          Ztr <- sapply(ms, function(m) wz(oof[[m]][, i], tr)[tr])
          Zte <- sapply(ms, function(m) {
            v <- oof[[m]][tr, i]
            q <- quantile(v, c(0.005, 0.995), names = FALSE)
            w <- pmin(pmax(v, q[1]), q[2])
            (pmin(pmax(lte[[m]][, i], q[1]), q[2]) - mean(w)) / sd(w)
          })
          Ztr <- matrix(Ztr, ncol = length(ms))
          Zte <- matrix(Zte, ncol = length(ms))
          D <- data.frame(t = tt[tr], e = ee[tr], Ztr)
          fc <- coxph(Surv(t, e) ~ ., data = D, ties = "efron")
          b <- coef(fc)
          put("Late PH horizon-tuned", sn, H[i], te,
              risk_h(as.vector(Ztr %*% b), as.vector(Zte %*% b), tr, H[i]))
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
       notes = notes, blocks = rbindlist(blocks),
       lams = rbindlist(lams)[, `:=`(frame = tk$fr, cohort = tk$co,
                                     seed = tk$seed)])
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
  bf <- file.path(OD, sprintf("arch_preds_blockstd_%s_%s.csv", fr,
                              COH[[co]]$code))
  chk(sprintf("%s %s: stage 18b predictions found", fr, co), basename(bf),
      if (file.exists(bf)) "PASS" else "FAIL")
  X <- list()
  for (m in unique(unlist(COH[[co]]$sets))) {
    cols <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cols])
    storage.mode(Xm) <- "double"
    X[[m]] <- Xm
  }
  fo <- folds[frame == fr & cohort == co]
  META[[paste(fr, co)]] <- list(fr = fr, co = co, sub = sub, bf = bf,
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
old <- list.files(OD, pattern = "^arch_preds_htuned_.*\\.csv$")
chk("previous stage 18d prediction files (archive them to v1)",
    length(old), if (length(old) == 0) "PASS" else "FAIL")
cat(sprintf("  %d tasks (frame x cohort x seed), 5 folds each\n",
            length(tasks)))
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "htune_checks.csv"))
  cat("\nSTOP: a check failed. see: htune_checks.csv\n")
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
  c(big = tm(cx(Xs)), small = tm(cx(B$ehr)))
}
# a horizon-tuned path costs about one cv.glmnet: wide fits are the early
# models and the 18b check on sets with the ECG, the block grid and the
# late ECG models; narrow fits are the rest
fits_per_fold <- function(tk) {
  mods <- unique(unlist(tk$sets))
  se <- sum(sapply(tk$sets, function(s) "ecg" %in% s))
  ns <- length(tk$sets)
  G <- 5^(length(tk$sets[[ns]]) - 1)
  c(big = 2 * se + G + 4 * ("ecg" %in% mods),
    small = 1 + 2 * (ns - se) + 3 + 4 * ("ctpa" %in% mods))
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
                    "(x%d): a fold %4.1f min\n"),
              kd, u[["big"]], nf[["big"]], u[["small"]], nf[["small"]],
              fold_s[kd] / 60))
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
PROG <- file.path(OD, "progress_18d.txt")
cat(sprintf("%s  started, %d tasks x 5 folds\n",
            format(Sys.time(), "%H:%M:%S"), length(tasks)), file = PROG)
for (i in seq_along(tasks)) tasks[[i]]$progress <- PROG
cat("  progress: one line per finished fold in", PROG, "\n")
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, ht_task, future.packages = PKG,
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
  same <- L[model == "check path identical", risk]
  chk(sprintf("%s: penalty path identical to cv.glmnet's", key),
      sprintf("%.0f%% of patient x set rows", 100 * mean(same == 1)),
      if (all(same == 1)) "PASS" else "FAIL")
  old <- fread(M$bf)
  cmp <- merge(L[model == "check 18b Early PH",
                 .(hadm_id, set, horizon, new = risk)],
               old[model == "Early PH", .(hadm_id, set, horizon,
                                          old = risk)],
               by = c("hadm_id", "set", "horizon"))
  dmax <- max(abs(cmp$new - cmp$old))
  chk(sprintf("%s: stage 18b Early PH reproduced (rows, max |diff|)", key),
      sprintf("%d, %.1e", nrow(cmp), dmax),
      if (nrow(cmp) == nrow(old[model == "Early PH"]) && dmax < 1e-6)
        "PASS" else "FAIL")
  L <- L[!grepl("^check ", model)]
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
LAM <- rbindlist(lapply(res, function(r) r$lams), use.names = TRUE)
eall <- LAM[grepl("horizon-tuned", model), .(share = mean(least)),
            by = model]
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
fwrite(wsum, file.path(RS, "htune_warnings.csv"))
chk("warnings from model fitting (see htune_warnings.csv)",
    sum(wsum$count), "INFO")
chk("glmnet tie-handling notices (harmless; Breslow ties throughout)",
    sum(sapply(res, function(r) r$notes)), "INFO")
ck <- rbindlist(CHK)
if (any(ck$status == "FAIL")) {
  fwrite(ck, file.path(RS, "htune_checks.csv"))
  cat("\nSTOP: a check failed; no predictions written.",
      "Outputs: htune_checks.csv\n")
  quit(status = 1)
}

cat("\n=== WRITE ===\n")
for (key in names(META)) {
  M <- META[[key]]
  fwrite(NEW[[key]], file.path(OD, sprintf("arch_preds_htuned_%s_%s.csv",
                                           M$fr, M$code)))
}
BLK <- rbindlist(lapply(res, function(r) r$blocks))
bsum <- BLK[, .(folds = .N), by = .(frame, cohort, horizon, choice)]
bsum[, share := folds / sum(folds), by = .(frame, cohort, horizon)]
setorder(bsum, frame, cohort, horizon, -folds)
fwrite(bsum, file.path(RS, "htune_blocks.csv"))
lsum <- LAM[, .(folds = .N, median_log10_lambda = median(log10_lambda),
                q25 = quantile(log10_lambda, 0.25),
                q75 = quantile(log10_lambda, 0.75),
                share_least_penalised = mean(least),
                share_most_penalised = mean(most)),
            by = .(frame, cohort, model, set, horizon)]
setorder(lsum, frame, cohort, model, set, horizon)
fwrite(lsum, file.path(RS, "htune_lambda.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "htune_checks.csv"))
cat("\n=== BLOCK PENALTIES CHOSEN AT EACH HORIZON (most common) ===\n")
print(bsum[, .SD[1], by = .(frame, cohort, horizon)], row.names = FALSE)
cat("\n=== PENALTY CHOSEN (median log10 lambda), early models ===\n")
print(dcast(lsum[grepl("^Early PH", model)],
            frame + cohort + model + set ~ horizon,
            value.var = "median_log10_lambda")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2)
                       else x)], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (arch_preds_htuned_*.csv, patient-level, local only)\n")
for (f in c("htune_checks.csv", "htune_blocks.csv", "htune_lambda.csv",
            "htune_warnings.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "htune_log.txt"), " (terminal log)\n")
cat("\nOutputs: htune_checks.csv, htune_blocks.csv, htune_lambda.csv,",
    "htune_warnings.csv\n")