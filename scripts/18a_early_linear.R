# 18a_early_linear.R
# Early and late linear fusion arms, with and without time-varying effects
# (PREREGISTRATION.md section 4; PREREGISTRATION_ADDENDUM_2.md, A10).
#
# Every linear arm except the registered block-standardised early arms
# (18b_early_blockstd.R). Only out-of-fold predictions are written;
# 18c_arch_eval.R evaluates all architectures together.
#
# Frame A (death at 10 and 30 days) and Frame B (death at 10, 30, 90, 180
# and 365 days), each in the EHR+ECG cohort and the three-modality cohort
# (sets EHR+ECG, EHR+CTPA, EHR+ECG+CTPA), with the whitelist and folds of
# 04_frames.py. In each training fold, features are median-imputed,
# non-binary features winsorised at their training 0.5th and 99.5th
# percentiles (D9) and all standardised. Every model is a ridge glmnet with
# standardize = FALSE (D10), an 80-value penalty path down to 1e-7 of the
# largest (D11) and lambda from inner 3-fold subject-grouped CV. PH arms are
# ridge-Cox, with a Breslow baseline on centred linear predictors so that it
# cannot overflow (D9); TV arms are a ridge logistic model per horizon.
#
# Arms: EHR alone; early plain concatenation, without block scaling (D10);
# early PH with block penalties on the full set (EHR factor 1, others 1, 4,
# 16, 64 or 256, chosen by inner CV; IPF, Boulesteix et al. 2017); early PH
# with the EHR unpenalised (Herrmann et al. 2021); late nested fusion (a
# ridge-Cox score per modality, out of fold inside the training fold, then
# a Cox stack for PH or a logistic stack per horizon for TV).
#
# Checks:  missing predictions (FAIL); risks of exactly 0 or 1, EHR-alone
#          AUROC more than 0.02 from 06a, over 20% of penalty choices at
#          the least-penalised end of the path (WARN).
# Notes:   stops if arch_preds_*.csv from an earlier run remain. Tasks
#          (frame, cohort, seed) run in parallel with future.apply after a
#          timing probe that pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_linear_*.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: data/processed/fusion_arch/arch_preds_linear_*.csv (patient
#          level, not in the repo); results/early_tuning.csv,
#          results/early_lambda.csv, results/early_checks.csv,
#          results/early_warnings.csv, results/early_log.txt
# Usage:   Rscript scripts/18a_early_linear.R [--probe-only]

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
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
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
arch_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  out <- list()
  tune <- list()
  edges <- list()
  put <- function(model, set, h, idx, val) {
    k <- paste(model, set, h, sep = "|")
    if (is.null(out[[k]])) out[[k]] <<- rep(NA_real_, n)
    out[[k]][idx] <<- val
  }
  rec <- function(model, set, h, cv) {
    K <- length(cv$lambda)
    i <- which(cv$lambda == cv$lambda.min)[1]
    edges[[length(edges) + 1]] <<- data.table(
      model = model, set = set, horizon = h, nlam = K,
      least = i == K, most = i == 1)
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
      tt <- tk$t
      ee <- tk$e
      cox_cv <- function(Xm, rows, fid, pf = NULL) {
        if (is.null(pf)) pf <- rep(1, ncol(Xm))
        cv.glmnet(Xm[rows, , drop = FALSE], Surv(tt[rows], ee[rows]),
                  family = "cox", alpha = 0, foldid = fid, nlambda = NLAM,
                  lambda.min.ratio = LMR, type.measure = "deviance",
                  penalty.factor = pf, standardize = FALSE)
      }
      lpred <- function(cv, Xm, rows) {
        as.vector(predict(cv, Xm[rows, , drop = FALSE], s = "lambda.min"))
      }
      ph_risk <- function(lp_tr, lp_te) {
        c0 <- max(lp_tr)
        sapply(H, function(h) {
          hs <- breslow(tt[tr], ee[tr], lp_tr, h + 0.5)
          if (hs == 0) return(rep(0, length(lp_te)))
          -expm1(-hs * exp(lp_te - c0))
        })
      }
      ph_model <- function(Xm, pf = NULL) {
        cv <- cox_cv(Xm, tr, inner, pf)
        list(risk = matrix(ph_risk(lpred(cv, Xm, tr), lpred(cv, Xm, te)),
                           ncol = length(H)),
             score = -min(cv$cvm), lp_te = lpred(cv, Xm, te), cv = cv)
      }
      tv_model <- function(Xm, model, set) {
        matrix(sapply(H, function(h) {
          y <- as.numeric(ee == 1 & tt <= h + 0.5)
          cv <- suppressWarnings(cv.glmnet(
            Xm[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0,
            foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
            type.measure = "deviance", standardize = FALSE))
          rec(model, set, h, cv)
          as.vector(predict(cv, Xm[te, , drop = FALSE], s = "lambda.min",
                            type = "response"))
        }), ncol = length(H))
      }
      # EHR alone
      fe <- ph_model(B$ehr)
      rec("EHR alone PH", "EHR", NA_real_, fe$cv)
      tve <- tv_model(B$ehr, "EHR alone TV", "EHR")
      for (i in seq_along(H)) {
        put("EHR alone PH", "EHR", H[i], te, fe$risk[, i])
        put("EHR alone TV", "EHR", H[i], te, tve[, i])
      }
      # early, plain concatenation, per set
      for (sn in names(tk$sets)) {
        Xp <- do.call(cbind, B[tk$sets[[sn]]])
        fp <- ph_model(Xp)
        rec("Early PH plain concatenation", sn, NA_real_, fp$cv)
        tp <- tv_model(Xp, "Early TV plain concatenation", sn)
        for (i in seq_along(H)) {
          put("Early PH plain concatenation", sn, H[i], te, fp$risk[, i])
          put("Early TV plain concatenation", sn, H[i], te, tp[, i])
        }
      }
      # early, block-aware, full set
      full <- tk$sets[[length(tk$sets)]]
      fsn <- names(tk$sets)[length(tk$sets)]
      Xs <- do.call(cbind, B[full])
      blk <- rep(full, sapply(B[full], ncol))
      oth <- setdiff(full, "ehr")
      grid <- as.matrix(expand.grid(rep(list(IPF), length(oth))))
      colnames(grid) <- oth
      best <- NULL
      for (g in seq_len(nrow(grid))) {
        pf <- ifelse(blk == "ehr", 1, grid[g, match(blk, oth)])
        f <- ph_model(Xs, pf)
        if (is.null(best) || f$score > best$score) {
          best <- f
          best$g <- grid[g, ]
        }
      }
      rec("Early PH block penalties", fsn, NA_real_, best$cv)
      tune[[length(tune) + 1]] <- data.table(
        frame = tk$fr, cohort = tk$co, seed = tk$seed, fold = k,
        choice = paste(sprintf("%s x%g", toupper(oth), best$g),
                       collapse = ", "))
      f0 <- ph_model(Xs, ifelse(blk == "ehr", 0, 1))
      rec("Early PH EHR unpenalised", fsn, NA_real_, f0$cv)
      for (i in seq_along(H)) {
        put("Early PH block penalties", fsn, H[i], te, best$risk[, i])
        put("Early PH EHR unpenalised", fsn, H[i], te, f0$risk[, i])
      }
      # late, nested: inner out-of-fold modality scores
      mods <- unique(unlist(tk$sets))
      oof <- list()
      lte <- list()
      for (m in mods) {
        v <- rep(NA_real_, n)
        for (j in 1:3) {
          a <- tr[inner != j]
          b <- tr[inner == j]
          cv <- cox_cv(B[[m]], a, grouped(a, tk$seed * 1000 + k * 10 + j, 3))
          rec("Late, modality model (inner)", toupper(m), NA_real_, cv)
          v[b] <- lpred(cv, B[[m]], b)
        }
        oof[[m]] <- v
        if (m == "ehr") {
          lte[[m]] <- fe$lp_te
        } else {
          fm <- ph_model(B[[m]])
          rec("Late, modality model", toupper(m), NA_real_, fm$cv)
          lte[[m]] <- fm$lp_te
        }
      }
      for (sn in names(tk$sets)) {
        ms <- tk$sets[[sn]]
        Ztr <- sapply(ms, function(m) wz(oof[[m]], tr)[tr])
        Zte <- sapply(ms, function(m) {
          q <- quantile(oof[[m]][tr], c(0.005, 0.995), names = FALSE)
          w <- pmin(pmax(oof[[m]][tr], q[1]), q[2])
          (pmin(pmax(lte[[m]], q[1]), q[2]) - mean(w)) / sd(w)
        })
        Ztr <- matrix(Ztr, ncol = length(ms))
        Zte <- matrix(Zte, ncol = length(ms))
        D <- data.frame(t = tt[tr], e = ee[tr], Ztr)
        fc <- coxph(Surv(t, e) ~ ., data = D, ties = "efron")
        b <- coef(fc)
        lp_tr <- as.vector(Ztr %*% b)
        lp_te <- as.vector(Zte %*% b)
        rk <- matrix(ph_risk(lp_tr, lp_te), ncol = length(H))
        for (i in seq_along(H)) {
          h <- H[i]
          y <- as.numeric(ee[tr] == 1 & tt[tr] <= h + 0.5)
          g <- suppressWarnings(glm.fit(cbind(1, Ztr), y,
                                        family = binomial()))
          cf <- g$coefficients
          cf[is.na(cf)] <- 0
          put("Late PH", sn, h, te, rk[, i])
          put("Late TV", sn, h, te,
              plogis(as.vector(cbind(1, Zte) %*% cf)))
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
  list(fr = tk$fr, co = tk$co, seed = tk$seed, out = out,
       tune = rbindlist(tune), warn = warn, notes = notes,
       edges = rbindlist(edges)[, `:=`(frame = tk$fr, cohort = tk$co)])
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
old <- list.files(OD, pattern = "^arch_preds_.*\\.csv$")
chk("previous prediction files in fusion_arch (archive them to v1)",
    length(old), if (length(old) == 0) "PASS" else "FAIL")
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in names(COH)) {
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  mods <- unique(unlist(COH[[co]]$sets))
  X <- list()
  for (m in mods) {
    cols <- wl[modality == m, column]
    cols <- intersect(cols, names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cols])
    storage.mode(Xm) <- "double"
    X[[m]] <- Xm
    chk(sprintf("%s %s: %s features / rows without a record", fr, co,
                toupper(m)),
        sprintf("%d / %d", ncol(Xm), sum(rowSums(!is.na(Xm)) == 0)),
        if (ncol(Xm) > 0) "PASS" else "FAIL")
  }
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  META[[paste(fr, co)]] <- list(fr = fr, co = co, sub = sub,
                                code = COH[[co]]$code)
  for (s in seeds) {
    fs <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(fs)) stop("fold missing: ", fr, " ", co, " seed ", s)
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = fs, subj = sub$subject_id,
      t = sub$time_death, e = as.integer(sub$event_death == 1),
      X = X, sets = COH[[co]]$sets, hz = HZ[[fr]])
  }
}
cat(sprintf("  %d tasks (frame x cohort x seed), 5 folds each\n",
            length(tasks)))
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "early_checks.csv"))
  cat("\nSTOP: a check failed. see: early_checks.csv\n")
  quit(status = 1)
}

cat("\n=== TIMING PROBE: the slowest kinds of fit, per frame x cohort ===\n")
# With the longer penalty path (D11), ridge-Cox fits on the wide feature
# sets are far slower where features outnumber deaths. The probe times
# one wide Cox fit (the full set), one narrow Cox fit (EHR) and one wide
# logistic fit on fold 0 of each frame x cohort, in parallel, and counts
# how many of each a fold needs.
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
  y <- as.numeric(tk$e == 1 & tk$t <= min(tk$hz) + 0.5)
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
  c(big = tm(cx(Xs)), small = tm(cx(B$ehr)),
    bin = tm(suppressWarnings(cv.glmnet(
      Xs[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0,
      foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
      type.measure = "deviance", standardize = FALSE))))
}
# fits per fold in this stage: wide (sets with the ECG, the block grid,
# EHR unpenalised, the late ECG models), narrow, and logistic
fits_per_fold <- function(tk) {
  mods <- unique(unlist(tk$sets))
  se <- sum(sapply(tk$sets, function(s) "ecg" %in% s))
  ns <- length(tk$sets)
  G <- 5^(length(tk$sets[[ns]]) - 1)
  c(big = se + G + 1 + 4 * ("ecg" %in% mods),
    small = 1 + (ns - se) + 3 + 4 * ("ctpa" %in% mods),
    bin = length(tk$hz) * (1 + ns))
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
                    "(x%d), logistic %4.1f s (x%d): a fold %4.1f min\n"),
              kd, u[["big"]], nf[["big"]], u[["small"]], nf[["small"]],
              u[["bin"]], nf[["bin"]], fold_s[kd] / 60))
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
PROG <- file.path(OD, "progress_18a.txt")
cat(sprintf("%s  started, %d tasks x 5 folds\n",
            format(Sys.time(), "%H:%M:%S"), length(tasks)), file = PROG)
for (i in seq_along(tasks)) tasks[[i]]$progress <- PROG
cat("  progress: one line per finished fold in", PROG, "\n")
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, arch_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)

# ---------------- assemble ----------------
cat("\n=== ASSEMBLE ===\n")
TUNE <- rbindlist(lapply(res, function(r) r$tune))
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
  nna <- sum(is.na(L$risk))
  rg <- range(L$risk, na.rm = TRUE)
  nsat <- sum(L$risk <= 0 | L$risk >= 1, na.rm = TRUE)
  chk(sprintf("%s: missing predictions", key), nna,
      if (nna == 0) "PASS" else "FAIL")
  chk(sprintf("%s: risks of exactly 0 or 1 (range)", key),
      sprintf("%d (%.2e to %.6f)", nsat, rg[1], rg[2]),
      if (nsat == 0) "PASS" else "WARN")
  hh <- if (M$fr == "A") 30 else c(30, 365)
  l1 <- fread(file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv",
                                              M$fr, M$code)),
              select = c("hadm_id", sprintf("ehr.cox.death.r%d", hh)))
  l1 <- l1[, lapply(.SD, mean), by = hadm_id]
  for (h in hh) {
    y <- as.integer(M$sub$event_death == 1 & M$sub$time_death <= h + 0.5)
    auc <- function(s) {
      r <- frank(s)
      n1 <- sum(y)
      (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1))
    }
    a1 <- auc(L[model == "EHR alone PH" & horizon == h][
      match(M$sub$hadm_id, hadm_id), risk])
    a0 <- auc(l1[match(M$sub$hadm_id, hadm_id)][[
      sprintf("ehr.cox.death.r%d", h)]])
    chk(sprintf("%s: EHR alone AUROC at %d d, this stage vs stage 6a",
                key, h), sprintf("%.4f vs %.4f", a1, a0),
        if (abs(a1 - a0) <= 0.02) "PASS" else "WARN")
  }
  fwrite(L, file.path(OD, sprintf("arch_preds_linear_%s_%s.csv", M$fr,
                                  M$code)))
}
EDG <- rbindlist(lapply(res, function(r) r$edges))
esum <- EDG[, .(fits = .N, median_path_length = median(nlam),
                share_least_penalised = mean(least),
                share_most_penalised = mean(most)),
            by = .(frame, cohort, model, set, horizon)]
setorder(esum, frame, cohort, model, set, horizon)
fwrite(esum, file.path(RS, "early_lambda.csv"))
eall <- EDG[, .(share = mean(least)), by = model]
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
fwrite(wsum, file.path(RS, "early_warnings.csv"))
chk("warnings from model fitting (see early_warnings.csv)",
    sum(wsum$count), "INFO")
chk("glmnet tie-handling notices (harmless; Breslow ties throughout)",
    sum(sapply(res, function(r) r$notes)), "INFO")
tsum <- TUNE[, .(folds = .N), by = .(frame, cohort, choice)]
tsum[, share := folds / sum(folds), by = .(frame, cohort)]
setorder(tsum, frame, cohort, -folds)
fwrite(tsum, file.path(RS, "early_tuning.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "early_checks.csv"))
cat("\n=== BLOCK PENALTIES CHOSEN (share of folds) ===\n")
print(tsum, row.names = FALSE)
cat("\n=== PENALTY CHOICES AT EITHER END OF THE PATH ===\n")
print(EDG[, .(fits = .N, least = round(mean(least), 3),
              most = round(mean(most), 3)), by = model], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (arch_preds_linear_*.csv, patient-level, local only)\n")
for (f in c("early_tuning.csv", "early_lambda.csv", "early_checks.csv",
            "early_warnings.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "early_log.txt"), " (terminal log)\n")
cat("\nOutputs: early_tuning.csv, early_lambda.csv, early_checks.csv,",
    "early_warnings.csv\n")
