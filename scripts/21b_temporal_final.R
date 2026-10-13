# 21b_temporal_final.R
# Temporal validation of the final models (A35): every level-1 model and the
# super learner are trained on the earlier years only and tested, untouched,
# on the later years.
#
# --frame=B is the discharge model, v7-core (the four EHR+ families with the
# same families' ECG and CTPA models; within 0.0016 AUROC of v7 in 18r).
# --frame=P is the at-diagnosis model (A32): the same four families.
# --outcomes=death,comp,cv selects outcomes (default all three; CV needs the
# death scores, so death is always fitted). Cohorts are EHR+ECG and
# three-modality.
#
# Development is admissions with approx_year_hi up to 2013
# (followup_365.csv; the year groups of 21a) and the test set is 2014
# onwards, reported pooled and as 2014-2016 and 2017 onwards. Test
# admissions of patients with any development admission are left out.
# Development folds: 5, grouped by patient, one split per seed (--seeds=K,
# default 3).
#
# Level 1 is run_task of 06q with training rows restricted to development:
# folds 1-5 give out-of-fold development scores, and fold 0 trains on all
# development rows and scores the test rows. Models are ridge-Cox retuned,
# spline Cox (EHR, D16), RSF (competing-risks forest for CV) and oblique
# RSF, Brier-tuned with tune_expand.R. Level 2 uses the candidates and super
# learner of 18r, with weights and recalibration fitted on the development
# out-of-fold scores and applied to the test rows.
#
# On the test rows: AUROC, IPA, log-likelihood skill, observed/expected and
# calibration slope at 10, 30, 90, 180 and 365 days (CV also among patients
# who did not die first), with a 2,000-sample patient-clustered bootstrap.
# Contrasts are the temporal model minus the main cross-validated model on
# the same test patients, SL minus SL-EHR and, at diagnosis, SL minus sPESI
# and PESI-10 (death). The script stops on missing inputs or years, a
# patient on both sides, missing scores or risks outside [0, 1].
#
# Inputs:  data/processed/ (frame_B.csv, discharge/, presentation/,
#          feature_whitelist.csv, super_learner/{slP,optB2}_preds_*.csv),
#          data/inputs/followup_365.csv, scripts/tune_expand.R
# Outputs: results/tv2_<FRAME>_{summary,contrasts,metrics,drift,
#          checks}.csv, results/tv2_<FRAME>_log.txt
# Usage:   Rscript scripts/21b_temporal_final.R --frame=B|P [--seeds=K]
#          [--outcomes=death,comp,cv] [--workers=N] [--probe-only]
# Notes:   timings and finished tasks are saved, so a rerun skips both.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(mgcv)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf not installed")
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
geta <- function(k, ok, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args,
                                           value = TRUE))
  if (!length(v) && !is.null(def)) return(def)
  if (length(v) != 1 || !all(strsplit(v, ",")[[1]] %in% ok)) {
    stop(sprintf("give --%s= from %s", k, paste(ok, collapse = ", ")))
  }
  v
}
FR <- geta("frame", c("P", "B"))
OCS <- strsplit(geta("outcomes", c("death", "comp", "cv"),
                     "death,comp,cv"), ",")[[1]]
OCS <- unique(c("death", OCS))
NSEED <- as.integer(sub("^--seeds=", "", c(grep("^--seeds=", args,
                                                value = TRUE), "3")[1]))
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
SLD <- file.path(PR, "super_learner")
TD <- file.path(PR, "temporal2")
RS <- file.path(P2, "results")
CKPT <- file.path(TD, paste0(FR, "_ckpt"))
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
TAG <- paste0("tv2_", FR)
source(file.path(P2, "scripts", "tune_expand.R"))
DEV_MAX <- 2013
SEEDS <- c(42L, 7L, 13L, 1L, 2L, 3L)[seq_len(min(6, NSEED))]
HZ <- c(10, 30, 90, 180, 365)
TCUTS <- c(1, 2, 4, 7, 14, 30, 60, 120, 240)
TMAX <- 366
NLAM <- 80L
LMR <- 1e-7
KNOTS <- 5L
MINU <- 10L
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5, 15, 40, 80, 160)
MTRYF <- c(0.5, 1, 3)
LEAF <- c(5, 20, 60)
OMTRYF <- c(0.5, 1, 2, 4)
BATCH <- 50L
TOL <- 0.001
MAXEXT <- 3L
NFOLD <- 5L
NBOOT <- 2000L
BSEED <- 42L
FAMS <- c("ridge-Cox retuned", "spline Cox", "RSF, Brier-tuned",
          "oblique RSF, Brier-tuned")
SHORT <- c("ridge-Cox retuned" = "cox2", "spline Cox" = "spline",
           "RSF, Brier-tuned" = "rsfb",
           "oblique RSF, Brier-tuned" = "orsfb")
SETS <- list("EHR+ECG" = list(code = "ehrecg",
                              sets = list("EHR+ECG" = c("ehr", "ecg"))),
             "3-mod" = list(code = "3mod",
                            sets = list("EHR+ECG" = c("ehr", "ecg"),
                                        "EHR+CTPA" = c("ehr", "ctpa"),
                                        "EHR+ECG+CTPA" =
                                          c("ehr", "ecg", "ctpa"))))
PKG <- c("survival", "glmnet", "mgcv", "data.table", "randomForestSRC",
         "aorsf")
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
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, paste0(TAG, "_checks.csv")))
    cat(sprintf("\nSTOP: a check failed. see: %s_checks.csv\n", TAG))
    quit(status = 1)
  }
}

# ---------------- level 1: one task = outcome x cohort x modality x seed ----
# stage 6q's run_task; the only change is the training rows: development
# rows (fold > 0) outside fold k. Fold 0 = the test rows.
run_task <- function(tk, folds_to_run = 0:5, probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
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
    sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
  }
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  out <- list(cox2 = rep(NA_real_, n), spline = rep(NA_real_, n),
              rsfb = rep(NA_real_, n), orsfb = rep(NA_real_, n))
  prm <- list()
  fails <- list()
  for (k in folds_to_run) {
    tr <- which(tk$fold != k & tk$fold > 0)
    te <- which(tk$fold == k)
    Xs <- prep(tk$X, tr)
    colnames(Xs) <- sprintf("f%d", seq_len(ncol(Xs)))
    # ridge-Cox retuned, path extended if the smallest penalty is chosen
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    lmr <- tk$lmr
    for (pass in 1:3) {
      cv <- cv.glmnet(Xs[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                      family = "cox", alpha = 0, foldid = inner,
                      nlambda = tk$nlam, lambda.min.ratio = lmr,
                      type.measure = "deviance",
                      penalty.factor = rep(1, ncol(Xs)),
                      standardize = FALSE)
      at_end <- which(cv$lambda == cv$lambda.min) == length(cv$lambda)
      if (!at_end || pass == 3) break
      lmr <- lmr / 100
    }
    out$cox2[te] <- as.vector(predict(cv, Xs[te, , drop = FALSE],
                                      s = "lambda.min"))
    # spline Cox (D14's fallback rule); linear terms ridge-penalised (D16);
    # fitted for the EHR only (ECG and CTPA use ridge-Cox, as in v4)
    nu <- apply(Xs[tr, , drop = FALSE], 2, function(v) length(unique(v)))
    if (!tk$spline) nu[] <- 0L
    sm <- colnames(Xs)[nu >= tk$minu]
    li <- setdiff(colnames(Xs), sm)
    dsp <- data.frame(t = tt[tr], e = ee[tr], Xs[tr, sm, drop = FALSE])
    dnw <- data.frame(Xs[te, sm, drop = FALSE])
    pp <- NULL
    if (length(li)) {
      dsp$L <- Xs[tr, li, drop = FALSE]
      dnw$L <- Xs[te, li, drop = FALSE]
      pp <- list(L = list(diag(length(li))))
    }
    sk <- 0L
    for (kk in if (tk$spline) unique(c(tk$knots, 4L, 3L)) else integer(0)) {
      fml <- as.formula(paste("t ~", paste(c(sprintf(
        "s(%s, k = %d, bs = \"cr\")", sm, kk), if (length(li)) "L"),
        collapse = " + ")))
      g <- tryCatch(mgcv::gam(fml, family = mgcv::cox.ph(), weights = e,
                              data = dsp, method = "REML", select = TRUE,
                              optimizer = "efs", paraPen = pp),
                    error = function(err) conditionMessage(err))
      if (is.character(g)) {
        fails[[length(fails) + 1]] <- data.table(
          fold = k, family = "spline Cox", setting = sprintf("k = %d", kk),
          message = substr(g, 1, 150))
        next
      }
      out$spline[te] <- as.vector(predict(g, newdata = dnw,
                                          type = "link"))
      sk <- kk
      break
    }
    if (sk == 0L) out$spline[te] <- out$cox2[te]
    if (!tk$spline) sk <- NA_integer_
    Ytr <- matrix(sapply(H, function(h) as.numeric(ee[tr] == 1 &
                                                      tt[tr] <= h + 0.5)),
                  ncol = length(H))
    # RSF, Brier-tuned, grid extended past its edges
    X <- tk$X
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- med[na[, 2]]
    colnames(X) <- sprintf("f%d", seq_len(ncol(X)))
    ev <- sort(unique(tt[tr][ee[tr] == 1]))
    at_h <- sapply(H, function(h) {
      v <- ev[ev <= h + 0.5]
      if (length(v)) max(v) else NA_real_
    })
    pick <- ev[unique(round(seq(1, length(ev),
                                length.out = min(tk$ntime, length(ev)))))]
    nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
    st_tr <- if (is.null(tk$cr)) ee[tr] else tk$cr[tr]
    dtr <- data.frame(time = tt[tr], status = st_tr, X[tr, , drop = FALSE])
    dte <- as.data.frame(X[te, , drop = FALSE])
    p <- ncol(X)
    gi <- 0L
    rsf_eval <- function(par) {
      gi <<- gi + 1L
      f <- rfsrc(Surv(time, status) ~ ., data = dtr, ntree = tk$ntree,
                 nodesize = par$node, mtry = par$mtry, nsplit = tk$nsplit,
                 ntime = nt, splitrule = "logrank", save.memory = TRUE,
                 seed = -as.integer(tk$seed * 1000 + k * 100 + gi))
      ro <- sapply(H, function(h) {
        i <- findInterval(h + 0.5, f$time.interest)
        if (i == 0) return(rep(0, length(tr)))
        # competing risks (CV): out-of-bag CIF of CV readmission (cause 1)
        if (is.null(tk$cr)) 1 - f$survival.oob[, i] else f$cif.oob[, i, 1]
      })
      br <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2))
      pp <- predict(f, newdata = dte)$predicted
      sc <- if (is.null(tk$cr)) pp else pp[, 1]
      rm(f)
      invisible(gc(verbose = FALSE))
      list(score = br, pred = sc)
    }
    rr <- tune_expand(
      list(node = tk$nodes,
           mtry = unique(pmin(p, ceiling(tk$mtryf * sqrt(p))))),
      rsf_eval,
      list(node = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) v * 2, lo = 1,
                       hi = floor(length(tr) / 4)),
           mtry = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) ceiling(v * 1.75), lo = 1, hi = p)),
      tol = tk$tol, max_ext = tk$maxext)
    if (!is.null(rr$fit)) out$rsfb[te] <- rr$fit$pred
    # oblique RSF, Brier-tuned, grid extended past its edges
    ps <- ncol(Xs)
    otr <- data.frame(time = tt[tr], status = ee[tr], Xs[tr, , drop = FALSE])
    ote <- as.data.frame(Xs[te, , drop = FALSE])
    nb <- max(1L, tk$ntree %/% tk$batch)
    oi <- 0L
    orsf_eval <- function(par) {
      oi <<- oi + 1L
      ro <- list()
      ms <- list()
      for (bb in seq_len(nb)) {
        f <- aorsf::orsf(otr, Surv(time, status) ~ ., n_tree = tk$batch,
                         mtry = par$mtry, leaf_min_obs = par$leaf,
                         split_min_obs = 2L * par$leaf,
                         oobag_pred_type = "mort",
                         tree_seeds = as.integer(tk$seed * 10000 +
                                                   k * 1000 + oi * 20 + bb),
                         n_thread = 1)
        ro[[bb]] <- matrix(predict(f, oobag = TRUE, pred_type = "risk",
                                   pred_horizon = H + 0.5,
                                   boundary_checks = FALSE),
                           ncol = length(H))
        ms[[bb]] <- as.vector(predict(f, new_data = ote,
                                      pred_type = "mort"))
        rm(f)
      }
      invisible(gc(verbose = FALSE))
      br <- sum(colMeans((Ytr - Reduce(`+`, ro) / nb)^2, na.rm = TRUE))
      list(score = br, pred = rowMeans(do.call(cbind, ms)))
    }
    ro <- tune_expand(
      list(leaf = tk$leaf,
           mtry = unique(pmax(1, pmin(ps, ceiling(tk$omtryf * sqrt(ps)))))),
      orsf_eval,
      list(leaf = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) v * 2, lo = 1,
                       hi = floor(length(tr) / 8)),
           mtry = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) ceiling(v * 1.75), lo = 1,
                       hi = ps)),
      tol = tk$tol, max_ext = tk$maxext)
    if (!is.null(ro$fit)) out$orsfb[te] <- ro$fit$pred
    prm[[length(prm) + 1]] <- data.table(
      fold = k, ridge_path_passes = pass, spline_k = sk,
      rsf_fails = rr$nfail, orsf_fails = ro$nfail,
      rsf_edge = if (is.null(rr$edge)) NA else any(rr$edge),
      orsf_edge = if (is.null(ro$edge)) NA else any(ro$edge))
  }
  res <- list(oc = tk$oc, co = tk$co, mod = tk$mod, seed = tk$seed,
              out = out, prm = rbindlist(prm),
              fails = rbindlist(fails, fill = TRUE),
              sec = proc.time()[["elapsed"]] - t_task)
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, tk$ckpt)
  }
  res
}

# ---------------- level 2: stage 18r's stacks and super learner ----------
nscore <- function(x, ref) {
  qnorm((findInterval(x, sort(ref)) + 0.5) / (length(ref) + 1))
}
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
stack_ph <- function(Ztr, Zte, ttr, etr, cuts) {
  D <- data.frame(t = ttr, e = etr, Ztr)
  b <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
  b[is.na(b)] <- 0
  lt <- as.vector(Ztr %*% b)
  lv <- as.vector(Zte %*% b)
  H <- bres_cuts(ttr, etr, lt, cuts)
  sapply(seq_along(cuts), function(i) -expm1(-H[i] * exp(lv - max(lt))))
}
stack_tv <- function(Ztr, Zte, ttr, etr, cuts) {
  sapply(cuts, function(cu) {
    y <- as.numeric(etr == 1 & ttr <= cu)
    g <- suppressWarnings(glm.fit(cbind(1, Ztr), y, family = binomial()))
    cf <- g$coefficients
    cf[is.na(cf)] <- 0
    plogis(as.vector(cbind(1, Zte) %*% cf))
  })
}
stack_tvc <- function(Ztr, Zte, ttr, etr, cuts, tcuts, tmax) {
  lb <- c(0, tcuts)
  ub <- c(tcuts, tmax)
  g <- log(c(0.5, sqrt(lb[-1] * ub[-1])))
  ep <- survSplit(Surv(t, e) ~ ., cut = tcuts, episode = "ep",
                  data = data.frame(t = ttr, e = etr, id = seq_along(ttr)))
  Zr <- Ztr[ep$id, , drop = FALSE]
  X <- cbind(Zr, Zr * g[ep$ep])
  b <- coef(coxph(Surv(ep$tstart, ep$t, ep$e) ~ X, ties = "efron"))
  b[is.na(b)] <- 0
  q <- ncol(Ztr)
  l0 <- as.vector(Ztr %*% b[seq_len(q)])
  l1 <- as.vector(Ztr %*% b[q + seq_len(q)])
  n0 <- as.vector(Zte %*% b[seq_len(q)])
  n1 <- as.vector(Zte %*% b[q + seq_len(q)])
  o <- order(ttr)
  ts <- ttr[o]
  es <- etr[o]
  ut <- unique(ts[es == 1])
  d <- tabulate(match(ts[es == 1], ut), length(ut))
  pos <- match(ut, ts)
  pint <- findInterval(ut, lb, left.open = TRUE)
  inc <- numeric(length(ut))
  cp <- numeric(length(lb))
  for (p in unique(pint)) {
    lp <- l0[o] + l1[o] * g[p]
    cp[p] <- max(lp)
    rs <- rev(cumsum(rev(exp(lp - cp[p]))))
    sel <- which(pint == p)
    inc[sel] <- d[sel] / rs[pos[sel]]
  }
  sapply(cuts, function(cu) {
    Hx <- rep(0, length(n0))
    for (p in unique(pint)) {
      sel <- which(pint == p & ut <= cu)
      if (length(sel)) {
        Hx <- Hx + sum(inc[sel]) * exp(n0 + n1 * g[p] - cp[p])
      }
    }
    -expm1(-Hx)
  })
}
km_risk <- function(ttr, etr, cuts, nte) {
  sf <- summary(survfit(Surv(ttr, etr) ~ 1), times = cuts, extend = TRUE)
  matrix(rep(1 - sf$surv, each = nte), nrow = nte)
}
proj_simplex <- function(v) {
  u <- sort(v, decreasing = TRUE)
  cs <- cumsum(u)
  r <- max(which(u - (cs - 1) / seq_along(u) > 0))
  pmax(v - (cs[r] - 1) / r, 0)
}
sl_weights <- function(X, y) {
  J <- ncol(X)
  a <- rep(1 / J, J)
  XtX <- crossprod(X)
  Xty <- as.vector(crossprod(X, y))
  L <- max(eigen(XtX, symmetric = TRUE, only.values = TRUE)$values)
  for (it in seq_len(20000)) {
    an <- proj_simplex(a - (as.vector(XtX %*% a) - Xty) / L)
    if (max(abs(an - a)) < 1e-12) {
      a <- an
      break
    }
    a <- an
  }
  a
}
candidates <- function(tk, a, b, cuts) {
  S <- tk$S
  sets <- tk$sets
  tt <- tk$t
  ee <- tk$e
  Z <- function(SS, fam, ms, rows) {
    matrix(sapply(ms, function(m) nscore(SS[[fam]][[m]][rows],
                                         SS[[fam]][[m]][a])),
           ncol = length(ms))
  }
  out <- list()
  if (!tk$cv) {
    out[["EHR|Kaplan-Meier"]] <- km_risk(tt[a], ee[a], cuts, length(b))
    for (fam in FAMS) {
      out[[paste0("EHR|EHR alone, ", fam)]] <-
        stack_ph(Z(S, fam, "ehr", a), Z(S, fam, "ehr", b), tt[a], ee[a],
                 cuts)
    }
    for (sn in names(sets)) {
      ms <- sets[[sn]]
      for (fam in FAMS) {
        za <- Z(S, fam, ms, a)
        zb <- Z(S, fam, ms, b)
        out[[sprintf("%s|%s PH", sn, fam)]] <-
          stack_ph(za, zb, tt[a], ee[a], cuts)
        out[[sprintf("%s|%s TV", sn, fam)]] <-
          stack_tv(za, zb, tt[a], ee[a], cuts)
        out[[sprintf("%s|%s TVC", sn, fam)]] <-
          stack_tvc(za, zb, tt[a], ee[a], cuts, tcuts = tk$tcuts,
                    tmax = tk$tmax)
      }
    }
  } else {
    D <- tk$D
    inc <- sapply(cuts, function(cu) mean(ee[a] == 1 & tt[a] <= cu))
    out[["EHR|marginal incidence"]] <- matrix(rep(inc, each = length(b)),
                                              nrow = length(b))
    for (fam in FAMS) {
      za <- cbind(Z(S, fam, "ehr", a), Z(D, fam, "ehr", a))
      zb <- cbind(Z(S, fam, "ehr", b), Z(D, fam, "ehr", b))
      out[[paste0("EHR|EHR alone, ", fam, " TV")]] <-
        stack_tv(za, zb, tt[a], ee[a], cuts)
    }
    for (sn in names(sets)) {
      ms <- sets[[sn]]
      for (fam in FAMS) {
        za <- cbind(Z(S, fam, ms, a), Z(D, fam, ms, a))
        zb <- cbind(Z(S, fam, ms, b), Z(D, fam, ms, b))
        out[[sprintf("%s|%s TV", sn, fam)]] <-
          stack_tv(za, zb, tt[a], ee[a], cuts)
      }
    }
  }
  lapply(out, function(m) matrix(m, nrow = length(b)))
}
# one super learner: fitted on the development rows (dev), applied to the
# test rows (te); weights from inner grouped 3-fold predictions on dev
sl_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  cuts <- H + 0.5
  sets <- tk$sets
  dev <- which(tk$dev)
  te <- which(!tk$dev)
  us <- unique(tk$subj[dev])
  set.seed(tk$seed * 100)
  fu <- sample(rep_len(1:3, length(us)))
  inner <- fu[match(tk$subj[dev], us)]
  out <- list()
  wts <- list()
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  PI <- NULL
  for (j in 1:3) {
    a <- dev[inner != j]
    b <- dev[inner == j]
    cj <- candidates(tk, a, b, cuts)
    if (is.null(PI)) {
      PI <- lapply(cj, function(m) matrix(NA_real_, length(dev), length(H)))
    }
    for (nm in names(cj)) PI[[nm]][inner == j, ] <- cj[[nm]]
  }
  CT <- candidates(tk, dev, te, cuts)
  Y <- matrix(sapply(cuts, function(cu) {
    as.numeric(tk$e[dev] == 1 & tk$t[dev] <= cu)
  }), ncol = length(H))
  fit_sl <- function(cn, arm, sn) {
    X <- do.call(cbind, lapply(cn, function(nm) as.vector(PI[[nm]])))
    w <- sl_weights(X, as.vector(Y))
    ens <- Reduce(`+`, Map(function(x, nm) x * CT[[nm]], w, cn))
    ein <- Reduce(`+`, Map(function(x, nm) x * PI[[nm]], w, cn))
    rc <- sapply(seq_along(H), function(i) {
      g <- tryCatch(suppressWarnings(glm.fit(
        cbind(1, qlogis(clip(ein[, i]))), Y[, i],
        family = binomial())), error = function(e) NULL)
      cf <- if (is.null(g)) c(NA, NA) else g$coefficients
      if (any(!is.finite(cf))) return(ens[, i])
      plogis(cf[1] + cf[2] * qlogis(clip(ens[, i])))
    })
    out[[paste(paste0(arm, ", recalibrated"), sn, sep = "|")]] <<-
      matrix(rc, ncol = length(H))
    wts[[length(wts) + 1]] <<- data.table(
      model = arm, set = sn, candidate = sub("^[^|]*\\|", "", cn),
      weight = w)
  }
  ce <- names(CT)[startsWith(names(CT), "EHR|")]
  fit_sl(ce, "Super learner-EHR", "EHR")
  for (sn in names(sets)) {
    ca <- c(ce, names(CT)[startsWith(names(CT), paste0(sn, "|"))])
    fit_sl(ca, "Super learner", sn)
  }
  list(job = tk$job, seed = tk$seed, out = out,
       wts = rbindlist(wts)[, `:=`(job = tk$job, seed = tk$seed)])
}

# ---------------- evaluation helpers ----------------
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(J$hz)) {
    y <- J$Y[idx, hi]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    keep <- if (is.null(J$DF)) NULL else J$DF[idx, hi] == 0
    for (k in names(J$P)) {
      raw <- J$P[[k]][idx, hi]
      key <- paste(k, J$hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, raw)
      if (!is.null(keep)) {
        v[paste(key, "auc_excl", sep = "|")] <- auc(y[keep], raw[keep])
      }
      if (k %in% J$scores) next
      p <- pmin(pmax(raw, 1e-6), 1 - 1e-6)
      ok <- pr > 0 && pr < 1
      v[paste(key, "ipa", sep = "|")] <- if (ok) {
        1 - mean((y - p)^2) / (pr * (1 - pr))
      } else NA_real_
      v[paste(key, "ll", sep = "|")] <- if (ok) {
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      } else NA_real_
      v[paste(key, "oe", sep = "|")] <- pr / mean(p)
    }
  }
  v
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by), length(J$rows_by),
                                     replace = TRUE)], use.names = FALSE)
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat(sprintf("=== BUILD (frame %s, outcomes %s, %d seeds) ===\n", FR,
            paste(OCS, collapse = ", "), length(SEEDS)))
DIR <- if (FR == "P") PP else file.path(PR, "discharge")
fframe <- if (FR == "P") file.path(PP, "frame_P_outcomes.csv") else
  file.path(PR, "frame_B.csv")
fm <- if (FR == "P") {
  c(ehr = "ehrP.csv", ecg = "ecgP.csv", ctpa = "ctpaP.csv")
} else {
  c(ehr = "ehrplusB.csv", ecg = "ecgB.csv", ctpa = "ctpaB.csv")
}
need <- c(fframe, file.path(DIR, fm), file.path(IN, "followup_365.csv"),
          file.path(PR, "feature_whitelist.csv"),
          file.path(P2, "scripts", "tune_expand.R"))
if (FR == "P") need <- c(need, file.path(PP, "scores_P.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- unique(fread(fframe), by = "hadm_id")
if (FR == "B") fp[, in_three_mod_p := in_three_mod]
fu <- unique(fread(file.path(IN, "followup_365.csv")), by = "hadm_id")
fp[, year := as.numeric(fu$approx_year_hi[match(hadm_id, fu$hadm_id)])]
chk("admissions with an approximate year",
    sprintf("%.1f%% of %d", 100 * mean(is.finite(fp$year)), nrow(fp)),
    if (all(is.finite(fp$year))) "PASS" else "FAIL")
stop_if_fail()
dev_subj <- unique(fp[year <= DEV_MAX, subject_id])
fp[, part := ifelse(year <= DEV_MAX, "development",
                    ifelse(subject_id %in% dev_subj, "excluded",
                           ifelse(year <= 2016, "2014-2016",
                                  "2017 onwards")))]
chk("test admissions left out (patient also in development)",
    sum(fp$part == "excluded"), "INFO")
fp <- fp[part != "excluded"]
both <- length(intersect(fp[part == "development", subject_id],
                         fp[part != "development", subject_id]))
chk("patients in both development and test", both,
    if (both == 0) "PASS" else "FAIL")
stop_if_fail()
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mat <- function(f, rows, keep = NULL) {
  D <- unique(fread(file.path(DIR, f)), by = "hadm_id")
  cc <- setdiff(names(D), c("hadm_id", "subject_id"))
  if (!is.null(keep)) cc <- intersect(keep, cc)
  X <- as.matrix(D[match(fp$hadm_id[rows], hadm_id), ..cc])
  storage.mode(X) <- "double"
  list(X = X, miss = sum(is.na(match(fp$hadm_id[rows], D$hadm_id))))
}
OUTD <- list(
  death = list(t = "time_death", e = function(d) d$event_death == 1),
  comp = list(t = "time_first", e = function(d) d$event_comp == 1),
  cv = list(t = "time_first", e = function(d) d$event_cr == 1))
DR <- list()
tasks <- list()
META <- list()
for (co in names(SETS)) {
  C <- SETS[[co]]
  rows <- if (co == "EHR+ECG") rep(TRUE, nrow(fp)) else
    fp$in_three_mod_p == 1
  sub <- fp[rows]
  mods <- unique(unlist(C$sets))
  MX <- list()
  for (m in mods) {
    mm <- mat(fm[[m]], rows, if (m == "ctpa") wl[modality == "ctpa",
                                                 column] else NULL)
    chk(sprintf("%s %s: features (admissions missing)", co, toupper(m)),
        sprintf("%d (%d)", ncol(mm$X), mm$miss),
        if (ncol(mm$X) > 0 && mm$miss == 0) "PASS" else "FAIL")
    MX[[m]] <- mm$X
  }
  stop_if_fail()
  isdev <- sub$part == "development"
  fmap <- list()
  for (s in SEEDS) {
    us <- unique(sub$subject_id[isdev])
    set.seed(s)
    fu_ <- sample(rep_len(seq_len(NFOLD), length(us)))
    f <- integer(nrow(sub))
    f[isdev] <- fu_[match(sub$subject_id[isdev], us)]
    fmap[[as.character(s)]] <- f
  }
  for (oc in OCS) {
    tt <- sub[[OUTD[[oc]]$t]]
    ee <- as.integer(OUTD[[oc]]$e(sub))
    for (pt in c("development", "2014-2016", "2017 onwards")) {
      r <- sub$part == pt
      DR[[length(DR) + 1]] <- data.table(
        cohort = co, outcome = oc, period = pt, admissions = sum(r),
        events_30 = sum(ee[r] == 1 & tt[r] <= 30.5),
        events_365 = sum(ee[r] == 1 & tt[r] <= 365.5),
        km_risk_365 = round(1 - summary(survfit(Surv(tt[r], ee[r]) ~ 1),
                                        times = 365.5,
                                        extend = TRUE)$surv, 4),
        year_min = min(sub$year[r]), year_max = max(sub$year[r]))
    }
    for (s in SEEDS) for (m in mods) {
      tasks[[length(tasks) + 1]] <- list(
        oc = oc, co = co, mod = m, seed = s, X = MX[[m]], t = tt, e = ee,
        cr = if (oc == "cv") as.integer(sub$event_cr) else NULL,
        subj = sub$subject_id, fold = fmap[[as.character(s)]], hz = HZ,
        nlam = NLAM, lmr = LMR, knots = KNOTS, minu = MINU, ntree = NTREE,
        nsplit = NSPLIT, ntime = NTIME, nodes = NODES, mtryf = MTRYF,
        leaf = LEAF, omtryf = OMTRYF, batch = BATCH, tol = TOL,
        maxext = MAXEXT, spline = (m == "ehr"),
        ckpt = file.path(CKPT, sprintf("%s_%s_%s_%d.rds", oc, C$code, m,
                                       s)))
    }
  }
  META[[co]] <- list(co = co, code = C$code, sub = sub, sets = C$sets,
                     mods = mods, fmap = fmap)
  cat(sprintf("  %-7s development %d / test %d (2014-2016 %d, 2017+ %d)\n",
              co, sum(isdev), sum(!isdev), sum(sub$part == "2014-2016"),
              sum(sub$part == "2017 onwards")))
}
drift <- rbindlist(DR)
small <- function(v) ifelse(v >= 1 & v <= 9, "<10", as.character(v))
drift[, `:=`(events_30 = small(events_30), events_365 = small(events_365))]
fwrite(drift, file.path(RS, paste0(TAG, "_drift.csv")))
cat("\n=== EVENTS BY PERIOD ===\n")
print(drift, row.names = FALSE)
kind <- sapply(tasks, function(t) paste(t$oc, t$co, t$mod))
# minutes per fold for each outcome x cohort x modality: from the saved
# plan, with a probe only for kinds not timed before
kc <- if (file.exists(PLANF)) readRDS(PLANF)$kc else numeric(0)
newk <- setdiff(unique(kind), names(kc))
if (!length(newk)) {
  cat("\n=== RESUMING (every task kind timed before); probe skipped ===\n")
} else {
  cat(sprintf(paste("\n=== TIMING PROBE (%d task kinds, real data, one",
                    "development fold, parallel) ===\n"), length(newk)))
  reps <- tasks[match(newk, kind)]
  names(reps) <- newk
  plan(multisession, workers = min(workers, length(reps)))
  pres <- future_lapply(reps, run_task, folds_to_run = 1L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  for (kd in names(pres)) {
    cat(sprintf("  %-26s %6.1f min per fold\n", kd, pres[[kd]]$sec / 60))
  }
  kc <- c(kc, sapply(pres, function(r) r$sec))
}
costs <- unname(kc[kind]) * (NFOLD + 1)
cat(sprintf("  %d tasks on %d workers: about %.1f hours in all\n",
            length(tasks), workers,
            max(sum(costs) / workers, max(costs)) * 1.1 / 3600))
saveRDS(list(kc = kc, when = format(Sys.time(), "%Y-%m-%d %H:%M")), PLANF)
if (probe_only) {
  cat("--probe-only: stopping here (timings kept, no models saved).\n")
  quit(status = 0)
}
ok_ck <- function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(r) && all(sapply(r$out, function(v) all(is.finite(v))))
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

# ---------------- level-1 checks ----------------
cat("\n=== LEVEL 1 CHECKS ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
nmiss <- 0L
nbad <- 0L
P1 <- list()
for (i in seq_along(tasks)) {
  tk <- tasks[[i]]
  r <- res[[i]]
  for (nm in names(r$out)) {
    nmiss <- nmiss + sum(!is.finite(r$out[[nm]]))
    if (nm == "spline" && !tk$spline) next
    v <- sapply(1:NFOLD, function(k) {
      ii <- which(tk$fold == k)
      cidx(tk$t[ii], tk$e[ii] == 1, r$out[[nm]][ii])
    })
    if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
  }
  P1[[i]] <- cbind(data.table(outcome = tk$oc, cohort = tk$co,
                              modality = tk$mod, seed = tk$seed), r$prm)
}
P1 <- rbindlist(P1, fill = TRUE)
chk("missing level-1 scores", nmiss, if (nmiss == 0) "PASS" else "FAIL")
chk("development scores failing the per-fold rule (D2)", nbad,
    if (nbad == 0) "PASS" else "WARN")
nfo <- sum(P1$rsf_fails) + sum(P1$orsf_fails)
chk("forest fits that failed (best remaining one used)", nfo,
    if (nfo == 0) "PASS" else "WARN")
pe <- P1[modality == "ehr"]
chk("EHR spline Cox fits: k = 5 / smaller / replaced by ridge-Cox",
    sprintf("%d / %d / %d", sum(pe$spline_k == KNOTS),
            sum(pe$spline_k %in% c(3L, 4L)), sum(pe$spline_k == 0L)),
    if (all(pe$spline_k == KNOTS)) "PASS" else "WARN")
chk("forest tunings still at a grid edge (share of fits)", sprintf(
  "%.0f%%", 100 * mean(c(P1$rsf_edge, P1$orsf_edge), na.rm = TRUE)),
  "INFO")
stop_if_fail()

# ---------------- super learners ----------------
cat("\n=== SUPER LEARNERS (development -> test) ===\n")
getr <- function(oc, co, m, s) {
  i <- which(sapply(tasks, function(t) t$oc == oc && t$co == co &&
                      t$mod == m && t$seed == s))
  res[[i]]$out
}
sl_tasks <- list()
for (co in names(META)) {
  M <- META[[co]]
  for (oc in OCS) {
    for (s in SEEDS) {
      build <- function(o2) {
        SS <- setNames(lapply(FAMS, function(z) list()), FAMS)
        for (fam in FAMS) for (m in M$mods) {
          sh <- SHORT[[fam]]
          if (sh == "spline" && m != "ehr") sh <- "cox2"
          SS[[fam]][[m]] <- getr(o2, co, m, s)[[sh]]
        }
        SS
      }
      sub <- M$sub
      sl_tasks[[length(sl_tasks) + 1]] <- list(
        job = paste(oc, co), seed = s, subj = sub$subject_id,
        dev = sub$part == "development", t = sub[[OUTD[[oc]]$t]],
        e = as.integer(OUTD[[oc]]$e(sub)), S = build(oc),
        D = if (oc == "cv") build("death") else NULL, cv = oc == "cv",
        sets = M$sets, hz = HZ, tcuts = TCUTS, tmax = TMAX)
    }
  }
}
plan(multisession, workers = workers)
slr <- future_lapply(sl_tasks, sl_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)

# ---------------- assemble and evaluate ----------------
cat("\n=== ASSEMBLE ===\n")
R <- function(m) paste0(m, ", recalibrated")
INT <- function(oc, code) {
  if (FR == "P" && oc == "death") {
    list(f = file.path(SLD, sprintf("slP_preds_%s.csv", code)),
         m = "Super learner P, recalibrated")
  } else {
    list(f = file.path(SLD, sprintf("optB2_preds_%s_%s_%s.csv", FR, oc,
                                    code)),
         m = "Super learner, recalibrated")
  }
}
SCO <- if (FR == "P") unique(fread(file.path(PP, "scores_P.csv")),
                             by = "hadm_id") else NULL
JOBS <- list()
LONG <- list()
for (co in names(META)) {
  M <- META[[co]]
  sub <- M$sub
  te <- which(sub$part != "development")
  for (oc in OCS) {
    job <- paste(oc, co)
    rv <- Filter(function(r) r$job == job, slr)
    P <- list()
    for (k in names(rv[[1]]$out)) {
      P[[k]] <- Reduce(`+`, lapply(rv, function(r) r$out[[k]])) /
        length(rv)
    }
    bad <- sum(sapply(P, function(m) sum(!is.finite(m))))
    rg <- range(unlist(P))
    chk(sprintf("%s: test risks missing / range", job),
        sprintf("%d / %.5f to %.5f", bad, rg[1], rg[2]),
        if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
    for (k in names(P)) {
      z <- strsplit(k, "|", fixed = TRUE)[[1]]
      LONG[[length(LONG) + 1]] <- data.table(
        hadm_id = rep(sub$hadm_id[te], length(HZ)), outcome = oc,
        cohort = co, model = z[1], set = z[2],
        horizon = rep(HZ, each = length(te)), risk = as.vector(P[[k]]))
    }
    it <- INT(oc, M$code)
    ni <- NA
    if (file.exists(it$f)) {
      V <- fread(it$f)[model == it$m]
      ni <- 0L
      for (sn in names(M$sets)) {
        x <- V[set == sn]
        mm <- sapply(HZ, function(h) {
          y <- x[horizon == h]
          y$risk[match(sub$hadm_id[te], y$hadm_id)]
        })
        ni <- ni + sum(is.na(mm))
        P[[paste("Internal CV model", sn, sep = "|")]] <- mm
      }
    }
    chk(sprintf("%s: main cross-validated predictions (missing)", job),
        if (is.na(ni)) "file not found" else ni,
        if (!is.na(ni) && ni == 0) "PASS" else "WARN")
    scores <- character(0)
    if (!is.null(SCO) && oc == "death") {
      sc <- SCO[match(sub$hadm_id[te], hadm_id)]
      P[["sPESI|score"]] <- matrix(as.numeric(sc$spesi), length(te),
                                   length(HZ))
      P[["PESI-10|score"]] <- matrix(as.numeric(sc$pesi10), length(te),
                                     length(HZ))
      scores <- c("sPESI|score", "PESI-10|score")
    }
    tt <- sub[[OUTD[[oc]]$t]][te]
    ee <- as.integer(OUTD[[oc]]$e(sub))[te]
    Y <- matrix(sapply(HZ, function(h) as.integer(ee == 1 &
                                                    tt <= h + 0.5)),
                ncol = length(HZ))
    DF <- if (oc == "cv") {
      matrix(sapply(HZ, function(h) as.integer(
        sub$event_cr[te] == 2 & sub$time_first[te] <= h + 0.5)),
        ncol = length(HZ))
    } else NULL
    for (pt in c("test, pooled", "2014-2016", "2017 onwards")) {
      rr <- if (pt == "test, pooled") seq_along(te) else
        which(sub$part[te] == pt)
      if (length(rr) < 50) next
      JOBS[[paste(job, pt)]] <- list(
        oc = oc, co = co, period = pt, hz = HZ, scores = scores,
        P = lapply(P, function(m) m[rr, , drop = FALSE]),
        Y = Y[rr, , drop = FALSE],
        DF = if (is.null(DF)) NULL else DF[rr, , drop = FALSE],
        sets = M$sets,
        rows_by = split(seq_along(rr), match(sub$subject_id[te][rr],
                                             unique(sub$subject_id[te][rr]))))
    }
  }
}
stop_if_fail()
fwrite(rbindlist(LONG), file.path(TD, paste0(TAG, "_preds.csv")))
cat(sprintf("\n=== BOOTSTRAP (%d resamples, clustered by patient) ===\n",
            NBOOT))
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
plan(multisession, workers = workers)
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  BL <- future_lapply(seq_len(NBOOT), boot_rep, J = J,
                      future.seed = BSEED, future.packages = PKG)
  B <- do.call(rbind, lapply(BL, function(v) v[names(pt)]))
  colnames(B) <- names(pt)
  cat(sprintf("  %s, %s, %s done\n", J$oc, J$co, J$period))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      outcome = J$oc, cohort = J$co, period = J$period, model = z[1],
      set = z[2], horizon = as.numeric(z[3]), measure = z[4],
      estimate = pt[[k]], lo = q[1], hi = q[2])
  }
  for (k in setdiff(names(J$P), J$scores)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        outcome = J$oc, cohort = J$co, period = J$period, model = z[1],
        set = z[2], horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  pairs <- list()
  for (sn in names(J$sets)) {
    pairs[[length(pairs) + 1]] <- c(R("Super learner"), sn,
                                    "Internal CV model", sn)
    pairs[[length(pairs) + 1]] <- c(R("Super learner"), sn,
                                    R("Super learner-EHR"), "EHR")
    for (sc in J$scores) {
      z <- strsplit(sc, "|", fixed = TRUE)[[1]]
      pairs[[length(pairs) + 1]] <- c(R("Super learner"), sn, z[1], z[2])
    }
  }
  for (cd in pairs) for (h in J$hz) {
    for (m0 in c("auc", "ipa", if (J$oc == "cv") "auc_excl")) {
      ka <- paste(cd[1], cd[2], h, m0, sep = "|")
      kb <- paste(cd[3], cd[4], h, m0, sep = "|")
      if (!(ka %in% colnames(B)) || !(kb %in% colnames(B))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        outcome = J$oc, cohort = J$co, period = J$period, set = cd[2],
        contrast = sprintf("%s - %s", cd[1], cd[3]), horizon = h,
        measure = m0, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
fwrite(mt, file.path(RS, paste0(TAG, "_metrics.csv")))
fwrite(cn, file.path(RS, paste0(TAG, "_contrasts.csv")))
prim <- (mt$cohort == "EHR+ECG" & mt$set == "EHR+ECG") |
  (mt$cohort == "3-mod" & mt$set == "EHR+ECG+CTPA")
smy <- mt[prim & model %in% c(R("Super learner"), "Internal CV model") &
            measure %in% c("auc", "auc_excl", "ipa", "oe", "slope")]
fwrite(smy, file.path(RS, paste0(TAG, "_summary.csv")))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, paste0(TAG, "_checks.csv")))
cat("\n=== TEMPORAL vs INTERNAL (AUROC, test patients, pooled) ===\n")
x <- smy[period == "test, pooled" & measure == "auc"]
print(dcast(x, outcome + cohort + model ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE)
cat("\n=== CALIBRATION IN THE TEST YEARS (observed / expected) ===\n")
x <- smy[period == "test, pooled" & measure == "oe" &
           model == R("Super learner")]
print(dcast(x, outcome + cohort ~ horizon, value.var = "estimate")[
  , lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)],
  row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (f in c("summary", "contrasts", "metrics", "drift", "checks")) {
  cat("  ", file.path(RS, sprintf("%s_%s.csv", TAG, f)), "\n")
}
cat("  ", file.path(RS, paste0(TAG, "_log.txt")),
    " (terminal log)\n")
cat("  ", file.path(TD, paste0(TAG, "_preds.csv")), " (local only)\n")
cat(sprintf(paste("\nOutputs: %s_summary.csv, %s_contrasts.csv,",
                  "%s_drift.csv, %s_checks.csv\n"), TAG, TAG, TAG, TAG))