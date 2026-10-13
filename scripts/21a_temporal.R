# 21a_temporal.R
# Temporal validation of the final pipeline (A24): everything is developed
# on the earlier admission years and tested, untouched, on the later ones.
#
# Frame B, death to 365 days after discharge. MIMIC shifts dates but keeps
# each patient's 3-year anchor window, so an index admission's approximate
# year is the window midpoint (approx_year_hi minus 1, from frame_B.csv or
# followup_365.csv). Development is 2013 or earlier; the test years, 2014
# onwards, are reported pooled and as 2014-2016 and 2017 onwards.
#
# The four families that carry Super learner v4's weight are tuned inside
# the development data only: ridge-Cox retuned (D11 path, deviance) for
# every modality, spline Cox for the EHR (ECG and CTPA use ridge-Cox), and
# Brier-tuned RSF (06c grid, out-of-bag Brier, as 06g) and oblique RSF (06g
# grid and rule) for every modality. Version 4 (D14, D15): a forest tuning
# candidate that fails is skipped, and a spline Cox fit that fails is
# retried with k = 4, then 3, then replaced by the ridge-Cox score for that
# fold. Development level-1 scores are 5-fold out of fold (seed 42). Test
# scores are the mean of the five development-fold models (cross-fitting),
# because a full refit puts the oblique forest on another scale: its
# predicted mortality grows with the number of distinct event times in
# training.
#
# Level 2 is as 18h (Kaplan-Meier, EHR alone per family, PH, TV and TVC
# stacks), with super learner weights from inner 3-fold cross-validation on
# the development scores and per-horizon logistic recalibration, all fitted
# on development only. AUROC, IPA, log-likelihood skill, observed/expected
# and calibration slope at 10, 30, 90, 180 and 365 days are reported for the
# test years pooled and by period, and for the development years by
# cross-validation. A 2,000-sample bootstrap of the test years (seed 42)
# gives contrasts of the recalibrated super learner minus the late-fusion
# analogue (ridge-Cox retuned PH stack), minus EHR alone (spline Cox) and
# minus the uncalibrated super learner. The script stops if fewer than 90%
# of admissions get a period, a score is missing, or the D2 per-fold rule
# fails on development scores.
#
# Inputs:  data/processed/frame_B.csv, feature_whitelist.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv,
#          followup_365.csv
# Outputs: results/temporal_{metrics,contrasts,drift,weights,checks,
#          failures}.csv, results/temporal_log.txt
# Usage:   Rscript scripts/21a_temporal.R [--probe-only] [--workers=N]
# Notes:   finished level-1 tasks are checkpointed and reused on a rerun.

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
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OD <- file.path(PR, "temporal")
CKPT <- file.path(OD, "temporal_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
TCUTS <- c(1, 2, 4, 7, 14, 30, 60, 120, 240)
TMAX <- 366
SPLIT_YEAR <- 2014L
PERIODS <- list("2014-2016" = c(2014, 2016), "2017 onwards" = c(2017, 9999))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
FAMS <- c("ridge-Cox retuned", "spline Cox", "RSF, Brier-tuned",
          "oblique RSF, Brier-tuned")
COH <- list("EHR+ECG" = list(code = "ehrecg",
                             sets = list("EHR+ECG" = c("ehr", "ecg"))),
            "3-mod" = list(code = "3mod",
                           sets = list("EHR+ECG" = c("ehr", "ecg"),
                                       "EHR+CTPA" = c("ehr", "ctpa"),
                                       "EHR+ECG+CTPA" =
                                         c("ehr", "ecg", "ctpa"))))
NFOLD <- 5L
SEED <- 42L
NLAM <- 80L
LMR <- 1e-7
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5L, 15L, 40L, 80L, 160L)
MTRYF <- c(0.5, 1, 3)
LEAF <- c(5L, 20L, 60L)
OMTRYF <- c(0.5, 1, 2, 4)
ONTREE <- 500L
BATCH <- 50L
KNOTS <- 5L
MINU <- 10L
NBOOT <- 2000L
PKG <- c("survival", "glmnet", "mgcv", "data.table", "randomForestSRC",
         "aorsf")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 220)
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
    fwrite(ck, file.path(RS, "temporal_checks.csv"))
    cat("\nSTOP: a check failed. see: temporal_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- level 1: one task = cohort, modality, fold ----------------
# each development-fold model scores its held-out development patients and
# every test-year patient (D15)
l1_task <- function(tk) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  tr <- tk$tr
  te <- tk$te
  tt <- tk$t
  ee <- tk$e
  H <- tk$hz
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
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  out <- list()
  nfail <- 0L
  fails <- list()
  note_fail <- function(fam, cfg, msg) {
    nfail <<- nfail + 1L
    fails[[length(fails) + 1]] <<- data.table::data.table(
      family = fam, setting = cfg, message = substr(msg, 1, 150))
  }
  spline_k <- NA_integer_
  Xs <- prep(tk$X, tr)
  colnames(Xs) <- sprintf("f%d", seq_len(ncol(Xs)))
  # ridge-Cox retuned
  inner <- grouped(tr, tk$seed, 3)
  cv <- cv.glmnet(Xs[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                  family = "cox", alpha = 0, foldid = inner,
                  nlambda = tk$nlam, lambda.min.ratio = tk$lmr,
                  type.measure = "deviance",
                  penalty.factor = rep(1, ncol(Xs)), standardize = FALSE)
  out[["ridge-Cox retuned"]] <- as.vector(predict(
    cv, Xs[te, , drop = FALSE], s = "lambda.min"))
  # spline Cox (EHR only)
  # (D14: a failed fit is retried with smaller bases, then replaced by
  # the ridge-Cox score for that fold)
  if (tk$mod == "ehr") {
    nu <- apply(Xs[tr, , drop = FALSE], 2, function(v) length(unique(v)))
    sm <- colnames(Xs)[nu >= tk$minu]
    li <- setdiff(colnames(Xs), sm)
    for (kk in unique(c(tk$knots, 4L, 3L))) {
      fml <- as.formula(paste("t ~", paste(c(sprintf(
        "s(%s, k = %d, bs = \"cr\")", sm, kk), li), collapse = " + ")))
      g <- tryCatch(mgcv::gam(fml, family = mgcv::cox.ph(), weights = e,
                              data = data.frame(t = tt[tr], e = ee[tr],
                                                Xs[tr, , drop = FALSE]),
                              method = "REML", select = TRUE),
                    error = function(err) conditionMessage(err))
      if (is.character(g)) {
        note_fail("spline Cox", sprintf("k = %d", kk), g)
        next
      }
      out[["spline Cox"]] <- as.vector(predict(
        g, newdata = data.frame(Xs[te, , drop = FALSE]), type = "link"))
      spline_k <- kk
      break
    }
    if (is.null(out[["spline Cox"]])) {
      out[["spline Cox"]] <- out[["ridge-Cox retuned"]]
      spline_k <- 0L
    }
  }
  Ytr <- matrix(sapply(H, function(h) as.numeric(ee[tr] == 1 &
                                                    tt[tr] <= h + 0.5)),
                ncol = length(H))
  # RSF, Brier-tuned (stage 6c data handling, stage 6g rule)
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
  dtr <- data.frame(time = tt[tr], status = ee[tr], X[tr, , drop = FALSE])
  dte <- as.data.frame(X[te, , drop = FALSE])
  p <- ncol(X)
  grid <- expand.grid(nodesize = tk$nodes,
                      mtry = unique(pmin(p, ceiling(tk$mtryf * sqrt(p)))))
  best <- Inf
  for (gi in seq_len(nrow(grid))) {
    f <- tryCatch(rfsrc(Surv(time, status) ~ ., data = dtr,
                        ntree = tk$ntree, nodesize = grid$nodesize[gi],
                        mtry = grid$mtry[gi], nsplit = tk$nsplit,
                        ntime = nt, splitrule = "logrank",
                        save.memory = TRUE,
                        seed = -as.integer(tk$seed * 100 + gi)),
                  error = function(err) conditionMessage(err))
    if (is.character(f)) {
      note_fail("RSF, Brier-tuned", sprintf("node %d, mtry %d",
                                            grid$nodesize[gi],
                                            grid$mtry[gi]), f)
      next
    }
    ro <- sapply(H, function(h) {
      i <- findInterval(h + 0.5, f$time.interest)
      if (i == 0) rep(0, length(tr)) else 1 - f$survival.oob[, i]
    })
    br <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2))
    if (is.finite(br) && br < best) {
      best <- br
      out[["RSF, Brier-tuned"]] <- predict(f, newdata = dte)$predicted
    }
    rm(f)
    invisible(gc(verbose = FALSE))
  }
  # oblique RSF, Brier-tuned (stage 6g)
  ps <- ncol(Xs)
  otr <- data.frame(time = tt[tr], status = ee[tr], Xs[tr, , drop = FALSE])
  ote <- as.data.frame(Xs[te, , drop = FALSE])
  og <- expand.grid(leaf = tk$leaf,
                    mtry = unique(pmax(1L, pmin(ps, ceiling(tk$omtryf *
                                                              sqrt(ps))))))
  nb <- max(1L, tk$ontree %/% tk$batch)
  best <- Inf
  for (gi in seq_len(nrow(og))) {
    ro <- list()
    ms <- list()
    for (bb in seq_len(nb)) {
      f <- tryCatch(aorsf::orsf(
        otr, Surv(time, status) ~ ., n_tree = min(tk$ontree, tk$batch),
        mtry = og$mtry[gi], leaf_min_obs = og$leaf[gi],
        split_min_obs = 2L * og$leaf[gi], oobag_pred_type = "mort",
        tree_seeds = as.integer(tk$seed * 1000 + gi * 20 + bb),
        n_thread = 1), error = function(err) conditionMessage(err))
      if (is.character(f)) {
        note_fail("oblique RSF, Brier-tuned",
                  sprintf("leaf %d, mtry %d, batch %d", og$leaf[gi],
                          og$mtry[gi], bb), f)
        break
      }
      ro[[bb]] <- matrix(predict(f, oobag = TRUE, pred_type = "risk",
                                 pred_horizon = H + 0.5,
                                 boundary_checks = FALSE), ncol = length(H))
      ms[[bb]] <- as.vector(predict(f, new_data = ote, pred_type = "mort"))
      rm(f)
      invisible(gc(verbose = FALSE))
    }
    if (length(ms) < nb) next
    br <- sum(colMeans((Ytr - Reduce(`+`, ro) / nb)^2, na.rm = TRUE))
    if (is.finite(br) && br < best) {
      best <- br
      out[["oblique RSF, Brier-tuned"]] <- rowMeans(do.call(cbind, ms))
    }
  }
  need <- c("ridge-Cox retuned", "RSF, Brier-tuned",
            "oblique RSF, Brier-tuned", if (tk$mod == "ehr") "spline Cox")
  res <- list(co = tk$co, mod = tk$mod, fold = tk$fold, te = te,
              scores = out, nfail = nfail, ok = all(need %in% names(out)),
              fails = data.table::rbindlist(fails), spline_k = spline_k,
              min = (proc.time()[["elapsed"]] - t_task) / 60)
  if (!is.null(tk$ckpt)) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, tk$ckpt)
  }
  res
}

# ---------------- level 2 (as stage 18h) ----------------
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
    Hx <- 0
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
candidates <- function(S, sets, a, b, tt, ee, cuts) {
  Z <- function(fam, ms, rows) {
    matrix(sapply(ms, function(m) nscore(S[[fam]][[m]][rows],
                                         S[[fam]][[m]][a])),
           ncol = length(ms))
  }
  out <- list()
  out[["EHR|Kaplan-Meier"]] <- km_risk(tt[a], ee[a], cuts, length(b))
  for (fam in FAMS) {
    out[[paste0("EHR|EHR alone, ", fam)]] <-
      stack_ph(Z(fam, "ehr", a), Z(fam, "ehr", b), tt[a], ee[a], cuts)
  }
  for (sn in names(sets)) {
    ms <- sets[[sn]]
    for (fam in FAMS) {
      za <- Z(fam, ms, a)
      zb <- Z(fam, ms, b)
      out[[sprintf("%s|%s PH", sn, fam)]] <-
        stack_ph(za, zb, tt[a], ee[a], cuts)
      out[[sprintf("%s|%s TV", sn, fam)]] <-
        stack_tv(za, zb, tt[a], ee[a], cuts)
      out[[sprintf("%s|%s TVC", sn, fam)]] <-
        stack_tvc(za, zb, tt[a], ee[a], cuts, TCUTS, TMAX)
    }
  }
  lapply(out, function(m) matrix(m, nrow = length(b)))
}
metrics <- function(P, Y, idx, hz) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(hz)) {
    y <- Y[idx, hi]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    for (k in names(P)) {
      p <- pmin(pmax(P[[k]][idx, hi], 1e-6), 1 - 1e-6)
      key <- paste(k, hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, p)
      v[paste(key, "ipa", sep = "|")] <- 1 - mean((y - p)^2) /
        (pr * (1 - pr))
      v[paste(key, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      v[paste(key, "oe", sep = "|")] <- pr / mean(p)
    }
  }
  v
}
slope <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  tryCatch(unname(coef(glm(y ~ qlogis(p), family = binomial))[2]),
           error = function(e) NA_real_)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fb <- fread(file.path(PR, "frame_B.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- lapply(MATS, function(f) unique(fread(file.path(IN, f)),
                                        by = "hadm_id"))
if (!("approx_year_hi" %in% names(fb))) {
  fu <- fread(file.path(IN, "followup_365.csv"))
  fb[, approx_year_hi := fu$approx_year_hi[match(hadm_id, fu$hadm_id)]]
}
fb[, year := suppressWarnings(as.numeric(approx_year_hi)) - 1]
na_y <- mean(!is.finite(fb$year))
chk("index admissions with an approximate year",
    sprintf("%.1f%% of %d", 100 * (1 - na_y), nrow(fb)),
    if (na_y <= 0.10) "PASS" else "FAIL")
stop_if_fail()
fb[, period := ifelse(year < SPLIT_YEAR, "development",
                      ifelse(year <= PERIODS[[1]][2], names(PERIODS)[1],
                             names(PERIODS)[2]))]
DRIFT <- list()
tasks <- list()
META <- list()
for (co in names(COH)) {
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  sub <- sub[is.finite(year)]
  tt <- sub$time_death
  ee <- as.integer(sub$event_death == 1)
  dev <- which(sub$period == "development")
  tst <- which(sub$period != "development")
  for (pn in c("development", names(PERIODS))) {
    ii <- which(sub$period == pn)
    if (!length(ii)) next
    km <- summary(survfit(Surv(tt[ii], ee[ii]) ~ 1), times = 365.5,
                  extend = TRUE)
    DRIFT[[length(DRIFT) + 1]] <- data.table(
      cohort = co, period = pn, n = length(ii),
      deaths_10 = sum(ee[ii] == 1 & tt[ii] <= 10.5),
      deaths_365 = sum(ee[ii] == 1 & tt[ii] <= 365.5),
      km_risk_365 = 1 - km$surv, year_min = min(sub$year[ii]),
      year_max = max(sub$year[ii]))
  }
  d365 <- c(sum(ee[dev] == 1 & tt[dev] <= 365.5),
            sum(ee[tst] == 1 & tt[tst] <= 365.5))
  chk(sprintf("%s: development / test patients (deaths by 365 days)", co),
      sprintf("%d (%d) / %d (%d)", length(dev), d365[1], length(tst),
              d365[2]), if (min(d365) >= 100) "PASS" else "WARN")
  set.seed(SEED)
  ud <- unique(sub$subject_id[dev])
  fdev <- sample(rep_len(seq_len(NFOLD), length(ud)))[
    match(sub$subject_id[dev], ud)]
  mods <- unique(unlist(COH[[co]]$sets))
  for (m in mods) {
    cc <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cc])
    storage.mode(Xm) <- "double"
    for (k in seq_len(NFOLD)) {
      tr <- dev[fdev != k]
      te <- c(dev[fdev == k], tst)
      tasks[[length(tasks) + 1]] <- list(
        co = co, mod = m, fold = k, X = Xm, t = tt, e = ee,
        subj = sub$subject_id, tr = tr, te = te, hz = HZ, seed = SEED + k,
        nlam = NLAM, lmr = LMR, knots = KNOTS, minu = MINU, ntree = NTREE,
        nsplit = NSPLIT, ntime = NTIME, nodes = NODES, mtryf = MTRYF,
        leaf = LEAF, omtryf = OMTRYF, ontree = ONTREE, batch = BATCH,
        ckpt = file.path(CKPT, sprintf("%s_%s_cf%d.rds", COH[[co]]$code,
                                       m, k)))
    }
  }
  META[[co]] <- list(co = co, code = COH[[co]]$code, sub = sub, t = tt,
                     e = ee, dev = dev, tst = tst, fdev = fdev, mods = mods,
                     sets = COH[[co]]$sets)
}
drift <- rbindlist(DRIFT)
print(drift[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                     else x)], row.names = FALSE)
stop_if_fail()
ok_ck <- function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(FALSE)
  if (!is.null(r$ok)) return(isTRUE(r$ok))
  r$nfail == 0
}
todo <- which(!sapply(tasks, function(t) ok_ck(t$ckpt)))
cat(sprintf("\n  %d level-1 tasks, %d to fit\n", length(tasks), length(todo)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
if (length(todo)) {
  est <- 0
  for (key in unique(sapply(tasks[todo], function(t) paste(t$co, t$mod)))) {
    tk <- Filter(function(t) paste(t$co, t$mod) == key, tasks)[[1]]
    tr <- tk$tr
    X <- tk$X
    for (j in seq_len(ncol(X))) {
      X[is.na(X[, j]), j] <- median(X[tr, j], na.rm = TRUE)
    }
    X[!is.finite(X)] <- 0
    # the oblique forest refuses columns that are constant in its data
    # (the fits themselves drop them in preprocessing)
    vary <- apply(X[tr, , drop = FALSE], 2, function(v) {
      length(unique(v)) > 1
    })
    X <- X[, vary, drop = FALSE]
    colnames(X) <- sprintf("f%d", seq_len(ncol(X)))
    d <- data.frame(time = tk$t[tr], status = tk$e[tr], X[tr, ])
    t0 <- Sys.time()
    invisible(rfsrc(Surv(time, status) ~ ., data = d, ntree = NTREE,
                    nodesize = 15, mtry = ceiling(sqrt(ncol(X))),
                    nsplit = NSPLIT, ntime = NTIME, save.memory = TRUE,
                    seed = -1L))
    t1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time()
    invisible(aorsf::orsf(d, Surv(time, status) ~ ., n_tree = BATCH,
                          mtry = ceiling(sqrt(ncol(X))), leaf_min_obs = 20,
                          split_min_obs = 40, oobag_pred_type = "mort",
                          tree_seeds = 1L, n_thread = 1))
    t2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n_k <- sum(sapply(tasks[todo], function(t) {
      paste(t$co, t$mod) == key
    }))
    per <- 15 * t1 + 12 * (ONTREE %/% BATCH) * t2 * 1.3 +
      if (tk$mod == "ehr") 30 else 5
    est <- est + per * n_k
    cat(sprintf("  %-14s RSF fit %5.1f s, oblique batch %5.1f s: about %.0f",
                key, t1, t2, per / 60), "min per task\n")
  }
  cat(sprintf("  %d tasks on %d workers: about %.1f hours\n", length(todo),
              workers, est / workers / 3600 * 1.1))
}
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
if (length(todo)) {
  plan(multisession, workers = workers)
  ord <- todo[order(-sapply(tasks[todo], function(t) {
    ncol(t$X) * length(t$tr)
  }))]
  invisible(future_lapply(tasks[ord], l1_task, future.packages = PKG,
                          future.seed = TRUE, future.chunk.size = 1L))
  plan(sequential)
}
res <- lapply(tasks, function(t) readRDS(t$ckpt))
# (D14) failures and substitutions
FL <- rbindlist(lapply(res, function(r) {
  if (is.null(r$fails) || !nrow(r$fails)) return(NULL)
  cbind(data.table(cohort = r$co, modality = r$mod,
                   fold = as.character(r$fold)), r$fails)
}), fill = TRUE)
if (!nrow(FL)) {
  FL <- data.table(cohort = character(0), modality = character(0),
                   fold = character(0), family = character(0),
                   setting = character(0), message = character(0))
}
fwrite(FL, file.path(RS, "temporal_failures.csv"))
nmiss_fam <- sum(!sapply(res, function(r) {
  if (!is.null(r$ok)) isTRUE(r$ok) else r$nfail == 0
}))
chk("level-1 tasks with a family that has no scores", nmiss_fam,
    if (nmiss_fam == 0) "PASS" else "FAIL")
nfc <- nrow(FL[family != "spline Cox"])
chk("forest tuning candidates that failed (best remaining one used)",
    nfc, if (nfc == 0) "PASS" else "WARN")
sk <- sapply(Filter(function(r) r$mod == "ehr", res), function(r) {
  if (is.null(r$spline_k)) KNOTS else r$spline_k
})
chk("spline Cox fits: k = 5 / smaller basis / replaced by ridge-Cox",
    sprintf("%d / %d / %d", sum(sk == KNOTS), sum(sk %in% c(3L, 4L)),
            sum(sk == 0L)), if (all(sk == KNOTS)) "PASS" else "WARN")
if (nrow(FL)) {
  cat("\n  failures (also in temporal_failures.csv):\n")
  print(FL, row.names = FALSE)
}
stop_if_fail()

# ---------------- level 2 and evaluation ----------------
cat("\n=== LEVEL 2 AND EVALUATION ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
MET <- list()
CON <- list()
WTS <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
KEEP <- c("Super learner, recalibrated", "Super learner",
          "SL candidate: ridge-Cox retuned PH",
          "SL candidate: EHR alone, spline Cox",
          "SL candidate: spline Cox TV",
          "SL candidate: oblique RSF, Brier-tuned TV",
          "SL candidate: RSF, Brier-tuned TV",
          "SL candidate: Kaplan-Meier")
for (co in names(META)) {
  M <- META[[co]]
  n <- nrow(M$sub)
  S <- setNames(lapply(FAMS, function(z) list()), FAMS)
  SUM <- S
  CNT <- S
  for (fam in FAMS) for (m in M$mods) {
    S[[fam]][[m]] <- rep(NA_real_, n)
    SUM[[fam]][[m]] <- rep(0, n)
    CNT[[fam]][[m]] <- rep(0L, n)
  }
  # development patients: out-of-fold; test patients: mean of the five
  # fold models (D15)
  for (r in Filter(function(r) r$co == co, res)) {
    isdev <- r$te %in% M$dev
    for (fam in names(r$scores)) {
      v <- r$scores[[fam]]
      S[[fam]][[r$mod]][r$te[isdev]] <- v[isdev]
      it <- r$te[!isdev]
      SUM[[fam]][[r$mod]][it] <- SUM[[fam]][[r$mod]][it] + v[!isdev]
      CNT[[fam]][[r$mod]][it] <- CNT[[fam]][[r$mod]][it] + 1L
    }
  }
  for (fam in FAMS) for (m in M$mods) {
    ok <- CNT[[fam]][[m]] == NFOLD
    S[[fam]][[m]][M$tst] <- ifelse(ok[M$tst], SUM[[fam]][[m]][M$tst] /
                                     NFOLD, NA_real_)
  }
  for (m in setdiff(M$mods, "ehr")) {
    S[["spline Cox"]][[m]] <- S[["ridge-Cox retuned"]][[m]]
  }
  nb <- sum(sapply(S, function(z) sum(sapply(z, function(v) {
    sum(!is.finite(v[c(M$dev, M$tst)]))
  }))))
  chk(sprintf("%s: missing level-1 scores", co), nb,
      if (nb == 0) "PASS" else "FAIL")
  stop_if_fail()
  nbad <- 0L
  for (fam in FAMS) for (m in M$mods) {
    v <- sapply(seq_len(NFOLD), function(k) {
      ii <- M$dev[M$fdev == k]
      cidx(M$t[ii], M$e[ii], S[[fam]][[m]][ii])
    })
    if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
  }
  chk(sprintf("%s: development scores failing the per-fold rule (D2)", co),
      nbad, if (nbad == 0) "PASS" else "FAIL")
  stop_if_fail()
  cuts <- HZ + 0.5
  a_rows <- M$dev
  b_rows <- M$tst
  set.seed(SEED)
  inner <- sample(rep_len(1:3, length(a_rows)))
  PI <- NULL
  for (j in 1:3) {
    cj <- candidates(S, M$sets, a_rows[inner != j], a_rows[inner == j],
                     M$t, M$e, cuts)
    if (is.null(PI)) {
      PI <- lapply(cj, function(z) matrix(NA_real_, length(a_rows),
                                          length(HZ)))
    }
    for (nm in names(cj)) PI[[nm]][inner == j, ] <- cj[[nm]]
  }
  CT <- candidates(S, M$sets, a_rows, b_rows, M$t, M$e, cuts)
  Ydev <- sapply(cuts, function(cu) as.numeric(M$e[a_rows] == 1 &
                                                 M$t[a_rows] <= cu))
  Ytst <- sapply(cuts, function(cu) as.numeric(M$e[b_rows] == 1 &
                                                 M$t[b_rows] <= cu))
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  for (sn in names(M$sets)) {
    cn <- names(CT)[startsWith(names(CT), "EHR|") |
                      startsWith(names(CT), paste0(sn, "|"))]
    X <- do.call(cbind, lapply(cn, function(nm) as.vector(PI[[nm]])))
    a <- sl_weights(X, as.vector(Ydev))
    ens <- Reduce(`+`, Map(function(w, nm) w * CT[[nm]], a, cn))
    ein <- Reduce(`+`, Map(function(w, nm) w * PI[[nm]], a, cn))
    rc <- sapply(seq_along(HZ), function(i) {
      g <- tryCatch(suppressWarnings(glm.fit(
        cbind(1, qlogis(clip(ein[, i]))), Ydev[, i], family = binomial())),
        error = function(e) NULL)
      cf <- if (is.null(g)) c(NA, NA) else g$coefficients
      if (any(!is.finite(cf))) return(ens[, i])
      plogis(cf[1] + cf[2] * qlogis(clip(ens[, i])))
    })
    WTS[[length(WTS) + 1]] <- data.table(
      cohort = co, set = sn, candidate = sub("^[^|]*\\|", "", cn),
      weight = a)
    PT <- list("Super learner, recalibrated" = matrix(rc, ncol = length(HZ)),
               "Super learner" = ens)
    PD <- list("Super learner" = ein)
    for (nm in cn) {
      lab <- paste("SL candidate:", sub("^[^|]*\\|", "", nm))
      PT[[lab]] <- CT[[nm]]
      PD[[lab]] <- PI[[nm]]
    }
    PT <- PT[intersect(KEEP, names(PT))]
    PD <- PD[intersect(KEEP, names(PD))]
    fwrite(rbindlist(lapply(names(PT), function(k) {
      data.table(hadm_id = M$sub$hadm_id[b_rows],
                 period = M$sub$period[b_rows], model = k, set = sn,
                 horizon = rep(HZ, each = length(b_rows)),
                 risk = as.vector(PT[[k]]))
    })), file.path(OD, sprintf("temporal_preds_%s_%s.csv", M$code,
                               gsub("[^A-Za-z]+", "", sn))))
    # internal (development, cross-validated) and temporal (test) metrics
    vi <- metrics(PD, Ydev, seq_along(a_rows), HZ)
    groups <- c(list(test = seq_along(b_rows)),
                lapply(setNames(names(PERIODS), names(PERIODS)),
                       function(pn) which(M$sub$period[b_rows] == pn)))
    for (gn in names(groups)) {
      idx <- groups[[gn]]
      if (length(idx) < 50) next
      pt <- metrics(PT, Ytst, idx, HZ)
      B <- do.call(rbind, lapply(seq_len(NBOOT), function(i) {
        metrics(PT, Ytst, sample(idx, length(idx), replace = TRUE), HZ)
      }))
      for (k in names(pt)) {
        z <- strsplit(k, "|", fixed = TRUE)[[1]]
        q <- q2(B[, k])
        MET[[length(MET) + 1]] <- data.table(
          cohort = co, set = sn, data = if (gn == "test")
            "test years (pooled)" else paste("test", gn),
          model = z[1], horizon = as.numeric(z[2]), measure = z[3],
          estimate = pt[[k]], lo = q[1], hi = q[2])
      }
      for (k in names(PT)) for (hi in seq_along(HZ)) {
        MET[[length(MET) + 1]] <- data.table(
          cohort = co, set = sn, data = if (gn == "test")
            "test years (pooled)" else paste("test", gn), model = k,
          horizon = HZ[hi], measure = "slope",
          estimate = slope(Ytst[idx, hi], PT[[k]][idx, hi]), lo = NA_real_,
          hi = NA_real_)
      }
      if (gn != "test") next
      ref <- list(c("Super learner, recalibrated",
                    "SL candidate: ridge-Cox retuned PH"),
                  c("Super learner, recalibrated",
                    "SL candidate: EHR alone, spline Cox"),
                  c("Super learner, recalibrated", "Super learner"),
                  c("SL candidate: ridge-Cox retuned PH",
                    "SL candidate: EHR alone, spline Cox"))
      for (pp in ref) for (h in HZ) for (ms in c("auc", "ipa", "ll")) {
        ka <- paste(pp[1], h, ms, sep = "|")
        kb <- paste(pp[2], h, ms, sep = "|")
        if (!(ka %in% names(pt)) || !(kb %in% names(pt))) next
        q <- q2(B[, ka] - B[, kb])
        CON[[length(CON) + 1]] <- data.table(
          cohort = co, set = sn, contrast = sprintf("%s - %s", pp[1], pp[2]),
          horizon = h, measure = ms, estimate = pt[[ka]] - pt[[kb]],
          lo = q[1], hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
      }
    }
    for (k in names(vi)) {
      z <- strsplit(k, "|", fixed = TRUE)[[1]]
      MET[[length(MET) + 1]] <- data.table(
        cohort = co, set = sn, data = "development (internal, CV)",
        model = z[1], horizon = as.numeric(z[2]), measure = z[3],
        estimate = vi[[k]], lo = NA_real_, hi = NA_real_)
    }
    cat(sprintf("  %s %s done\n", co, sn))
  }
}
mt <- rbindlist(MET)
cn <- rbindlist(CON)
wt <- rbindlist(WTS)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "temporal_metrics.csv"))
fwrite(cn, file.path(RS, "temporal_contrasts.csv"))
fwrite(drift, file.path(RS, "temporal_drift.csv"))
fwrite(wt, file.path(RS, "temporal_weights.csv"))
fwrite(ck, file.path(RS, "temporal_checks.csv"))
cat("\n=== AUROC: internal (development, CV) vs test years ===\n")
x <- mt[measure == "auc" & model %in% c("Super learner, recalibrated",
                                        "Super learner",
                                        "SL candidate: ridge-Cox retuned PH")]
print(dcast(x, cohort + set + model + data ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE)
cat("\n=== CONTRASTS, test years pooled, AUROC ===\n")
print(cn[measure == "auc", .(cohort, set, contrast, horizon,
                             est = round(estimate, 4), lo = round(lo, 4),
                             hi = round(hi, 4), excludes_zero)],
      row.names = FALSE, nrows = 500)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (temporal_preds_*.csv and temporal_ckpt/, local only)\n")
for (f in c("temporal_metrics.csv", "temporal_contrasts.csv",
            "temporal_drift.csv", "temporal_weights.csv",
            "temporal_checks.csv", "temporal_failures.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "temporal_log.txt"), " (terminal log)\n")
cat("\nOutputs: temporal_metrics.csv, temporal_contrasts.csv,",
    "temporal_drift.csv, temporal_weights.csv, temporal_checks.csv,",
    "temporal_failures.csv\n")