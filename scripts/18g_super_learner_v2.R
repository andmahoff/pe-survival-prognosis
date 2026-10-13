# 18g_super_learner_v2.R
# Super learner with a larger library, and the effect of retuning the
# level-1 ridge-Cox (PREREGISTRATION_ADDENDUM_2.md, A16).
#
# Extends the A14 super learner of 18f_super_learner.R, whose eleven
# candidates are a subset of this library. Layer 1 is the out-of-fold death
# scores per seed and outer fold, as rank-based normal scores on the rows
# each stack is fitted to (as S4, 11c_eval_s4.R), from five families:
# ridge-Cox (06a_level1_linear.py; Harrell's C over five penalties),
# ridge-Cox retuned (06d_level1_extra.R; D11 path, deviance tuning), RSF
# (06b_level1_forests.R), and oblique RSF (aorsf) and CoxBoost (both
# 06d_level1_extra.R). The candidates are Kaplan-Meier; the EHR score alone
# for each family; for each set and family, the Cox stack (PH), the
# per-horizon logistic stacks (TV) and the log-time-varying Cox stack of A13
# (TVC); and A14's Cox stack on ridge-Cox and RSF scores together.
#
# The weights are fitted as in 18f (non-negative, summing to one, minimising
# the squared error of inner out-of-fold risks over the horizons). The
# models are Super learner v2 (one weight vector) and Super learner v2,
# best single.
#
# Evaluation: AUROC, IPA, log-likelihood skill score, observed/expected and
# calibration slope (point estimate), with 2,000-sample bootstrap intervals
# (seed 42) clustered by patient in Frame A. References are A14's super
# learner and 18a's Late PH and Late TV. Contrasts: Super learner v2 minus
# A14's super learner, Late PH and its own best single; ridge-Cox retuned PH
# minus ridge-Cox PH; oblique RSF minus RSF (PH and TV); CoxBoost PH minus
# ridge-Cox retuned PH. The layer-1 comparison (A15) is each modality's
# ridge-Cox risks, retuned minus 06a.
#
# Checks:  inputs found and aligned; A14's eleven candidates reproduced to
#          1e-9; no missing risks; risks in [0, 1] (all FAIL).
# Notes:   a timing probe runs the largest task and pauses 5 s. Tasks then
#          run in parallel as in 18f, so the inner folds match A14's.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          level1/oof_linear_*.csv, oof_forest_*.csv, oof_extra_*.csv,
#          super_learner/sl_preds_*.csv, fusion_arch/arch_preds_linear_*.csv
# Outputs: data/processed/super_learner/sl2_preds_*.csv (patient level, not
#          in the repo); results/sl2_metrics.csv, results/sl2_contrasts.csv,
#          results/sl2_weights.csv, results/sl2_checks.csv,
#          results/sl2_log.txt
# Usage:   Rscript scripts/18g_super_learner_v2.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
L1 <- file.path(PR, "level1")
FA <- file.path(PR, "fusion_arch")
RS <- file.path(P2, "results")
OD <- file.path(PR, "super_learner")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
TCUTS <- list(A = c(1, 2, 4, 7, 14), B = c(1, 2, 4, 7, 14, 30, 60, 120, 240))
TMAX <- c(A = 31, B = 366)
FAMS <- c("ridge-Cox", "ridge-Cox retuned", "RSF", "oblique RSF",
          "CoxBoost")
NBOOT <- 2000L
SEED <- 42L
PKG <- c("survival", "data.table")
options(width = 220)
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-64s %s\n", status, item, value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "sl2_checks.csv"))
    cat("\nSTOP: a check failed. see: sl2_checks.csv\n")
    quit(status = 1)
  }
}
COH <- list("EHR+ECG" = list(code = "ehrecg",
                             sets = list("EHR+ECG" = c("ehr", "ecg"))),
            "3-mod" = list(code = "3mod",
                           sets = list("EHR+ECG" = c("ehr", "ecg"),
                                       "EHR+CTPA" = c("ehr", "ctpa"),
                                       "EHR+ECG+CTPA" =
                                         c("ehr", "ecg", "ctpa"))))

# ---------------- stacks (sent to workers) ----------------
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
# every candidate, fitted on rows a and predicting rows b
candidates <- function(S, sets, a, b, tt, ee, cuts, tcuts, tmax) {
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
        stack_tvc(za, zb, tt[a], ee[a], cuts, tcuts, tmax)
    }
    out[[sprintf("%s|both families PH", sn)]] <-
      stack_ph(cbind(Z("ridge-Cox", ms, a), Z("RSF", ms, a)),
               cbind(Z("ridge-Cox", ms, b), Z("RSF", ms, b)),
               tt[a], ee[a], cuts)
  }
  lapply(out, function(m) matrix(m, nrow = length(b)))
}

# ---------------- one task = frame x cohort x seed ----------------
sl_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  cuts <- H + 0.5
  sets <- tk$sets
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  out <- list()
  wts <- list()
  put <- function(model, set, idx, M) {
    k <- paste(model, set, sep = "|")
    if (is.null(out[[k]])) out[[k]] <<- matrix(NA_real_, n, length(H))
    out[[k]][idx, ] <<- M
  }
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    PI <- NULL
    for (j in 1:3) {
      a <- tr[inner != j]
      b <- tr[inner == j]
      cj <- candidates(tk$S, sets, a, b, tk$t, tk$e, cuts, tk$tcuts,
                       tk$tmax)
      if (is.null(PI)) {
        PI <- lapply(cj, function(m) matrix(NA_real_, length(tr), length(H)))
      }
      for (nm in names(cj)) PI[[nm]][inner == j, ] <- cj[[nm]]
    }
    CT <- candidates(tk$S, sets, tr, te, tk$t, tk$e, cuts, tk$tcuts,
                     tk$tmax)
    Y <- matrix(sapply(cuts, function(cu) {
      as.numeric(tk$e[tr] == 1 & tk$t[tr] <= cu)
    }), ncol = length(H))
    for (nm in names(CT)) {
      z <- strsplit(nm, "|", fixed = TRUE)[[1]]
      put(paste("SL candidate:", z[2]), z[1], te, CT[[nm]])
    }
    for (sn in names(sets)) {
      cn <- names(CT)[startsWith(names(CT), "EHR|") |
                        startsWith(names(CT), paste0(sn, "|"))]
      X <- do.call(cbind, lapply(cn, function(nm) as.vector(PI[[nm]])))
      a <- sl_weights(X, as.vector(Y))
      loss <- colMeans((X - as.vector(Y))^2)
      best <- which.min(loss)
      put("Super learner v2", sn, te,
          Reduce(`+`, Map(function(w, nm) w * CT[[nm]], a, cn)))
      put("Super learner v2, best single", sn, te, CT[[cn[best]]])
      wts[[length(wts) + 1]] <- data.table(
        set = sn, fold = k, candidate = sub("^[^|]*\\|", "", cn),
        weight = a, best = seq_along(cn) == best, inner_loss = loss)
    }
  }
  list(fr = tk$fr, co = tk$co, seed = tk$seed, out = out,
       wts = rbindlist(wts)[, `:=`(frame = tk$fr, cohort = tk$co,
                                   seed = tk$seed)])
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
    for (k in names(J$P)) {
      p <- pmin(pmax(J$P[[k]][idx, hi], 1e-6), 1 - 1e-6)
      key <- paste(k, J$hz[hi], sep = "|")
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
boot_rep <- function(i, J) {
  if (J$fr == "A") {
    idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                       length(J$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(nrow(J$Y), nrow(J$Y), replace = TRUE)
  }
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in names(COH)) {
  code <- COH[[co]]$code
  key <- paste(fr, co)
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  n <- nrow(sub)
  hz <- HZ[[fr]]
  mods <- unique(unlist(COH[[co]]$sets))
  ff <- file.path(L1, sprintf("oof_%s_%s_%s.csv",
                              c("linear", "forest", "extra"), fr, code))
  chk(sprintf("%s: stage 6a, 6b and 6d files found", key),
      paste(basename(ff), collapse = ", "),
      if (all(file.exists(ff))) "PASS" else "FAIL")
  stop_if_fail()
  LT <- fread(ff[1], select = c("seed", "hadm_id",
                                sprintf("%s.cox.death.lp", mods),
                                as.vector(outer(mods, hz, function(m, h)
                                  sprintf("%s.cox.death.r%d", m, h)))))
  FT <- fread(ff[2], select = c("seed", "hadm_id",
                                sprintf("%s.rsf.death.score", mods)))
  XT <- fread(ff[3])
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(LT$seed))
  ok <- identical(seeds, sort(unique(FT$seed))) &&
    identical(seeds, sort(unique(XT$seed))) &&
    identical(seeds, sort(unique(fo$seed))) &&
    nrow(LT) == n * length(seeds) && nrow(FT) == n * length(seeds) &&
    nrow(XT) == n * length(seeds)
  chk(sprintf("%s: seeds and rows aligned (6a, 6b, 6d, folds)", key),
      sprintf("%d seeds, %d x %d rows", length(seeds), n, length(seeds)),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  for (s in seeds) {
    a <- LT[LT$seed == s][match(sub$hadm_id, hadm_id)]
    b <- FT[FT$seed == s][match(sub$hadm_id, hadm_id)]
    x <- XT[XT$seed == s][match(sub$hadm_id, hadm_id)]
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    S <- setNames(lapply(FAMS, function(z) list()), FAMS)
    for (m in mods) {
      S[["ridge-Cox"]][[m]] <- a[[sprintf("%s.cox.death.lp", m)]]
      S[["ridge-Cox retuned"]][[m]] <- x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF"]][[m]] <- b[[sprintf("%s.rsf.death.score", m)]]
      S[["oblique RSF"]][[m]] <- x[[sprintf("%s.orsf.death.score", m)]]
      S[["CoxBoost"]][[m]] <- x[[sprintf("%s.cboost.death.lp", m)]]
    }
    miss <- sum(is.na(f)) + sum(sapply(S, function(z)
      sum(sapply(z, function(v) sum(!is.finite(v))))))
    if (miss > 0) {
      chk(sprintf("%s seed %s: missing scores or folds", key, s), miss,
          "FAIL")
      stop_if_fail()
    }
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = f, subj = sub$subject_id,
      t = sub$time_death, e = as.integer(sub$event_death == 1), S = S,
      sets = COH[[co]]$sets, hz = hz, tcuts = TCUTS[[fr]],
      tmax = TMAX[[fr]])
  }
  # layer-1 ridge-Cox risks, stage 6a and stage 6d, seed-averaged
  r6 <- LT[, lapply(.SD, mean), by = hadm_id, .SDcols = !c("seed")][
    match(sub$hadm_id, hadm_id)]
  rx <- XT[, lapply(.SD, mean), by = hadm_id,
           .SDcols = as.vector(outer(mods, hz, function(m, h)
             sprintf("%s.cox2.death.r%d", m, h)))][
               match(sub$hadm_id, hadm_id)]
  L1R <- list()
  for (m in mods) {
    L1R[[paste("Layer 1 ridge-Cox", toupper(m), sep = "|")]] <-
      sapply(hz, function(h) r6[[sprintf("%s.cox.death.r%d", m, h)]])
    L1R[[paste("Layer 1 ridge-Cox retuned", toupper(m), sep = "|")]] <-
      sapply(hz, function(h) rx[[sprintf("%s.cox2.death.r%d", m, h)]])
  }
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub,
                      sets = COH[[co]]$sets, L1R = L1R)
}
cat(sprintf("  %d tasks (frame x cohort x seed)\n", length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) length(t$t) *
                          length(t$sets) * length(t$hz)))
t0 <- Sys.time()
invisible(sl_task(tasks[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  largest task %.0f s; %d tasks on %d workers: at most",
                  "about %.1f min\n"), tp, length(tasks), workers,
            tp * ceiling(length(tasks) / workers) * 1.2 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, sl_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE ===\n")
JOBS <- list()
A14 <- c("Kaplan-Meier", "EHR alone, ridge-Cox", "EHR alone, RSF",
         "ridge-Cox PH", "ridge-Cox TV", "ridge-Cox TVC", "RSF PH",
         "RSF TV", "RSF TVC", "both families PH")
for (key in names(META)) {
  M <- META[[key]]
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  P <- list()
  for (k in names(rr[[1]]$out)) {
    P[[k]] <- Reduce(`+`, lapply(rr, function(r) r$out[[k]])) / length(rr)
  }
  bad <- sum(sapply(P, function(m) sum(!is.finite(m))))
  rg <- range(unlist(P))
  chk(sprintf("%s: missing risks / range", key),
      sprintf("%d / %.5f to %.5f", bad, rg[1], rg[2]),
      if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
  hz <- HZ[[M$fr]]
  LONG <- rbindlist(lapply(names(P), function(k) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = rep(M$sub$hadm_id, length(hz)), model = z[1],
               set = z[2], horizon = rep(hz, each = nrow(M$sub)),
               risk = as.vector(P[[k]]))
  }))
  fwrite(LONG, file.path(OD, sprintf("sl2_preds_%s_%s.csv", M$fr,
                                     M$code)))
  # A14's candidates reproduced, and A14's super learner as a reference
  sf <- file.path(OD, sprintf("sl_preds_%s_%s.csv", M$fr, M$code))
  if (!file.exists(sf)) {
    chk(sprintf("%s: A14 predictions (%s)", key, basename(sf)), "missing",
        "FAIL")
  } else {
    S1 <- fread(sf)
    dm <- 0
    nm <- 0L
    for (cand in A14) {
      md <- paste("SL candidate:", cand)
      for (sn in unique(S1[model == md, set])) {
        k <- paste(md, sn, sep = "|")
        if (is.null(P[[k]])) next
        for (hi in seq_along(hz)) {
          x <- S1[model == md & set == sn & horizon == hz[hi]]
          dm <- max(dm, abs(x$risk[match(M$sub$hadm_id, x$hadm_id)] -
                              P[[k]][, hi]))
          nm <- nm + 1L
        }
      }
    }
    chk(sprintf("%s: A14 candidates reproduced (n, max |diff|)", key),
        sprintf("%d, %.1e", nm, dm),
        if (nm > 0 && dm <= 1e-9) "PASS" else "FAIL")
    for (sn in names(M$sets)) {
      x <- S1[model == "Super learner" & set == sn]
      P[[paste("Super learner (A14)", sn, sep = "|")]] <-
        matrix(sapply(hz, function(h) {
          y <- x[horizon == h]
          y$risk[match(M$sub$hadm_id, y$hadm_id)]
        }), ncol = length(hz))
    }
  }
  f18 <- file.path(FA, sprintf("arch_preds_linear_%s_%s.csv", M$fr,
                               M$code))
  if (file.exists(f18)) {
    R <- fread(f18)[model %in% c("Late PH", "Late TV")]
    for (mn in c("Late PH", "Late TV")) for (sn in names(M$sets)) {
      x <- R[model == mn & set == sn]
      if (!nrow(x)) next
      P[[paste(mn, sn, sep = "|")]] <- matrix(sapply(hz, function(h) {
        y <- x[horizon == h]
        y$risk[match(M$sub$hadm_id, y$hadm_id)]
      }), ncol = length(hz))
    }
  } else {
    chk(sprintf("%s: stage 18a references", key), "missing", "WARN")
  }
  for (k in names(M$L1R)) P[[k]] <- matrix(M$L1R[[k]], ncol = length(hz))
  Y <- sapply(hz, function(h) as.integer(M$sub$event_death == 1 &
                                           M$sub$time_death <= h + 0.5))
  JOBS[[key]] <- list(fr = M$fr, co = M$co, hz = hz, P = P,
                      Y = matrix(Y, ncol = length(hz)),
                      rows_by = split(seq_len(nrow(M$sub)),
                                      match(M$sub$subject_id,
                                            unique(M$sub$subject_id))),
                      sets = M$sets,
                      mods = unique(unlist(M$sets)))
}
stop_if_fail()
W <- rbindlist(lapply(res, function(r) r$wts))
wsum <- W[, .(folds = .N, mean_weight = mean(weight),
              share_nonzero = mean(weight > 1e-6),
              share_best_single = mean(best)),
          by = .(frame, cohort, set, candidate)]
setorder(wsum, frame, cohort, set, -mean_weight)
fwrite(wsum, file.path(RS, "sl2_weights.csv"))

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples) ===\n", NBOOT))
CONS <- function(sn) {
  list(c("Super learner v2", "Super learner (A14)", sn, sn),
       c("Super learner v2", "Late PH", sn, sn),
       c("Super learner v2", "Super learner v2, best single", sn, sn),
       c("SL candidate: ridge-Cox retuned PH", "SL candidate: ridge-Cox PH",
         sn, sn),
       c("SL candidate: oblique RSF PH", "SL candidate: RSF PH", sn, sn),
       c("SL candidate: oblique RSF TV", "SL candidate: RSF TV", sn, sn),
       c("SL candidate: CoxBoost PH", "SL candidate: ridge-Cox retuned PH",
         sn, sn))
}
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  Frame %s %s done\n", J$fr, J$co))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in names(J$P)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, model = z[1], set = z[2],
        horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  cds <- unlist(lapply(names(J$sets), CONS), recursive = FALSE)
  for (m in J$mods) {
    cds[[length(cds) + 1]] <- c("Layer 1 ridge-Cox retuned",
                                "Layer 1 ridge-Cox", toupper(m),
                                toupper(m))
  }
  for (cd in cds) {
    for (h in J$hz) for (ms in c("auc", "ipa", "ll")) {
      ka <- paste(cd[1], cd[3], h, ms, sep = "|")
      kb <- paste(cd[2], cd[4], h, ms, sep = "|")
      if (!(ka %in% names(pt)) || !(kb %in% names(pt))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, set = cd[3],
        contrast = sprintf("%s - %s", cd[1], cd[2]), horizon = h,
        measure = ms, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "sl2_metrics.csv"))
fwrite(cn, file.path(RS, "sl2_contrasts.csv"))
fwrite(ck, file.path(RS, "sl2_checks.csv"))

cat("\n=== MEAN WEIGHTS ABOVE 0.05, Frame B, EHR+ECG cohort ===\n")
print(wsum[frame == "B" & cohort == "EHR+ECG" & mean_weight > 0.05,
           .(candidate, mean_weight = round(mean_weight, 3),
             share_best_single = round(share_best_single, 2))],
      row.names = FALSE)
cat("\n=== CONTRASTS, AUROC ===\n")
print(cn[measure == "auc", .(frame, cohort, set, contrast, horizon,
                             est = round(estimate, 4), lo = round(lo, 4),
                             hi = round(hi, 4), excludes_zero)],
      row.names = FALSE, nrows = 2000)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (sl2_preds_*.csv, patient-level, local only)\n")
for (f in c("sl2_metrics.csv", "sl2_contrasts.csv", "sl2_weights.csv",
            "sl2_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "sl2_log.txt"), " (terminal log)\n")
cat("\nOutputs: sl2_metrics.csv, sl2_contrasts.csv, sl2_weights.csv,",
    "sl2_checks.csv\n")