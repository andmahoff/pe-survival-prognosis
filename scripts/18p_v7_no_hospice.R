# 18p_v7_no_hospice.R
# Super learner v7nh: how much the hospice destination indicator adds to
# the at-discharge model (PREREGISTRATION_ADDENDUM_2.md A31).
#
# Builds on 18o, in Frame B only. Four families are added, each labelled
# "(EHR+ no hospice)": ridge-Cox retuned, spline Cox, RSF Brier-tuned and
# oblique RSF Brier-tuned. Their EHR score comes from the no-hospice run
# of 6j (EHR+ without the hospice indicator); their ECG and CTPA scores are
# the same family's as in v4. Library v7nh is v6 plus these four families,
# fitted, weighted and recalibrated exactly as v7 (v6 plus the EHR+
# families); every other library is as in 18o. Evaluation is as in 18o
# (2,000-sample bootstrap, seed 42): on every patient, and with hospice
# discharges removed.
#
# Contrasts (all recalibrated): v7 minus v7nh (what the hospice indicator
# adds); v7nh minus v6 (what EHR+ adds without it); v7nh minus Late PH;
# each EHR+ candidate minus its no-hospice counterpart (EHR alone and PH
# stacks). The analysis is descriptive: no model is replaced, and v7 stays
# the final at-discharge model. sl7nh_summary.csv holds the AUROC and IPA
# of the first two contrasts in the primary cohort (Frame B, EHR+ECG).
#
# Checks: inputs found and aligned; 18o's super learners v2 to v7 and
# A16's candidates reproduced (difference at most 1e-9); no missing risks;
# risks in [0, 1]. A failed check stops the script. Tasks run in the
# workers, as in A14 and A16, so the random inner folds match.
#
# Inputs:  under data/processed/: frame_A.csv, frame_B.csv, folds.csv,
#          level1/oof_{linear,forest,extra,proper,enet,ehr79,ehrplus}_*.csv,
#          level1/oof_ehrplus_nh_B_*.csv, spline/spline_ckpt/*.rds,
#          super_learner/sl_preds_*.csv, sl2_preds_*.csv, sl7_preds_*.csv,
#          fusion_arch/arch_preds_linear_*.csv, blockopt_oof_*.csv;
#          data/inputs/followup_365.csv
# Outputs: results/sl7nh_metrics.csv, results/sl7nh_contrasts.csv,
#          results/sl7nh_weights.csv, results/sl7nh_summary.csv,
#          results/sl7nh_checks.csv, results/sl7nh_log.txt;
#          data/processed/super_learner/sl7nh_preds_B_*.csv (patient
#          level, not in the repo)
# Usage:   Rscript scripts/18p_v7_no_hospice.R [--probe-only]
#          [--workers=N]   (--probe-only stops after the timing probe)

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
SPD <- file.path(PR, "spline", "spline_ckpt")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
TCUTS <- list(A = c(1, 2, 4, 7, 14), B = c(1, 2, 4, 7, 14, 30, 60, 120, 240))
TMAX <- c(A = 31, B = 366)
LIB2 <- c("ridge-Cox", "ridge-Cox retuned", "RSF", "oblique RSF",
          "CoxBoost")
LIB3 <- c(LIB2, "spline Cox")
FAMS <- c(LIB3, "RSF, Brier-tuned", "oblique RSF, wide grid",
          "oblique RSF, Brier-tuned", "elastic-net Cox", "CoxBoost 500")
BFN <- "Block forest, optimised"
FAMS79 <- c("ridge-Cox retuned (EHR 79)", "spline Cox (EHR 79)",
            "RSF, Brier-tuned (EHR 79)", "oblique RSF, Brier-tuned (EHR 79)")
FAMSP <- c("ridge-Cox retuned (EHR+)", "spline Cox (EHR+)",
           "RSF, Brier-tuned (EHR+)", "oblique RSF, Brier-tuned (EHR+)")
FAMSNH <- c("ridge-Cox retuned (EHR+ no hospice)",
            "spline Cox (EHR+ no hospice)",
            "RSF, Brier-tuned (EHR+ no hospice)",
            "oblique RSF, Brier-tuned (EHR+ no hospice)")
FAMS_ALL <- c(FAMS, FAMS79, FAMSP, FAMSNH)
LIBS <- list(v2 = LIB2, v3 = LIB3, v4 = FAMS, v5 = c(FAMS, BFN),
             v6 = c(FAMS, BFN, FAMS79), v7 = c(FAMS, BFN, FAMS79, FAMSP),
             v7nh = c(FAMS, BFN, FAMS79, FAMSNH))
PH_LIBS <- character(0)
NBOOT <- 2000L
SEED <- 42L
PKG <- c("survival", "data.table")
options(width = 220)
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
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
    fwrite(ck, file.path(RS, "sl7nh_checks.csv"))
    cat("\nSTOP: a check failed. see: sl7nh_checks.csv\n")
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
candidates <- function(S, sets, a, b, tt, ee, cuts, tcuts, tmax,
                       BF = NULL) {
  Z <- function(fam, ms, rows) {
    matrix(sapply(ms, function(m) nscore(S[[fam]][[m]][rows],
                                         S[[fam]][[m]][a])),
           ncol = length(ms))
  }
  out <- list()
  out[["EHR|Kaplan-Meier"]] <- km_risk(tt[a], ee[a], cuts, length(b))
  for (fam in FAMS_ALL) {
    out[[paste0("EHR|EHR alone, ", fam)]] <-
      stack_ph(Z(fam, "ehr", a), Z(fam, "ehr", b), tt[a], ee[a], cuts)
  }
  for (sn in names(sets)) {
    ms <- sets[[sn]]
    for (fam in FAMS_ALL) {
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
    # A27: the optimised block forest's own out-of-fold risks
    if (!is.null(BF[[sn]])) {
      out[[sprintf("%s|%s", sn, BFN)]] <- BF[[sn]][b, , drop = FALSE]
    }
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
                       tk$tmax, tk$BF)
      if (is.null(PI)) {
        PI <- lapply(cj, function(m) matrix(NA_real_, length(tr), length(H)))
      }
      for (nm in names(cj)) PI[[nm]][inner == j, ] <- cj[[nm]]
    }
    CT <- candidates(tk$S, sets, tr, te, tk$t, tk$e, cuts, tk$tcuts,
                     tk$tmax, tk$BF)
    Y <- matrix(sapply(cuts, function(cu) {
      as.numeric(tk$e[tr] == 1 & tk$t[tr] <= cu)
    }), ncol = length(H))
    for (nm in names(CT)) {
      z <- strsplit(nm, "|", fixed = TRUE)[[1]]
      put(paste("SL candidate:", z[2]), z[1], te, CT[[nm]])
    }
    fam_of <- function(lab) {
      f <- sub("^EHR alone, ", "", lab)
      sub(" (PH|TV|TVC)$", "", f)
    }
    clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
    for (sn in names(sets)) {
      ca <- names(CT)[startsWith(names(CT), "EHR|") |
                        startsWith(names(CT), paste0(sn, "|"))]
      lab <- sub("^[^|]*\\|", "", ca)
      fm <- fam_of(lab)
      for (lv in names(LIBS)) {
        cn <- ca[fm %in% c(LIBS[[lv]], "Kaplan-Meier", "both families")]
        X <- do.call(cbind, lapply(cn, function(nm) as.vector(PI[[nm]])))
        a <- sl_weights(X, as.vector(Y))
        loss <- colMeans((X - as.vector(Y))^2)
        ens <- Reduce(`+`, Map(function(w, nm) w * CT[[nm]], a, cn))
        arm <- sprintf("Super learner %s", lv)
        put(arm, sn, te, ens)
        if (lv != "v2") {
          ein <- Reduce(`+`, Map(function(w, nm) w * PI[[nm]], a, cn))
          rc <- sapply(seq_along(H), function(i) {
            g <- tryCatch(suppressWarnings(glm.fit(
              cbind(1, qlogis(clip(ein[, i]))), Y[, i],
              family = binomial())), error = function(e) NULL)
            cf <- if (is.null(g)) c(NA, NA) else g$coefficients
            if (any(!is.finite(cf))) return(ens[, i])
            plogis(cf[1] + cf[2] * qlogis(clip(ens[, i])))
          })
          put(paste0(arm, ", recalibrated"), sn, te,
              matrix(rc, ncol = length(H)))
        }
        if (lv == "v4") {
          put("Super learner v4, best single", sn, te,
              CT[[cn[which.min(loss)]]])
        }
        if (lv != "v2") {
          wts[[length(wts) + 1]] <- data.table(
            library = lv, weights = "pooled", horizon = NA_real_,
            set = sn, fold = k,
            candidate = sub("^[^|]*\\|", "", cn), weight = a,
            best = seq_along(cn) == which.min(loss), inner_loss = loss)
        }
        # A28: weights fitted separately at each horizon
        if (lv %in% PH_LIBS) {
          eh <- matrix(NA_real_, length(te), length(H))
          ih <- matrix(NA_real_, length(tr), length(H))
          for (i in seq_along(H)) {
            Xi <- sapply(cn, function(nm) PI[[nm]][, i])
            ai <- sl_weights(Xi, Y[, i])
            eh[, i] <- as.vector(sapply(cn, function(nm) CT[[nm]][, i]) %*%
                                   ai)
            ih[, i] <- as.vector(Xi %*% ai)
            wts[[length(wts) + 1]] <- data.table(
              library = lv, weights = "per-horizon", horizon = H[i],
              set = sn, fold = k, candidate = sub("^[^|]*\\|", "", cn),
              weight = ai, best = NA, inner_loss = NA_real_)
          }
          arh <- sprintf("Super learner %s, per-horizon weights", lv)
          put(arh, sn, te, eh)
          rh <- sapply(seq_along(H), function(i) {
            g <- tryCatch(suppressWarnings(glm.fit(
              cbind(1, qlogis(clip(ih[, i]))), Y[, i],
              family = binomial())), error = function(e) NULL)
            cf <- if (is.null(g)) c(NA, NA) else g$coefficients
            if (any(!is.finite(cf))) return(eh[, i])
            plogis(cf[1] + cf[2] * qlogis(clip(eh[, i])))
          })
          put(paste0(arh, ", recalibrated"), sn, te,
              matrix(rh, ncol = length(H)))
        }
      }
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
FU <- unique(fread(file.path(P2, "data", "inputs", "followup_365.csv")),
             by = "hadm_id")
folds <- fread(file.path(PR, "folds.csv"))
tasks <- list()
META <- list()
for (fr in "B") for (co in names(COH)) {
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
  fpp <- file.path(L1, sprintf("oof_%s_%s_%s.csv", c("proper", "enet"),
                               fr, code))
  chk(sprintf("%s: stage 6g and 6h files found", key),
      paste(basename(fpp), collapse = ", "),
      if (all(file.exists(fpp))) "PASS" else "FAIL")
  stop_if_fail()
  GT <- fread(fpp[1])
  ET <- fread(fpp[2])
  fbo <- file.path(FA, sprintf("blockopt_oof_%s_%s.csv", fr, code))
  chk(sprintf("%s: stage 18l file found", key), basename(fbo),
      if (file.exists(fbo)) "PASS" else "FAIL")
  stop_if_fail()
  BT <- fread(fbo)
  f79 <- file.path(L1, sprintf("oof_ehr79_%s_%s.csv", fr, code))
  chk(sprintf("%s: stage 6i file found", key), basename(f79),
      if (file.exists(f79)) "PASS" else "FAIL")
  stop_if_fail()
  T79 <- fread(f79)
  fpl <- file.path(L1, sprintf("oof_ehrplus_%s_%s.csv", fr, code))
  chk(sprintf("%s: stage 6j file found", key), basename(fpl),
      if (file.exists(fpl)) "PASS" else "FAIL")
  stop_if_fail()
  TPL <- fread(fpl)
  fnh <- file.path(L1, sprintf("oof_ehrplus_nh_%s_%s.csv", fr, code))
  chk(sprintf("%s: stage 6j no-hospice file found", key), basename(fnh),
      if (file.exists(fnh)) "PASS" else "FAIL")
  stop_if_fail()
  TNH <- fread(fnh)
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(LT$seed))
  ok <- identical(seeds, sort(unique(FT$seed))) &&
    identical(seeds, sort(unique(XT$seed))) &&
    identical(seeds, sort(unique(fo$seed))) &&
    nrow(LT) == n * length(seeds) && nrow(FT) == n * length(seeds) &&
    nrow(XT) == n * length(seeds) && nrow(GT) == n * length(seeds) &&
    nrow(ET) == n * length(seeds) &&
    identical(seeds, sort(unique(GT$seed))) &&
    identical(seeds, sort(unique(ET$seed)))
  chk(sprintf("%s: seeds and rows aligned (6a-6h, folds)", key),
      sprintf("%d seeds, %d x %d rows", length(seeds), n, length(seeds)),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  for (s in seeds) {
    a <- LT[LT$seed == s][match(sub$hadm_id, hadm_id)]
    b <- FT[FT$seed == s][match(sub$hadm_id, hadm_id)]
    x <- XT[XT$seed == s][match(sub$hadm_id, hadm_id)]
    gg <- GT[GT$seed == s][match(sub$hadm_id, hadm_id)]
    en <- ET[ET$seed == s][match(sub$hadm_id, hadm_id)]
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    S <- setNames(lapply(FAMS_ALL, function(z) list()), FAMS_ALL)
    for (m in mods) {
      S[["ridge-Cox"]][[m]] <- a[[sprintf("%s.cox.death.lp", m)]]
      S[["ridge-Cox retuned"]][[m]] <- x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF"]][[m]] <- b[[sprintf("%s.rsf.death.score", m)]]
      S[["oblique RSF"]][[m]] <- x[[sprintf("%s.orsf.death.score", m)]]
      S[["CoxBoost"]][[m]] <- x[[sprintf("%s.cboost.death.lp", m)]]
      S[["spline Cox"]][[m]] <- x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF, Brier-tuned"]][[m]] <- gg[[sprintf("%s.rsfb.death.score", m)]]
      S[["oblique RSF, wide grid"]][[m]] <-
        gg[[sprintf("%s.orsfc.death.score", m)]]
      S[["oblique RSF, Brier-tuned"]][[m]] <-
        gg[[sprintf("%s.orsfb.death.score", m)]]
      S[["elastic-net Cox"]][[m]] <- en[[sprintf("%s.enet.death.lp", m)]]
      S[["CoxBoost 500"]][[m]] <- en[[sprintf("%s.cb500.death.lp", m)]]
      # A28: the EHR 79 families keep their own ECG and CTPA scores
      S[["ridge-Cox retuned (EHR 79)"]][[m]] <-
        x[[sprintf("%s.cox2.death.lp", m)]]
      S[["spline Cox (EHR 79)"]][[m]] <- x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF, Brier-tuned (EHR 79)"]][[m]] <-
        gg[[sprintf("%s.rsfb.death.score", m)]]
      S[["oblique RSF, Brier-tuned (EHR 79)"]][[m]] <-
        gg[[sprintf("%s.orsfb.death.score", m)]]
      # A30: the EHR+ families keep their own ECG and CTPA scores
      S[["ridge-Cox retuned (EHR+)"]][[m]] <-
        x[[sprintf("%s.cox2.death.lp", m)]]
      S[["spline Cox (EHR+)"]][[m]] <- x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF, Brier-tuned (EHR+)"]][[m]] <-
        gg[[sprintf("%s.rsfb.death.score", m)]]
      S[["oblique RSF, Brier-tuned (EHR+)"]][[m]] <-
        gg[[sprintf("%s.orsfb.death.score", m)]]
      # A31: the no-hospice families keep their own ECG and CTPA scores
      S[["ridge-Cox retuned (EHR+ no hospice)"]][[m]] <-
        x[[sprintf("%s.cox2.death.lp", m)]]
      S[["spline Cox (EHR+ no hospice)"]][[m]] <-
        x[[sprintf("%s.cox2.death.lp", m)]]
      S[["RSF, Brier-tuned (EHR+ no hospice)"]][[m]] <-
        gg[[sprintf("%s.rsfb.death.score", m)]]
      S[["oblique RSF, Brier-tuned (EHR+ no hospice)"]][[m]] <-
        gg[[sprintf("%s.orsfb.death.score", m)]]
    }
    y79 <- T79[T79$seed == s][match(sub$hadm_id, hadm_id)]
    if (any(y79$fold != f, na.rm = TRUE) || anyNA(y79$fold)) {
      chk(sprintf("%s seed %s: stage 6i folds match", key, s), "no", "FAIL")
      stop_if_fail()
    }
    S[["ridge-Cox retuned (EHR 79)"]][["ehr"]] <- y79$ehr79.cox2.death.lp
    S[["spline Cox (EHR 79)"]][["ehr"]] <- y79$ehr79.spline.death.lp
    S[["RSF, Brier-tuned (EHR 79)"]][["ehr"]] <- y79$ehr79.rsfb.death.score
    S[["oblique RSF, Brier-tuned (EHR 79)"]][["ehr"]] <-
      y79$ehr79.orsfb.death.score
    ypl <- TPL[TPL$seed == s][match(sub$hadm_id, hadm_id)]
    if (any(ypl$fold != f, na.rm = TRUE) || anyNA(ypl$fold)) {
      chk(sprintf("%s seed %s: stage 6j folds match", key, s), "no", "FAIL")
      stop_if_fail()
    }
    S[["ridge-Cox retuned (EHR+)"]][["ehr"]] <- ypl$ehrp.cox2.death.lp
    S[["spline Cox (EHR+)"]][["ehr"]] <- ypl$ehrp.spline.death.lp
    S[["RSF, Brier-tuned (EHR+)"]][["ehr"]] <- ypl$ehrp.rsfb.death.score
    S[["oblique RSF, Brier-tuned (EHR+)"]][["ehr"]] <-
      ypl$ehrp.orsfb.death.score
    ynh <- TNH[TNH$seed == s][match(sub$hadm_id, hadm_id)]
    if (any(ynh$fold != f, na.rm = TRUE) || anyNA(ynh$fold)) {
      chk(sprintf("%s seed %s: no-hospice folds match", key, s), "no",
          "FAIL")
      stop_if_fail()
    }
    S[["ridge-Cox retuned (EHR+ no hospice)"]][["ehr"]] <-
      ynh$ehrp.cox2.death.lp
    S[["spline Cox (EHR+ no hospice)"]][["ehr"]] <- ynh$ehrp.spline.death.lp
    S[["RSF, Brier-tuned (EHR+ no hospice)"]][["ehr"]] <-
      ynh$ehrp.rsfb.death.score
    S[["oblique RSF, Brier-tuned (EHR+ no hospice)"]][["ehr"]] <-
      ynh$ehrp.orsfb.death.score
    spf <- file.path(SPD, sprintf("%s_%s_%d.rds", fr, code, s))
    sp <- if (file.exists(spf)) readRDS(spf)$lp else NULL
    if (is.null(sp) || length(sp) != n || any(!is.finite(sp))) {
      chk(sprintf("%s seed %s: spline Cox scores (stage 19b)", key, s),
          basename(spf), "FAIL")
      stop_if_fail()
    }
    S[["spline Cox"]][["ehr"]] <- sp
    miss <- sum(is.na(f)) + sum(sapply(S, function(z)
      sum(sapply(z, function(v) sum(!is.finite(v))))))
    if (miss > 0) {
      chk(sprintf("%s seed %s: missing scores or folds", key, s), miss,
          "FAIL")
      stop_if_fail()
    }
    BF <- list()
    for (sn in names(COH[[co]]$sets)) {
      x <- BT[BT$seed == s & BT$set == sn][match(sub$hadm_id, hadm_id)]
      BF[[sn]] <- matrix(sapply(hz, function(h) x[[sprintf("r%d", h)]]),
                         ncol = length(hz))
      bad <- sum(!is.finite(BF[[sn]])) + sum(x$fold != f, na.rm = TRUE) +
        sum(is.na(x$fold))
      if (bad > 0) {
        chk(sprintf("%s seed %s %s: block forest risks and folds", key, s,
                    sn), bad, "FAIL")
        stop_if_fail()
      }
    }
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = f, subj = sub$subject_id,
      t = sub$time_death, e = as.integer(sub$event_death == 1), S = S,
      sets = COH[[co]]$sets, hz = hz, tcuts = TCUTS[[fr]],
      tmax = TMAX[[fr]], BF = BF)
  }
  # layer-1 ridge-Cox risks, stage 6a and stage 6d, seed-averaged
  r6 <- LT[, lapply(.SD, mean), by = hadm_id, .SDcols = !c("seed")][
    match(sub$hadm_id, hadm_id)]
  rx <- XT[, lapply(.SD, mean), by = hadm_id,
           .SDcols = as.vector(outer(mods, hz, function(m, h)
             sprintf("%s.cox2.death.r%d", m, h)))][
               match(sub$hadm_id, hadm_id)]
  avg_cols <- function(D, pat) {
    cc <- as.vector(outer(mods, hz, function(m, h) sprintf(pat, m, h)))
    D[, lapply(.SD, mean), by = hadm_id, .SDcols = cc][
      match(sub$hadm_id, hadm_id)]
  }
  ga <- lapply(c("rsfc", "rsfb", "orsfc9", "orsfc", "orsfb"), function(cg) {
    avg_cols(GT, paste0("%s.", cg, ".death.r%d"))
  })
  names(ga) <- c("rsfc", "rsfb", "orsfc9", "orsfc", "orsfb")
  ea <- avg_cols(ET, "%s.enet.death.r%d")
  L1N <- list("Layer 1 elastic-net Cox" = list(ea, "enet"),
              "Layer 1 RSF wide, C-tuned" = list(ga$rsfc, "rsfc"),
              "Layer 1 RSF wide, Brier-tuned" = list(ga$rsfb, "rsfb"),
              "Layer 1 oblique RSF, C-tuned" = list(ga$orsfc9, "orsfc9"),
              "Layer 1 oblique RSF wide, C-tuned" = list(ga$orsfc, "orsfc"),
              "Layer 1 oblique RSF wide, Brier-tuned" =
                list(ga$orsfb, "orsfb"))
  L1R <- list()
  for (m in mods) {
    L1R[[paste("Layer 1 ridge-Cox", toupper(m), sep = "|")]] <-
      sapply(hz, function(h) r6[[sprintf("%s.cox.death.r%d", m, h)]])
    L1R[[paste("Layer 1 ridge-Cox retuned", toupper(m), sep = "|")]] <-
      sapply(hz, function(h) rx[[sprintf("%s.cox2.death.r%d", m, h)]])
    for (nm in names(L1N)) {
      D <- L1N[[nm]][[1]]
      tg <- L1N[[nm]][[2]]
      L1R[[paste(nm, toupper(m), sep = "|")]] <-
        sapply(hz, function(h) D[[sprintf("%s.%s.death.r%d", m, tg, h)]])
    }
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
  fwrite(LONG, file.path(OD, sprintf("sl7nh_preds_%s_%s.csv", M$fr,
                                     M$code)))
  # stage 18o's super learners reproduced exactly
  s4 <- file.path(OD, sprintf("sl7_preds_%s_%s.csv", M$fr, M$code))
  if (!file.exists(s4)) {
    chk(sprintf("%s: stage 18o predictions", key), basename(s4), "FAIL")
  } else {
    S4 <- fread(s4)
    d4 <- 0
    n4 <- 0L
    for (md in c("Super learner v2", "Super learner v3",
                 "Super learner v3, recalibrated", "Super learner v4",
                 "Super learner v4, recalibrated",
                 "Super learner v4, best single", "Super learner v5",
                 "Super learner v5, recalibrated", "Super learner v6",
                 "Super learner v6, recalibrated", "Super learner v7",
                 "Super learner v7, recalibrated")) {
      for (sn in names(M$sets)) {
        k <- paste(md, sn, sep = "|")
        x <- S4[model == md & set == sn]
        if (is.null(P[[k]]) || !nrow(x)) next
        for (hi in seq_along(hz)) {
          y <- x[horizon == hz[hi]]
          d4 <- max(d4, abs(y$risk[match(M$sub$hadm_id, y$hadm_id)] -
                              P[[k]][, hi]))
          n4 <- n4 + 1L
        }
      }
    }
    chk(sprintf("%s: stage 18o super learners reproduced (n, max |diff|)",
                key), sprintf("%d, %.1e", n4, d4),
        if (n4 > 0 && d4 <= 1e-9) "PASS" else "FAIL")
  }
  # A16's candidates reproduced; A16's and A14's super learners as
  # references
  sf <- file.path(OD, sprintf("sl2_preds_%s_%s.csv", M$fr, M$code))
  s1 <- file.path(OD, sprintf("sl_preds_%s_%s.csv", M$fr, M$code))
  if (!file.exists(sf) || !file.exists(s1)) {
    chk(sprintf("%s: A16 and A14 predictions", key),
        paste(basename(c(sf, s1)), collapse = ", "), "FAIL")
  } else {
    S2 <- fread(sf)
    dm <- 0
    nm <- 0L
    for (md in grep("^SL candidate:", unique(S2$model), value = TRUE)) {
      for (sn in unique(S2[model == md, set])) {
        k <- paste(md, sn, sep = "|")
        if (is.null(P[[k]])) next
        for (hi in seq_along(hz)) {
          x <- S2[model == md & set == sn & horizon == hz[hi]]
          dm <- max(dm, abs(x$risk[match(M$sub$hadm_id, x$hadm_id)] -
                              P[[k]][, hi]))
          nm <- nm + 1L
        }
      }
    }
    chk(sprintf("%s: A16 candidates reproduced (n, max |diff|)", key),
        sprintf("%d, %.1e", nm, dm),
        if (nm > 0 && dm <= 1e-9) "PASS" else "FAIL")
    dv <- 0
    for (sn in names(M$sets)) {
      x <- S2[model == "Super learner v2" & set == sn]
      for (hi in seq_along(hz)) {
        y <- x[horizon == hz[hi]]
        dv <- max(dv, abs(y$risk[match(M$sub$hadm_id, y$hadm_id)] -
                            P[[paste("Super learner v2", sn,
                                     sep = "|")]][, hi]))
      }
    }
    chk(sprintf("%s: Super learner v2 reproduced (max |diff|)", key),
        sprintf("%.1e", dv), if (dv <= 1e-9) "PASS" else "FAIL")
    S1 <- fread(s1)
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
                      mods = unique(unlist(M$sets)), evaluation = "all")
  # the same predictions, evaluated without hospice discharges
  dl <- toupper(trimws(as.character(
    FU$discharge_location[match(M$sub$hadm_id, FU$hadm_id)])))
  hsp <- !is.na(dl) & dl == "HOSPICE"
  kr <- which(!hsp)
  chk(sprintf("%s: hospice discharges removed in the second evaluation",
              key), sprintf("%d of %d", sum(hsp), length(hsp)), "INFO")
  JOBS[[paste(key, "no hospice")]] <- list(
    fr = M$fr, co = M$co, hz = hz,
    P = lapply(P, function(m) m[kr, , drop = FALSE]),
    Y = matrix(Y, ncol = length(hz))[kr, , drop = FALSE],
    rows_by = split(seq_along(kr), match(M$sub$subject_id[kr],
                                         unique(M$sub$subject_id[kr]))),
    sets = M$sets, mods = unique(unlist(M$sets)),
    evaluation = "hospice excluded")
}
stop_if_fail()
W <- rbindlist(lapply(res, function(r) r$wts))
wsum <- W[, .(folds = .N, mean_weight = mean(weight),
              share_nonzero = mean(weight > 1e-6),
              share_best_single = mean(best)),
          by = .(library, weights, horizon, frame, cohort, set, candidate)]
setorder(wsum, library, weights, horizon, frame, cohort, set, -mean_weight)
fwrite(wsum, file.path(RS, "sl7nh_weights.csv"))

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples) ===\n", NBOOT))
CONS <- function(sn) {
  C <- function(a, b) c(a, b, sn, sn)
  K <- function(x) paste("SL candidate:", x)
  R <- function(lv) sprintf("Super learner %s, recalibrated", lv)
  out <- list(C(R("v7"), R("v7nh")), C(R("v7nh"), R("v6")),
              C(R("v7nh"), "Late PH"))
  for (f in FAMSNH) {
    fp <- sub(" no hospice\\)$", ")", f)
    out[[length(out) + 1]] <- C(K(paste0("EHR alone, ", fp)),
                                K(paste0("EHR alone, ", f)))
    out[[length(out) + 1]] <- C(K(paste(fp, "PH")), K(paste(f, "PH")))
  }
  out
}
L1PAIRS <- list()
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cds <- unlist(lapply(names(J$sets), CONS), recursive = FALSE)
  cds <- lapply(cds, function(cd) {
    if (grepl("^SL candidate: EHR alone|^SL candidate: Kaplan", cd[1])) {
      cd[3] <- "EHR"
    }
    if (grepl("^SL candidate: EHR alone|^SL candidate: Kaplan", cd[2])) {
      cd[4] <- "EHR"
    }
    cd
  })
  cds <- cds[!duplicated(sapply(cds, paste, collapse = "|"))]
  for (m in J$mods) for (pp in L1PAIRS) {
    cds[[length(cds) + 1]] <- c(pp[1], pp[2], toupper(m), toupper(m))
  }
  need <- unique(c(unlist(lapply(cds, function(cd) {
    c(paste(cd[1], cd[3], sep = "|"), paste(cd[2], cd[4], sep = "|"))
  })), grep("^(Super learner|Late)", names(J$P), value = TRUE)))
  JB <- J
  JB$P <- J$P[intersect(names(J$P), need)]
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = JB,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  Frame %s %s (%s) done\n", J$fr, J$co, J$evaluation))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- if (k %in% colnames(B)) q2(B[, k]) else c(NA_real_, NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, evaluation = J$evaluation,
      model = z[1], set = z[2],
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
        frame = J$fr, cohort = J$co, evaluation = J$evaluation,
        model = z[1], set = z[2],
        horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  for (cd in cds) {
    for (h in J$hz) for (ms in c("auc", "ipa", "ll")) {
      ka <- paste(cd[1], cd[3], h, ms, sep = "|")
      kb <- paste(cd[2], cd[4], h, ms, sep = "|")
      if (!(ka %in% colnames(B)) || !(kb %in% colnames(B))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, evaluation = J$evaluation,
        set = cd[3],
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
fwrite(mt, file.path(RS, "sl7nh_metrics.csv"))
fwrite(cn, file.path(RS, "sl7nh_contrasts.csv"))

# ---------------- summary: what the hospice indicator adds ----------------
keyc <- c(paste("Super learner v7, recalibrated -",
                "Super learner v7nh, recalibrated"),
          paste("Super learner v7nh, recalibrated -",
                "Super learner v6, recalibrated"))
smy <- cn[frame == "B" & cohort == "EHR+ECG" & set == "EHR+ECG" &
            measure %in% c("auc", "ipa") & contrast %in% keyc,
          .(evaluation, contrast, measure, horizon, estimate, lo, hi,
            excludes_zero)]
fwrite(smy, file.path(RS, "sl7nh_summary.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "sl7nh_checks.csv"))

cat("\n=== MEAN WEIGHTS ABOVE 0.05, v7nh ===\n")
print(wsum[library == "v7nh" & mean_weight > 0.05,
           .(cohort, set, candidate, mean_weight = round(mean_weight, 3))],
      row.names = FALSE, nrows = 500)
cat("\n=== WHAT THE HOSPICE INDICATOR ADDS (primary cohort) ===\n")
print(smy[, .(evaluation, contrast = sub(", recalibrated", "", contrast),
              measure, horizon, est = round(estimate, 4),
              lo = round(lo, 4), hi = round(hi, 4), excludes_zero)],
      row.names = FALSE, nrows = 200)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (sl7nh_preds_B_*.csv, patient-level, local only)\n")
for (f in c("sl7nh_metrics.csv", "sl7nh_contrasts.csv",
            "sl7nh_weights.csv", "sl7nh_summary.csv", "sl7nh_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "sl7nh_log.txt"), " (terminal log)\n")
cat("\nOutputs: sl7nh_summary.csv, sl7nh_contrasts.csv, sl7nh_metrics.csv,",
    "sl7nh_checks.csv\n")