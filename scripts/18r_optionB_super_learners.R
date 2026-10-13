# 18r_optionB_super_learners.R
# Super learners for the composite and CV-readmission outcomes in both
# pipelines (A33, A34), with the v7-core check on death.
#
# Level 1 uses the out-of-fold scores of 06q (composite, CV readmission and
# EHR refits without the hospice indicator) and the death scores of 06p,
# 06j, 06d and 06g, in four families: ridge-Cox retuned, spline Cox (EHR
# only; ECG and CTPA use ridge-Cox), RSF and oblique RSF (Brier-tuned).
#
# Composite candidates are Kaplan-Meier, EHR alone and the PH, TV and TVC
# stacks per family and set (as 18q). CV readmission, with death competing,
# uses per-horizon logistic stacks on "CV readmission first by the horizon"
# (death first = 0) with each modality's CV and death scores as inputs. The
# super learner (SL) pools weights over all candidates of a set and SL-EHR
# over the EHR-alone ones; both are recalibrated per horizon.
#
# A34 adds component-based composite candidates (inputs: each modality's
# death and CV scores) to form SL+, which replaces SL for the composite if
# it is significantly better on AUROC at two or more horizons and worse at
# none in the primary cohort (EHR+ECG). It also rebuilds the discharge
# composite and CV models from no-hospice EHR scores, contrasts them with
# the full models and re-evaluates every discharge super learner without
# hospice discharges. Version 3 (A36): at-diagnosis inputs rebuilt, so the
# check that SL and SL-EHR reproduce version 1 (to 1e-9) covers discharge
# only.
#
# v7-core check (A33): the discharge death super learner over the four EHR+
# families is representative of Super learner v7 (18o) if it is within 0.005
# AUROC of v7 at every horizon in the primary cohort. Metrics are AUROC,
# IPA, log-likelihood skill, observed/expected and calibration slope at 10,
# 30, 90, 180 and 365 days, with a 2,000-sample patient-clustered bootstrap
# (seed 42). The CV secondary analysis (as in the dissertation) is AUROC
# among patients who did not die first (auc_excl). Contrasts: SL minus
# SL-EHR, SL+ minus SL, SL minus SL without hospice, v7 minus v7-core. The
# script stops on missing inputs, misaligned folds or risks outside [0, 1].
#
# Inputs:  data/processed/frame_B.csv, presentation/frame_P_outcomes.csv,
#          level1/oof_*.csv, super_learner/{sl7,optB}_preds_B_*.csv (all
#          under data/processed/), data/inputs/followup_365.csv
# Outputs: results/optB2_sl_{summary,contrasts,metrics,weights,checks}.csv,
#          results/optB2_{decision,core_check}.csv, optB2_sl_log.txt
# Usage:   Rscript scripts/18r_optionB_super_learners.R [--probe-only]
#          [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
OD <- file.path(PR, "super_learner")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
TCUTS <- c(1, 2, 4, 7, 14, 30, 60, 120, 240)
TMAX <- 366
CORE_TOL <- 0.005
FAMS <- c("ridge-Cox retuned", "spline Cox", "RSF, Brier-tuned",
          "oblique RSF, Brier-tuned")
SETS <- list("EHR+ECG" = list(code = "ehrecg",
                              sets = list("EHR+ECG" = c("ehr", "ecg"))),
             "3-mod" = list(code = "3mod",
                            sets = list("EHR+ECG" = c("ehr", "ecg"),
                                        "EHR+CTPA" = c("ehr", "ctpa"),
                                        "EHR+ECG+CTPA" =
                                          c("ehr", "ecg", "ctpa"))))
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
    fwrite(ck, file.path(RS, "optB2_sl_checks.csv"))
    cat("\nSTOP: a check failed. see: optB2_sl_checks.csv\n")
    quit(status = 1)
  }
}

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
    # a vector from the start, so a cut before the first event gives zero
    # risk for every patient (stage 18o started from a scalar 0, which is
    # numerically identical whenever events precede the cut)
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
# every candidate, fitted on rows a and predicting rows b. Survival
# outcomes (composite, death): as stage 18q, plus (A34, composite only)
# component-based stacks on the death and CV scores. CV: per-horizon
# logistic stacks on CV and death scores.
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
    # A34: the composite from its components' scores
    if (!is.null(tk$C1)) {
      for (fam in FAMS) {
        za <- cbind(Z(tk$C1, fam, "ehr", a), Z(tk$C2, fam, "ehr", a))
        zb <- cbind(Z(tk$C1, fam, "ehr", b), Z(tk$C2, fam, "ehr", b))
        out[[paste0("EHR|EHR alone, ", fam, " (components)")]] <-
          stack_ph(za, zb, tt[a], ee[a], cuts)
      }
      for (sn in names(sets)) {
        ms <- sets[[sn]]
        for (fam in FAMS) {
          za <- cbind(Z(tk$C1, fam, ms, a), Z(tk$C2, fam, ms, a))
          zb <- cbind(Z(tk$C1, fam, ms, b), Z(tk$C2, fam, ms, b))
          out[[sprintf("%s|%s PH (components)", sn, fam)]] <-
            stack_ph(za, zb, tt[a], ee[a], cuts)
          out[[sprintf("%s|%s TV (components)", sn, fam)]] <-
            stack_tv(za, zb, tt[a], ee[a], cuts)
        }
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

# ---------------- one task = frame x outcome x cohort x seed ----------------
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
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  fit_sl <- function(cn, PI, CT, Y, arm, sn, k, te) {
    X <- do.call(cbind, lapply(cn, function(nm) as.vector(PI[[nm]])))
    a <- sl_weights(X, as.vector(Y))
    ens <- Reduce(`+`, Map(function(w, nm) w * CT[[nm]], a, cn))
    ein <- Reduce(`+`, Map(function(w, nm) w * PI[[nm]], a, cn))
    put(arm, sn, te, ens)
    rc <- sapply(seq_along(H), function(i) {
      g <- tryCatch(suppressWarnings(glm.fit(
        cbind(1, qlogis(clip(ein[, i]))), Y[, i],
        family = binomial())), error = function(e) NULL)
      cf <- if (is.null(g)) c(NA, NA) else g$coefficients
      if (any(!is.finite(cf))) return(ens[, i])
      plogis(cf[1] + cf[2] * qlogis(clip(ens[, i])))
    })
    put(paste0(arm, ", recalibrated"), sn, te, matrix(rc, ncol = length(H)))
    wts[[length(wts) + 1]] <<- data.table(
      model = arm, set = sn, fold = k,
      candidate = sub("^[^|]*\\|", "", cn), weight = a)
  }
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    PI <- NULL
    for (j in 1:3) {
      a <- tr[inner != j]
      b <- tr[inner == j]
      cj <- candidates(tk, a, b, cuts)
      if (is.null(PI)) {
        PI <- lapply(cj, function(m) matrix(NA_real_, length(tr), length(H)))
      }
      for (nm in names(cj)) PI[[nm]][inner == j, ] <- cj[[nm]]
    }
    CT <- candidates(tk, tr, te, cuts)
    Y <- matrix(sapply(cuts, function(cu) {
      as.numeric(tk$e[tr] == 1 & tk$t[tr] <= cu)
    }), ncol = length(H))
    for (nm in names(CT)) {
      z <- strsplit(nm, "|", fixed = TRUE)[[1]]
      put(paste("SL candidate:", z[2]), z[1], te, CT[[nm]])
    }
    comp <- grepl("(components)", names(CT), fixed = TRUE)
    ce <- names(CT)[startsWith(names(CT), "EHR|") & !comp]
    fit_sl(ce, PI, CT, Y, "Super learner-EHR", "EHR", k, te)
    for (sn in names(sets)) {
      cs <- names(CT)[startsWith(names(CT), paste0(sn, "|"))]
      ca <- c(ce, setdiff(cs, names(CT)[comp]))
      fit_sl(ca, PI, CT, Y, "Super learner", sn, k, te)
      if (any(comp)) {
        cp <- c(ca, names(CT)[comp & (startsWith(names(CT), "EHR|") |
                                        startsWith(names(CT),
                                                   paste0(sn, "|")))])
        fit_sl(cp, PI, CT, Y, "Super learner+", sn, k, te)
      }
    }
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
      p <- pmin(pmax(raw, 1e-6), 1 - 1e-6)
      v[paste(key, "auc", sep = "|")] <- auc(y, raw)
      v[paste(key, "ipa", sep = "|")] <- 1 - mean((y - p)^2) /
        (pr * (1 - pr))
      v[paste(key, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      v[paste(key, "oe", sep = "|")] <- pr / mean(p)
      # CV secondary: among patients not dying first by this horizon
      if (!is.null(keep)) {
        v[paste(key, "auc_excl", sep = "|")] <- auc(y[keep], raw[keep])
      }
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
cat("=== BUILD ===\n")
LJOBS <- list(c("P", "comp"), c("P", "cv"), c("B", "comp"), c("B", "cv"),
              c("B", "death"), c("B", "comp_nh"), c("B", "cv_nh"))
lf <- function(fr, key, code) {
  file.path(L1, switch(key,
                       comp = sprintf("oof_%scomp_%s.csv", fr, code),
                       cv = sprintf("oof_%scv_%s.csv", fr, code),
                       comp_nh = sprintf("oof_Bcomp_nh_%s.csv", code),
                       cv_nh = sprintf("oof_Bcv_nh_%s.csv", code),
                       Pdeath = sprintf("oof_P_%s.csv", code),
                       ehrplusB = sprintf("oof_ehrplus_B_%s.csv", code),
                       ehrplus_nhB = sprintf("oof_ehrplus_nh_B_%s.csv",
                                             code),
                       extraB = sprintf("oof_extra_B_%s.csv", code),
                       properB = sprintf("oof_proper_B_%s.csv", code)))
}
# where each family's score lives: (file key, column). kind: comp, cv,
# death, and at discharge comp_nh, cv_nh, death_nh (the EHR without the
# hospice indicator; ECG and CTPA as the full versions)
colfor <- function(fr, kind, m, fam) {
  short <- c("ridge-Cox retuned" = "cox2", "spline Cox" = "spline",
             "RSF, Brier-tuned" = "rsfb",
             "oblique RSF, Brier-tuned" = "orsfb")[[fam]]
  if (short == "spline" && m != "ehr") short <- "cox2"
  suf <- if (short %in% c("cox2", "spline")) "lp" else "score"
  nh <- grepl("_nh$", kind)
  base <- sub("_nh$", "", kind)
  if (base %in% c("comp", "cv")) {
    key <- if (nh && m == "ehr") kind else base
    return(c(key, sprintf("%s.%s.%s.%s", m, short, base, suf)))
  }
  if (fr == "P") return(c("Pdeath", sprintf("%s.%s.death.%s", m, short,
                                            suf)))
  if (m == "ehr") {
    return(c(if (nh) "ehrplus_nhB" else "ehrplusB",
             sprintf("ehrp.%s.death.%s", short, suf)))
  }
  c(if (short %in% c("cox2", "spline")) "extraB" else "properB",
    sprintf("%s.%s.death.%s", m, if (short == "spline") "cox2" else short,
            suf))
}
FRM <- list(P = file.path(PP, "frame_P_outcomes.csv"),
            B = file.path(PR, "frame_B.csv"))
keys_for <- function(fr, kind) {
  nh <- grepl("_nh$", kind)
  base <- sub("_nh$", "", kind)
  dk <- if (fr == "P") "Pdeath" else c(if (nh) "ehrplus_nhB" else
    "ehrplusB", "extraB", "properB")
  switch(base,
         death = dk,
         comp = unique(c(if (nh) "comp_nh", "comp", if (nh) "cv_nh", "cv",
                         dk)),
         cv = unique(c(if (nh) "cv_nh", "cv", dk)))
}
need <- c(unlist(FRM), file.path(IN, "followup_365.csv"))
for (jb in LJOBS) for (co in names(SETS)) {
  code <- SETS[[co]]$code
  need <- c(need, sapply(keys_for(jb[1], jb[2]), function(k)
    lf(jb[1], k, code)))
  if (jb[1] == "B" && jb[2] %in% c("comp", "cv", "death")) {
    need <- c(need, file.path(OD, sprintf("optB_preds_%s_%s_%s.csv",
                                          jb[1], jb[2], code)))
  }
}
need <- c(need, file.path(OD, sprintf("sl7_preds_B_%s.csv",
                                      c("ehrecg", "3mod"))))
need <- unique(need)
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else
  sprintf("all %d", length(need)), if (length(miss) == 0) "PASS" else
    "FAIL")
stop_if_fail()
FU <- unique(fread(file.path(IN, "followup_365.csv")), by = "hadm_id")
tasks <- list()
META <- list()
CACHE <- list()
getf <- function(path) {
  if (is.null(CACHE[[path]])) CACHE[[path]] <<- fread(path)
  CACHE[[path]]
}
for (jb in LJOBS) {
  fr <- jb[1]
  kind <- jb[2]
  base <- sub("_nh$", "", kind)
  nh <- grepl("_nh$", kind)
  F0 <- getf(FRM[[fr]])
  if (fr == "B" && !("in_three_mod_p" %in% names(F0))) {
    F0[, in_three_mod_p := in_three_mod]
  }
  for (co in names(SETS)) {
    code <- SETS[[co]]$code
    sets <- SETS[[co]]$sets
    mods <- unique(unlist(sets))
    sub <- if (co == "EHR+ECG") F0 else F0[in_three_mod_p == 1]
    n <- nrow(sub)
    tt <- if (base == "death") sub$time_death else sub$time_first
    ee <- switch(base, comp = as.integer(sub$event_comp == 1),
                 cv = as.integer(sub$event_cr == 1),
                 death = as.integer(sub$event_death == 1))
    job <- paste(fr, kind, co)
    main <- lf(fr, keys_for(fr, kind)[1], code)
    MT <- getf(main)
    seeds <- sort(unique(MT$seed))
    for (s in seeds) {
      mt <- MT[MT$seed == s][match(sub$hadm_id, hadm_id)]
      f <- mt$fold
      if (anyNA(f)) {
        chk(sprintf("%s seed %s: folds complete", job, s), sum(is.na(f)),
            "FAIL")
        stop_if_fail()
      }
      build <- function(k2) {
        SS <- setNames(lapply(FAMS, function(z) list()), FAMS)
        for (fam in FAMS) for (m in mods) {
          cf <- colfor(fr, k2, m, fam)
          D0 <- getf(lf(fr, cf[1], code))
          x <- D0[D0$seed == s][match(sub$hadm_id, hadm_id)]
          if (!(cf[2] %in% names(x)) || anyNA(x$hadm_id)) {
            chk(sprintf("%s seed %s: %s", job, s, cf[2]), "missing",
                "FAIL")
            stop_if_fail()
          }
          if ("fold" %in% names(x) && any(x$fold != f)) {
            chk(sprintf("%s seed %s: folds of %s match", job, s,
                        basename(lf(fr, cf[1], code))), "no", "FAIL")
            stop_if_fail()
          }
          SS[[fam]][[m]] <- x[[cf[2]]]
        }
        SS
      }
      dkind <- if (nh) "death_nh" else "death"
      S <- build(kind)
      D <- if (base == "cv") build(dkind) else NULL
      C1 <- if (base == "comp") build(dkind) else NULL
      C2 <- if (base == "comp") build(if (nh) "cv_nh" else "cv") else NULL
      bad <- sum(sapply(c(S, D, C1, C2), function(z)
        sum(sapply(z, function(v) sum(!is.finite(v))))))
      if (bad > 0) {
        chk(sprintf("%s seed %s: missing scores", job, s), bad, "FAIL")
        stop_if_fail()
      }
      tasks[[length(tasks) + 1]] <- list(
        job = job, seed = s, fold = f, subj = sub$subject_id, t = tt,
        e = ee, S = S, D = D, C1 = C1, C2 = C2, cv = (base == "cv"),
        sets = sets, hz = HZ, tcuts = TCUTS, tmax = TMAX)
    }
    DF <- if (base == "cv") {
      matrix(sapply(HZ, function(h) as.integer(sub$event_cr == 2 &
                                                 sub$time_first <= h + 0.5)),
             ncol = length(HZ))
    } else NULL
    dl <- toupper(trimws(as.character(
      FU$discharge_location[match(sub$hadm_id, FU$hadm_id)])))
    META[[job]] <- list(fr = fr, kind = kind, base = base, nh = nh,
                        co = co, code = code, sub = sub, sets = sets,
                        t = tt, e = ee, DF = DF,
                        hosp = !is.na(dl) & dl == "HOSPICE")
    cat(sprintf("  %-24s n=%d  events by 30/365 days %d/%d  seeds %d\n",
                job, n, sum(ee == 1 & tt <= 30.5), sum(ee == 1 & tt <= 365.5),
                length(seeds)))
  }
}
cat(sprintf("  %d tasks (frame x outcome x cohort x seed)\n",
            length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) length(t$t) * length(t$sets) *
                          (1 + !is.null(t$C1))))
t0 <- Sys.time()
invisible(sl_task(tasks[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  largest task %.0f s; %d tasks on %d workers: at most",
                  "about %.1f min, plus the bootstrap\n"), tp,
            length(tasks), workers,
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
plan(sequential)

# ---------------- assemble ----------------
cat("\n=== ASSEMBLE ===\n")
R <- function(m) paste0(m, ", recalibrated")
PALL <- list()
for (job in names(META)) {
  M <- META[[job]]
  n <- nrow(M$sub)
  rv <- Filter(function(r) r$job == job, res)
  P <- list()
  for (k in names(rv[[1]]$out)) {
    P[[k]] <- Reduce(`+`, lapply(rv, function(r) r$out[[k]])) / length(rv)
  }
  bad <- sum(sapply(P, function(m) sum(!is.finite(m))))
  rg <- range(unlist(P))
  chk(sprintf("%s: missing risks / range", job),
      sprintf("%d / %.5f to %.5f", bad, rg[1], rg[2]),
      if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
  LONG <- rbindlist(lapply(names(P), function(k) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = rep(M$sub$hadm_id, length(HZ)), model = z[1],
               set = z[2], horizon = rep(HZ, each = n),
               risk = as.vector(P[[k]]))
  }))
  fwrite(LONG, file.path(OD, sprintf("optB2_preds_%s_%s_%s.csv", M$fr,
                                     M$kind, M$code)))
  # version 1 reproduced exactly (discharge only; the at-diagnosis inputs
  # were rebuilt in A36)
  if (!M$nh && M$fr == "B") {
    V1 <- fread(file.path(OD, sprintf("optB_preds_%s_%s_%s.csv", M$fr,
                                      M$kind, M$code)))
    d1 <- 0
    n1 <- 0L
    for (k in grep("^Super learner(-EHR)?, recalibrated\\|", names(P),
                   value = TRUE)) {
      z <- strsplit(k, "|", fixed = TRUE)[[1]]
      x <- V1[model == z[1] & set == z[2]]
      if (!nrow(x)) next
      for (hi in seq_along(HZ)) {
        y <- x[horizon == HZ[hi]]
        d1 <- max(d1, abs(y$risk[match(M$sub$hadm_id, y$hadm_id)] -
                            P[[k]][, hi]))
        n1 <- n1 + 1L
      }
    }
    chk(sprintf("%s: first version reproduced (n, max |diff|)", job),
        sprintf("%d, %.1e", n1, d1),
        if (n1 > 0 && d1 <= 1e-9) "PASS" else "FAIL")
  }
  if (M$base == "death") {
    V <- fread(file.path(OD, sprintf("sl7_preds_B_%s.csv", M$code)))
    V <- V[model == "Super learner v7, recalibrated"]
    nv <- 0L
    for (sn in names(M$sets)) {
      x <- V[set == sn]
      mm <- sapply(HZ, function(h) {
        y <- x[horizon == h]
        y$risk[match(M$sub$hadm_id, y$hadm_id)]
      })
      nv <- nv + sum(is.na(mm))
      P[[paste("Super learner v7, recalibrated", sn, sep = "|")]] <- mm
    }
    chk(sprintf("%s: Super learner v7 predictions present (missing)", job),
        nv, if (nv == 0) "PASS" else "FAIL")
  }
  PALL[[job]] <- P
}
stop_if_fail()
# the no-hospice super learners join their full counterparts' jobs, so the
# contrast is paired on the same patients and resamples
for (job in names(META)) {
  M <- META[[job]]
  if (!M$nh) next
  tgt <- paste(M$fr, M$base, M$co)
  for (sn in c(names(M$sets), "EHR")) {
    for (mn in c("Super learner", "Super learner-EHR")) {
      k <- paste(R(mn), sn, sep = "|")
      if (is.null(PALL[[job]][[k]])) next
      PALL[[tgt]][[paste(R(paste(mn, "(no hospice)")), sn,
                         sep = "|")]] <- PALL[[job]][[k]]
    }
  }
}
JOBS <- list()
for (job in names(META)) {
  M <- META[[job]]
  if (M$nh) next
  n <- nrow(M$sub)
  P <- PALL[[job]]
  Y <- matrix(sapply(HZ, function(h) as.integer(M$e == 1 &
                                                  M$t <= h + 0.5)),
              ncol = length(HZ))
  keepP <- grep("^(Super learner|SL candidate)", names(P), value = TRUE)
  mk <- function(rows, ev) {
    list(job = job, M = M, hz = HZ, evaluation = ev,
         P = lapply(P[keepP], function(m) m[rows, , drop = FALSE]),
         Y = Y[rows, , drop = FALSE],
         DF = if (is.null(M$DF)) NULL else M$DF[rows, , drop = FALSE],
         rows_by = split(seq_along(rows),
                         match(M$sub$subject_id[rows],
                               unique(M$sub$subject_id[rows]))))
  }
  JOBS[[paste(job, "all")]] <- mk(seq_len(n), "all")
  if (M$fr == "B") {
    kr <- which(!M$hosp)
    chk(sprintf("%s: hospice discharges removed in the second evaluation",
                job), sprintf("%d of %d", sum(M$hosp), n), "INFO")
    JOBS[[paste(job, "no hospice")]] <- mk(kr, "hospice excluded")
  }
}
W <- rbindlist(lapply(res, function(r) r$wts))
wsum <- W[, .(folds = .N, mean_weight = mean(weight),
              share_nonzero = mean(weight > 1e-6)),
          by = .(job, model, set, candidate)]
setorder(wsum, job, model, set, -mean_weight)
fwrite(wsum, file.path(RS, "optB2_sl_weights.csv"))

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples, clustered by patient) ===\n",
            NBOOT))
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
plan(multisession, workers = workers)
for (J in JOBS) {
  M <- J$M
  JB <- J
  JB$P <- J$P[grep("^Super learner", names(J$P), value = TRUE)]
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  BL <- future_lapply(seq_len(NBOOT), boot_rep, J = JB,
                      future.seed = SEED, future.packages = PKG)
  bn <- names(stat_vec(JB, seq_len(nrow(J$Y))))
  B <- do.call(rbind, lapply(BL, function(v) v[bn]))
  colnames(B) <- bn
  cat(sprintf("  %s (%s) done\n", J$job, J$evaluation))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- if (k %in% colnames(B)) q2(B[, k]) else c(NA_real_, NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      frame = M$fr, outcome = M$base, cohort = M$co,
      evaluation = J$evaluation, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in grep("^Super learner", names(J$P), value = TRUE)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        frame = M$fr, outcome = M$base, cohort = M$co,
        evaluation = J$evaluation, model = z[1], set = z[2],
        horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  pairs <- list()
  for (sn in names(M$sets)) {
    pairs[[length(pairs) + 1]] <- c(R("Super learner"), sn,
                                    R("Super learner-EHR"), "EHR")
    if (M$base == "comp") {
      pairs[[length(pairs) + 1]] <- c(R("Super learner+"), sn,
                                      R("Super learner"), sn)
    }
    if (M$base == "death") {
      pairs[[length(pairs) + 1]] <- c(R("Super learner v7"), sn,
                                      R("Super learner"), sn)
    }
    if (M$fr == "B" && M$base %in% c("comp", "cv")) {
      pairs[[length(pairs) + 1]] <- c(R("Super learner"), sn,
                                      R("Super learner (no hospice)"), sn)
    }
  }
  ms <- c("auc", "ipa", if (M$base == "cv") "auc_excl")
  for (cd in pairs) for (h in J$hz) for (m0 in ms) {
    ka <- paste(cd[1], cd[2], h, m0, sep = "|")
    kb <- paste(cd[3], cd[4], h, m0, sep = "|")
    if (!(ka %in% colnames(B)) || !(kb %in% colnames(B))) next
    q <- q2(B[, ka] - B[, kb])
    CON[[length(CON) + 1]] <- data.table(
      frame = M$fr, outcome = M$base, cohort = M$co,
      evaluation = J$evaluation, set = cd[2],
      contrast = sprintf("%s - %s", cd[1], cd[3]), horizon = h,
      measure = m0, estimate = pt[[ka]] - pt[[kb]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0)
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
fwrite(mt, file.path(RS, "optB2_sl_metrics.csv"))
fwrite(cn, file.path(RS, "optB2_sl_contrasts.csv"))

# ---------------- decisions, checks and summary ----------------
cc <- cn[outcome == "death" & cohort == "EHR+ECG" & set == "EHR+ECG" &
           evaluation == "all" & measure == "auc" & grepl("v7", contrast)]
rep_ok <- nrow(cc) == length(HZ) && all(abs(cc$estimate) <= CORE_TOL)
core <- cc[, .(horizon, v7_minus_core = estimate, lo, hi)]
core[, representative := rep_ok]
fwrite(core, file.path(RS, "optB2_core_check.csv"))
chk("v7-core within 0.005 AUROC of v7 at every horizon (primary cohort)",
    sprintf("max |diff| %.4f", max(abs(cc$estimate))),
    if (rep_ok) "PASS" else "WARN")
DEC <- list()
for (fr in c("P", "B")) {
  x <- cn[frame == fr & outcome == "comp" & cohort == "EHR+ECG" &
            set == "EHR+ECG" & evaluation == "all" & measure == "auc" &
            grepl("^Super learner\\+", contrast)]
  better <- sum(x$excludes_zero & x$estimate > 0)
  worse <- sum(x$excludes_zero & x$estimate < 0)
  adopt <- nrow(x) == length(HZ) && better >= 2 && worse == 0
  DEC[[fr]] <- data.table(
    pipeline = if (fr == "P") "at diagnosis" else "at discharge",
    contrast = "SL+ (components) minus SL, composite, AUROC",
    better = better, worse = worse, adopt = adopt,
    final_composite_model = if (adopt) R("Super learner+") else
      R("Super learner"))
  chk(sprintf("A34 rule, composite %s", DEC[[fr]]$pipeline),
      DEC[[fr]]$final_composite_model, "INFO")
}
fwrite(rbindlist(DEC), file.path(RS, "optB2_decision.csv"))
prim <- (mt$cohort == "EHR+ECG" & mt$set %in% c("EHR+ECG", "EHR")) |
  (mt$cohort == "3-mod" & mt$set %in% c("EHR+ECG+CTPA", "EHR"))
smy <- mt[prim & grepl("^Super learner", model) &
            measure %in% c("auc", "auc_excl", "ipa", "oe", "slope")]
fwrite(smy, file.path(RS, "optB2_sl_summary.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "optB2_sl_checks.csv"))
cat("\n=== SUPER LEARNER AUROC (primary sets, every patient) ===\n")
x <- smy[measure %in% c("auc", "auc_excl") & evaluation == "all"]
print(dcast(x, frame + outcome + cohort + model + measure ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE, nrows = 300)
cat("\n=== A34 DECISIONS ===\n")
print(rbindlist(DEC), row.names = FALSE)
cat("\n=== v7-CORE CHECK ===\n")
print(core, row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (optB2_preds_*.csv, patient-level, local only)\n")
for (f in c("optB2_sl_summary.csv", "optB2_sl_contrasts.csv",
            "optB2_decision.csv", "optB2_sl_metrics.csv",
            "optB2_core_check.csv", "optB2_sl_weights.csv",
            "optB2_sl_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "optB2_sl_log.txt"), " (terminal log)\n")
cat("\nOutputs: optB2_sl_summary.csv, optB2_sl_contrasts.csv,",
    "optB2_decision.csv, optB2_sl_checks.csv\n")