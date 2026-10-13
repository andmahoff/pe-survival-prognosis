# 18q_super_learner_P.R
# Super learner P, the final at-diagnosis death model, compared with sPESI,
# PESI-10, the EHR alone and Frame A's super learner v4
# (PREREGISTRATION_ADDENDUM_2.md A32; D17, D18).
#
# Layer 1 uses the out-of-fold scores of 6p for the EHR at presentation,
# the first ECG and the presentation CTPA, in four families: ridge-Cox
# retuned, spline Cox, RSF Brier-tuned and oblique RSF Brier-tuned (the
# ECG and CTPA spline Cox candidate uses ridge-Cox scores, as in v4).
# Candidates (Kaplan-Meier, EHR alone, and PH, TV and TVC stacks for each
# set and family), weights and per-horizon recalibration are as in 18h.
# Four super learners are fitted: P (every candidate of the set), P-EHR
# (the EHR-alone candidates only), and Psens and Psens-EHR, which use the
# EHR scores of 6p's --sens run. Horizons are 10, 30, 90, 180 and 365 days
# from admission. Cohorts: EHR+ECG (every Frame P admission) and the
# three-modality cohort with a presentation CTPA, which adds the sets
# EHR+CTPA and EHR+ECG+CTPA.
#
# Benchmarks: sPESI and PESI-10 (4p; sensitivity versions go with Psens)
# and Frame A's super learner v4, recalibrated (10 and 30 days). Measures
# are as in 18h (AUROC only for the scores), plus the model's specificity
# at sPESI's sensitivity (sPESI of 1 or more is high risk; threshold
# matched within each resample) minus sPESI's specificity. Bootstrap:
# 2,000 samples, seed 42, clustered by patient, on every patient and (D17)
# on patients with presentation vital signs only.
#
# Contrasts: P minus sPESI and PESI-10 (AUROC, specificity); P minus P-EHR
# and Frame A v4 minus P at 10 and 30 days (AUROC, IPA); Psens minus the
# sensitivity sPESI and PESI-10, and minus Psens-EHR.
# Checks (a failure stops the script): inputs found and aligned; 6p's main
# and sensitivity folds match; sPESI, PESI-10 and Frame A v4 risks present;
# risks in [0, 1].
#
# Inputs:  data/processed/presentation/frame_P.csv, scores_P.csv;
#          data/processed/level1/oof_P_*.csv, oof_Psens_*.csv;
#          data/processed/super_learner/sl5_preds_A_*.csv
# Outputs: results/slP_metrics.csv, results/slP_contrasts.csv,
#          results/slP_weights.csv, results/slP_summary.csv,
#          results/slP_checks.csv, results/slP_log.txt;
#          data/processed/super_learner/slP_preds_*.csv (patient level,
#          not in the repo)
# Usage:   Rscript scripts/18q_super_learner_P.R [--probe-only]
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
PP <- file.path(PR, "presentation")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
OD <- file.path(PR, "super_learner")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
TCUTS <- c(1, 2, 4, 7, 14, 30, 60, 120, 240)
TMAX <- 366
FAMS <- c("ridge-Cox retuned", "spline Cox", "RSF, Brier-tuned",
          "oblique RSF, Brier-tuned")
COL <- list(ehr = c("ridge-Cox retuned" = "ehr.cox2.death.lp",
                    "spline Cox" = "ehr.spline.death.lp",
                    "RSF, Brier-tuned" = "ehr.rsfb.death.score",
                    "oblique RSF, Brier-tuned" = "ehr.orsfb.death.score"),
            ecg = c("ridge-Cox retuned" = "ecg.cox2.death.lp",
                    "spline Cox" = "ecg.cox2.death.lp",
                    "RSF, Brier-tuned" = "ecg.rsfb.death.score",
                    "oblique RSF, Brier-tuned" = "ecg.orsfb.death.score"),
            ctpa = c("ridge-Cox retuned" = "ctpa.cox2.death.lp",
                     "spline Cox" = "ctpa.cox2.death.lp",
                     "RSF, Brier-tuned" = "ctpa.rsfb.death.score",
                     "oblique RSF, Brier-tuned" =
                       "ctpa.orsfb.death.score"))
COH <- list("EHR+ECG" = list(code = "ehrecg",
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
    fwrite(ck, file.path(RS, "slP_checks.csv"))
    cat("\nSTOP: a check failed. see: slP_checks.csv\n")
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
  }
  lapply(out, function(m) matrix(m, nrow = length(b)))
}

# ---------------- one task = cohort x seed x version ----------------
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
      put(paste0("SL candidate", tk$sfx, ": ", z[2]), z[1], te, CT[[nm]])
    }
    ce <- names(CT)[startsWith(names(CT), "EHR|")]
    fit_sl(ce, PI, CT, Y, paste0("Super learner ", tk$ver, "-EHR"), "EHR",
           k, te)
    for (sn in names(sets)) {
      ca <- c(ce, names(CT)[startsWith(names(CT), paste0(sn, "|"))])
      fit_sl(ca, PI, CT, Y, paste0("Super learner ", tk$ver), sn, k, te)
    }
  }
  list(co = tk$co, seed = tk$seed, ver = tk$ver, out = out,
       wts = rbindlist(wts)[, `:=`(cohort = tk$co, seed = tk$seed)])
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
      raw <- J$P[[k]][idx, hi]
      if (anyNA(raw)) next
      key <- paste(k, J$hz[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, raw)
      if (k %in% J$scores) next
      p <- pmin(pmax(raw, 1e-6), 1 - 1e-6)
      v[paste(key, "ipa", sep = "|")] <- 1 - mean((y - p)^2) /
        (pr * (1 - pr))
      v[paste(key, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
      v[paste(key, "oe", sep = "|")] <- pr / mean(p)
    }
    # specificity at sPESI's sensitivity (sPESI 1 or more = high risk):
    # the model's threshold gives the same sensitivity in this resample
    # (every key is always written, NA when a resample has no events or
    # no survivors, so all resamples return the same columns in order)
    ok <- sum(y) > 0 && sum(y) < length(y)
    for (sp in names(J$spec)) {
      sc <- J$P[[J$spec[[sp]]$score]][idx, hi]
      se <- if (ok) mean(sc[y == 1] >= 1) else NA_real_
      sps <- if (ok) mean(sc[y == 0] < 1) else NA_real_
      for (k in J$spec[[sp]]$models) {
        val <- NA_real_
        if (ok) {
          p <- J$P[[k]][idx, hi]
          pe <- sort(p[y == 1], decreasing = TRUE)
          thr <- pe[max(1, ceiling(se * length(pe)))]
          val <- mean(p[y == 0] < thr) - sps
        }
        v[paste(k, J$hz[hi], paste("specdiff vs", sp), sep = "|")] <- val
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
need <- c(file.path(PP, c("frame_P.csv", "scores_P.csv")),
          file.path(L1, sprintf("oof_%s_%s.csv",
                                rep(c("P", "Psens"), each = 2),
                                c("ehrecg", "3mod"))),
          file.path(OD, sprintf("sl5_preds_A_%s.csv", c("ehrecg", "3mod"))))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- fread(file.path(PP, "frame_P.csv"))
SCO <- unique(fread(file.path(PP, "scores_P.csv")), by = "hadm_id")
tasks <- list()
META <- list()
for (co in names(COH)) {
  C <- COH[[co]]
  sub <- if (co == "EHR+ECG") fp else fp[in_three_mod_p == 1]
  n <- nrow(sub)
  mods <- unique(unlist(C$sets))
  T1 <- fread(file.path(L1, sprintf("oof_P_%s.csv", C$code)))
  TS <- fread(file.path(L1, sprintf("oof_Psens_%s.csv", C$code)))
  seeds <- sort(unique(T1$seed))
  ok <- identical(seeds, sort(unique(TS$seed))) &&
    nrow(T1) == n * length(seeds) && nrow(TS) == n * length(seeds)
  chk(sprintf("%s: stage 6p seeds and rows aligned (main, sens)", co),
      sprintf("%d seeds, %d x %d rows", length(seeds), n, length(seeds)),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  for (s in seeds) {
    a <- T1[T1$seed == s][match(sub$hadm_id, hadm_id)]
    b <- TS[TS$seed == s][match(sub$hadm_id, hadm_id)]
    if (anyNA(a$fold) || any(a$fold != b$fold)) {
      chk(sprintf("%s seed %s: main and sens folds match", co, s), "no",
          "FAIL")
      stop_if_fail()
    }
    S <- setNames(lapply(FAMS, function(z) list()), FAMS)
    for (fam in FAMS) for (m in mods) S[[fam]][[m]] <- a[[COL[[m]][[fam]]]]
    Ss <- S
    for (fam in FAMS) Ss[[fam]][["ehr"]] <- b[[COL$ehr[[fam]]]]
    bad <- sum(sapply(c(S, Ss), function(z)
      sum(sapply(z, function(v) sum(!is.finite(v))))))
    if (bad > 0) {
      chk(sprintf("%s seed %s: missing scores", co, s), bad, "FAIL")
      stop_if_fail()
    }
    for (ver in c("P", "Psens")) {
      tasks[[length(tasks) + 1]] <- list(
        co = co, seed = s, ver = ver,
        sfx = if (ver == "P") "" else " (sens)", fold = a$fold,
        subj = sub$subject_id, t = sub$time_death,
        e = as.integer(sub$event_death == 1),
        S = if (ver == "P") S else Ss, sets = C$sets, hz = HZ,
        tcuts = TCUTS, tmax = TMAX)
    }
  }
  META[[co]] <- list(co = co, code = C$code, sub = sub, sets = C$sets)
  cat(sprintf("  %-7s n=%d deaths by 365 days=%d seeds=%d\n", co, n,
              sum(sub$event_death), length(seeds)))
}
cat(sprintf("  %d tasks (cohort x seed x version)\n", length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) length(t$t) * length(t$sets)))
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

# ---------------- assemble ----------------
cat("\n=== ASSEMBLE ===\n")
R <- function(v) sprintf("Super learner %s, recalibrated", v)
JOBS <- list()
for (co in names(META)) {
  M <- META[[co]]
  n <- nrow(M$sub)
  P <- list()
  for (ver in c("P", "Psens")) {
    rv <- Filter(function(r) r$co == co && r$ver == ver, res)
    for (k in names(rv[[1]]$out)) {
      P[[k]] <- Reduce(`+`, lapply(rv, function(r) r$out[[k]])) /
        length(rv)
    }
  }
  bad <- sum(sapply(P, function(m) sum(!is.finite(m))))
  rg <- range(unlist(P))
  chk(sprintf("%s: missing risks / range", co),
      sprintf("%d / %.5f to %.5f", bad, rg[1], rg[2]),
      if (bad == 0 && rg[1] >= 0 && rg[2] <= 1) "PASS" else "FAIL")
  LONG <- rbindlist(lapply(names(P), function(k) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = rep(M$sub$hadm_id, length(HZ)), model = z[1],
               set = z[2], horizon = rep(HZ, each = n),
               risk = as.vector(P[[k]]))
  }))
  fwrite(LONG, file.path(OD, sprintf("slP_preds_%s.csv", M$code)))
  sc <- SCO[match(M$sub$hadm_id, hadm_id)]
  scm <- function(v) matrix(as.numeric(v), n, length(HZ))
  P[["sPESI|score"]] <- scm(sc$spesi)
  P[["PESI-10|score"]] <- scm(sc$pesi10)
  P[["sPESI (sens)|score"]] <- scm(sc$spesi_sens)
  P[["PESI-10 (sens)|score"]] <- scm(sc$pesi10_sens)
  ns <- sum(sapply(P[grep("\\|score$", names(P))], function(m)
    sum(!is.finite(m))))
  chk(sprintf("%s: sPESI and PESI-10 present (missing)", co), ns,
      if (ns == 0) "PASS" else "FAIL")
  A <- fread(file.path(OD, sprintf("sl5_preds_A_%s.csv", M$code)))
  A <- A[model == "Super learner v4, recalibrated"]
  na <- 0L
  for (sn in names(M$sets)) {
    x <- A[set == sn]
    if (!nrow(x)) next
    mm <- sapply(HZ, function(h) {
      y <- x[horizon == h]
      if (!nrow(y)) return(rep(NA_real_, n))
      y$risk[match(M$sub$hadm_id, y$hadm_id)]
    })
    na <- na + sum(is.na(mm[, HZ <= 30]))
    P[[paste("Frame A v4, recalibrated", sn, sep = "|")]] <- mm
  }
  chk(sprintf("%s: Frame A v4 risks at 10 and 30 days (missing)", co), na,
      if (na == 0) "PASS" else "FAIL")
  Y <- sapply(HZ, function(h) as.integer(M$sub$event_death == 1 &
                                           M$sub$time_death <= h + 0.5))
  Y <- matrix(Y, ncol = length(HZ))
  spec <- list(
    "sPESI" = list(score = "sPESI|score",
                   models = paste(R("P"), names(M$sets), sep = "|")),
    "sPESI (sens)" = list(score = "sPESI (sens)|score",
                          models = paste(R("Psens"), names(M$sets),
                                         sep = "|")))
  scores <- grep("\\|score$", names(P), value = TRUE)
  mk_job <- function(rows, ev) {
    list(co = co, hz = HZ, P = lapply(P, function(m) m[rows, ,
                                                         drop = FALSE]),
         Y = Y[rows, , drop = FALSE],
         rows_by = split(seq_along(rows),
                         match(M$sub$subject_id[rows],
                               unique(M$sub$subject_id[rows]))),
         sets = M$sets, scores = scores, spec = spec, evaluation = ev)
  }
  JOBS[[paste(co, "all")]] <- mk_job(seq_len(n), "all")
  kv <- which(M$sub$has_pres_vitals == 1)
  chk(sprintf("%s: patients with presentation vital signs (D17)", co),
      sprintf("%d of %d", length(kv), n), "INFO")
  JOBS[[paste(co, "vitals")]] <- mk_job(kv,
                                        "with presentation vital signs")
}
stop_if_fail()
plan(sequential)
W <- rbindlist(lapply(res, function(r) r$wts))
wsum <- W[, .(folds = .N, mean_weight = mean(weight),
              share_nonzero = mean(weight > 1e-6)),
          by = .(cohort, model, set, candidate)]
setorder(wsum, cohort, model, set, -mean_weight)
fwrite(wsum, file.path(RS, "slP_weights.csv"))

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples, clustered by patient) ===\n",
            NBOOT))
CONS <- function(sn) {
  list(c(R("P"), sn, "sPESI", "score"),
       c(R("P"), sn, "PESI-10", "score"),
       c(R("P"), sn, R("P-EHR"), "EHR"),
       c("Frame A v4, recalibrated", sn, R("P"), sn),
       c(R("Psens"), sn, "sPESI (sens)", "score"),
       c(R("Psens"), sn, "PESI-10 (sens)", "score"),
       c(R("Psens"), sn, R("Psens-EHR"), "EHR"))
}
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
plan(multisession, workers = workers)
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(nrow(J$Y)))
  BL <- future_lapply(seq_len(NBOOT), boot_rep, J = J,
                      future.seed = SEED, future.packages = PKG)
  B <- do.call(rbind, lapply(BL, function(v) v[names(pt)]))
  colnames(B) <- names(pt)
  cat(sprintf("  %s (%s) done\n", J$co, J$evaluation))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- if (k %in% colnames(B)) q2(B[, k]) else c(NA_real_, NA_real_)
    MET[[length(MET) + 1]] <- data.table(
      cohort = J$co, evaluation = J$evaluation, model = z[1], set = z[2],
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in setdiff(names(J$P), J$scores)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      if (anyNA(J$P[[k]][, hi])) next
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        cohort = J$co, evaluation = J$evaluation, model = z[1], set = z[2],
        horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  cds <- unlist(lapply(names(J$sets), CONS), recursive = FALSE)
  for (cd in cds) {
    for (h in J$hz) for (ms in c("auc", "ipa")) {
      ka <- paste(cd[1], cd[2], h, ms, sep = "|")
      kb <- paste(cd[3], cd[4], h, ms, sep = "|")
      if (!(ka %in% colnames(B)) || !(kb %in% colnames(B))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        cohort = J$co, evaluation = J$evaluation, set = cd[2],
        contrast = sprintf("%s - %s", cd[1], cd[3]), horizon = h,
        measure = ms, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
  for (k in grep("\\|specdiff vs ", names(pt), value = TRUE)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    CON[[length(CON) + 1]] <- data.table(
      cohort = J$co, evaluation = J$evaluation, set = z[2],
      contrast = sprintf("%s - %s", z[1], sub("specdiff vs ", "", z[4])),
      horizon = as.numeric(z[3]),
      measure = "specificity at sPESI's sensitivity",
      estimate = pt[[k]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0)
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
fwrite(mt, file.path(RS, "slP_metrics.csv"))
fwrite(cn, file.path(RS, "slP_contrasts.csv"))

# ---------------- summary ----------------
pick <- mt[evaluation == "all" & measure == "auc" &
             ((cohort == "EHR+ECG" & set %in% c("EHR+ECG", "EHR",
                                                  "score")) |
                (cohort == "3-mod" & set %in% c("EHR+ECG+CTPA", "EHR",
                                                  "score"))) &
             model %in% c(R("P"), R("P-EHR"), "sPESI", "PESI-10",
                          "Frame A v4, recalibrated")]
smy <- rbind(
  pick[, .(cohort, set, horizon, item = paste("AUROC:", model),
           estimate, lo, hi, excludes_zero = NA)],
  cn[evaluation == "all" &
       ((cohort == "EHR+ECG" & set == "EHR+ECG") |
          (cohort == "3-mod" & set == "EHR+ECG+CTPA")) &
       !grepl("sens", contrast),
     .(cohort, set, horizon, item = paste0(measure, ": ", contrast),
       estimate, lo, hi, excludes_zero)])
setorder(smy, cohort, horizon, item)
fwrite(smy, file.path(RS, "slP_summary.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "slP_checks.csv"))
cat("\n=== THE AT-DIAGNOSIS MODEL AGAINST sPESI AND PESI-10 (AUROC) ===\n")
print(dcast(smy[grepl("^AUROC", item)], cohort + item ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE)
cat("\n=== CONTRASTS (every patient) ===\n")
print(smy[!grepl("^AUROC", item),
          .(cohort, horizon, item, est = round(estimate, 3),
            lo = round(lo, 3), hi = round(hi, 3), excludes_zero)],
      row.names = FALSE, nrows = 500)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (slP_preds_*.csv, patient-level, local only)\n")
for (f in c("slP_summary.csv", "slP_metrics.csv", "slP_contrasts.csv",
            "slP_weights.csv", "slP_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "slP_log.txt"), " (terminal log)\n")
cat("\nOutputs: slP_summary.csv, slP_contrasts.csv, slP_metrics.csv,",
    "slP_weights.csv, slP_checks.csv\n")