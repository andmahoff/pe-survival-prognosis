# 19b_spline_cox.R
# Penalised-spline Cox model for the EHR (A19): is the forests' EHR
# advantage smooth non-linearity that an additive Cox model can capture?
#
# Level 1 gives EHR death scores for every frame, cohort and seed, out of
# fold over the folds of 04, with the preprocessing of 18a (median
# imputation, winsorising, standardisation on the training fold). The model
# is mgcv::gam(family = cox.ph). Every feature with at least 10 distinct
# values in the training fold enters as a cubic regression spline (s(x, k =
# 5)) and the others enter linearly. Smoothness is chosen by REML with
# select = TRUE, so a whole term can be shrunk to zero. The outputs are the
# linear predictor and the risk at every horizon (1 - S(h + 0.5)), at 10 and
# 30 days in Frame A and 10 to 365 days in Frame B. The scores are compared
# with the EHR scores of ridge-Cox retuned and oblique RSF (06d) and RSF
# (06b), each through the same stacks.
#
# Stacks are fitted per seed and outer fold on rank-based normal scores of
# the training rows (as S4 and A14): a Cox stack (PH) and per-horizon
# logistic stacks (TV) for EHR alone, EHR + ECG and, in the three-modality
# cohort, EHR + ECG + CTPA. The ECG and CTPA scores are 06d's retuned
# ridge-Cox. Measures are AUROC, IPA, log-likelihood skill and calibration
# slope, and the spline model's own risks are compared with ridge-Cox
# retuned's. Contrasts, for EHR alone and fused: spline minus ridge-Cox
# retuned, oblique RSF minus spline and RSF minus spline, with a
# 2,000-sample bootstrap (seed 42), patient-clustered in Frame A.
#
# Checks: every fit succeeded and every patient scored (stop); fits short of
# full REML convergence (warning above 10%); the D2 per-fold rule (stop);
# inputs aligned (stop). Every finished level-1 task is saved in
# spline_ckpt/, so a rerun fits only the missing ones.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_extra_*.csv,
#          level1/oof_forest_*.csv; data/inputs/ehr28.csv
# Outputs: results/spline_{metrics,contrasts,checks}.csv,
#          results/spline_log.txt, data/processed/spline/
#          spline_scores_*.csv and spline_ckpt/ (patient level, not in
#          the repo)
# Usage:   Rscript scripts/19b_spline_cox.R [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
OD <- file.path(PR, "spline")
CKPT <- file.path(OD, "spline_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
KNOTS <- 5L
MINU <- 10L
NBOOT <- 2000L
SEED <- 42L
FAMS <- c("ridge-Cox retuned", "spline Cox", "RSF", "oblique RSF")
PKG <- c("survival", "mgcv", "data.table")
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
    fwrite(ck, file.path(RS, "spline_checks.csv"))
    cat("\nSTOP: a check failed. see: spline_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- level 1: additive Cox model ----------------
spline_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                        probe = FALSE) {
  n <- length(tk$t)
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
  lp <- rep(NA_real_, n)
  risk <- matrix(NA_real_, n, length(tk$hz))
  info <- list()
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- prep(tk$X, tr)
    colnames(X) <- sprintf("f%d", seq_len(ncol(X)))
    nu <- apply(X[tr, , drop = FALSE], 2, function(v) length(unique(v)))
    sm <- colnames(X)[nu >= tk$minu]
    li <- setdiff(colnames(X), sm)
    rhs <- c(sprintf("s(%s, k = %d, bs = \"cr\")", sm, tk$knots), li)
    fml <- as.formula(paste("t ~", paste(rhs, collapse = " + ")))
    dtr <- data.frame(t = tk$t[tr], e = tk$e[tr], X[tr, , drop = FALSE])
    dte <- data.frame(X[te, , drop = FALSE])
    fit <- tryCatch(gam(fml, family = cox.ph(), weights = e, data = dtr,
                        method = "REML", select = TRUE),
                    error = function(err) conditionMessage(err))
    if (is.character(fit)) {
      info[[length(info) + 1]] <- data.table(
        fold = k, ok = FALSE, conv = "error", msg = substr(fit, 1, 120),
        edf = NA_real_, n_smooth = length(sm))
      next
    }
    lp[te] <- as.vector(predict(fit, newdata = dte, type = "link"))
    for (hi in seq_along(tk$hz)) {
      nd <- dte
      nd$t <- tk$hz[hi] + tk$off
      risk[te, hi] <- 1 - as.vector(predict(fit, newdata = nd,
                                            type = "response"))
    }
    cv <- fit$outer.info$conv
    info[[length(info) + 1]] <- data.table(
      fold = k, ok = TRUE, conv = if (is.null(cv)) "" else as.character(cv),
      msg = "", edf = sum(fit$edf), n_smooth = length(sm))
    rm(fit)
    invisible(gc(verbose = FALSE))
  }
  out <- list(fr = tk$fr, co = tk$co, seed = tk$seed, lp = lp,
              risk = risk, info = rbindlist(info))
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, tk$ckpt)
  }
  out
}

# ---------------- level 2: stacks ----------------
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
stack_fit <- function(S, ms, tr, te, tt, ee, cuts) {
  Ztr <- matrix(sapply(ms, function(m) nscore(S[[m]][tr], S[[m]][tr])),
                ncol = length(ms))
  Zte <- matrix(sapply(ms, function(m) nscore(S[[m]][te], S[[m]][tr])),
                ncol = length(ms))
  D <- data.frame(t = tt[tr], e = ee[tr], Ztr)
  b <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
  b[is.na(b)] <- 0
  lt <- as.vector(Ztr %*% b)
  lv <- as.vector(Zte %*% b)
  H <- bres_cuts(tt[tr], ee[tr], lt, cuts)
  ph <- sapply(seq_along(cuts), function(i) -expm1(-H[i] *
                                                     exp(lv - max(lt))))
  tv <- sapply(cuts, function(cu) {
    y <- as.numeric(ee[tr] == 1 & tt[tr] <= cu)
    g <- suppressWarnings(glm.fit(cbind(1, Ztr), y, family = binomial()))
    cf <- g$coefficients
    cf[is.na(cf)] <- 0
    plogis(as.vector(cbind(1, Zte) %*% cf))
  })
  list(PH = matrix(ph, nrow = length(te)), TV = matrix(tv,
                                                       nrow = length(te)))
}
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
    }
  }
  v
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                     length(J$rows_by), replace = TRUE)],
                use.names = FALSE)
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat(sprintf("=== BUILD (mgcv %s) ===\n", as.character(packageVersion("mgcv"))))
fr_d <- list(A = fread(file.path(PR, "frame_A.csv")),
             B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
c28 <- intersect(wl[modality == "ehr", column], names(e28))
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  key <- paste(fr, co)
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fr_d[[fr]] else fr_d[[fr]][in_three_mod == 1]
  X <- as.matrix(e28[match(sub$hadm_id, hadm_id), ..c28])
  storage.mode(X) <- "double"
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  fx <- file.path(L1, sprintf("oof_extra_%s_%s.csv", fr, code))
  ff <- file.path(L1, sprintf("oof_forest_%s_%s.csv", fr, code))
  ok <- file.exists(fx) && file.exists(ff)
  chk(sprintf("%s: stage 6d and 6b files found", key),
      paste(basename(c(fx, ff)), collapse = ", "), if (ok) "PASS" else "FAIL")
  stop_if_fail()
  fmap <- list()
  for (s in seeds) {
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(f)) {
      chk(sprintf("%s seed %s: folds complete", key, s), sum(is.na(f)),
          "FAIL")
      stop_if_fail()
    }
    fmap[[as.character(s)]] <- as.integer(f)
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, X = X, t = sub$time_death,
      e = as.integer(sub$event_death == 1), fold = as.integer(f),
      hz = HZ[[fr]], off = 0.5, knots = KNOTS, minu = MINU,
      ckpt = file.path(CKPT, sprintf("%s_%s_%d.rds", fr, code, s)))
  }
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub,
                      seeds = seeds, fmap = fmap, fx = fx, ff = ff,
                      mods = c("ehr", "ecg", if (co == "3-mod") "ctpa"))
}
cat(sprintf("  %d level-1 tasks (frame x cohort x seed), %d EHR features\n",
            length(tasks), length(c28)))
done <- which(sapply(tasks, function(t) {
  if (!file.exists(t$ckpt)) return(FALSE)
  r <- tryCatch(readRDS(t$ckpt), error = function(e) NULL)
  !is.null(r) && all(is.finite(r$lp)) && all(r$info$ok)
}))
todo <- setdiff(seq_along(tasks), done)

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
if (length(todo)) {
  big <- todo[which.max(sapply(tasks[todo], function(t) length(t$t)))]
  t0 <- Sys.time()
  invisible(spline_task(tasks[[big]], folds_to_run = 0L, probe = TRUE))
  tp <- as.numeric(difftime(Sys.time(), t0, units = "secs")) * 5
  cat(sprintf(paste("  largest task about %.0f s; %d of %d tasks to fit on",
                    "%d workers: about %.1f min, plus the bootstrap\n"),
              tp, length(todo), length(tasks), workers,
              tp * ceiling(length(todo) / workers) * 1.2 / 60))
} else {
  cat("  every level-1 task already saved\n")
}
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
if (length(todo)) {
  ord <- todo[order(-sapply(tasks[todo], function(t) length(t$t)))]
  invisible(future_lapply(tasks[ord], spline_task, future.packages = PKG,
                          future.seed = TRUE, future.chunk.size = 1L))
}
res <- lapply(tasks, function(t) readRDS(t$ckpt))
INFO <- rbindlist(lapply(res, function(r) r$info[, `:=`(
  frame = r$fr, cohort = r$co, seed = r$seed)]))
nfail <- sum(!INFO$ok)
chk("spline fits that failed (first message)",
    sprintf("%d%s", nfail, if (nfail) paste0(" (", INFO$msg[!INFO$ok][1],
                                              ")") else ""),
    if (nfail == 0) "PASS" else "FAIL")
nc <- INFO$ok & !grepl("full convergence", INFO$conv) & nzchar(INFO$conv)
chk("spline fits short of full REML convergence (first status)",
    sprintf("%d of %d%s", sum(nc), nrow(INFO),
            if (any(nc)) paste0(" (", INFO$conv[nc][1], ")") else ""),
    if (mean(nc) <= 0.10) "PASS" else "WARN")
chk("spline terms per fit / median total edf",
    sprintf("%d / %.1f", max(INFO$n_smooth), median(INFO$edf, na.rm = TRUE)),
    "INFO")
stop_if_fail()

# ---------------- stacks and evaluation ----------------
cat("\n=== STACKS ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
JOBS <- list()
for (key in names(META)) {
  M <- META[[key]]
  sub <- M$sub
  n <- nrow(sub)
  hz <- HZ[[M$fr]]
  tt <- sub$time_death
  ee <- as.integer(sub$event_death == 1)
  XT <- fread(M$fx)
  FT <- fread(M$ff, select = c("seed", "hadm_id", "ehr.rsf.death.score"))
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  SETS <- list("EHR" = "ehr", "EHR+ECG" = c("ehr", "ecg"))
  if (M$co == "3-mod") SETS[["EHR+ECG+CTPA"]] <- c("ehr", "ecg", "ctpa")
  P <- list()
  nbad <- 0L
  dv <- c()
  for (s in M$seeds) {
    x <- XT[XT$seed == s][match(sub$hadm_id, hadm_id)]
    b <- FT[FT$seed == s][match(sub$hadm_id, hadm_id)]
    r <- Filter(function(z) z$seed == s, rr)[[1]]
    EH <- list("ridge-Cox retuned" = x[["ehr.cox2.death.lp"]],
               "spline Cox" = r$lp,
               "RSF" = b[["ehr.rsf.death.score"]],
               "oblique RSF" = x[["ehr.orsf.death.score"]])
    f <- M$fmap[[as.character(s)]]
    for (k in sort(unique(f))) {
      ii <- which(f == k)
      dv <- c(dv, cidx(tt[ii], ee[ii], r$lp[ii]))
    }
    for (fam in FAMS) {
      S <- list(ehr = EH[[fam]], ecg = x[["ecg.cox2.death.lp"]])
      if (M$co == "3-mod") S$ctpa <- x[["ctpa.cox2.death.lp"]]
      if (any(!is.finite(unlist(S)))) {
        nbad <- nbad + 1L
        next
      }
      for (sn in names(SETS)) {
        for (form in c("PH", "TV")) {
          kk <- paste(fam, sn, form, sep = " | ")
          if (is.null(P[[kk]])) P[[kk]] <- matrix(0, n, length(hz))
        }
        for (k in sort(unique(f))) {
          tr <- which(f != k)
          te <- which(f == k)
          sf <- stack_fit(S, SETS[[sn]], tr, te, tt, ee, hz + 0.5)
          for (form in c("PH", "TV")) {
            kk <- paste(fam, sn, form, sep = " | ")
            P[[kk]][te, ] <- P[[kk]][te, ] + sf[[form]] / length(M$seeds)
          }
        }
      }
    }
  }
  chk(sprintf("%s: seeds with missing scores", key), nbad,
      if (nbad == 0) "PASS" else "FAIL")
  chk(sprintf("%s: spline per-fold C below 0.40 (D2)", key),
      sprintf("mean %.3f", mean(dv, na.rm = TRUE)),
      if (mean(dv, na.rm = TRUE) >= 0.40) "PASS" else "FAIL")
  stop_if_fail()
  rk <- Reduce(`+`, lapply(rr, function(r) r$risk)) / length(rr)
  P[["spline Cox | EHR | own risks"]] <- rk
  rx <- XT[, lapply(.SD, mean), by = hadm_id,
           .SDcols = sprintf("ehr.cox2.death.r%d", hz)][
             match(sub$hadm_id, hadm_id)]
  P[["ridge-Cox retuned | EHR | own risks"]] <-
    as.matrix(rx[, sprintf("ehr.cox2.death.r%d", hz), with = FALSE])
  fwrite(data.table(hadm_id = sub$hadm_id,
                    spline_lp = Reduce(`+`, lapply(rr, function(r) r$lp)) /
                      length(rr)),
         file.path(OD, sprintf("spline_scores_%s_%s.csv", M$fr, M$code)))
  Y <- sapply(hz, function(h) as.integer(ee == 1 & tt <= h + 0.5))
  JOBS[[key]] <- list(fr = M$fr, co = M$co, hz = hz, P = P, n = n,
                      Y = matrix(Y, ncol = length(hz)),
                      rows_by = split(seq_len(n),
                                      match(sub$subject_id,
                                            unique(sub$subject_id))),
                      sets = names(SETS))
}

cat(sprintf("\n=== BOOTSTRAP (%d resamples) ===\n", NBOOT))
PAIRS <- list(c("spline Cox", "ridge-Cox retuned"),
              c("oblique RSF", "spline Cox"), c("RSF", "spline Cox"))
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
MET <- list()
CON <- list()
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(J$n))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  Frame %s %s done\n", J$fr, J$co))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, family = trimws(z[1]),
      set = trimws(z[2]), form = trimws(z[3]), horizon = as.numeric(z[4]),
      measure = z[5], estimate = pt[[k]], lo = q[1], hi = q[2])
  }
  for (k in names(J$P)) {
    z <- strsplit(k, " | ", fixed = TRUE)[[1]]
    for (hi in seq_along(J$hz)) {
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, family = z[1], set = z[2],
        form = z[3], horizon = J$hz[hi], measure = "slope", estimate = sl,
        lo = NA_real_, hi = NA_real_)
    }
  }
  forms <- list()
  for (sn in J$sets) for (form in c("PH", "TV")) {
    forms[[length(forms) + 1]] <- c(sn, form)
  }
  forms[[length(forms) + 1]] <- c("EHR", "own risks")
  for (sf in forms) for (pp in PAIRS) {
    for (h in J$hz) for (ms in c("auc", "ipa", "ll")) {
      ka <- paste(paste(pp[1], sf[1], sf[2], sep = " | "), h, ms, sep = "|")
      kb <- paste(paste(pp[2], sf[1], sf[2], sep = " | "), h, ms, sep = "|")
      if (!(ka %in% names(pt)) || !(kb %in% names(pt))) next
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, set = sf[1], form = sf[2],
        contrast = sprintf("%s - %s", pp[1], pp[2]), horizon = h,
        measure = ms, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "spline_metrics.csv"))
fwrite(cn, file.path(RS, "spline_contrasts.csv"))
fwrite(ck, file.path(RS, "spline_checks.csv"))
cat("\n=== EHR ALONE (PH STACK), AUROC ===\n")
print(dcast(mt[set == "EHR" & form == "PH" & measure == "auc"],
            frame + cohort + horizon ~ family, value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                       else x)], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (spline_scores_*.csv and spline_ckpt/, local only)\n")
for (f in c("spline_metrics.csv", "spline_contrasts.csv",
            "spline_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "spline_log.txt"), " (terminal log)\n")
cat("\nOutputs: spline_metrics.csv, spline_contrasts.csv, spline_checks.csv\n")