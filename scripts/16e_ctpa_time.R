# 16e_ctpa_time.R
# The CTPA's time pattern measured like for like across datasets and feature
# sets, and the tuning of its level-1 ridge-Cox
# (PREREGISTRATION_ADDENDUM_2.md, A15).
#
# Scores: ridge-Cox death scores, out of fold over the folds of 04 and 04b
# and every seed, fitted as 18a_early_linear.R fits the EHR alone, and
# seed-averaged. They cover Frame B (three-modality cohort, discharge clock)
# and Frame M (MIMIC, admission clock) with CTPA 46, CTPA 38 and EHR; Frame
# I (INSPECT, CTPA scan clock) with CTPA 38 and EHR; and Frame A
# (three-modality cohort) with CTPA 46, for part 2 only. CTPA 46 is all
# MIMIC CTPA features, including the device block; CTPA 38 is the features
# shared with INSPECT.
#
# Part 1, time pattern (Frames B, M, I; winsorised, standardised scores):
# per-SD hazard ratios for death in days 0-30 and 31-365 and their ratio, as
# in S5, for the CTPA alone and with the EHR score (both with early and late
# terms); incident AUROC as in A9 (share of patients at risk with a lower
# score at each death, kernel-smoothed on log time, bandwidth 0.4) and
# cumulative IPCW AUROC at 10, 30, 90, 180 and 365 days, with the CTPA-EHR
# gap and its change from 10 to 365 days; contrasts of the early/late ratio,
# 46 over 38 features in B and M (paired) and M over I with 38 features
# (independent). Part 2, CTPA tuning (Frames A and B, CTPA 46): these
# ridge-Cox risks against those of 06a (tuned on Harrell's C over five
# penalties) by AUROC, IPA, log-likelihood skill score, observed/expected
# and calibration slope. Intervals: 2,000-sample bootstrap, seed 42, scores
# held fixed (patient-clustered where a patient has several admissions).
#
# Checks:  stops if inputs or folds are missing or Frame B's EHR scores
#          differ from the EHR-alone PH risks of 18a by over 1e-6; warns if
#          over 20% of penalty choices are at the least-penalised end.
# Notes:   all scoring tasks run in the parallel workers, as in 18a, so that
#          the random inner folds match it.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, frame_M.csv,
#          frame_I.csv, folds.csv, folds_IM.csv, feature_whitelist.csv,
#          feature_whitelist_IM.csv, level1/oof_linear_A_3mod.csv,
#          level1/oof_linear_B_3mod.csv,
#          fusion_arch/arch_preds_linear_B_3mod.csv; data/inputs/ehr28.csv,
#          ctpa46.csv, inspect_ehr.csv, inspect/inspect_imp_features_v2.csv
# Outputs: results/ctpa_hr.csv, ctpa_incident.csv, ctpa_cumulative.csv,
#          ctpa_contrasts.csv, ctpa_retune.csv, ctpa_checks.csv,
#          ctpa_log.txt
# Usage:   Rscript scripts/16e_ctpa_time.R [--probe-only]

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
OD <- file.path(PR, "ctpa_time")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ5 <- c(10, 30, 90, 180, 365)
NLAM <- 80L
LMR <- 1e-7
BW <- 0.4
NBOOT <- 2000L
SEED <- 42L
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
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "ctpa_checks.csv"))
    cat("\nSTOP: a check failed. see: ctpa_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- scoring task (as stage 18a's EHR alone) ----------------
score_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
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
  lp <- rep(NA_real_, n)
  risk <- matrix(NA_real_, n, length(H))
  edges <- list()
  tt <- tk$t
  ee <- tk$e
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    X <- prep(tk$X, tr)
    cv <- cv.glmnet(X[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                    family = "cox", alpha = 0, foldid = inner,
                    nlambda = NLAM, lambda.min.ratio = LMR,
                    type.measure = "deviance",
                    penalty.factor = rep(1, ncol(X)), standardize = FALSE)
    K <- length(cv$lambda)
    i <- which(cv$lambda == cv$lambda.min)[1]
    edges[[length(edges) + 1]] <- data.table(least = i == K, most = i == 1)
    l_tr <- as.vector(predict(cv, X[tr, , drop = FALSE], s = "lambda.min"))
    l_te <- as.vector(predict(cv, X[te, , drop = FALSE], s = "lambda.min"))
    lp[te] <- l_te
    c0 <- max(l_tr)
    for (hi in seq_along(H)) {
      hs <- breslow(tt[tr], ee[tr], l_tr, H[hi] + tk$off)
      risk[te, hi] <- if (hs == 0) 0 else -expm1(-hs * exp(l_te - c0))
    }
  }
  list(unit = tk$unit, seed = tk$seed, lp = lp, risk = risk,
       edges = rbindlist(edges))
}

# ---------------- part 1 statistics ----------------
wz <- function(v) {
  q <- quantile(v, c(0.005, 0.995), names = FALSE)
  w <- pmin(pmax(v, q[1]), q[2])
  (w - mean(w)) / sd(w)
}
# early/late Cox models for every feature set of a frame, with patient
# weights w (1 for the point estimate, bootstrap counts otherwise)
hr_stats <- function(F, w) {
  v <- numeric(0)
  wr <- w[F$ep$id]
  for (fs in names(F$Z)) {
    zc <- F$Z[[fs]][F$ep$id]
    ze <- F$Z[["EHR"]][F$ep$id]
    if (fs == "EHR") next
    X1 <- cbind(zc * F$early, zc * (1 - F$early))
    X2 <- cbind(X1, ze * F$early, ze * (1 - F$early))
    for (mm in c("alone", "with EHR")) {
      X <- if (mm == "alone") X1 else X2
      k <- wr > 0
      b <- tryCatch(coef(coxph(Surv(F$ep$tstart[k], F$ep$t[k], F$ep$e[k]) ~
                                 X[k, , drop = FALSE], weights = wr[k],
                               ties = "efron")),
                    error = function(e) rep(NA_real_, ncol(X)))
      key <- paste(fs, mm, sep = "|")
      v[paste(key, "log HR early", sep = "|")] <- b[1]
      v[paste(key, "log HR late", sep = "|")] <- b[2]
      v[paste(key, "log ratio", sep = "|")] <- b[1] - b[2]
    }
  }
  v
}
inc_stats <- function(F, w) {
  v <- numeric(0)
  den <- as.vector(F$A %*% w)
  wd <- w[F$d]
  ok <- den > 0
  for (sc in names(F$M)) {
    ci <- as.vector(F$M[[sc]] %*% w) / pmax(den, 1e-12)
    for (h in HZ5) {
      kw <- wd * F$K[, as.character(h)] * ok
      v[paste(sc, h, sep = "|")] <- sum(kw * ci) / sum(kw)
    }
  }
  v
}
frame_stats <- function(F, w) c(hr_stats(F, w), inc_stats(F, w))
boot_frame <- function(i, Fx) {
  idx <- unlist(Fx$rows_by[sample.int(length(Fx$rows_by),
                                      length(Fx$rows_by), replace = TRUE)],
                use.names = FALSE)
  frame_stats(Fx, tabulate(idx, Fx$n))
}
cum_auc <- function(t, e, s, h) {
  ut <- sort(unique(t))
  nr <- length(t) - findInterval(ut, sort(t), left.open = TRUE)
  nc <- tabulate(match(t[e == 0], ut), length(ut))
  G <- cumprod(1 - nc / nr)
  g_at <- function(x, left) c(1, G)[findInterval(x, ut, left.open = left) + 1]
  case <- e == 1 & t <= h
  ctrl <- t > h
  wc <- 1 / g_at(t[case], TRUE)
  rc <- sort(s[ctrl])
  lo <- findInterval(s[case], rc, left.open = TRUE)
  up <- findInterval(s[case], rc)
  sum(wc * (lo + 0.5 * (up - lo)) / length(rc)) / sum(wc)
}

# ---------------- part 2 statistics ----------------
bin_stats <- function(R, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(R$hz)) {
    y <- R$Y[idx, hi]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    for (k in names(R$P)) {
      p <- pmin(pmax(R$P[[k]][idx, hi], 1e-6), 1 - 1e-6)
      key <- paste(k, R$hz[hi], sep = "|")
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
boot_bin <- function(i, R) {
  idx <- unlist(R$rows_by[sample.int(length(R$rows_by),
                                     length(R$rows_by), replace = TRUE)],
                use.names = FALSE)
  bin_stats(R, idx)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PR, c("frame_A.csv", "frame_B.csv", "frame_M.csv",
                          "frame_I.csv", "folds.csv", "folds_IM.csv",
                          "feature_whitelist.csv",
                          "feature_whitelist_IM.csv",
                          "level1/oof_linear_A_3mod.csv",
                          "level1/oof_linear_B_3mod.csv",
                          "fusion_arch/arch_preds_linear_B_3mod.csv")),
          file.path(IN, c("ehr28.csv", "ctpa46.csv", "inspect_ehr.csv",
                          "inspect/inspect_imp_features_v2.csv")))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
wl <- fread(file.path(PR, "feature_whitelist.csv"))
wli <- fread(file.path(PR, "feature_whitelist_IM.csv"))
me <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
mc <- unique(fread(file.path(IN, "ctpa46.csv")), by = "hadm_id")
ie <- unique(fread(file.path(IN, "inspect_ehr.csv")), by = "person_id")
ic <- unique(fread(file.path(IN, "inspect", "inspect_imp_features_v2.csv")),
             by = "person_id")
c46 <- intersect(wl[modality == "ctpa", column], names(mc))
c38 <- intersect(wli[modality == "ctpa", column], names(mc))
e28 <- intersect(wl[modality == "ehr", column], names(me))
e28i <- intersect(wli[modality == "ehr", column], names(ie))
chk("features: CTPA 46 / CTPA 38 / EHR (MIMIC) / EHR (INSPECT)",
    sprintf("%d / %d / %d / %d", length(c46), length(c38), length(e28),
            length(e28i)), "INFO")
mat <- function(dt, key, ids, cols) {
  m <- as.matrix(dt[match(ids, dt[[key]]), ..cols])
  storage.mode(m) <- "double"
  m
}
fo <- fread(file.path(PR, "folds.csv"))
foi <- fread(file.path(PR, "folds_IM.csv"))
fA <- fread(file.path(PR, "frame_A.csv"))[in_three_mod == 1]
fB <- fread(file.path(PR, "frame_B.csv"))[in_three_mod == 1]
fM <- fread(file.path(PR, "frame_M.csv"))
fI <- fread(file.path(PR, "frame_I.csv"))
FR <- list(
  A = list(id = fA$hadm_id, subj = fA$subject_id, t = fA$time_death,
           e = as.integer(fA$event_death == 1), off = 0.5, hz = c(10, 30),
           cut = 30.5,
           folds = function(s) fo[frame == "A" & cohort == "3-mod" &
                                    seed == s][match(fA$hadm_id, hadm_id),
                                               fold],
           seeds = sort(unique(fo[frame == "A" & cohort == "3-mod", seed]))),
  B = list(id = fB$hadm_id, subj = fB$subject_id, t = fB$time_death,
           e = as.integer(fB$event_death == 1), off = 0.5, hz = HZ5,
           cut = 30.5,
           folds = function(s) fo[frame == "B" & cohort == "3-mod" &
                                    seed == s][match(fB$hadm_id, hadm_id),
                                               fold],
           seeds = sort(unique(fo[frame == "B" & cohort == "3-mod", seed]))),
  M = list(id = fM$hadm_id, subj = fM$hadm_id, t = fM$time_death,
           e = as.integer(fM$event_death == 1), off = 0.5, hz = HZ5,
           cut = 30.5,
           folds = function(s) foi[frame == "M" & seed == s][
             match(fM$hadm_id, id), fold],
           seeds = sort(unique(foi[frame == "M", seed]))),
  I = list(id = fI$person_id, subj = fI$person_id,
           t = pmax(fI$time_death, 1e-3),
           e = as.integer(fI$event_death == 1), off = 0, hz = HZ5,
           cut = 30,
           folds = function(s) foi[frame == "I" & seed == s][
             match(fI$person_id, id), fold],
           seeds = sort(unique(foi[frame == "I", seed]))))
UNITS <- list(
  "A|CTPA 46" = mat(mc, "hadm_id", fA$hadm_id, c46),
  "B|CTPA 46" = mat(mc, "hadm_id", fB$hadm_id, c46),
  "B|CTPA 38" = mat(mc, "hadm_id", fB$hadm_id, c38),
  "B|EHR" = mat(me, "hadm_id", fB$hadm_id, e28),
  "M|CTPA 46" = mat(mc, "hadm_id", fM$hadm_id, c46),
  "M|CTPA 38" = mat(mc, "hadm_id", fM$hadm_id, c38),
  "M|EHR" = mat(me, "hadm_id", fM$hadm_id, e28i),
  "I|CTPA 38" = mat(ic, "person_id", fI$person_id, c38),
  "I|EHR" = mat(ie, "person_id", fI$person_id, e28i))
tasks <- list()
for (u in names(UNITS)) {
  f <- strsplit(u, "|", fixed = TRUE)[[1]][1]
  Fr <- FR[[f]]
  nrec <- sum(rowSums(!is.na(UNITS[[u]])) == 0)
  chk(sprintf("%s: patients / features / rows with no record", u),
      sprintf("%d / %d / %d", nrow(UNITS[[u]]), ncol(UNITS[[u]]), nrec),
      "INFO")
  for (s in Fr$seeds) {
    fs <- Fr$folds(s)
    if (anyNA(fs)) {
      chk(sprintf("%s seed %s: folds complete", u, s), sum(is.na(fs)),
          "FAIL")
      stop_if_fail()
    }
    tasks[[length(tasks) + 1]] <- list(
      unit = u, seed = s, fold = fs, subj = Fr$subj, t = Fr$t, e = Fr$e,
      X = UNITS[[u]], hz = Fr$hz, off = Fr$off)
  }
}
cat(sprintf("  %d scoring tasks (unit x seed)\n", length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) nrow(t$X) * ncol(t$X)))
t0 <- Sys.time()
invisible(score_task(tasks[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  largest scoring task %.0f s; %d tasks on %d workers:",
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
res <- future_lapply(tasks, score_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)

# ---------------- seed-averaged scores and checks ----------------
cat("\n=== SCORES ===\n")
SC <- list()
for (u in names(UNITS)) {
  rr <- Filter(function(r) r$unit == u, res)
  SC[[u]] <- list(
    lp = Reduce(`+`, lapply(rr, function(r) r$lp)) / length(rr),
    risk = Reduce(`+`, lapply(rr, function(r) r$risk)) / length(rr))
  ed <- rbindlist(lapply(rr, function(r) r$edges))
  chk(sprintf("%s: penalty choices at the least / most penalised end", u),
      sprintf("%.1f%% / %.1f%%", 100 * mean(ed$least), 100 * mean(ed$most)),
      if (mean(ed$least) <= 0.20) "PASS" else "WARN")
}
a18 <- fread(file.path(PR, "fusion_arch", "arch_preds_linear_B_3mod.csv"))[
  model == "EHR alone PH" & set == "EHR"]
dm <- max(sapply(seq_along(HZ5), function(hi) {
  x <- a18[horizon == HZ5[hi]]
  max(abs(x$risk[match(fB$hadm_id, x$hadm_id)] - SC[["B|EHR"]]$risk[, hi]))
}))
chk("Frame B EHR risks reproduce stage 18a's EHR alone PH (max |diff|)",
    sprintf("%.1e", dm), if (is.finite(dm) && dm < 1e-6) "PASS" else "FAIL")
stop_if_fail()
fwrite(rbindlist(lapply(names(SC), function(u) {
  z <- strsplit(u, "|", fixed = TRUE)[[1]]
  data.table(frame = z[1], features = z[2], id = FR[[z[1]]]$id,
             score = SC[[u]]$lp)
})), file.path(OD, "ctpa_scores.csv"))

# ---------------- part 1 ----------------
cat(sprintf("\n=== PART 1: TIME PATTERN (%d resamples) ===\n", NBOOT))
PF <- list()
for (f in c("B", "M", "I")) {
  Fr <- FR[[f]]
  us <- grep(paste0("^", f, "\\|"), names(SC), value = TRUE)
  Z <- setNames(lapply(us, function(u) wz(SC[[u]]$lp)),
                sub("^.\\|", "", us))
  ep <- survSplit(Surv(t, e) ~ ., cut = Fr$cut, episode = "per",
                  data = data.frame(t = Fr$t, e = Fr$e,
                                    id = seq_along(Fr$t)))
  d <- which(Fr$e == 1 & Fr$t <= 365 + Fr$off)
  A <- outer(Fr$t[d], Fr$t, "<") * 1
  M <- lapply(Z, function(z) {
    A * (outer(z[d], z, ">") + 0.5 * outer(z[d], z, "=="))
  })
  K <- sapply(HZ5, function(h) exp(-0.5 * ((log(Fr$t[d]) - log(h)) / BW)^2))
  K <- matrix(K, ncol = length(HZ5), dimnames = list(NULL, HZ5))
  PF[[f]] <- list(n = length(Fr$t), Z = Z, early = as.numeric(ep$per == 1),
                  ep = list(tstart = ep$tstart, t = ep$t, e = ep$e,
                            id = ep$id),
                  d = d, A = A, M = M, K = K,
                  rows_by = split(seq_along(Fr$t),
                                  match(Fr$subj, unique(Fr$subj))))
}
P1 <- list()
for (f in names(PF)) {
  pt <- frame_stats(PF[[f]], rep(1, PF[[f]]$n))
  if (f == "B") {
    t0 <- Sys.time()
    invisible(lapply(1:5, boot_frame, Fx = PF[[f]]))
    cat(sprintf("  bootstrap: %.2f s per resample in Frame B\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5))
  }
  Bm <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_frame,
                                     Fx = PF[[f]], future.seed = SEED,
                                     future.packages = PKG))
  P1[[f]] <- list(pt = pt, B = Bm)
  cat(sprintf("  Frame %s done\n", f))
}
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
HRT <- list()
INC <- list()
for (f in names(P1)) {
  pt <- P1[[f]]$pt
  Bm <- P1[[f]]$B
  for (k in grep("log ", names(pt), value = TRUE)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(Bm[, k])
    HRT[[length(HRT) + 1]] <- data.table(
      frame = f, features = z[1], model = z[2],
      quantity = sub("log ", "", z[3]), estimate = exp(pt[[k]]),
      lo = exp(q[1]), hi = exp(q[2]))
  }
  fs <- setdiff(names(PF[[f]]$Z), "EHR")
  for (sc in names(PF[[f]]$Z)) {
    for (h in HZ5) {
      k <- paste(sc, h, sep = "|")
      q <- q2(Bm[, k])
      INC[[length(INC) + 1]] <- data.table(
        frame = f, score = sc, horizon = as.character(h),
        estimate = pt[[k]], lo = q[1], hi = q[2])
    }
    a <- paste(sc, 365, sep = "|")
    b <- paste(sc, 10, sep = "|")
    q <- q2(Bm[, a] - Bm[, b])
    INC[[length(INC) + 1]] <- data.table(
      frame = f, score = sc, horizon = "change, 365 minus 10",
      estimate = pt[[a]] - pt[[b]], lo = q[1], hi = q[2])
  }
  for (c1 in fs) {
    for (h in HZ5) {
      a <- paste("EHR", h, sep = "|")
      b <- paste(c1, h, sep = "|")
      q <- q2(Bm[, a] - Bm[, b])
      INC[[length(INC) + 1]] <- data.table(
        frame = f, score = sprintf("gap, EHR minus %s", c1),
        horizon = as.character(h), estimate = pt[[a]] - pt[[b]],
        lo = q[1], hi = q[2])
    }
    gd <- function(X, h) X[[paste("EHR", h, sep = "|")]] -
      X[[paste(c1, h, sep = "|")]]
    est <- gd(as.list(pt), 365) - gd(as.list(pt), 10)
    bv <- (Bm[, paste("EHR", 365, sep = "|")] -
             Bm[, paste(c1, 365, sep = "|")]) -
      (Bm[, paste("EHR", 10, sep = "|")] - Bm[, paste(c1, 10, sep = "|")])
    q <- q2(bv)
    INC[[length(INC) + 1]] <- data.table(
      frame = f, score = sprintf("gap, EHR minus %s", c1),
      horizon = "change, 365 minus 10", estimate = est, lo = q[1],
      hi = q[2])
  }
}
CTR <- list()
lr <- function(f, fs, mm) paste(fs, mm, "log ratio", sep = "|")
for (mm in c("alone", "with EHR")) {
  for (f in c("B", "M")) {
    a <- lr(f, "CTPA 46", mm)
    b <- lr(f, "CTPA 38", mm)
    est <- P1[[f]]$pt[[a]] - P1[[f]]$pt[[b]]
    q <- q2(P1[[f]]$B[, a] - P1[[f]]$B[, b])
    CTR[[length(CTR) + 1]] <- data.table(
      contrast = sprintf("Frame %s: 46 over 38 features", f), model = mm,
      ratio_of_ratios = exp(est), lo = exp(q[1]), hi = exp(q[2]),
      excludes_one = q[1] > 0 | q[2] < 0)
  }
  a <- lr("M", "CTPA 38", mm)
  est <- P1$M$pt[[a]] - P1$I$pt[[a]]
  q <- q2(P1$M$B[, a] - P1$I$B[, a])
  CTR[[length(CTR) + 1]] <- data.table(
    contrast = "MIMIC (Frame M) over INSPECT, 38 features", model = mm,
    ratio_of_ratios = exp(est), lo = exp(q[1]), hi = exp(q[2]),
    excludes_one = q[1] > 0 | q[2] < 0)
}
CUM <- list()
for (f in names(PF)) {
  Fr <- FR[[f]]
  for (sc in names(PF[[f]]$Z)) for (h in HZ5) {
    CUM[[length(CUM) + 1]] <- data.table(
      frame = f, score = sc, horizon = h,
      cumulative_auc = cum_auc(Fr$t, Fr$e, PF[[f]]$Z[[sc]], h + Fr$off))
  }
}

# ---------------- part 2 ----------------
cat(sprintf("\n=== PART 2: CTPA TUNING (%d resamples) ===\n", NBOOT))
RT <- list()
for (f in c("A", "B")) {
  Fr <- FR[[f]]
  sub <- if (f == "A") fA else fB
  cols <- sprintf("ctpa.cox.death.r%d", Fr$hz)
  L6 <- fread(file.path(PR, "level1", sprintf("oof_linear_%s_3mod.csv", f)),
              select = c("hadm_id", cols))
  L6 <- L6[, lapply(.SD, mean), by = hadm_id][match(sub$hadm_id, hadm_id)]
  P6 <- as.matrix(L6[, ..cols])
  chk(sprintf("Frame %s: stage 6a CTPA risks found / missing", f),
      sprintf("%d / %d", nrow(P6), sum(is.na(P6))),
      if (!anyNA(P6)) "PASS" else "FAIL")
  stop_if_fail()
  R <- list(hz = Fr$hz,
            P = list("this stage (D11 path, deviance)" =
                       SC[[paste(f, "CTPA 46", sep = "|")]]$risk,
                     "stage 6a (C, five penalties)" = P6),
            Y = matrix(sapply(Fr$hz, function(h) {
              as.integer(Fr$e == 1 & Fr$t <= h + Fr$off)
            }), ncol = length(Fr$hz)),
            rows_by = split(seq_along(Fr$t),
                            match(Fr$subj, unique(Fr$subj))))
  pt <- bin_stats(R, seq_along(Fr$t))
  Bm <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_bin, R = R,
                                     future.seed = SEED,
                                     future.packages = PKG))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(Bm[, k])
    RT[[length(RT) + 1]] <- data.table(
      frame = f, model = z[1], horizon = as.numeric(z[2]),
      measure = z[3], estimate = pt[[k]], lo = q[1], hi = q[2])
  }
  for (hi in seq_along(Fr$hz)) for (ms in c("auc", "ipa", "ll")) {
    a <- paste(names(R$P)[1], Fr$hz[hi], ms, sep = "|")
    b <- paste(names(R$P)[2], Fr$hz[hi], ms, sep = "|")
    q <- q2(Bm[, a] - Bm[, b])
    RT[[length(RT) + 1]] <- data.table(
      frame = f, model = "difference, this stage minus stage 6a",
      horizon = Fr$hz[hi], measure = ms, estimate = pt[[a]] - pt[[b]],
      lo = q[1], hi = q[2])
  }
  for (k in names(R$P)) for (hi in seq_along(Fr$hz)) {
    p <- pmin(pmax(R$P[[k]][, hi], 1e-6), 1 - 1e-6)
    sl <- tryCatch(unname(coef(glm(R$Y[, hi] ~ qlogis(p),
                                   family = binomial))[2]),
                   error = function(e) NA_real_)
    RT[[length(RT) + 1]] <- data.table(
      frame = f, model = k, horizon = Fr$hz[hi], measure = "slope",
      estimate = sl, lo = NA_real_, hi = NA_real_)
  }
  cat(sprintf("  Frame %s done\n", f))
}
plan(sequential)

# ---------------- write ----------------
hr <- rbindlist(HRT)
inc <- rbindlist(INC)
ctr <- rbindlist(CTR)
cum <- rbindlist(CUM)
rt <- rbindlist(RT)
ck <- rbindlist(CHK)
fwrite(hr, file.path(RS, "ctpa_hr.csv"))
fwrite(inc, file.path(RS, "ctpa_incident.csv"))
fwrite(cum, file.path(RS, "ctpa_cumulative.csv"))
fwrite(ctr, file.path(RS, "ctpa_contrasts.csv"))
fwrite(rt, file.path(RS, "ctpa_retune.csv"))
fwrite(ck, file.path(RS, "ctpa_checks.csv"))
cat("\n=== EARLY / LATE HAZARD RATIOS (per SD) ===\n")
print(dcast(hr, frame + features + model ~ quantity,
            value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                       else x)], row.names = FALSE)
cat("\n=== RATIO CONTRASTS ===\n")
print(ctr[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)],
      row.names = FALSE)
cat("\n=== INCIDENT AUROC ===\n")
print(dcast(inc[!grepl("^gap", score)], frame + score ~ horizon,
            value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3)
                       else x)], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(OD, "ctpa_scores.csv"), " (patient-level, local only)\n")
for (f in c("ctpa_hr.csv", "ctpa_incident.csv", "ctpa_cumulative.csv",
            "ctpa_contrasts.csv", "ctpa_retune.csv", "ctpa_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "ctpa_log.txt"), " (terminal log)\n")
cat("\nOutputs: ctpa_hr.csv, ctpa_incident.csv, ctpa_cumulative.csv,",
    "ctpa_contrasts.csv, ctpa_retune.csv, ctpa_checks.csv\n")