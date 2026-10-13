# 16b_transfer_s6.R
# INSPECT <-> MIMIC transfer read through the Debray framework, with a
# missing-data sensitivity analysis (PREREGISTRATION.md S6; DEVIATIONS.md
# D6; PREREGISTRATION_ADDENDUM_2.md A7, A8).
#
# Version 2 (D12): D11 penalty path, source-percentile winsorising and a
# record of penalty choices at the ends of the path. The first run is
# archived in results/s6_v1; the script stops if it is still in results/.
#
# Frames: INSPECT (I), 3,300 PE-positive persons, clock at the index CTPA,
# death within 366 days, real censoring; MIMIC-IV (M), 1,490 first
# admissions with a CTPA, clock at admission, complete follow-up. Features:
# the 28 EHR and 38 CTPA features shared by both. Models: ridge-Cox per
# modality (glmnet, alpha 0, lambda by cross-validated C, 80 penalties down
# to 1e-7 of the largest) and a Cox stack on winsorised modality scores
# (EHR, CTPA, EHR + CTPA) fitted on out-of-fold scores, with Breslow risks
# at 30 and 365 days and six seeds averaged. Imputation follows A8, then
# non-binary features are winsorised at the source's 0.5th and 99.5th
# percentiles for source and target. Internal: 5-fold CV within each dataset
# (folds from 04b). Transfer: trained on one dataset, applied to the other.
#
# Step 1, relatedness: membership-model AUROC (ridge logistic, out of fold),
# feature SMDs and missingness, shift in the transferred score, Kaplan-Meier
# mortality. Step 2, performance in the target at 30 and 365 days: IPCW AUC,
# Brier and IPA with influence-function intervals and paired contrasts
# (riskRegression::Score), Uno's C, KM observed / expected, calibration
# slope, decile tables. Step 3, updating of the fused model by 5-fold CV in
# the target: baseline re-estimated (score as offset), baseline and slope
# re-estimated, full refit. Censoring bias (INSPECT): observed risk and
# AUROC when censored patients are dropped or counted as survivors, against
# KM and IPCW. Missing data (A8): source medians (main), plus missing
# indicators, regression imputation fitted on the source without the
# outcome, or the target's own medians; every transfer runs under all four.
#
# Checks:  warns if over 20% of penalty choices hit the least-penalised end.
#
# Inputs:  data/processed/frame_I.csv, frame_M.csv, folds_IM.csv,
#          feature_whitelist_IM.csv; data/inputs/inspect_ehr.csv,
#          inspect/inspect_imp_features_v2.csv, ehr28.csv, ctpa46.csv
# Outputs: results/s6_performance.csv, s6_contrasts.csv, s6_relatedness.csv,
#          s6_features.csv, s6_censoring.csv, s6_deciles.csv, s6_lambda.csv,
#          s6_checks.csv, s6_log.txt
# Usage:   Rscript scripts/16b_transfer_s6.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(prodlim)
  library(riskRegression)
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
HZ <- c(30, 365)
STRATS <- c("median (source)", "median + missing indicators",
            "regression imputation (source)", "median (target's own)")
MAIN <- STRATS[1]
STACKS <- c("EHR", "CTPA", "EHR + CTPA")
NLAM <- 80L
LMR <- 1e-7
PKG <- c("survival", "glmnet", "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-58s %s\n", status, item, value))
}
auc1 <- function(y, s) {
  y <- as.logical(y)
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---------------- modelling functions (sent to workers) ----------------
impute <- function(Xs, Xt, strategy) {
  meds <- apply(Xs, 2, median, na.rm = TRUE)
  meds[!is.finite(meds)] <- 0
  fill <- function(X, m) {
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- m[na[, 2]]
    X
  }
  if (strategy == "median (source)") {
    return(list(s = fill(Xs, meds), t = fill(Xt, meds)))
  }
  if (strategy == "median (target's own)") {
    mt <- apply(Xt, 2, median, na.rm = TRUE)
    mt[!is.finite(mt)] <- meds[!is.finite(mt)]
    return(list(s = fill(Xs, meds), t = fill(Xt, mt)))
  }
  miss <- which(colSums(is.na(Xs)) > 0 | colSums(is.na(Xt)) > 0)
  if (strategy == "median + missing indicators") {
    Is <- 1 * is.na(Xs[, miss, drop = FALSE])
    It <- 1 * is.na(Xt[, miss, drop = FALSE])
    colnames(Is) <- colnames(It) <- paste0("missing_", colnames(Xs)[miss])
    return(list(s = cbind(fill(Xs, meds), Is), t = cbind(fill(Xt, meds), It)))
  }
  S0 <- fill(Xs, meds)
  T0 <- fill(Xt, meds)
  Ss <- S0
  Tt <- T0
  for (j in miss) {
    obs <- !is.na(Xs[, j])
    if (sum(obs) < 30) next
    b <- lm.fit(cbind(1, S0[obs, -j, drop = FALSE]), Xs[obs, j])$coefficients
    b[is.na(b)] <- 0
    ms <- is.na(Xs[, j])
    mt <- is.na(Xt[, j])
    if (any(ms)) Ss[ms, j] <- cbind(1, S0[ms, -j, drop = FALSE]) %*% b
    if (any(mt)) Tt[mt, j] <- cbind(1, T0[mt, -j, drop = FALSE]) %*% b
  }
  list(s = Ss, t = Tt)
}
# non-binary features winsorised at the source's 0.5th and 99.5th
# percentiles, source and target alike (D12)
winsorise <- function(imp) {
  for (j in seq_len(ncol(imp$s))) {
    if (length(unique(imp$s[, j])) > 2) {
      q <- quantile(imp$s[, j], c(0.005, 0.995), names = FALSE)
      imp$s[, j] <- pmin(pmax(imp$s[, j], q[1]), q[2])
      imp$t[, j] <- pmin(pmax(imp$t[, j], q[1]), q[2])
    }
  }
  imp
}
fit_predict <- function(Xtr, ttr, etr, ftr, Xte, cuts, seed, groups) {
  make_wz <- function(v) {
    q <- quantile(v, c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    list(q = q, mu = mean(w), sd = sd(w))
  }
  apply_wz <- function(p, v) (pmin(pmax(v, p$q[1]), p$q[2]) - p$mu) / p$sd
  set.seed(seed)
  lp_oof <- list()
  lp_te <- list()
  edges <- list()
  rec <- function(m, fit, cv) {
    K <- length(cv$lambda)
    i <- which(cv$lambda == cv$lambda.min)[1]
    edges[[length(edges) + 1]] <<- data.table::data.table(
      group = m, fit = fit, nlam = K, least = i == K, most = i == 1)
  }
  for (m in names(groups)) {
    cols <- intersect(groups[[m]], colnames(Xtr))
    A <- Xtr[, cols, drop = FALSE]
    B <- Xte[, cols, drop = FALSE]
    oof <- numeric(nrow(A))
    for (k in sort(unique(ftr))) {
      tr <- ftr != k
      inner <- sample(rep_len(1:3, sum(tr)))
      cv <- glmnet::cv.glmnet(A[tr, , drop = FALSE],
                              survival::Surv(ttr[tr], etr[tr]),
                              family = "cox", alpha = 0, foldid = inner,
                              type.measure = "C", nlambda = NLAM,
                              lambda.min.ratio = LMR)
      rec(m, "inner", cv)
      oof[!tr] <- as.vector(predict(cv, A[!tr, , drop = FALSE],
                                    s = "lambda.min", type = "link"))
    }
    cvf <- glmnet::cv.glmnet(A, survival::Surv(ttr, etr), family = "cox",
                             alpha = 0, foldid = ftr, type.measure = "C",
                             nlambda = NLAM, lambda.min.ratio = LMR)
    rec(m, "full", cvf)
    lp_oof[[m]] <- oof
    lp_te[[m]] <- as.vector(predict(cvf, B, s = "lambda.min",
                                    type = "link"))
  }
  stacks <- list("EHR" = "EHR", "CTPA" = "CTPA",
                 "EHR + CTPA" = c("EHR", "CTPA"))
  out <- list()
  for (sn in names(stacks)) {
    ms <- stacks[[sn]]
    P <- lapply(ms, function(m) make_wz(lp_oof[[m]]))
    Ztr <- matrix(sapply(seq_along(ms), function(i)
      apply_wz(P[[i]], lp_oof[[ms[i]]])), ncol = length(ms))
    Zte <- matrix(sapply(seq_along(ms), function(i)
      apply_wz(P[[i]], lp_te[[ms[i]]])), ncol = length(ms))
    f <- survival::coxph(survival::Surv(ttr, etr) ~ Ztr, ties = "efron")
    bh <- survival::basehaz(f, centered = FALSE)
    lp2 <- as.vector(Zte %*% coef(f))
    R <- sapply(cuts, function(ct) {
      i <- findInterval(ct, bh$time)
      H0 <- if (i == 0) 0 else bh$hazard[i]
      1 - exp(-H0 * exp(lp2))
    })
    out[[sn]] <- list(lp = lp2, risk = matrix(R, ncol = length(cuts)))
  }
  list(stacks = out, edges = data.table::rbindlist(edges))
}
run_task <- function(tk) {
  res <- tryCatch({
    imp <- winsorise(impute(tk$Xs, tk$Xt, tk$strategy))
    Xte <- if (tk$type == "transfer") rbind(imp$t, imp$s) else imp$t
    fp <- fit_predict(imp$s, tk$ts, tk$es, tk$ftr, Xte, tk$cuts, tk$seed,
                      tk$groups)
    nt <- nrow(imp$t)
    tgt <- lapply(fp$stacks, function(z) list(
      lp = z$lp[seq_len(nt)],
      risk = z$risk[seq_len(nt), , drop = FALSE]))
    src_lp <- if (tk$type == "transfer") {
      fp$stacks[["EHR + CTPA"]]$lp[-seq_len(nt)]
    } else NULL
    list(ok = TRUE, tgt = tgt, src_lp = src_lp, edges = fp$edges, msg = "")
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  c(list(type = tk$type, key = tk$key, seed = tk$seed,
         strategy = tk$strategy, rows = tk$rows), res)
}

# ---------------- evaluation functions (main process) ----------------
breslow_H <- function(tt, ee, w, cut) {
  o <- order(tt)
  tt <- tt[o]
  ee <- ee[o]
  w <- w[o]
  rs <- rev(cumsum(rev(w)))
  ut <- unique(tt[ee == 1 & tt <= cut])
  if (!length(ut)) return(0)
  d <- sapply(ut, function(u) sum(ee[tt == u]))
  sum(d / rs[match(ut, tt)])
}
eval_target <- function(tt, ee, cut, ymax, R) {
  df <- data.frame(time = tt, status = as.integer(ee))
  sc <- Score(R, formula = Hist(time, status) ~ 1, data = df, times = cut,
              metrics = c("auc", "brier"), summary = "ipa",
              cens.model = "km", null.model = TRUE, conf.int = TRUE,
              contrasts = TRUE)
  au <- as.data.table(sc$AUC$score)
  br <- as.data.table(sc$Brier$score)
  sf <- summary(survfit(Surv(time, status) ~ 1, data = df), times = cut,
                extend = TRUE)
  obs <- 1 - sf$surv
  th <- pmin(tt, cut)
  eh <- as.integer(ee == 1 & tt <= cut)
  rows <- list()
  for (m in names(R)) {
    r <- pmin(pmax(R[[m]], 1e-6), 1 - 1e-6)
    cll <- log(-log(1 - r))
    f <- tryCatch(coxph(Surv(th, eh) ~ cll), error = function(e) NULL)
    sl <- if (is.null(f)) rep(NA_real_, 3) else {
      b <- coef(f)[[1]]
      s <- sqrt(vcov(f)[1, 1])
      c(b, b - 1.96 * s, b + 1.96 * s)
    }
    uc <- tryCatch(as.numeric(concordance(
      Surv(tt, ee) ~ r, reverse = TRUE, timewt = "n/G2",
      ymax = ymax)$concordance), error = function(e) NA_real_)
    a <- au[as.character(model) == m]
    b <- br[as.character(model) == m]
    rows[[m]] <- data.table(
      model = m, n = length(tt), events_by_h = sum(eh), km_risk = obs,
      mean_pred = mean(r), oe = obs / mean(r),
      oe_lo = (1 - sf$upper) / mean(r), oe_hi = (1 - sf$lower) / mean(r),
      slope = sl[1], slope_lo = sl[2], slope_hi = sl[3], uno_c = uc,
      auc = a$AUC, auc_lo = a$lower, auc_hi = a$upper, brier = b$Brier,
      brier_lo = b$lower, brier_hi = b$upper, ipa = b$IPA)
  }
  list(perf = rbindlist(rows), ca = as.data.table(sc$AUC$contrasts),
       cb = as.data.table(sc$Brier$contrasts))
}
dec_tab <- function(tt, ee, cut, r) {
  g <- ceiling(10 * frank(r, ties.method = "first") / length(r))
  rbindlist(lapply(1:10, function(k) {
    ii <- g == k
    sf <- summary(survfit(Surv(tt[ii], ee[ii]) ~ 1), times = cut,
                  extend = TRUE)
    data.table(decile = k, n = sum(ii), mean_pred = mean(r[ii]),
               km_risk = 1 - sf$surv, km_lo = 1 - sf$upper,
               km_hi = 1 - sf$lower)
  }))
}
member_auc <- function(X, y) {
  set.seed(42)
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 5,
                  keep = TRUE, type.measure = "auc", nlambda = NLAM,
                  lambda.min.ratio = LMR)
  auc1(y == 1, cv$fit.preval[, match(cv$lambda.min, cv$lambda)])
}

# ---------------- data ----------------
cat("=== DATA ===\n")
if (file.exists(file.path(RS, "s6_performance.csv"))) {
  cat("STOP: results/s6_performance.csv exists. Archive the first S6 run",
      "to results/s6_v1 before rerunning (D12).\n")
  quit(status = 1)
}
fi <- fread(file.path(PR, "frame_I.csv"))
fm <- fread(file.path(PR, "frame_M.csv"))
wl <- fread(file.path(PR, "feature_whitelist_IM.csv"))
fo <- fread(file.path(PR, "folds_IM.csv"))
ehr_c <- wl[modality == "ehr", column]
ct_c <- wl[modality == "ctpa", column]
ie <- fread(file.path(IN, "inspect_ehr.csv"))
ic <- fread(file.path(IN, "inspect", "inspect_imp_features_v2.csv"))
me <- fread(file.path(IN, "ehr28.csv"))
mc <- fread(file.path(IN, "ctpa46.csv"))
chk("INSPECT EHR: extra rows for the same person",
    sum(duplicated(ie$person_id)),
    if (anyDuplicated(ie$person_id)) "WARN" else "PASS")
chk("INSPECT CTPA: extra rows for the same person",
    sum(duplicated(ic$person_id)),
    if (anyDuplicated(ic$person_id)) "WARN" else "PASS")
ie <- unique(ie, by = "person_id")
ic <- unique(ic, by = "person_id")
me <- unique(me, by = "hadm_id")
mc <- unique(mc, by = "hadm_id")
num <- function(dt, key, ids, cols) {
  m <- as.matrix(dt[match(ids, dt[[key]]), ..cols])
  storage.mode(m) <- "double"
  m
}
XI <- cbind(num(ie, "person_id", fi$person_id, ehr_c),
            num(ic, "person_id", fi$person_id, ct_c))
XM <- cbind(num(me, "hadm_id", fm$hadm_id, ehr_c),
            num(mc, "hadm_id", fm$hadm_id, ct_c))
colnames(XI) <- colnames(XM) <- c(ehr_c, ct_c)
DS <- list(
  I = list(name = "INSPECT", X = XI, t = fi$time_death,
           e = as.integer(fi$event_death), id = fi$person_id, off = 0,
           ctau = 0.001),
  M = list(name = "MIMIC", X = XM, t = fm$time_death,
           e = as.integer(fm$event_death), id = fm$hadm_id, off = 0.5,
           ctau = 1))
for (f in c("I", "M")) {
  D <- DS[[f]]
  chk(sprintf("%s: patients / deaths within 366 d", D$name),
      sprintf("%d / %d", nrow(D$X), sum(D$e)), "INFO")
  chk(sprintf("%s: patients with no feature row at all", D$name),
      sum(rowSums(!is.na(D$X)) == 0),
      if (any(rowSums(!is.na(D$X)) == 0)) "WARN" else "PASS")
  chk(sprintf("%s: cells missing, EHR / CTPA", D$name),
      sprintf("%.3f / %.3f", mean(is.na(D$X[, ehr_c])),
              mean(is.na(D$X[, ct_c]))), "INFO")
}
seeds <- sort(unique(fo$seed))
fold_of <- function(f, s) {
  x <- fo[fo$frame == f & fo$seed == s]
  x$fold[match(DS[[f]]$id, x$id)]
}
for (f in c("I", "M")) for (s in seeds) {
  if (anyNA(fold_of(f, s))) stop("fold missing: ", f, " seed ", s)
}
GROUPS <- list(EHR = c(ehr_c, paste0("missing_", ehr_c)),
               CTPA = c(ct_c, paste0("missing_", ct_c)))

# ---------------- tasks ----------------
tasks <- list()
for (f in c("I", "M")) for (s in seeds) {
  fold <- fold_of(f, s)
  for (k in 0:4) {
    tr <- which(fold != k)
    te <- which(fold == k)
    set.seed(s * 100 + k)
    tasks[[length(tasks) + 1]] <- list(
      type = "internal", key = f, seed = s, strategy = MAIN, rows = te,
      Xs = DS[[f]]$X[tr, , drop = FALSE], ts = pmax(DS[[f]]$t[tr], 1e-3),
      es = DS[[f]]$e[tr], ftr = sample(rep_len(1:5, length(tr))),
      Xt = DS[[f]]$X[te, , drop = FALSE], cuts = HZ + DS[[f]]$off,
      groups = GROUPS)
  }
}
for (tgt in c("I", "M")) {
  src <- setdiff(c("I", "M"), tgt)
  for (s in seeds) for (st in STRATS) {
    tasks[[length(tasks) + 1]] <- list(
      type = "transfer", key = tgt, seed = s, strategy = st,
      rows = seq_len(nrow(DS[[tgt]]$X)), Xs = DS[[src]]$X,
      ts = pmax(DS[[src]]$t, 1e-3), es = DS[[src]]$e,
      ftr = fold_of(src, s) + 1L, Xt = DS[[tgt]]$X,
      cuts = HZ + DS[[tgt]]$off, groups = GROUPS)
  }
}
cat(sprintf("  %d tasks (%d internal, %d transfer)\n", length(tasks),
            sum(sapply(tasks, function(t) t$type == "internal")),
            sum(sapply(tasks, function(t) t$type == "transfer"))))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
i1 <- which(sapply(tasks, function(t) t$type == "internal"))[1]
i2 <- which(sapply(tasks, function(t) t$type == "transfer" &&
                     t$key == "M"))[1]
t0 <- Sys.time()
p1 <- run_task(tasks[[i1]])
ta <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
t0 <- Sys.time()
p2 <- run_task(tasks[[i2]])
tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (!p1$ok || !p2$ok) {
  cat("STOP: probe task failed:", p1$msg, p2$msg, "\n")
  quit(status = 1)
}
ni <- sum(sapply(tasks, function(t) t$type == "internal"))
nt_ <- length(tasks) - ni
proj <- (ta * ni + tb * nt_) / workers * 1.2
cat(sprintf("  internal task %.1f s, transfer task %.1f s; projected",
            ta, tb), sprintf("about %.1f min\n", proj / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, run_task, future.packages = PKG,
                     future.seed = TRUE)
plan(sequential)
nfail <- sum(!sapply(res, function(r) r$ok))
chk("tasks that failed", nfail, if (nfail == 0) "PASS" else "FAIL")
if (nfail > 0) {
  for (r in Filter(function(r) !r$ok, res)) cat("   ", r$msg, "\n")
  fwrite(rbindlist(CHK), file.path(RS, "s6_checks.csv"))
  quit(status = 1)
}

# ---------------- penalty choices at the ends of the path ----------------
EDG <- rbindlist(lapply(res, function(r) {
  r$edges[, `:=`(type = r$type, target = r$key, strategy = r$strategy)]
}))
esum <- EDG[, .(fits = .N, median_path_length = median(nlam),
                share_least_penalised = mean(least),
                share_most_penalised = mean(most)),
            by = .(type, target, strategy, group, fit)]
fwrite(esum, file.path(RS, "s6_lambda.csv"))
for (g in sort(unique(EDG$group))) {
  sh <- EDG[group == g, mean(least)]
  chk(sprintf("%s models: choices at the least-penalised end", g),
      sprintf("%.1f%%", 100 * sh), if (sh <= 0.20) "PASS" else "WARN")
}

# ---------------- assemble predictions ----------------
INT <- list()
TRS <- list()
LPT <- list()
LPS <- list()
for (f in c("I", "M")) {
  n <- nrow(DS[[f]]$X)
  INT[[f]] <- list()
  for (sn in STACKS) {
    arr <- array(NA_real_, c(n, length(HZ), length(seeds)))
    for (r in Filter(function(r) r$type == "internal" && r$key == f, res)) {
      arr[r$rows, , match(r$seed, seeds)] <- r$tgt[[sn]]$risk
    }
    INT[[f]][[sn]] <- apply(arr, c(1, 2), mean)
  }
  TRS[[f]] <- list()
  for (st in STRATS) {
    TRS[[f]][[st]] <- list()
    rr <- Filter(function(r) r$type == "transfer" && r$key == f &&
                   r$strategy == st, res)
    for (sn in STACKS) {
      TRS[[f]][[st]][[sn]] <- Reduce(`+`, lapply(rr, function(r)
        r$tgt[[sn]]$risk)) / length(rr)
    }
    if (st == MAIN) {
      LPT[[f]] <- setNames(lapply(rr, function(r)
        r$tgt[["EHR + CTPA"]]$lp), sapply(rr, function(r) r$seed))
      LPS[[f]] <- Reduce(`+`, lapply(rr, function(r) r$src_lp)) / length(rr)
    }
  }
}
rng <- range(unlist(lapply(c(INT, TRS), function(z) unlist(z))))
chk("risk range across all predictions", sprintf("%.5f to %.5f", rng[1],
                                                 rng[2]),
    if (rng[1] >= 0 && rng[2] <= 1) "PASS" else "FAIL")

# ---------------- step 3: updating in the target ----------------
UPD <- list()
for (tgt in c("I", "M")) {
  D <- DS[[tgt]]
  n <- nrow(D$X)
  U0 <- array(NA_real_, c(n, length(HZ), length(seeds)))
  U1 <- U0
  for (si in seq_along(seeds)) {
    lp <- LPT[[tgt]][[as.character(seeds[si])]]
    fold <- fold_of(tgt, seeds[si])
    for (k in 0:4) {
      tr <- fold != k
      te <- fold == k
      b <- coef(coxph(Surv(D$t[tr], D$e[tr]) ~ lp[tr]))[[1]]
      for (hi in seq_along(HZ)) {
        cut <- HZ[hi] + D$off
        H0a <- breslow_H(D$t[tr], D$e[tr], exp(lp[tr]), cut)
        H0b <- breslow_H(D$t[tr], D$e[tr], exp(b * lp[tr]), cut)
        U0[te, hi, si] <- 1 - exp(-H0a * exp(lp[te]))
        U1[te, hi, si] <- 1 - exp(-H0b * exp(b * lp[te]))
      }
    }
  }
  UPD[[tgt]] <- list(base = apply(U0, c(1, 2), mean),
                     slope = apply(U1, c(1, 2), mean))
}

# ---------------- step 2: performance, and censoring bias ----------------
cat("\n=== EVALUATION ===\n")
PERF <- list()
CON <- list()
DEC <- list()
CEN <- list()
for (tgt in c("I", "M")) {
  D <- DS[[tgt]]
  src <- setdiff(c("I", "M"), tgt)
  dirl <- sprintf("%s -> %s", DS[[src]]$name, D$name)
  for (hi in seq_along(HZ)) {
    cut <- HZ[hi] + D$off
    ymax <- HZ[hi] + D$ctau
    for (sn in STACKS) {
      R <- list()
      R[["internal (refit in target)"]] <- INT[[tgt]][[sn]][, hi]
      for (st in STRATS) {
        R[[paste0("transfer: ", st)]] <- TRS[[tgt]][[st]][[sn]][, hi]
      }
      if (sn == "EHR + CTPA") {
        R[["transfer + baseline update"]] <- UPD[[tgt]]$base[, hi]
        R[["transfer + baseline and slope update"]] <-
          UPD[[tgt]]$slope[, hi]
      }
      ev <- eval_target(D$t, D$e, cut, ymax, R)
      lab <- data.table(target = D$name, direction = dirl,
                        horizon = HZ[hi], modality = sn)
      PERF[[length(PERF) + 1]] <- cbind(lab, ev$perf)
      if (nrow(ev$ca)) {
        CON[[length(CON) + 1]] <- cbind(lab, measure = "AUC", ev$ca)
      }
      if (nrow(ev$cb)) {
        CON[[length(CON) + 1]] <- cbind(lab, measure = "Brier", ev$cb)
      }
      if (sn == "EHR + CTPA") {
        for (m in c("internal (refit in target)",
                    paste0("transfer: ", MAIN),
                    "transfer + baseline update")) {
          DEC[[length(DEC) + 1]] <- cbind(lab, model = m,
                                          dec_tab(D$t, D$e, cut, R[[m]]))
        }
        if (tgt == "I") {
          died <- D$e == 1 & D$t <= cut
          known <- !(D$e == 0 & D$t < cut)
          for (m in c("internal (refit in target)",
                      paste0("transfer: ", MAIN))) {
            pr <- ev$perf[model == m]
            CEN[[length(CEN) + 1]] <- data.table(
              horizon = HZ[hi], model = m,
              censored_before_h = sum(!known),
              km_risk = pr$km_risk,
              risk_dropping_censored = mean(died[known]),
              risk_censored_as_alive = mean(died),
              auc_ipcw = pr$auc,
              auc_dropping_censored = auc1(died[known], R[[m]][known]),
              auc_censored_as_alive = auc1(died, R[[m]]))
          }
        }
      }
    }
  }
  cat(sprintf("  %s done\n", dirl))
}

# ---------------- step 1: relatedness ----------------
cat("\n=== RELATEDNESS ===\n")
REL <- list()
Xp <- rbind(XM, XI)
yp <- c(rep(0L, nrow(XM)), rep(1L, nrow(XI)))
imp_p <- impute(Xp, Xp, "median (source)")$s
imp_pi <- impute(Xp, Xp, "median + missing indicators")$s
mems <- list(
  "all features" = imp_p,
  "all features + missing indicators" = imp_pi,
  "EHR features only" = imp_p[, ehr_c, drop = FALSE],
  "CTPA features only" = imp_p[, ct_c, drop = FALSE])
for (mn in names(mems)) {
  REL[[length(REL) + 1]] <- data.table(
    part = "membership model (INSPECT vs MIMIC), out-of-fold AUROC",
    item = mn, value = member_auc(mems[[mn]], yp))
}
for (tgt in c("I", "M")) {
  src <- setdiff(c("I", "M"), tgt)
  ls_ <- LPS[[tgt]]
  lt <- Reduce(`+`, LPT[[tgt]]) / length(LPT[[tgt]])
  dirl <- sprintf("%s -> %s", DS[[src]]$name, DS[[tgt]]$name)
  REL[[length(REL) + 1]] <- data.table(
    part = sprintf("fused score, %s", dirl),
    item = c("mean in source", "SD in source", "mean in target",
             "SD in target", "mean difference in source SDs",
             "SD ratio target / source"),
    value = c(mean(ls_), sd(ls_), mean(lt), sd(lt),
              (mean(lt) - mean(ls_)) / sd(ls_), sd(lt) / sd(ls_)))
}
for (f in c("I", "M")) {
  D <- DS[[f]]
  sf <- summary(survfit(Surv(D$t, D$e) ~ 1), times = HZ + D$off,
                extend = TRUE)
  REL[[length(REL) + 1]] <- data.table(
    part = sprintf("%s Kaplan-Meier mortality", D$name),
    item = sprintf("%d days", HZ), value = 1 - sf$surv)
}
rel <- rbindlist(REL)
print(rel, row.names = FALSE)
FEAT <- rbindlist(lapply(colnames(XM), function(cn) {
  a <- XM[, cn]
  b <- XI[, cn]
  s <- sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
  data.table(feature = cn, modality = if (cn %in% ehr_c) "EHR" else "CTPA",
             mean_mimic = mean(a, na.rm = TRUE),
             mean_inspect = mean(b, na.rm = TRUE),
             smd_inspect_minus_mimic = (mean(b, na.rm = TRUE) -
                                          mean(a, na.rm = TRUE)) / s,
             missing_mimic = mean(is.na(a)),
             missing_inspect = mean(is.na(b)))
}))
setorder(FEAT, -smd_inspect_minus_mimic)

# ---------------- write ----------------
pf <- rbindlist(PERF, fill = TRUE)
cn <- rbindlist(CON, fill = TRUE)
fwrite(pf, file.path(RS, "s6_performance.csv"))
fwrite(cn, file.path(RS, "s6_contrasts.csv"))
fwrite(rel, file.path(RS, "s6_relatedness.csv"))
fwrite(FEAT, file.path(RS, "s6_features.csv"))
fwrite(rbindlist(CEN), file.path(RS, "s6_censoring.csv"))
fwrite(rbindlist(DEC), file.path(RS, "s6_deciles.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "s6_checks.csv"))

cat("\n=== FUSED MODEL (EHR + CTPA), 365 DAYS ===\n")
print(pf[modality == "EHR + CTPA" & horizon == 365,
         .(direction, model, auc = round(auc, 3), uno_c = round(uno_c, 3),
           oe = round(oe, 2), slope = round(slope, 2),
           ipa = round(ipa, 3))], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("s6_performance.csv", "s6_contrasts.csv", "s6_relatedness.csv",
            "s6_features.csv", "s6_censoring.csv", "s6_deciles.csv",
            "s6_lambda.csv", "s6_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "s6_log.txt"), " (terminal log)\n")
cat("\nOutputs: s6_performance.csv, s6_contrasts.csv, s6_relatedness.csv,",
    "s6_features.csv, s6_censoring.csv, s6_deciles.csv, s6_lambda.csv,",
    "s6_checks.csv\n")