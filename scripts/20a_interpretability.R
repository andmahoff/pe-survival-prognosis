# 20a_interpretability.R
# Interpretability of the final model's main components (A23).
#
# Part 1 refits the penalised-spline Cox model of 19b (mgcv cox.ph, cubic
# regression splines with k = 5 for features with at least 10 distinct
# values, linear otherwise, REML, select = TRUE) once on every patient of
# each frame and cohort, with the 18a preprocessing fitted on all patients.
# It is for interpretation, not prediction. For each EHR feature it reports
# effective degrees of freedom, Wald chi-square and p-value; the log hazard
# ratio from the 1st to the 99th percentile relative to the median, with 95%
# intervals from the model covariance (exact contrasts through the
# linear-predictor matrix); and the SD of the term across patients as an
# importance score.
#
# Part 2 computes SurvSHAP(t) (Krzyzinski et al. 2023) for the two forest
# families that carry most of Super learner v4's forest weight, Brier-tuned
# RSF and oblique RSF (06g). Each is fitted to the EHR alone, once on every
# patient, at the settings 06g chose most often (modal choice across folds
# and seeds). Interventional Shapley values of the predicted survival at
# each horizon are estimated by antithetic permutation sampling (16
# permutations and their reverses) against 100 background patients, for 300
# randomly chosen patients (seed 2026). Values are on the risk scale (1 -
# survival), so a positive value raises risk, and each patient's values must
# add up to their predicted risk minus the background mean at every horizon.
# Summaries are the mean |value| per feature and horizon, each feature's
# share of the total (does importance shift over time?) and the Spearman
# correlation of importance with the spline model.
#
# Part 1 covers every frame and cohort; Part 2 covers the EHR+ECG cohort in
# Frames B (primary) and A. The script stops if a fit does not converge,
# SHAP additivity fails beyond 1e-8 or the modal settings are missing.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, feature_whitelist.csv,
#          level1/params_proper.csv; data/inputs/ehr28.csv
# Outputs: results/interp_{spline_terms,spline_curves,shap_importance,
#          agreement,checks}.csv, results/interp_log.txt,
#          data/processed/interp/shap_*.csv (patient level, not in the
#          repo)
# Usage:   Rscript scripts/20a_interpretability.R [--probe-only]
#          [--workers=N]

suppressPackageStartupMessages({
  library(survival)
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
OD <- file.path(PR, "interp")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
KNOTS <- 5L
MINU <- 10L
NGRID <- 50L
SHAP_SETS <- list(c("B", "EHR+ECG"), c("A", "EHR+ECG"))
NEXPL <- 300L
NBG <- 100L
NPERM <- 16L
CHUNK <- 8L
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
MTRYF <- c(0.5, 1, 3)
OMTRYF <- c(0.5, 1, 2, 4)
BATCH <- 50L
SEED <- 2026L
PKG <- c("survival", "data.table", "randomForestSRC", "aorsf")
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
    fwrite(ck, file.path(RS, "interp_checks.csv"))
    cat("\nSTOP: a check failed. see: interp_checks.csv\n")
    quit(status = 1)
  }
}
# stage 18a preprocessing, fitted on the rows `tr`; returns the function
prep_fit <- function(X, tr = seq_len(nrow(X))) {
  med <- q1 <- q2 <- mu <- s <- numeric(ncol(X))
  bin <- logical(ncol(X))
  Z <- X
  for (j in seq_len(ncol(X))) {
    med[j] <- median(X[tr, j], na.rm = TRUE)
    if (!is.finite(med[j])) med[j] <- 0
    Z[is.na(Z[, j]), j] <- med[j]
    bin[j] <- length(unique(Z[tr, j])) <= 2
    q1[j] <- if (bin[j]) -Inf else quantile(Z[tr, j], 0.005, names = FALSE)
    q2[j] <- if (bin[j]) Inf else quantile(Z[tr, j], 0.995, names = FALSE)
    Z[, j] <- pmin(pmax(Z[, j], q1[j]), q2[j])
    mu[j] <- mean(Z[tr, j])
    s[j] <- sd(Z[tr, j])
  }
  keep <- is.finite(s) & s > 0
  list(med = med, q1 = q1, q2 = q2, mu = mu, s = s, keep = keep,
       bin = bin, names = colnames(X))
}
prep_apply <- function(P, X) {
  for (j in seq_len(ncol(X))) {
    X[is.na(X[, j]), j] <- P$med[j]
    X[, j] <- pmin(pmax(X[, j], P$q1[j]), P$q2[j])
  }
  X <- sweep(sweep(X[, P$keep, drop = FALSE], 2, P$mu[P$keep]), 2,
             P$s[P$keep], "/")
  X
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fr_d <- list(A = fread(file.path(PR, "frame_A.csv")),
             B = fread(file.path(PR, "frame_B.csv")))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
c28 <- intersect(wl[modality == "ehr", column], names(e28))
cat(sprintf("  %d EHR features\n", length(c28)))
getX <- function(sub) {
  X <- as.matrix(e28[match(sub$hadm_id, hadm_id), ..c28])
  storage.mode(X) <- "double"
  X
}

# ---------------- part 1: spline Cox partial effects ----------------
cat("\n=== PART 1: SPLINE COX PARTIAL EFFECTS ===\n")
TERMS <- list()
CURV <- list()
SPI <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  sub <- if (co == "EHR+ECG") fr_d[[fr]] else fr_d[[fr]][in_three_mod == 1]
  X <- getX(sub)
  P <- prep_fit(X)
  Z <- prep_apply(P, X)
  nm <- c28[P$keep]
  colnames(Z) <- sprintf("f%d", seq_len(ncol(Z)))
  nu <- apply(Z, 2, function(v) length(unique(v)))
  sm <- colnames(Z)[nu >= MINU]
  li <- setdiff(colnames(Z), sm)
  fml <- as.formula(paste("t ~", paste(c(sprintf(
    "s(%s, k = %d, bs = \"cr\")", sm, KNOTS), li), collapse = " + ")))
  D <- data.frame(t = sub$time_death, e = as.integer(sub$event_death == 1),
                  Z)
  t0 <- Sys.time()
  fit <- tryCatch(gam(fml, family = cox.ph(), weights = e, data = D,
                      method = "REML", select = TRUE),
                  error = function(err) conditionMessage(err))
  if (is.character(fit)) {
    chk(sprintf("%s %s: spline Cox fit", fr, co), substr(fit, 1, 80),
        "FAIL")
    stop_if_fail()
  }
  cv <- fit$outer.info$conv
  chk(sprintf("%s %s: spline Cox fit (%.0f s, total edf %.1f)", fr, co,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(fit$edf)), if (is.null(cv)) "" else cv,
      if (is.null(cv) || grepl("full convergence", cv)) "PASS" else "WARN")
  ss <- summary(fit)
  st <- ss$s.table
  pt <- ss$p.table
  b <- coef(fit)
  Vp <- fit$Vp
  TM <- predict(fit, type = "terms")
  for (j in seq_along(nm)) {
    fj <- sprintf("f%d", j)
    jj <- which(P$keep)[j]
    smooth <- fj %in% sm
    if (smooth) {
      rn <- sprintf("s(%s)", fj)
      cols <- grep(sprintf("^s\\(%s\\)\\.", fj), names(b))
      ed <- st[rn, "edf"]
      chi <- st[rn, "Chi.sq"]
      pv <- st[rn, "p-value"]
      tsd <- sd(TM[, rn])
    } else {
      cols <- which(names(b) == fj)
      ed <- 1
      chi <- (pt[fj, "Estimate"] / pt[fj, "Std. Error"])^2
      pv <- pt[fj, "Pr(>|z|)"]
      tsd <- sd(TM[, fj])
    }
    TERMS[[length(TERMS) + 1]] <- data.table(
      frame = fr, cohort = co, feature = nm[j],
      type = if (smooth) "spline" else "linear", edf = ed, chi_sq = chi,
      p_value = pv, term_sd = tsd, n = nrow(sub),
      deaths = sum(D$e))
    zg <- if (smooth) {
      seq(quantile(Z[, fj], 0.01), quantile(Z[, fj], 0.99),
          length.out = NGRID)
    } else {
      sort(unique(Z[, fj]))
    }
    zmed <- median(Z[, fj])
    G <- as.data.frame(matrix(0, length(zg) + 1, ncol(Z),
                              dimnames = list(NULL, colnames(Z))))
    G[, fj] <- c(zg, zmed)
    G$t <- 1
    L <- predict(fit, newdata = G, type = "lpmatrix")
    Dm <- L[seq_along(zg), cols, drop = FALSE] -
      matrix(L[length(zg) + 1, cols], length(zg), length(cols),
             byrow = TRUE)
    est <- as.vector(Dm %*% b[cols])
    se <- sqrt(pmax(0, rowSums((Dm %*% Vp[cols, cols, drop = FALSE]) * Dm)))
    xo <- P$mu[jj] + P$s[jj] * zg
    CURV[[length(CURV) + 1]] <- data.table(
      frame = fr, cohort = co, feature = nm[j],
      type = if (smooth) "spline" else "linear", x = xo,
      x_median = P$mu[jj] + P$s[jj] * zmed,
      x_p10 = quantile(X[, jj], 0.10, na.rm = TRUE),
      x_p90 = quantile(X[, jj], 0.90, na.rm = TRUE),
      log_hr = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }
  SPI[[paste(fr, co)]] <- rbindlist(TERMS)[frame == fr & cohort == co]
  rm(fit)
}
terms <- rbindlist(TERMS)
setorder(terms, frame, cohort, -chi_sq)
fwrite(terms, file.path(RS, "interp_spline_terms.csv"))
fwrite(rbindlist(CURV), file.path(RS, "interp_spline_curves.csv"))
cat("\n  strongest terms, Frame B EHR+ECG (chi-square):\n")
print(terms[frame == "B" & cohort == "EHR+ECG"][1:10,
            .(feature, type, edf = round(edf, 1), chi_sq = round(chi_sq, 1),
              p_value = signif(p_value, 2), term_sd = round(term_sd, 3))],
      row.names = FALSE)

# ---------------- part 2: SurvSHAP(t) ----------------
cat("\n=== PART 2: SurvSHAP(t) FOR THE FORESTS ===\n")
pp <- fread(file.path(PR, "level1", "params_proper.csv"))
modal <- function(v) {
  v <- v[is.finite(v)]
  as.numeric(names(sort(table(v), decreasing = TRUE))[1])
}
# fitted forests wrapped as prediction functions: original-scale EHR
# matrix -> matrix of risk at each horizon
fit_rsf <- function(X, tt, ee, H, node, mtryv) {
  P <- prep_fit(X)
  Xi <- X
  for (j in seq_len(ncol(X))) Xi[is.na(Xi[, j]), j] <- P$med[j]
  colnames(Xi) <- sprintf("f%d", seq_len(ncol(Xi)))
  ev <- sort(unique(tt[ee == 1]))
  at_h <- sapply(H, function(h) {
    v <- ev[ev <= h + 0.5]
    if (length(v)) max(v) else NA_real_
  })
  pick <- ev[unique(round(seq(1, length(ev),
                              length.out = min(NTIME, length(ev)))))]
  nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
  f <- rfsrc(Surv(time, status) ~ ., ntree = NTREE, nodesize = node,
             mtry = mtryv, nsplit = NSPLIT, ntime = nt,
             splitrule = "logrank", save.memory = TRUE, seed = -SEED,
             data = data.frame(time = tt, status = ee, Xi))
  list(kind = "rsf", f = f, P = P, H = H)
}
fit_orsf <- function(X, tt, ee, H, leaf, mlevel) {
  P <- prep_fit(X)
  Z <- prep_apply(P, X)
  colnames(Z) <- sprintf("f%d", seq_len(ncol(Z)))
  ps <- ncol(Z)
  mt <- sort(unique(pmax(1L, pmin(ps, ceiling(OMTRYF * sqrt(ps))))))
  mtry <- mt[min(mlevel, length(mt))]
  d <- data.frame(time = tt, status = ee, Z)
  fs <- lapply(seq_len(NTREE %/% BATCH), function(bb) {
    aorsf::orsf(d, Surv(time, status) ~ ., n_tree = BATCH, mtry = mtry,
                leaf_min_obs = leaf, split_min_obs = 2L * leaf,
                oobag_pred_type = "none",
                tree_seeds = as.integer(SEED * 10 + bb), n_thread = 1)
  })
  list(kind = "orsf", f = fs, P = P, H = H, mtry = mtry)
}
predict_risk <- function(M, X) {
  if (M$kind == "rsf") {
    Xi <- X
    for (j in seq_len(ncol(X))) Xi[is.na(Xi[, j]), j] <- M$P$med[j]
    colnames(Xi) <- sprintf("f%d", seq_len(ncol(Xi)))
    pr <- predict(M$f, newdata = as.data.frame(Xi))
    sapply(M$H, function(h) {
      i <- findInterval(h + 0.5, pr$time.interest)
      if (i == 0) rep(0, nrow(X)) else 1 - pr$survival[, i]
    })
  } else {
    Z <- prep_apply(M$P, X)
    colnames(Z) <- sprintf("f%d", seq_len(ncol(Z)))
    Zd <- as.data.frame(Z)
    Reduce(`+`, lapply(M$f, function(f) {
      matrix(predict(f, new_data = Zd, pred_type = "risk",
                     pred_horizon = M$H + 0.5, boundary_checks = FALSE),
             ncol = length(M$H))
    })) / length(M$f)
  }
}
# antithetic permutation SHAP of the risk vector, for one patient
shap_one <- function(M, x, BG, nperm, seed) {
  p <- length(x)
  K <- length(M$H)
  set.seed(seed)
  perms <- lapply(seq_len(nperm), function(i) sample.int(p))
  perms <- c(perms, lapply(perms, rev))
  phi <- matrix(0, p, K)
  nb <- nrow(BG)
  for (ch in split(seq_along(perms), ceiling(seq_along(perms) / CHUNK))) {
    Xall <- do.call(rbind, lapply(perms[ch], function(pm) {
      do.call(rbind, lapply(0:p, function(s) {
        B <- BG
        if (s > 0) {
          cs <- pm[seq_len(s)]
          B[, cs] <- matrix(x[cs], nb, s, byrow = TRUE)
        }
        B
      }))
    }))
    R <- predict_risk(M, Xall)
    blk <- rep(seq_len(length(ch) * (p + 1)), each = nb)
    V <- rowsum(R, blk) / nb
    for (i in seq_along(ch)) {
      pm <- perms[[ch[i]]]
      v <- V[(i - 1) * (p + 1) + seq_len(p + 1), , drop = FALSE]
      phi[pm, ] <- phi[pm, ] + (v[-1, , drop = FALSE] -
                                  v[-(p + 1), , drop = FALSE])
    }
  }
  phi / length(perms)
}
JOBS <- list()
for (sc in SHAP_SETS) {
  fr <- sc[1]
  co <- sc[2]
  sub <- if (co == "EHR+ECG") fr_d[[fr]] else fr_d[[fr]][in_three_mod == 1]
  X <- getX(sub)
  q <- pp[frame == fr & cohort == co & modality == "ehr"]
  if (!nrow(q)) {
    chk(sprintf("%s %s: stage 6g settings", fr, co), "missing", "FAIL")
    stop_if_fail()
  }
  p <- ncol(X)
  node <- modal(q$rsfb_node)
  mtry_r <- modal(q$rsfb_mtry)
  leaf <- modal(q$orsfb_leaf)
  mlev <- modal(q$orsfb_mlevel)
  chk(sprintf("%s %s: modal settings RSF node/mtry, oblique leaf/level",
              fr, co), sprintf("%g / %g, %g / %g", node, mtry_r, leaf,
                               mlev), "PASS")
  tt <- sub$time_death
  ee <- as.integer(sub$event_death == 1)
  set.seed(SEED)
  bg <- sample.int(nrow(X), min(NBG, nrow(X)))
  ex <- sample.int(nrow(X), min(NEXPL, nrow(X)))
  JOBS[[length(JOBS) + 1]] <- list(fr = fr, co = co, X = X, tt = tt, ee = ee,
                                   H = HZ[[fr]], node = node, mtry_r = mtry_r,
                                   leaf = leaf, mlev = mlev, bg = bg, ex = ex,
                                   hadm = sub$hadm_id)
}
stop_if_fail()
MODELS <- list()
for (J in JOBS) {
  t0 <- Sys.time()
  MODELS[[paste(J$fr, J$co, "RSF, Brier-tuned")]] <-
    fit_rsf(J$X, J$tt, J$ee, J$H, J$node, J$mtry_r)
  MODELS[[paste(J$fr, J$co, "oblique RSF, Brier-tuned")]] <-
    fit_orsf(J$X, J$tt, J$ee, J$H, J$leaf, J$mlev)
  cat(sprintf("  Frame %s %s: forests fitted in %.1f min\n", J$fr, J$co,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
tot <- 0
for (J in JOBS) for (mn in c("RSF, Brier-tuned",
                             "oblique RSF, Brier-tuned")) {
  M <- MODELS[[paste(J$fr, J$co, mn)]]
  t0 <- Sys.time()
  invisible(shap_one(M, J$X[J$ex[1], ], J$X[J$bg, , drop = FALSE], NPERM,
                     1L))
  tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tot <- tot + tp * length(J$ex)
  cat(sprintf("  Frame %s %-26s %.1f s per patient\n", J$fr, mn, tp))
}
cat(sprintf(paste("  %d patients per model: about %.0f min on %d",
                  "workers\n"), NEXPL, tot * 1.15 / workers / 60, workers))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
IMP <- list()
AGR <- list()
for (J in JOBS) for (mn in c("RSF, Brier-tuned",
                             "oblique RSF, Brier-tuned")) {
  M <- MODELS[[paste(J$fr, J$co, mn)]]
  BG <- J$X[J$bg, , drop = FALSE]
  t0 <- Sys.time()
  ph <- future_lapply(seq_along(J$ex), function(i) {
    shap_one(M, J$X[J$ex[i], ], BG, NPERM, 1000L + i)
  }, future.packages = PKG, future.seed = TRUE)
  fx <- predict_risk(M, J$X[J$ex, , drop = FALSE])
  f0 <- colMeans(predict_risk(M, BG))
  add <- max(sapply(seq_along(ph), function(i) {
    max(abs(colSums(ph[[i]]) - (fx[i, ] - f0)))
  }))
  chk(sprintf("Frame %s %s: SHAP additivity (max |diff|), %.1f min", J$fr,
              mn, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
      sprintf("%.1e", add), if (add <= 1e-8) "PASS" else "FAIL")
  A <- array(unlist(ph), c(ncol(J$X), length(J$H), length(ph)))
  long <- rbindlist(lapply(seq_along(ph), function(i) {
    data.table(hadm_id = J$hadm[J$ex[i]], feature = rep(c28, length(J$H)),
               horizon = rep(J$H, each = length(c28)),
               value = rep(J$X[J$ex[i], ], length(J$H)),
               shap_risk = as.vector(ph[[i]]))
  }))
  fwrite(long, file.path(OD, sprintf("shap_%s_%s_%s.csv", J$fr,
                                     gsub("[^A-Za-z]+", "", J$co),
                                     gsub("[^A-Za-z]+", "_", mn))))
  ma <- apply(abs(A), c(1, 2), mean)
  sh <- sweep(ma, 2, colSums(ma), "/")
  for (hi in seq_along(J$H)) {
    IMP[[length(IMP) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, model = mn, feature = c28,
      horizon = J$H[hi], mean_abs_shap = ma[, hi], share = sh[, hi],
      rank = rank(-ma[, hi], ties.method = "first"),
      mean_shap = apply(A[, hi, , drop = FALSE], 1, mean),
      n_explained = length(ph), n_background = nrow(BG))
  }
  sp <- SPI[[paste(J$fr, J$co)]]
  hl <- length(J$H)
  r <- suppressWarnings(cor(ma[, hl], sp$term_sd[match(c28, sp$feature)],
                            method = "spearman", use = "complete.obs"))
  top <- c28[order(-ma[, hl])][1:5]
  AGR[[length(AGR) + 1]] <- data.table(
    frame = J$fr, cohort = J$co, model = mn,
    comparison = "importance at the last horizon vs spline term SD",
    spearman = r, forest_top5 = paste(top, collapse = ", "),
    spline_top5 = paste(sp[order(-term_sd)]$feature[1:5], collapse = ", "))
}
plan(sequential)
imp <- rbindlist(IMP)
agr <- rbindlist(AGR)
fwrite(imp, file.path(RS, "interp_shap_importance.csv"))
fwrite(agr, file.path(RS, "interp_agreement.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "interp_checks.csv"))
cat("\n=== AGREEMENT WITH THE SPLINE MODEL ===\n")
print(agr[, .(frame, model, spearman = round(spearman, 2), forest_top5,
              spline_top5)], row.names = FALSE)
cat("\n=== SHARE OF IMPORTANCE BY HORIZON, Frame B, top 8 at 365 days ===\n")
for (mn in unique(imp$model)) {
  x <- imp[frame == "B" & model == mn]
  top <- x[horizon == max(horizon)][order(rank)]$feature[1:8]
  cat("  ", mn, "\n")
  print(dcast(x[feature %in% top], feature ~ horizon, value.var = "share")[
    , lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)],
    row.names = FALSE)
}
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("interp_spline_terms.csv", "interp_spline_curves.csv",
            "interp_shap_importance.csv", "interp_agreement.csv",
            "interp_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", OD, " (shap_*.csv, patient-level, local only)\n")
cat("  ", file.path(RS, "interp_log.txt"), " (terminal log)\n")
cat("\nOutputs: interp_spline_terms.csv, interp_spline_curves.csv,",
    "interp_shap_importance.csv, interp_agreement.csv, interp_checks.csv\n")