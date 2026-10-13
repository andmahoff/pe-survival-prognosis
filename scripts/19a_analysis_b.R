# 19a_analysis_b.R
# Analysis B (A18): does the ECG's added value, and its front-loaded weight,
# survive a much stronger EHR model?
#
# Frame B only, because the extra labs run to discharge and could leak
# in-hospital deaths in Frame A. Death scores are ridge-Cox as in 18a and
# 06d (median imputation, winsorising, standardisation, glmnet alpha 0, D11
# penalty path, inner 3-fold subject-grouped deviance tuning, the folds of
# 04, six seeds) for three EHR versions: EHR 28 (the 28 registered features,
# reusing the 06d scores), EHR 79 (the 28 plus the 51 extra index-admission
# labs of ehr_extra.csv, with vancomycin level and blood-gas temperature
# excluded as recorded at the pull) and EHR 79 + MI (EHR 79 plus a
# missing-value indicator for every extra lab with any missing value,
# EHR+ECG cohort only). ECG and CTPA scores are 06d's retuned ridge-Cox.
#
# Scores are fused per seed and outer fold on rank-based normal scores of
# the training rows (as S4), with a Cox stack (PH, Breslow risks) and
# per-horizon logistic stacks (TV); risks are averaged over seeds. The
# EHR+ECG cohort has EHR alone and EHR + ECG for each EHR version; the
# three-modality cohort has EHR alone, + ECG, + CTPA and + ECG + CTPA (EHR
# 28 and 79).
#
# The primary's estimand (07_primary.R) is the ECG hazard ratio in days 0-30
# over days 31-365 (Cox model split at 30.5 on standardised seed-averaged
# scores), computed with the EHR 28 score and with the EHR 79 score, and the
# ratio of the two. Measures are AUROC, IPA and log-likelihood skill at 10,
# 30, 90, 180 and 365 days, and the calibration slope (point estimate).
# Contrasts: EHR 79 minus EHR 28 (alone and fused), EHR 79 + MI minus EHR
# 79, and the fusion gain of each modality set with each EHR version and the
# difference of those gains (does a stronger EHR leave less to add?).
# Intervals use a 2,000-sample patient-clustered bootstrap (seed 42).
#
# Checks: inputs found and aligned (stop); extra labs found (stop if fewer
# than 40); every score present (stop); the D2 per-fold rule (stop); penalty
# choices at the least-penalised end (warning above 20%).
#
# Inputs:  data/processed/frame_B.csv, folds.csv, feature_whitelist.csv,
#          level1/oof_extra_B_*.csv; data/inputs/ehr28.csv, ehr_extra.csv
# Outputs: results/ab_{metrics,contrasts,ecg_ratio,checks}.csv,
#          results/ab_log.txt, data/processed/analysis_b/ab_scores_*.csv
#          (patient level, not in the repo)
# Usage:   Rscript scripts/19a_analysis_b.R [--probe-only] [--workers=N]

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
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
OD <- file.path(PR, "analysis_b")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
CUT <- 30.5
NLAM <- 80L
LMR <- 1e-7
NBOOT <- 2000L
SEED <- 42L
EXCL <- c("vancomycin", "_50825_")
PKG <- c("survival", "glmnet", "data.table")
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
    fwrite(ck, file.path(RS, "ab_checks.csv"))
    cat("\nSTOP: a check failed. see: ab_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- level 1: ridge-Cox on an EHR version ----------------
ridge_task <- function(tk) {
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
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  lp <- rep(NA_real_, n)
  edges <- c()
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    X <- prep(tk$X, tr)
    cv <- cv.glmnet(X[tr, , drop = FALSE], Surv(tk$t[tr], tk$e[tr]),
                    family = "cox", alpha = 0, foldid = inner,
                    nlambda = NLAM, lambda.min.ratio = LMR,
                    type.measure = "deviance",
                    penalty.factor = rep(1, ncol(X)), standardize = FALSE)
    i <- which(cv$lambda == cv$lambda.min)[1]
    edges <- c(edges, i == length(cv$lambda))
    lp[te] <- as.vector(predict(cv, X[te, , drop = FALSE],
                                s = "lambda.min"))
  }
  list(co = tk$co, ver = tk$ver, seed = tk$seed, lp = lp, least = edges)
}

# ---------------- level 2: stacks (as S4 and A14) ----------------
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

# ---------------- evaluation ----------------
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(0)
  for (hi in seq_along(HZ)) {
    y <- J$Y[idx, hi]
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    for (k in names(J$P)) {
      p <- pmin(pmax(J$P[[k]][idx, hi], 1e-6), 1 - 1e-6)
      key <- paste(k, HZ[hi], sep = "|")
      v[paste(key, "auc", sep = "|")] <- auc(y, p)
      v[paste(key, "ipa", sep = "|")] <- 1 - mean((y - p)^2) /
        (pr * (1 - pr))
      v[paste(key, "ll", sep = "|")] <-
        1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
    }
  }
  v
}
ratio_vec <- function(R, w) {
  v <- numeric(0)
  k <- w[R$ep$id] > 0
  for (nm in names(R$Z)) {
    ze <- R$Z[[nm]][R$ep$id]
    zc <- R$ecg[R$ep$id]
    X <- cbind(ze * R$early, ze * (1 - R$early), zc * R$early,
               zc * (1 - R$early))
    b <- tryCatch(coef(coxph(Surv(R$ep$tstart[k], R$ep$t[k], R$ep$e[k]) ~
                               X[k, , drop = FALSE], weights = w[R$ep$id][k],
                             ties = "efron")),
                  error = function(e) rep(NA_real_, 4))
    v[paste(nm, "ECG log HR early", sep = "|")] <- b[3]
    v[paste(nm, "ECG log HR late", sep = "|")] <- b[4]
    v[paste(nm, "ECG log ratio", sep = "|")] <- b[3] - b[4]
    v[paste(nm, "EHR log ratio", sep = "|")] <- b[1] - b[2]
  }
  v
}
boot_rep <- function(i, J, R) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                     length(J$rows_by), replace = TRUE)],
                use.names = FALSE)
  out <- stat_vec(J, idx)
  if (!is.null(R)) out <- c(out, ratio_vec(R, tabulate(idx, J$n)))
  out
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PR, c("frame_B.csv", "folds.csv",
                          "feature_whitelist.csv",
                          "level1/oof_extra_B_ehrecg.csv",
                          "level1/oof_extra_B_3mod.csv")),
          file.path(IN, c("ehr28.csv", "ehr_extra.csv")))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
ex <- fread(file.path(IN, "ehr_extra.csv"))
chk("ehr_extra.csv keyed by hadm_id", "hadm_id" %in% names(ex),
    if ("hadm_id" %in% names(ex)) "PASS" else "FAIL")
stop_if_fail()
ex <- unique(ex, by = "hadm_id")
labs <- grep("^lab_", names(ex), value = TRUE)
drop <- labs[grepl(paste(EXCL, collapse = "|"), labs)]
labs <- setdiff(labs, drop)
chk("extra labs used / excluded", sprintf("%d / %s", length(labs),
                                          paste(drop, collapse = ", ")),
    if (length(labs) >= 40) "PASS" else "FAIL")
stop_if_fail()
c28 <- intersect(wl[modality == "ehr", column], names(e28))
tasks <- list()
META <- list()
for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  n <- nrow(sub)
  X28 <- as.matrix(e28[match(sub$hadm_id, hadm_id), ..c28])
  XL <- as.matrix(ex[match(sub$hadm_id, hadm_id), ..labs])
  storage.mode(X28) <- "double"
  storage.mode(XL) <- "double"
  cov <- colMeans(!is.na(XL))
  MI <- 1 * is.na(XL[, cov < 1, drop = FALSE])
  colnames(MI) <- paste0("mi_", colnames(MI))
  chk(sprintf("%s: patients / median extra-lab coverage / indicators",
              co), sprintf("%d / %.2f / %d", n, median(cov), ncol(MI)),
      "INFO")
  VERS <- list("EHR 79" = cbind(X28, XL))
  if (co == "EHR+ECG") VERS[["EHR 79 + MI"]] <- cbind(X28, XL, MI)
  fo <- folds[frame == "B" & cohort == co]
  seeds <- sort(unique(fo$seed))
  XT <- fread(file.path(L1, sprintf("oof_extra_B_%s.csv", code)))
  ok <- identical(seeds, sort(unique(XT$seed))) &&
    nrow(XT) == n * length(seeds)
  chk(sprintf("%s: stage 6d scores aligned with the folds", co),
      sprintf("%d seeds, %d rows", length(seeds), nrow(XT)),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  fmap <- list()
  for (s in seeds) {
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(f)) {
      chk(sprintf("%s seed %s: folds complete", co, s), sum(is.na(f)),
          "FAIL")
      stop_if_fail()
    }
    fmap[[as.character(s)]] <- as.integer(f)
    for (vn in names(VERS)) {
      tasks[[length(tasks) + 1]] <- list(
        co = co, ver = vn, seed = s, X = VERS[[vn]], t = sub$time_death,
        e = as.integer(sub$event_death == 1), subj = sub$subject_id,
        fold = as.integer(f))
    }
  }
  META[[co]] <- list(co = co, code = code, sub = sub, seeds = seeds,
                     fmap = fmap, XT = XT, vers = names(VERS))
}
cat(sprintf("  %d level-1 tasks (cohort x EHR version x seed)\n",
            length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) nrow(t$X) * ncol(t$X)))
t0 <- Sys.time()
invisible(ridge_task(tasks[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  largest task %.0f s; %d tasks on %d workers: about",
                  "%.1f min, plus stacks and bootstrap\n"), tp,
            length(tasks), workers,
            tp * ceiling(length(tasks) / workers) * 1.2 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, ridge_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1L)

# ---------------- stacks, per cohort ----------------
cat("\n=== STACKS ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
JOBS <- list()
RAT <- NULL
for (co in names(META)) {
  M <- META[[co]]
  sub <- M$sub
  n <- nrow(sub)
  tt <- sub$time_death
  ee <- as.integer(sub$event_death == 1)
  cuts <- HZ + 0.5
  least <- unlist(lapply(Filter(function(r) r$co == co, res),
                         function(r) r$least))
  chk(sprintf("%s: EHR 79 penalty at the least-penalised end", co),
      sprintf("%.1f%%", 100 * mean(least)),
      if (mean(least) <= 0.20) "PASS" else "WARN")
  SETS <- list(c("EHR 28"), c("EHR 79"), c("EHR 28", "ECG"),
               c("EHR 79", "ECG"))
  if ("EHR 79 + MI" %in% M$vers) {
    SETS <- c(SETS, list(c("EHR 79 + MI"), c("EHR 79 + MI", "ECG")))
  }
  if (co == "3-mod") {
    SETS <- c(SETS, list(c("EHR 28", "CTPA"), c("EHR 79", "CTPA"),
                         c("EHR 28", "ECG", "CTPA"),
                         c("EHR 79", "ECG", "CTPA")))
  }
  P <- list()
  SC <- list()
  for (s in M$seeds) {
    x <- M$XT[M$XT$seed == s][match(sub$hadm_id, hadm_id)]
    S <- list("EHR 28" = x[["ehr.cox2.death.lp"]],
              "ECG" = x[["ecg.cox2.death.lp"]])
    if (co == "3-mod") S[["CTPA"]] <- x[["ctpa.cox2.death.lp"]]
    for (r in Filter(function(r) r$co == co && r$seed == s, res)) {
      S[[r$ver]] <- r$lp
    }
    nb <- sum(sapply(S, function(v) sum(!is.finite(v))))
    if (nb > 0) {
      chk(sprintf("%s seed %s: missing scores", co, s), nb, "FAIL")
      stop_if_fail()
    }
    SC[[as.character(s)]] <- S
    f <- M$fmap[[as.character(s)]]
    for (st in SETS) {
      nm <- paste(st, collapse = " + ")
      for (form in c("PH", "TV")) {
        k <- paste(nm, form, sep = " | ")
        if (is.null(P[[k]])) P[[k]] <- matrix(0, n, length(HZ))
      }
      for (k in sort(unique(f))) {
        tr <- which(f != k)
        te <- which(f == k)
        sf <- stack_fit(S, st, tr, te, tt, ee, cuts)
        P[[paste(nm, "PH", sep = " | ")]][te, ] <-
          P[[paste(nm, "PH", sep = " | ")]][te, ] + sf$PH / length(M$seeds)
        P[[paste(nm, "TV", sep = " | ")]][te, ] <-
          P[[paste(nm, "TV", sep = " | ")]][te, ] + sf$TV / length(M$seeds)
      }
    }
  }
  nbad <- 0L
  for (vn in M$vers) {
    v <- c()
    for (s in M$seeds) {
      f <- M$fmap[[as.character(s)]]
      for (k in sort(unique(f))) {
        ii <- which(f == k)
        v <- c(v, cidx(tt[ii], ee[ii], SC[[as.character(s)]][[vn]][ii]))
      }
    }
    if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
  }
  chk(sprintf("%s: new EHR scores failing the per-fold rule (D2)", co),
      nbad, if (nbad == 0) "PASS" else "FAIL")
  stop_if_fail()
  avg <- lapply(names(SC[[1]]), function(v) {
    Reduce(`+`, lapply(SC, function(z) z[[v]])) / length(SC)
  })
  names(avg) <- names(SC[[1]])
  fwrite(cbind(data.table(hadm_id = sub$hadm_id), as.data.table(avg)),
         file.path(OD, sprintf("ab_scores_%s.csv", M$code)))
  Y <- sapply(HZ, function(h) as.integer(ee == 1 & tt <= h + 0.5))
  JOBS[[co]] <- list(co = co, n = n, P = P, Y = matrix(Y, ncol = length(HZ)),
                     rows_by = split(seq_len(n),
                                     match(sub$subject_id,
                                           unique(sub$subject_id))),
                     sets = sapply(SETS, paste, collapse = " + "))
  if (co == "EHR+ECG") {
    z <- function(v) (v - mean(v)) / sd(v)
    ep <- survSplit(Surv(t, e) ~ ., cut = CUT, episode = "per",
                    data = data.frame(t = tt, e = ee, id = seq_len(n)))
    vz <- c("EHR 28", M$vers)
    RAT <- list(Z = setNames(lapply(vz, function(v) z(avg[[v]])), vz),
                ecg = z(avg[["ECG"]]), early = as.numeric(ep$per == 1),
                ep = list(tstart = ep$tstart, t = ep$t, e = ep$e,
                          id = ep$id))
  }
}

# ---------------- bootstrap ----------------
cat(sprintf("\n=== BOOTSTRAP (%d resamples) ===\n", NBOOT))
gainpairs <- function(J) {
  out <- list()
  vs <- c("EHR 28", "EHR 79", "EHR 79 + MI")
  vs <- vs[vs %in% J$sets]
  add <- list("ECG" = "ECG", "CTPA" = "CTPA", "ECG + CTPA" = "ECG + CTPA")
  for (form in c("PH", "TV")) {
    K <- function(s) paste(s, form, sep = " | ")
    for (a in setdiff(vs, "EHR 28")) {
      out[[length(out) + 1]] <- list(
        name = sprintf("%s minus EHR 28, alone", a), form = form,
        plus = K(a), minus = K("EHR 28"))
    }
    if ("EHR 79 + MI" %in% vs) {
      out[[length(out) + 1]] <- list(
        name = "EHR 79 + MI minus EHR 79, alone", form = form,
        plus = K("EHR 79 + MI"), minus = K("EHR 79"))
    }
    for (ad in names(add)) {
      for (v in vs) {
        f <- paste(v, add[[ad]], sep = " + ")
        if (!(f %in% J$sets)) next
        out[[length(out) + 1]] <- list(
          name = sprintf("gain from %s with %s", ad, v), form = form,
          plus = K(f), minus = K(v))
      }
      for (v in setdiff(vs, "EHR 28")) {
        f1 <- paste(v, add[[ad]], sep = " + ")
        f0 <- paste("EHR 28", add[[ad]], sep = " + ")
        if (!(f1 %in% J$sets) || !(f0 %in% J$sets)) next
        out[[length(out) + 1]] <- list(
          name = sprintf("gain from %s: with %s minus with EHR 28", ad, v),
          form = form, plus = c(K(f1), K("EHR 28")),
          minus = c(K(v), K(f0)))
        out[[length(out) + 1]] <- list(
          name = sprintf("%s + %s minus EHR 28 + %s", v, ad, ad),
          form = form, plus = K(f1), minus = K(f0))
      }
    }
  }
  out
}
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
MET <- list()
CON <- list()
RT <- list()
for (J in JOBS) {
  R <- if (J$co == "EHR+ECG") RAT else NULL
  pt <- stat_vec(J, seq_len(J$n))
  if (!is.null(R)) pt <- c(pt, ratio_vec(R, rep(1, J$n)))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    R = R, future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  %s cohort done\n", J$co))
  for (k in grep("\\|(auc|ipa|ll)$", names(pt), value = TRUE)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      cohort = J$co, model = trimws(z[1]), form = trimws(z[2]),
      horizon = as.numeric(z[3]), measure = z[4], estimate = pt[[k]],
      lo = q[1], hi = q[2])
  }
  for (k in names(J$P)) {
    z <- strsplit(k, " | ", fixed = TRUE)[[1]]
    for (hi in seq_along(HZ)) {
      p <- pmin(pmax(J$P[[k]][, hi], 1e-6), 1 - 1e-6)
      sl <- tryCatch(unname(coef(glm(J$Y[, hi] ~ qlogis(p),
                                     family = binomial))[2]),
                     error = function(e) NA_real_)
      MET[[length(MET) + 1]] <- data.table(
        cohort = J$co, model = z[1], form = z[2], horizon = HZ[hi],
        measure = "slope", estimate = sl, lo = NA_real_, hi = NA_real_)
    }
  }
  for (g in gainpairs(J)) for (h in HZ) for (ms in c("auc", "ipa", "ll")) {
    kk <- function(x) paste(x, h, ms, sep = "|")
    if (!all(kk(c(g$plus, g$minus)) %in% names(pt))) next
    est <- sum(pt[kk(g$plus)]) - sum(pt[kk(g$minus)])
    bv <- rowSums(B[, kk(g$plus), drop = FALSE]) -
      rowSums(B[, kk(g$minus), drop = FALSE])
    q <- q2(bv)
    CON[[length(CON) + 1]] <- data.table(
      cohort = J$co, contrast = g$name, form = g$form, horizon = h,
      measure = ms, estimate = est, lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0)
  }
  if (!is.null(R)) {
    for (k in grep("log (HR|ratio)", names(pt), value = TRUE)) {
      z <- strsplit(k, "|", fixed = TRUE)[[1]]
      q <- q2(B[, k])
      RT[[length(RT) + 1]] <- data.table(
        ehr_version = z[1], quantity = sub("log ", "", z[2]),
        estimate = exp(pt[[k]]), lo = exp(q[1]), hi = exp(q[2]))
    }
    for (v in setdiff(names(R$Z), "EHR 28")) {
      a <- paste(v, "ECG log ratio", sep = "|")
      b <- paste("EHR 28", "ECG log ratio", sep = "|")
      q <- q2(B[, a] - B[, b])
      RT[[length(RT) + 1]] <- data.table(
        ehr_version = sprintf("%s over EHR 28", v),
        quantity = "ECG ratio of ratios", estimate = exp(pt[[a]] - pt[[b]]),
        lo = exp(q[1]), hi = exp(q[2]))
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
rt <- rbindlist(RT)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "ab_metrics.csv"))
fwrite(cn, file.path(RS, "ab_contrasts.csv"))
fwrite(rt, file.path(RS, "ab_ecg_ratio.csv"))
fwrite(ck, file.path(RS, "ab_checks.csv"))
cat("\n=== ECG EARLY / LATE HAZARD RATIOS (EHR+ECG cohort) ===\n")
print(rt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)],
      row.names = FALSE)
cat("\n=== KEY CONTRASTS, AUROC, PH STACK ===\n")
print(cn[measure == "auc" & form == "PH",
         .(cohort, contrast, horizon, est = round(estimate, 4),
           lo = round(lo, 4), hi = round(hi, 4), excludes_zero)],
      row.names = FALSE, nrows = 2000)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (ab_scores_*.csv, patient-level, local only)\n")
for (f in c("ab_metrics.csv", "ab_contrasts.csv", "ab_ecg_ratio.csv",
            "ab_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "ab_log.txt"), " (terminal log)\n")
cat("\nOutputs: ab_metrics.csv, ab_contrasts.csv, ab_ecg_ratio.csv,",
    "ab_checks.csv\n")