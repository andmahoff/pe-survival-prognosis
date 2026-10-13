# 11d_eval_s4_wide.R
# S4 sensitivity analysis with the wide-grid forests (DEVIATIONS.md D11,
# D13; stage 6c).
#
# Asks whether the S4 conclusions (forests against ridge-Cox, and the fused
# models) change when the forests are tuned over the wider grid. The stage
# 6a ridge-Cox and L2 logistic scores are unchanged. For every modality,
# the survival family (death, composite and CV readmission) is ridge-Cox,
# RSF (stage 6b), RSF wide (stage 6c) and the rank average of ridge-Cox and
# RSF wide. The binary family covers the death labels only (the D13 scope
# of stage 6c): L2 logistic, random forest (6b), random forest wide (6c)
# and the rank average with the wide forest. The fused death models are a
# Cox stack per seed and fold on rank-based normal scores of the modality
# scores, as in S4, for ridge-Cox, RSF, RSF wide and the wide rank average.
# Endpoints or labels missing from the stage 6c files are skipped.
#
# Measures are AUROC at every horizon, Harrell's C over the follow-up, IPA
# and log-likelihood skill score for every model with absolute risks, and
# calibration slope (point estimate). Contrasts are the wide forest minus
# the linear learner and the wide rank average minus the linear learner
# (the S4 contrasts with the wide forests), and the wide forest minus the
# registered forest (what the wider grid changed). Intervals come from
# 2,000 bootstrap samples (seed 42), clustered by patient in Frame A.
#
# Checks: inputs found and aligned (FAIL); the registered RSF and ridge-Cox
# AUROCs, single-modality and fused, reproduce those of stage 11c (FAIL
# beyond 1e-9; WARN if results/s4_metrics.csv is missing).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          data/processed/level1/oof_linear_*.csv, oof_forest_*.csv,
#          oof_forest_wide_*.csv, results/s4_metrics.csv (for one check)
# Outputs: results/s4w_metrics.csv, s4w_contrasts.csv, s4w_checks.csv,
#          s4w_log.txt (terminal log)
# Usage:   Rscript scripts/11d_eval_s4_wide.R [--probe-only]

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
RS <- file.path(P2, "results")
NBOOT <- 2000L
SEED <- 42L
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
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
    fwrite(ck, file.path(RS, "s4w_checks.csv"))
    cat("\nSTOP: a check failed. see: s4w_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- helpers ----------------
auc_fast <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- data.table::frank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cidx <- function(t, e, s) {
  out <- tryCatch(survival::concordance(
    survival::Surv(t, e) ~ s, reverse = TRUE)$concordance,
    error = function(err) NA_real_)
  as.numeric(out)
}
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
nrank <- function(x) data.table::frank(x) / length(x)
TGT <- c(death = "death", composite = "composite", cv = "cs_cv")
PAIRS <- list(
  survival = list(c("RSF wide", "ridge-Cox"),
                  c("rank average wide", "ridge-Cox"),
                  c("RSF wide", "RSF")),
  binary = list(c("random forest wide", "L2 logistic"),
                c("rank average wide", "L2 logistic"),
                c("random forest wide", "random forest")),
  fused = list(c("RSF wide", "ridge-Cox"),
               c("rank average wide", "ridge-Cox"),
               c("RSF wide", "RSF")))

# items for one data source; g(name) returns a column aligned to the
# cohort (NULL if absent); fz is a list (family -> n x horizons risks)
build_items <- function(g, mods, hz, fz) {
  AU <- list()
  CI <- list()
  RK <- list()
  for (m in mods) for (en in names(TGT)) {
    u <- toupper(m)
    lp <- g(sprintf("%s.cox.%s.lp", m, TGT[[en]]))
    rs <- g(sprintf("%s.rsf.%s.score", m, TGT[[en]]))
    ws <- g(sprintf("%s.rsf.%s.score", m, TGT[[en]]), wide = TRUE)
    if (is.null(ws)) next
    ss <- list("ridge-Cox" = lp, "RSF" = rs, "RSF wide" = ws,
               "rank average wide" = (nrank(lp) + nrank(ws)) / 2)
    for (lr in names(ss)) {
      CI[[paste(u, en, "survival", lr, sep = "|")]] <-
        list(s = ss[[lr]], en = en)
      for (h in hz) {
        AU[[paste(u, en, "survival", lr, h, sep = "|")]] <-
          list(s = ss[[lr]], y = paste(en, h, sep = "|"))
      }
    }
    if (en != "cv") for (h in hz) {
      pr <- list("ridge-Cox" = g(sprintf("%s.cox.%s.r%d", m, TGT[[en]], h)),
                 "RSF" = g(sprintf("%s.rsf.%s.r%d", m, TGT[[en]], h)),
                 "RSF wide" = g(sprintf("%s.rsf.%s.r%d", m, TGT[[en]], h),
                                wide = TRUE))
      for (lr in names(pr)) {
        RK[[paste(u, en, "survival", lr, h, sep = "|")]] <-
          list(p = pr[[lr]], y = paste(en, h, sep = "|"))
      }
    }
  }
  for (m in mods) for (h in hz) {
    u <- toupper(m)
    pw <- g(sprintf("%s.rf.death_h%d", m, h), wide = TRUE)
    if (is.null(pw)) next
    pl <- g(sprintf("%s.lr.death_h%d", m, h))
    pf <- g(sprintf("%s.rf.death_h%d", m, h))
    bb <- list("L2 logistic" = pl, "random forest" = pf,
               "random forest wide" = pw,
               "rank average wide" = (nrank(pl) + nrank(pw)) / 2)
    for (lr in names(bb)) {
      AU[[paste(u, "death", "binary", lr, h, sep = "|")]] <-
        list(s = bb[[lr]], y = paste("death", h, sep = "|"))
    }
    for (lr in c("L2 logistic", "random forest", "random forest wide")) {
      RK[[paste(u, "death", "binary", lr, h, sep = "|")]] <-
        list(p = bb[[lr]], y = paste("death", h, sep = "|"))
    }
  }
  un <- paste("fused", paste(toupper(mods), collapse = "+"))
  for (fam in names(fz)) for (i in seq_along(hz)) {
    k <- paste(un, "death", "fused", fam, hz[i], sep = "|")
    AU[[k]] <- list(s = fz[[fam]][, i], y = paste("death", hz[i],
                                                  sep = "|"))
    RK[[k]] <- list(p = fz[[fam]][, i], y = paste("death", hz[i],
                                                  sep = "|"))
  }
  list(AU = AU, CI = CI, RK = RK)
}
stat_vec <- function(J, idx) {
  v <- numeric(0)
  Y <- lapply(J$Y, function(y) y[idx])
  for (k in names(J$IT$AU)) {
    it <- J$IT$AU[[k]]
    v[paste(k, "auc", sep = "|")] <- auc_fast(Y[[it$y]], it$s[idx])
  }
  for (k in names(J$IT$CI)) {
    it <- J$IT$CI[[k]]
    te <- J$TE[[it$en]]
    v[paste(k, "NA", "C", sep = "|")] <-
      cidx(te$t[idx], te$e[idx], it$s[idx])
  }
  for (k in names(J$IT$RK)) {
    it <- J$IT$RK[[k]]
    y <- Y[[it$y]]
    p <- pmin(pmax(it$p[idx], 1e-6), 1 - 1e-6)
    pr <- mean(y)
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    v[paste(k, "ipa", sep = "|")] <- 1 - mean((y - p)^2) / (pr * (1 - pr))
    v[paste(k, "ll", sep = "|")] <-
      1 - (-mean(y * log(p) + (1 - y) * log(1 - p))) / l0
  }
  v
}
boot_rep <- function(i, J) {
  if (J$fr == "A") {
    idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                       length(J$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(J$n, J$n, replace = TRUE)
  }
  stat_vec(J, idx)
}
parse_keys <- function(k) {
  z <- tstrsplit(k, "|", fixed = TRUE)
  data.table(id = k, unit = z[[1]], endpoint = z[[2]], family = z[[3]],
             learner = z[[4]], horizon = suppressWarnings(
               as.numeric(z[[5]])), measure = z[[6]])
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
JOBS <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  key <- paste(fr, co)
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  n <- nrow(sub)
  hz <- HZ[[fr]]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  fl <- file.path(L1, sprintf("oof_linear_%s_%s.csv", fr, code))
  ff <- file.path(L1, sprintf("oof_forest_%s_%s.csv", fr, code))
  fw <- file.path(L1, sprintf("oof_forest_wide_%s_%s.csv", fr, code))
  ok <- all(file.exists(c(fl, ff, fw)))
  chk(sprintf("%s: stage 6a, 6b and 6c files found", key),
      paste(basename(c(fl, ff, fw)), collapse = ", "),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  LT <- fread(fl)
  FT <- fread(ff)
  WT <- fread(fw)
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(LT$seed))
  same <- identical(seeds, sort(unique(FT$seed))) &&
    identical(seeds, sort(unique(WT$seed))) &&
    identical(seeds, sort(unique(fo$seed))) &&
    nrow(LT) == n * length(seeds) && nrow(FT) == n * length(seeds) &&
    nrow(WT) == n * length(seeds)
  chk(sprintf("%s: seeds and rows aligned (6a, 6b, 6c, folds)", key),
      sprintf("%d seeds, %d x %d rows", length(seeds), n, length(seeds)),
      if (same) "PASS" else "FAIL")
  stop_if_fail()
  num <- function(D) setdiff(names(D), c("seed", "subject_id", "hadm_id",
                                         "fold"))
  avg <- function(D) {
    D[, lapply(.SD, mean), by = hadm_id, .SDcols = num(D)][
      match(sub$hadm_id, hadm_id)]
  }
  la <- avg(LT)
  fa <- avg(FT)
  wa <- avg(WT)
  getter <- function(A, B, W) {
    function(nm, wide = FALSE) {
      if (wide) return(if (nm %in% names(W)) W[[nm]] else NULL)
      if (nm %in% names(A)) A[[nm]] else B[[nm]]
    }
  }
  # fused death models, per seed and fold, as in S4
  td <- sub$time_death
  ed <- as.integer(sub$event_death == 1)
  FZ <- list()
  for (s in seeds) {
    a <- LT[LT$seed == s][match(sub$hadm_id, hadm_id)]
    b <- FT[FT$seed == s][match(sub$hadm_id, hadm_id)]
    w <- WT[WT$seed == s][match(sub$hadm_id, hadm_id)]
    fold <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    sc <- list("ridge-Cox" = list(), "RSF" = list(), "RSF wide" = list(),
               "rank average wide" = list())
    for (m in mods) {
      lp <- a[[sprintf("%s.cox.death.lp", m)]]
      sc[["ridge-Cox"]][[m]] <- lp
      sc[["RSF"]][[m]] <- b[[sprintf("%s.rsf.death.score", m)]]
      sc[["RSF wide"]][[m]] <- w[[sprintf("%s.rsf.death.score", m)]]
      sc[["rank average wide"]][[m]] <-
        (nrank(lp) + nrank(sc[["RSF wide"]][[m]])) / 2
    }
    FZ[[as.character(s)]] <- lapply(sc, function(S) {
      R <- matrix(NA_real_, n, length(hz))
      for (k in sort(unique(fold))) {
        tr <- which(fold != k)
        te <- which(fold == k)
        Ztr <- sapply(mods, function(m) nscore(S[[m]][tr], S[[m]][tr]))
        Zte <- sapply(mods, function(m) nscore(S[[m]][te], S[[m]][tr]))
        Ztr <- matrix(Ztr, ncol = length(mods))
        Zte <- matrix(Zte, ncol = length(mods))
        D <- data.frame(t = td[tr], e = ed[tr], Ztr)
        bb <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
        bb[is.na(bb)] <- 0
        lt <- as.vector(Ztr %*% bb)
        lv <- as.vector(Zte %*% bb)
        H <- bres_cuts(td[tr], ed[tr], lt, hz + 0.5)
        for (i in seq_along(hz)) {
          R[te, i] <- -expm1(-H[i] * exp(lv - max(lt)))
        }
      }
      R
    })
  }
  fz_avg <- lapply(names(FZ[[1]]), function(fam) {
    Reduce(`+`, lapply(FZ, function(x) x[[fam]])) / length(FZ)
  })
  names(fz_avg) <- names(FZ[[1]])
  IT <- build_items(getter(la, fa, wa), mods, hz, fz_avg)
  nbad <- sum(sapply(IT$AU, function(it) sum(!is.finite(it$s)))) +
    sum(sapply(IT$RK, function(it) sum(!is.finite(it$p))))
  chk(sprintf("%s: missing scores or risks", key), nbad,
      if (nbad == 0) "PASS" else "FAIL")
  stop_if_fail()
  ends <- unique(sapply(strsplit(names(IT$CI), "|", fixed = TRUE),
                        `[`, 2))
  nbin <- sum(grepl("\\|binary\\|", names(IT$AU)))
  chk(sprintf("%s: survival endpoints / binary items with wide forests",
              key), sprintf("%s / %d", paste(ends, collapse = ", "), nbin),
      "INFO")
  TE <- list(
    death = list(t = sub$time_death, e = sub$event_death == 1),
    composite = list(t = sub$time_first, e = sub$event_comp == 1),
    cv = list(t = sub$time_first, e = sub$event_cr == 1))
  Y <- list()
  for (en in names(TE)) for (h in hz) {
    Y[[paste(en, h, sep = "|")]] <-
      as.integer(TE[[en]]$e & TE[[en]]$t <= h + 0.5)
  }
  JOBS[[key]] <- list(fr = fr, co = co, n = n, IT = IT, Y = Y, TE = TE,
                      mods = mods, hz = hz,
                      rows_by = split(seq_len(n),
                                      match(sub$subject_id,
                                            unique(sub$subject_id))))
}

# ---------------- check: S4 reproduced ----------------
s4f <- file.path(RS, "s4_metrics.csv")
if (!file.exists(s4f)) {
  chk("S4 reproduction: results/s4_metrics.csv", "missing", "WARN")
} else {
  s4 <- fread(s4f)[measure == "auc"]
  dmax <- 0
  nm <- 0L
  for (J in JOBS) {
    pt <- stat_vec(J, seq_len(J$n))
    K <- parse_keys(names(pt))
    K[, est := pt[id]]
    K <- K[measure == "auc" & endpoint == "death" &
             learner %in% c("ridge-Cox", "RSF") &
             family %in% c("survival", "fused")]
    K[, `:=`(frame = J$fr, cohort = J$co)]
    m <- merge(K, s4, by = c("frame", "cohort", "unit", "endpoint",
                             "family", "learner", "horizon", "measure"))
    nm <- nm + nrow(m)
    if (nrow(m)) dmax <- max(dmax, abs(m$est - m$estimate))
  }
  chk("S4 reproduction: registered RSF and ridge-Cox AUROCs (n, max)",
      sprintf("%d, %.1e", nm, dmax),
      if (nm > 0 && dmax <= 1e-9) "PASS" else "FAIL")
  stop_if_fail()
}

# ---------------- timing probe ----------------
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
plan(multisession, workers = workers)
nb <- workers * 2L
proj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nb), boot_rep, J = J, future.seed = 1L,
                          future.packages = PKG))
  proj <- proj + as.numeric(difftime(Sys.time(), t0, units = "secs")) *
    NBOOT / nb
}
cat(sprintf("  projected: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- bootstrap, metrics, contrasts ----------------
MET <- list()
CON <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
for (J in JOBS) {
  cat(sprintf("=== Frame %s %s: %d bootstrap samples ===\n", J$fr, J$co,
              NBOOT))
  pt <- stat_vec(J, seq_len(J$n))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  K <- parse_keys(names(pt))
  K[, `:=`(estimate = pt[id], lo = NA_real_, hi = NA_real_)]
  for (i in seq_len(nrow(K))) {
    q <- q2(B[, K$id[i]])
    K$lo[i] <- q[1]
    K$hi[i] <- q[2]
  }
  for (k in names(J$IT$RK)) {
    it <- J$IT$RK[[k]]
    y <- J$Y[[it$y]]
    p <- pmin(pmax(it$p, 1e-6), 1 - 1e-6)
    sl <- tryCatch(unname(coef(glm(y ~ qlogis(p), family = binomial))[2]),
                   error = function(e) NA_real_)
    z <- parse_keys(paste(k, "slope", sep = "|"))
    z[, `:=`(estimate = sl, lo = NA_real_, hi = NA_real_)]
    K <- rbind(K, z)
  }
  MET[[length(MET) + 1]] <- K[, `:=`(frame = J$fr, cohort = J$co)]
  KK <- K[measure %in% c("auc", "C", "ipa", "ll")]
  for (i in seq_len(nrow(KK))) {
    a <- KK[i]
    for (pp in PAIRS[[a$family]]) {
      if (a$learner != pp[1]) next
      bk <- paste(a$unit, a$endpoint, a$family, pp[2],
                  ifelse(is.na(a$horizon), "NA", a$horizon), a$measure,
                  sep = "|")
      if (!(bk %in% names(pt))) next
      q <- q2(B[, a$id] - B[, bk])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, unit = a$unit, endpoint = a$endpoint,
        family = a$family, contrast = sprintf("%s - %s", pp[1], pp[2]),
        horizon = a$horizon, measure = a$measure,
        estimate = pt[[a$id]] - pt[[bk]], lo = q[1], hi = q[2],
        excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
setcolorder(mt, c("frame", "cohort"))
mt[, id := NULL]
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "s4w_metrics.csv"))
fwrite(cn, file.path(RS, "s4w_contrasts.csv"))
fwrite(ck, file.path(RS, "s4w_checks.csv"))

cat("\n=== DEATH: WIDE FOREST MINUS RIDGE-COX AND MINUS REGISTERED FOREST",
    "(AUROC) ===\n")
print(cn[endpoint == "death" & measure == "auc" &
           grepl("^(RSF|random forest) wide", contrast) &
           (horizon == 30 | horizon == 365),
         .(frame, cohort, unit, family, contrast, horizon,
           est = round(estimate, 4), lo = round(lo, 4), hi = round(hi, 4),
           excludes_zero)], row.names = FALSE, nrows = 2000)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("s4w_metrics.csv", "s4w_contrasts.csv", "s4w_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "s4w_log.txt"), " (terminal log)\n")
cat("\nOutputs: s4w_metrics.csv, s4w_contrasts.csv, s4w_checks.csv\n")