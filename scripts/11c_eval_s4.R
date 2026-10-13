# 11c_eval_s4.R
# S4, learner families (PREREGISTRATION.md section 6; DEVIATIONS.md D11).
#
# Compares ridge-Cox with random survival forests (RSF) and their rank
# average, and L2 logistic regression with random forests and their rank
# average, per modality and in late fusion built on each family. Neural
# learners are evaluated in stage 11e.
#
# Level 1 takes the seed-averaged out-of-fold scores of stages 6a and 6b
# for each frame, cohort and modality, at horizons of 10 and 30 days
# (Frame A) and 10, 30, 90, 180 and 365 days (Frame B). Endpoints are
# death, composite and CV readmission (first event: a death first is
# censored for the C-index and counts as no event for the labels, as in
# the cs_cv and cvu targets). The survival family is ridge-Cox, RSF and
# their rank average (mean of the two normalised ranks); the binary family
# is L2 logistic, random forest and their rank average, one model per
# horizon. Discrimination is AUROC at every horizon (with complete
# follow-up, the cumulative/dynamic AUC) and, for survival scores,
# Harrell's C over the frame's follow-up (Uno's C truncated there when
# follow-up is complete). Models with absolute risks (survival risks for
# death and composite, every binary model) also get IPA (event-rate
# reference, D4), log-likelihood skill score, observed/expected and
# calibration slope (point estimates).
#
# Fused models (death, each cohort's full modality set) are a Cox stack
# per family, seed and outer fold on rank-based normal scores of the
# modality scores, fitted on the training fold (the same transform for
# every family), with Breslow risks at every horizon, averaged over seeds.
#
# Contrasts are forest minus linear and rank average minus linear, within
# each family and for the fused models: AUROC and C, plus IPA and
# log-likelihood where both models have risks. Intervals come from 2,000
# bootstrap samples (seed 42), clustered by patient in Frame A (D4). Each
# AUROC contrast is also given as the mean and SD of the per-fold
# differences over seeds and folds (D2). The share of tuning choices at
# either end of the stage 6a and 6b grids is tabulated (D11).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          data/processed/level1/oof_linear_*.csv, oof_forest_*.csv,
#          params_linear.csv, params_forest.csv
# Outputs: results/s4_metrics.csv, s4_contrasts.csv, s4_perfold.csv,
#          s4_tuning.csv, s4_checks.csv, s4_log.txt (terminal log)
# Usage:   Rscript scripts/11c_eval_s4.R [--probe-only]

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
    fwrite(ck, file.path(RS, "s4_checks.csv"))
    cat("\nSTOP: a check failed. see: s4_checks.csv\n")
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
BASE <- c(survival = "ridge-Cox", binary = "L2 logistic",
          fused = "ridge-Cox")
TGT <- c(death = "death", composite = "composite", cv = "cs_cv")
LAB <- c(death = "death", composite = "composite", cv = "cvu")

# items for one data source; g(name) returns a column aligned to the
# cohort; fz is a list (family -> n x horizons risk matrix) or NULL
build_items <- function(g, mods, hz, fz) {
  AU <- list()
  CI <- list()
  RK <- list()
  for (m in mods) for (en in names(TGT)) {
    u <- toupper(m)
    lp <- g(sprintf("%s.cox.%s.lp", m, TGT[[en]]))
    sf <- g(sprintf("%s.rsf.%s.score", m, TGT[[en]]))
    ss <- list("ridge-Cox" = lp, "RSF" = sf,
               "rank average" = (nrank(lp) + nrank(sf)) / 2)
    for (lr in names(ss)) {
      CI[[paste(u, en, "survival", lr, sep = "|")]] <-
        list(s = ss[[lr]], en = en)
      for (h in hz) {
        AU[[paste(u, en, "survival", lr, h, sep = "|")]] <-
          list(s = ss[[lr]], y = paste(en, h, sep = "|"))
      }
    }
    for (h in hz) {
      if (en != "cv") {
        RK[[paste(u, en, "survival", "ridge-Cox", h, sep = "|")]] <-
          list(p = g(sprintf("%s.cox.%s.r%d", m, TGT[[en]], h)),
               y = paste(en, h, sep = "|"))
        RK[[paste(u, en, "survival", "RSF", h, sep = "|")]] <-
          list(p = g(sprintf("%s.rsf.%s.r%d", m, TGT[[en]], h)),
               y = paste(en, h, sep = "|"))
      }
      pl <- g(sprintf("%s.lr.%s_h%d", m, LAB[[en]], h))
      pf <- g(sprintf("%s.rf.%s_h%d", m, LAB[[en]], h))
      bb <- list("L2 logistic" = pl, "random forest" = pf,
                 "rank average" = (nrank(pl) + nrank(pf)) / 2)
      for (lr in names(bb)) {
        AU[[paste(u, en, "binary", lr, h, sep = "|")]] <-
          list(s = bb[[lr]], y = paste(en, h, sep = "|"))
      }
      for (lr in c("L2 logistic", "random forest")) {
        RK[[paste(u, en, "binary", lr, h, sep = "|")]] <-
          list(p = bb[[lr]], y = paste(en, h, sep = "|"))
      }
    }
  }
  if (!is.null(fz)) {
    un <- paste("fused", paste(toupper(mods), collapse = "+"))
    for (fam in names(fz)) for (i in seq_along(hz)) {
      k <- paste(un, "death", "fused", fam, hz[i], sep = "|")
      AU[[k]] <- list(s = fz[[fam]][, i], y = paste("death", hz[i],
                                                    sep = "|"))
      RK[[k]] <- list(p = fz[[fam]][, i], y = paste("death", hz[i],
                                                    sep = "|"))
    }
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
    v[paste(k, "oe", sep = "|")] <- pr / mean(p)
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

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
JOBS <- list()
PERFOLD <- list()
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
  ok <- file.exists(fl) && file.exists(ff)
  chk(sprintf("%s: stage 6a and 6b files found", key),
      paste(basename(c(fl, ff)), collapse = ", "),
      if (ok) "PASS" else "FAIL")
  stop_if_fail()
  need_l <- c()
  need_f <- c()
  for (m in mods) {
    for (en in names(TGT)) {
      need_l <- c(need_l, sprintf("%s.cox.%s.lp", m, TGT[[en]]),
                  sprintf("%s.lr.%s_h%d", m, LAB[[en]], hz))
      need_f <- c(need_f, sprintf("%s.rsf.%s.score", m, TGT[[en]]),
                  sprintf("%s.rf.%s_h%d", m, LAB[[en]], hz))
      if (en != "cv") {
        need_l <- c(need_l, sprintf("%s.cox.%s.r%d", m, TGT[[en]], hz))
        need_f <- c(need_f, sprintf("%s.rsf.%s.r%d", m, TGT[[en]], hz))
      }
    }
  }
  head_l <- names(fread(fl, nrows = 0))
  head_f <- names(fread(ff, nrows = 0))
  ml <- setdiff(c("seed", "hadm_id", need_l), head_l)
  mf <- setdiff(c("seed", "hadm_id", need_f), head_f)
  chk(sprintf("%s: columns needed / missing", key),
      sprintf("%d / %d%s", length(need_l) + length(need_f),
              length(ml) + length(mf),
              if (length(ml) + length(mf))
                paste0(" (e.g. ", paste(head(c(ml, mf), 3),
                                        collapse = ", "), ")") else ""),
      if (length(ml) + length(mf) == 0) "PASS" else "FAIL")
  stop_if_fail()
  LT <- fread(fl, select = c("seed", "hadm_id", need_l))
  FT <- fread(ff, select = c("seed", "hadm_id", need_f))
  fo <- folds[frame == fr & cohort == co]
  sl <- sort(unique(LT$seed))
  sf <- sort(unique(FT$seed))
  so <- sort(unique(fo$seed))
  chk(sprintf("%s: seeds in 6a / 6b / folds", key),
      sprintf("%s / %s / %s", paste(sl, collapse = ","),
              paste(sf, collapse = ","), paste(so, collapse = ",")),
      if (identical(sl, sf) && identical(sl, so)) "PASS" else "FAIL")
  rows_ok <- nrow(LT) == n * length(sl) && nrow(FT) == n * length(sf)
  chk(sprintf("%s: rows = admissions x seeds (6a, 6b)", key),
      sprintf("%d, %d = %d x %d", nrow(LT), nrow(FT), n, length(sl)),
      if (rows_ok) "PASS" else "FAIL")
  stop_if_fail()
  seeds <- sl
  # per-seed tables aligned to the cohort
  per_seed <- lapply(seeds, function(s) {
    a <- LT[LT$seed == s][match(sub$hadm_id, hadm_id)]
    b <- FT[FT$seed == s][match(sub$hadm_id, hadm_id)]
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    list(L = a, F = b, fold = f)
  })
  names(per_seed) <- as.character(seeds)
  la <- LT[, lapply(.SD, mean), by = hadm_id, .SDcols = need_l][
    match(sub$hadm_id, hadm_id)]
  fa <- FT[, lapply(.SD, mean), by = hadm_id, .SDcols = need_f][
    match(sub$hadm_id, hadm_id)]
  nmiss <- sum(is.na(la[, ..need_l])) + sum(is.na(fa[, ..need_f]))
  nfold <- sum(sapply(per_seed, function(z) sum(is.na(z$fold))))
  chk(sprintf("%s: missing scores after averaging / missing folds", key),
      sprintf("%d / %d", nmiss, nfold),
      if (nmiss == 0 && nfold == 0) "PASS" else "FAIL")
  stop_if_fail()
  getter <- function(A, B) {
    function(nm) if (nm %in% names(A)) A[[nm]] else B[[nm]]
  }
  # fused death models, per seed and fold, for each family
  td <- sub$time_death
  ed <- as.integer(sub$event_death == 1)
  FZ <- list()
  for (s in names(per_seed)) {
    z <- per_seed[[s]]
    g <- getter(z$L, z$F)
    sc <- list("ridge-Cox" = list(), "RSF" = list(),
               "rank average" = list())
    for (m in mods) {
      lp <- g(sprintf("%s.cox.death.lp", m))
      rf <- g(sprintf("%s.rsf.death.score", m))
      sc[["ridge-Cox"]][[m]] <- lp
      sc[["RSF"]][[m]] <- rf
      sc[["rank average"]][[m]] <- (nrank(lp) + nrank(rf)) / 2
    }
    FZ[[s]] <- lapply(sc, function(S) {
      R <- matrix(NA_real_, n, length(hz))
      for (k in sort(unique(z$fold))) {
        tr <- which(z$fold != k)
        te <- which(z$fold == k)
        Ztr <- sapply(mods, function(m) nscore(S[[m]][tr], S[[m]][tr]))
        Zte <- sapply(mods, function(m) nscore(S[[m]][te], S[[m]][tr]))
        Ztr <- matrix(Ztr, ncol = length(mods))
        Zte <- matrix(Zte, ncol = length(mods))
        D <- data.frame(t = td[tr], e = ed[tr], Ztr)
        b <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
        b[is.na(b)] <- 0
        lt <- as.vector(Ztr %*% b)
        lv <- as.vector(Zte %*% b)
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
  nfz <- sum(sapply(fz_avg, function(x) sum(!is.finite(x))))
  chk(sprintf("%s: fused risks missing", key), nfz,
      if (nfz == 0) "PASS" else "FAIL")
  stop_if_fail()
  IT <- build_items(getter(la, fa), mods, hz, fz_avg)
  Y <- list()
  TE <- list(
    death = list(t = sub$time_death, e = sub$event_death == 1),
    composite = list(t = sub$time_first, e = sub$event_comp == 1),
    cv = list(t = sub$time_first, e = sub$event_cr == 1))
  for (en in names(TE)) for (h in hz) {
    Y[[paste(en, h, sep = "|")]] <-
      as.integer(TE[[en]]$e & TE[[en]]$t <= h + 0.5)
  }
  # sanity: every death score ranks better than chance at the last horizon
  hlast <- max(hz)
  bad <- 0L
  for (k in names(IT$AU)) {
    it <- IT$AU[[k]]
    if (it$y == paste("death", hlast, sep = "|")) {
      a <- auc_fast(Y[[it$y]], it$s)
      if (!is.finite(a) || a <= 0.5) bad <- bad + 1L
    }
  }
  chk(sprintf("%s: death scores at or below chance at %d d", key, hlast),
      bad,
      if (bad == 0) "PASS" else "WARN")
  JOBS[[key]] <- list(fr = fr, co = co, n = n, IT = IT, Y = Y, TE = TE,
                      rows_by = split(seq_len(n),
                                      match(sub$subject_id,
                                            unique(sub$subject_id))))
  # D2 per-fold AUROCs from per-seed predictions
  pf <- list()
  for (s in names(per_seed)) {
    z <- per_seed[[s]]
    its <- build_items(getter(z$L, z$F), mods, hz, FZ[[s]])
    for (k in sort(unique(z$fold))) {
      te <- which(z$fold == k)
      a <- sapply(its$AU, function(it) auc_fast(Y[[it$y]][te], it$s[te]))
      pf[[length(pf) + 1]] <- data.table(seed = s, fold = k,
                                         item = names(a), auc = a)
    }
  }
  PERFOLD[[key]] <- rbindlist(pf)
  cat(sprintf("  %s: n=%d, %d seeds, %d AUROC, %d C and %d risk items\n",
              key, n, length(seeds), length(IT$AU), length(IT$CI),
              length(IT$RK)))
}

# ---------------- tuning edges (stages 6a and 6b) ----------------
cat("\n=== TUNING CHOICES AT THE GRID ENDS ===\n")
TUN <- list()
pl <- file.path(L1, "params_linear.csv")
if (file.exists(pl)) {
  P <- fread(pl)
  if (!("param" %in% names(P))) P[, param := NA_real_]
  mc <- intersect(c("modality", "mod"), names(P))[1]
  if (is.na(mc)) P[, modality := "all"] else setnames(P, mc, "modality")
  P[, lg := round(log10(as.numeric(param)), 6)]
  TUN[[1]] <- P[learner %in% c("cox", "lr"), .(
    choices = .N,
    share_least_regularised = mean(lg == 0),
    share_most_regularised = mean(lg == ifelse(learner == "cox", 4, -4)),
    share_larger_mtry = NA_real_),
    by = .(learner, modality)][, stage := "6a"]
} else {
  chk("stage 6a tuning file (params_linear.csv)", "not found", "INFO")
}
pf2 <- file.path(L1, "params_forest.csv")
if (file.exists(pf2)) {
  Q <- fread(pf2)
  Q[, big := mtry == max(mtry, na.rm = TRUE), by = .(modality)]
  TUN[[2]] <- Q[, .(
    choices = .N,
    share_least_regularised = mean(nodesize == 5, na.rm = TRUE),
    share_most_regularised = mean(nodesize == 40, na.rm = TRUE),
    share_larger_mtry = mean(big, na.rm = TRUE)),
    by = .(learner, modality)][, stage := "6b"]
} else {
  chk("stage 6b tuning file (params_forest.csv)", "not found", "INFO")
}
tun <- rbindlist(TUN, use.names = TRUE)
if (nrow(tun)) setcolorder(tun, "stage")
fwrite(tun, file.path(RS, "s4_tuning.csv"))
print(tun[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)],
      row.names = FALSE)

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
PFO <- list()
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
parse_keys <- function(k) {
  z <- tstrsplit(k, "|", fixed = TRUE)
  data.table(id = k, unit = z[[1]], endpoint = z[[2]], family = z[[3]],
             learner = z[[4]], horizon = suppressWarnings(
               as.numeric(z[[5]])), measure = z[[6]])
}
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
  # calibration slope, point estimate
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
  # contrasts: each learner minus its family's base learner
  KK <- K[measure %in% c("auc", "C", "ipa", "ll")]
  for (i in which(KK$learner != BASE[KK$family])) {
    a <- KK[i]
    bk <- paste(a$unit, a$endpoint, a$family, BASE[[a$family]],
                ifelse(is.na(a$horizon), "NA", a$horizon), a$measure,
                sep = "|")
    if (!(bk %in% names(pt))) next
    bv <- B[, a$id] - B[, bk]
    q <- q2(bv)
    CON[[length(CON) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, unit = a$unit, endpoint = a$endpoint,
      family = a$family,
      contrast = sprintf("%s - %s", a$learner, BASE[[a$family]]),
      horizon = a$horizon, measure = a$measure,
      estimate = pt[[a$id]] - pt[[bk]], lo = q[1], hi = q[2],
      excludes_zero = q[1] > 0 | q[2] < 0)
  }
  # D2 per-fold differences
  PF <- PERFOLD[[paste(J$fr, J$co)]]
  PF <- cbind(PF, parse_keys(paste(PF$item, "auc", sep = "|"))[, -1])
  bas <- PF[learner == BASE[family], .(seed, fold, unit, endpoint, family,
                                       horizon, base = auc)]
  oth <- PF[learner != BASE[family]]
  mm <- merge(oth, bas, by = c("seed", "fold", "unit", "endpoint",
                               "family", "horizon"))
  PFO[[length(PFO) + 1]] <- mm[, .(
    folds = .N, mean_diff = mean(auc - base, na.rm = TRUE),
    sd_diff = sd(auc - base, na.rm = TRUE)),
    by = .(unit, endpoint, family, learner, horizon)][
      , `:=`(frame = J$fr, cohort = J$co)]
}
plan(sequential)
mt <- rbindlist(MET)
setcolorder(mt, c("frame", "cohort"))
mt[, id := NULL]
cn <- rbindlist(CON)
pfo <- rbindlist(PFO)
pfo[, contrast := sprintf("%s - %s", learner, BASE[family])]
pfo <- pfo[, .(frame, cohort, unit, endpoint, family, contrast, horizon,
               folds, mean_diff, sd_diff)]
fwrite(mt, file.path(RS, "s4_metrics.csv"))
fwrite(cn, file.path(RS, "s4_contrasts.csv"))
fwrite(pfo, file.path(RS, "s4_perfold.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "s4_checks.csv"))

cat("\n=== DEATH: AUROC AT 30 DAYS AND THE LAST HORIZON ===\n")
show <- mt[endpoint == "death" & measure == "auc" &
             (horizon == 30 | (frame == "B" & horizon == 365))]
print(dcast(show, frame + cohort + unit + family + learner ~ horizon,
            value.var = "estimate")[
              , lapply(.SD, function(x) if (is.numeric(x)) round(x, 4)
                       else x)], row.names = FALSE)
cat("\n=== DEATH CONTRASTS (AUROC and C) ===\n")
print(cn[endpoint == "death" & measure %in% c("auc", "C"),
         .(frame, cohort, unit, family, contrast, horizon, measure,
           est = round(estimate, 4), lo = round(lo, 4), hi = round(hi, 4),
           excludes_zero)], row.names = FALSE, nrows = 2000)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("s4_metrics.csv", "s4_contrasts.csv", "s4_perfold.csv",
            "s4_tuning.csv", "s4_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "s4_log.txt"), " (terminal log)\n")
cat("\nOutputs: s4_metrics.csv, s4_contrasts.csv, s4_perfold.csv,",
    "s4_tuning.csv, s4_checks.csv\n")