# 15b_stability.R
# Stability of the predictions (Riley and Collins, Biometrical Journal 2023)
# and performance by subgroup (PREREGISTRATION_ADDENDUM.md, A3).
# Exploratory; intervals only.
#
# Part A, seed instability: the six seeds' out-of-fold Cox-stack and
# binary-stack death predictions (script 11b). Per patient, the mean
# absolute deviation of each seed's risk from the seed average, summarised
# overall and by risk decile, and the share of patients whose classification
# at 10% or 20% risk differs between seeds. Frame A at 30 days and Frame B
# at 30 and 365 days, in both cohorts.
#
# Part B, bootstrap instability of the fusion stage: the Cox stack and
# binary stack are refitted (including winsorising and standardising, D3) to
# the seed-averaged out-of-fold level-1 scores in the full cohort (the
# original model) and in 1,000 bootstrap samples (seed 42; patient-clustered
# in Frame A). Each bootstrap model predicts the original sample. Per
# patient: mean absolute prediction error (MAPE) against the original model,
# and classification instability at 10% and 20% (the share of bootstrap
# models that classify the patient differently). Per bootstrap model: ICI
# against the observed outcomes and net benefit at 5%, 10% and 20%.
# Instability-plot data are pooled percentiles of the bootstrap predictions
# in 20 bins of original risk (800 randomly chosen patients). Level-1 models
# are held fixed, so this measures the fusion stage only.
#
# Part C, subgroups: Frame A death at 30 days and Frame B death at 365 days,
# EHR+ECG cohort, Cox stack and binary stack. For age bands (<65, 65-79,
# 80+) and for sex if an EHR feature records it: n, events, observed and
# mean predicted risk, O/E, AUROC, calibration slope and ICI, with
# 1,000-sample bootstrap intervals (patient-clustered in Frame A).
#
# Notes:   a timing probe runs real bootstrap work for Parts B and C first.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv,
#          calibration/cal_preds_*.csv, level1/oof_linear_*.csv;
#          data/inputs/ehr28.csv
# Outputs: results/stab_seed.csv, stab_summary.csv, stab_bins.csv,
#          subgroups.csv, stab_log.txt
# Usage:   Rscript scripts/15b_stability.R [--probe-only]

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
RS <- file.path(P2, "results")
SDIR <- file.path(PR, "stability")
dir.create(SDIR, showWarnings = FALSE, recursive = TRUE)
NBOOT <- 1000L
NBOOT_SUB <- 1000L
SEED <- 42L
TH_CLS <- c(0.10, 0.20)
TH_NB <- c(0.05, 0.10, 0.20)
NSUB <- 800L
NKEEP <- 200L
HZ <- list(A = 30, B = c(30, 365))
MLAB <- c(coxstack = "Cox stack", binstack = "binary stack")
PKG <- c("survival", "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

avg <- function(path, cols) {
  x <- fread(path, select = c("hadm_id", cols))
  x[, lapply(.SD, mean), by = hadm_id, .SDcols = cols]
}
ici1 <- function(y, p) {
  f <- tryCatch(fitted(loess(y ~ p, span = 0.75, degree = 2)),
                error = function(e) rep(NA_real_, length(p)))
  mean(abs(f - p))
}
nb1 <- function(y, p, th) {
  n <- length(y)
  vapply(th, function(t) {
    tr <- p >= t
    (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / n
  }, numeric(1))
}

# ---------------- level-2 refit (Part B) ----------------
fit_stacks <- function(J, rows) {
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
  out <- list()
  Z <- sapply(J$mods, function(m) wz(J$lp[[m]], rows))
  D <- data.frame(tt = J$t, ee = J$e, Z)
  f <- tryCatch(coxph(Surv(tt, ee) ~ ., data = D[rows, ], ties = "efron"),
                error = function(e) NULL)
  if (!is.null(f)) {
    bh <- basehaz(f, centered = FALSE)
    lp <- as.vector(as.matrix(D[, J$mods, drop = FALSE]) %*%
                      coef(f)[J$mods])
    for (h in J$hs) {
      i <- findInterval(h + 0.5, bh$time)
      H0 <- if (i == 0) 0 else bh$hazard[i]
      out[[sprintf("coxstack|%d", h)]] <- 1 - exp(-H0 * exp(lp))
    }
  }
  for (h in J$hs) {
    L <- as.data.frame(sapply(J$mods, function(m)
      wz(qlogis(clip(J$lr[[sprintf("%s|%d", m, h)]])), rows)))
    y <- J$y[[as.character(h)]]
    g <- tryCatch(suppressWarnings(glm(
      y ~ ., family = binomial,
      data = cbind(y = y[rows], L[rows, , drop = FALSE]))),
      error = function(e) NULL)
    if (!is.null(g)) {
      out[[sprintf("binstack|%d", h)]] <-
        as.vector(predict(g, newdata = L, type = "response"))
    }
  }
  out
}
boot_chunk <- function(reps, J, keys, p_orig, th_cls, th_nb, sub_idx) {
  nbf <- function(y, p, th) {
    n <- length(y)
    vapply(th, function(t) {
      tr <- p >= t
      (sum(tr & y == 1) - sum(tr & y == 0) * t / (1 - t)) / n
    }, numeric(1))
  }
  acc <- list()
  for (k in keys) {
    acc[[k]] <- list(sad = numeric(J$n),
                     fl = matrix(0, J$n, length(th_cls)),
                     ici = numeric(0),
                     nb = matrix(numeric(0), 0, length(th_nb)),
                     sub = matrix(numeric(0), 0, length(sub_idx)),
                     nok = 0L)
  }
  for (r in reps) {
    if (J$fr == "A") {
      rows <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                          length(J$rows_by),
                                          replace = TRUE)],
                     use.names = FALSE)
    } else {
      rows <- sample.int(J$n, J$n, replace = TRUE)
    }
    pr <- fit_stacks(J, rows)
    for (k in keys) {
      pb <- pr[[k]]
      if (is.null(pb) || anyNA(pb)) next
      po <- p_orig[[k]]
      y <- J$y[[strsplit(k, "|", fixed = TRUE)[[1]][2]]]
      a <- acc[[k]]
      a$sad <- a$sad + abs(pb - po)
      for (j in seq_along(th_cls)) {
        a$fl[, j] <- a$fl[, j] + ((pb >= th_cls[j]) != (po >= th_cls[j]))
      }
      fit_l <- tryCatch(fitted(loess(y ~ pb, span = 0.75, degree = 2)),
                        error = function(e) rep(NA_real_, length(pb)))
      a$ici <- c(a$ici, mean(abs(fit_l - pb)))
      a$nb <- rbind(a$nb, nbf(y, pb, th_nb))
      a$sub <- rbind(a$sub, pb[sub_idx])
      a$nok <- a$nok + 1L
      acc[[k]] <- a
    }
  }
  acc
}
combine <- function(parts, keys) {
  out <- list()
  for (k in keys) {
    out[[k]] <- list(
      sad = Reduce(`+`, lapply(parts, function(p) p[[k]]$sad)),
      fl = Reduce(`+`, lapply(parts, function(p) p[[k]]$fl)),
      ici = unlist(lapply(parts, function(p) p[[k]]$ici)),
      nb = do.call(rbind, lapply(parts, function(p) p[[k]]$nb)),
      sub = do.call(rbind, lapply(parts, function(p) p[[k]]$sub)),
      nok = sum(sapply(parts, function(p) p[[k]]$nok)))
  }
  out
}

# ---------------- subgroups (Part C) ----------------
sub_vec <- function(S, idx) {
  met <- function(y, p) {
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    n1 <- sum(y)
    n0 <- length(y) - n1
    auc <- NA_real_
    if (n1 > 0 && n0 > 0) {
      r <- data.table::frank(p)
      auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }
    sl <- tryCatch(suppressWarnings(unname(glm.fit(
      cbind(1, qlogis(p)), y, family = binomial())$coefficients[2])),
      error = function(e) NA_real_)
    f <- tryCatch(fitted(loess(y ~ p, span = 0.75, degree = 2)),
                  error = function(e) rep(NA_real_, length(p)))
    c(oe = mean(y) / mean(p), auc = auc, slope = sl,
      ici = mean(abs(f - p)))
  }
  v <- numeric(0)
  y <- S$y[idx]
  for (m in names(S$P)) {
    p <- S$P[[m]][idx]
    for (gv in names(S$G)) {
      g <- S$G[[gv]][idx]
      for (lv in S$levels[[gv]]) {
        ii <- which(g == lv)
        r <- if (length(ii) < 20 || sum(y[ii]) < 5) {
          c(oe = NA, auc = NA, slope = NA, ici = NA)
        } else {
          met(y[ii], p[ii])
        }
        v[paste(m, gv, lv, names(r), sep = "|")] <- r
      }
    }
  }
  v
}
sub_boot <- function(i, S) {
  if (S$fr == "A") {
    idx <- unlist(S$rows_by[sample.int(length(S$rows_by),
                                       length(S$rows_by), replace = TRUE)],
                  use.names = FALSE)
  } else {
    idx <- sample.int(S$n, S$n, replace = TRUE)
  }
  sub_vec(S, idx)
}

# ---------------- cells ----------------
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
CELLS <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  Fd <- if (fr == "A") fa else fb
  CELLS[[paste(fr, co)]] <- list(
    fr = fr, co = co, code = if (co == "EHR+ECG") "ehrecg" else "3mod",
    sub = if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1],
    hs = HZ[[fr]], mods = c("ehr", "ecg", if (co == "3-mod") "ctpa"))
}
make_J <- function(C) {
  lpc <- sprintf("%s.cox.death.lp", C$mods)
  lrc <- as.vector(outer(C$mods, C$hs, function(m, h)
    sprintf("%s.lr.death_h%d", m, h)))
  lrk <- as.vector(outer(C$mods, C$hs, function(m, h)
    sprintf("%s|%d", m, h)))
  x <- avg(file.path(PR, "level1", sprintf("oof_linear_%s_%s.csv",
                                           C$fr, C$code)), c(lpc, lrc))
  x <- x[match(C$sub$hadm_id, hadm_id)]
  list(fr = C$fr, co = C$co, n = nrow(C$sub), hs = C$hs, mods = C$mods,
       t = C$sub$time_death, e = as.integer(C$sub$event_death == 1),
       lp = setNames(lapply(lpc, function(cc) x[[cc]]), C$mods),
       lr = setNames(lapply(lrc, function(cc) x[[cc]]), lrk),
       y = setNames(lapply(C$hs, function(h) as.integer(
         C$sub$event_death == 1 & C$sub$time_death <= h + 0.5)),
         as.character(C$hs)),
       rows_by = split(seq_len(nrow(C$sub)),
                       match(C$sub$subject_id,
                             unique(C$sub$subject_id))))
}
JS <- lapply(CELLS, make_J)

# subgroup data
ehrm <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
nm <- names(ehrm)
agec <- nm[tolower(nm) %in% c("age", "age_at_admit", "anchor_age")][1]
sexc <- nm[grepl(paste0("^(sex|gender|male|is_male|sex_male|gender_m|",
                        "female|is_female|gender_f)$"), tolower(nm))][1]
cat(sprintf("age column: %s; sex column: %s\n",
            ifelse(is.na(agec), "NOT FOUND", agec),
            ifelse(is.na(sexc), "NOT FOUND (sex skipped)", sexc)))
grp_sex <- function(v, nmv) {
  if (is.character(v) || is.factor(v)) {
    u <- toupper(substr(as.character(v), 1, 1))
    return(ifelse(u == "M", "male", ifelse(u == "F", "female",
                                           "sex missing")))
  }
  lv <- tolower(nmv)
  if (lv %in% c("male", "is_male", "sex_male", "gender_m")) {
    return(ifelse(is.na(v), "sex missing",
                  ifelse(v == 1, "male", "female")))
  }
  if (lv %in% c("female", "is_female", "gender_f")) {
    return(ifelse(is.na(v), "sex missing",
                  ifelse(v == 1, "female", "male")))
  }
  ifelse(is.na(v), "sex missing", paste0(nmv, " = ", v))
}
make_S <- function(fr, h) {
  C <- CELLS[[paste(fr, "EHR+ECG")]]
  x <- avg(file.path(PR, "calibration", sprintf("cal_preds_%s_ehrecg.csv",
                                                fr)),
           sprintf("death|%s|%d", names(MLAB), h))
  x <- x[match(C$sub$hadm_id, hadm_id)]
  e <- ehrm[match(C$sub$hadm_id, hadm_id)]
  G <- list()
  if (!is.na(agec)) {
    a <- e[[agec]]
    G[["age band"]] <- ifelse(is.na(a), "age missing",
                              ifelse(a < 65, "<65",
                                     ifelse(a < 80, "65-79", "80+")))
  }
  if (!is.na(sexc)) G[["sex"]] <- grp_sex(e[[sexc]], sexc)
  list(fr = fr, h = h, n = nrow(C$sub),
       y = as.integer(C$sub$event_death == 1 &
                        C$sub$time_death <= h + 0.5),
       P = setNames(lapply(names(MLAB), function(m)
         x[[sprintf("death|%s|%d", m, h)]]), MLAB),
       G = G, levels = lapply(G, function(g) sort(unique(g))),
       rows_by = split(seq_len(nrow(C$sub)),
                       match(C$sub$subject_id, unique(C$sub$subject_id))))
}
SS <- list(make_S("A", 30), make_S("B", 365))

# ---------------- timing probe ----------------
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
J1 <- JS[["A EHR+ECG"]]
po1 <- fit_stacks(J1, seq_len(J1$n))
t0 <- Sys.time()
invisible(boot_chunk(1:2, J1, names(po1), po1, TH_CLS, TH_NB, 1:10))
t_b <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 2
t0 <- Sys.time()
invisible(sub_boot(1, SS[[1]]))
t_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
proj <- (t_b * NBOOT * length(JS) + t_s * NBOOT_SUB * length(SS)) /
  workers * 1.3
cat(sprintf("  Part B %.2f s per bootstrap model set; Part C %.2f s per",
            t_b, t_s),
    sprintf("sample; projected total about %.1f min\n", proj / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)

# ---------------- Part A: seeds ----------------
cat("\n=== PART A: seed instability ===\n")
SA <- list()
for (C in CELLS) {
  pc <- file.path(PR, "calibration",
                  sprintf("cal_preds_%s_%s.csv", C$fr, C$code))
  cols <- as.vector(outer(names(MLAB), C$hs, function(m, h)
    sprintf("death|%s|%d", m, h)))
  x <- fread(pc, select = c("seed", "hadm_id", cols))
  seeds <- sort(unique(x$seed))
  for (cc in cols) {
    M <- sapply(seeds, function(s) {
      xs <- x[x$seed == s]
      xs[[cc]][match(C$sub$hadm_id, xs$hadm_id)]
    })
    pbar <- rowMeans(M)
    dev <- rowMeans(abs(M - pbar))
    pz <- strsplit(cc, "|", fixed = TRUE)[[1]]
    flip <- function(ii) sapply(TH_CLS, function(t) {
      r <- rowSums(M[ii, , drop = FALSE] >= t)
      mean(r > 0 & r < ncol(M))
    })
    dec <- ceiling(10 * frank(pbar, ties.method = "first") / length(pbar))
    groups <- c(list(all = seq_along(pbar)),
                setNames(lapply(1:10, function(g) which(dec == g)),
                         paste("risk decile", 1:10)))
    for (gn in names(groups)) {
      ii <- groups[[gn]]
      fl <- flip(ii)
      SA[[length(SA) + 1]] <- data.table(
        frame = C$fr, cohort = C$co, horizon = as.integer(pz[3]),
        model = MLAB[[pz[2]]], seeds = length(seeds), group = gn,
        n = length(ii), mean_risk = mean(pbar[ii]),
        mean_abs_dev = mean(dev[ii]),
        p95_abs_dev = quantile(dev[ii], 0.95, names = FALSE),
        max_abs_dev = max(dev[ii]), flip_at_10pct = fl[1],
        flip_at_20pct = fl[2])
    }
  }
}
sa <- rbindlist(SA)
fwrite(sa, file.path(RS, "stab_seed.csv"))

# ---------------- Part B: bootstrap ----------------
SB <- list()
BIN <- list()
for (J in JS) {
  cat(sprintf("=== PART B: Frame %s %s, %d bootstrap samples ===\n",
              J$fr, J$co, NBOOT))
  po <- fit_stacks(J, seq_len(J$n))
  keys <- names(po)
  set.seed(SEED)
  sub_idx <- sort(sample.int(J$n, min(NSUB, J$n)))
  chunks <- split(seq_len(NBOOT), ceiling(seq_len(NBOOT) / 20))
  parts <- future_lapply(chunks, boot_chunk, J = J, keys = keys,
                         p_orig = po, th_cls = TH_CLS, th_nb = TH_NB,
                         sub_idx = sub_idx, future.seed = SEED,
                         future.packages = PKG)
  A <- combine(parts, keys)
  for (k in keys) {
    pz <- strsplit(k, "|", fixed = TRUE)[[1]]
    y <- J$y[[pz[2]]]
    a <- A[[k]]
    mape <- a$sad / a$nok
    inst <- a$fl / a$nok
    nbo <- nb1(y, po[[k]], TH_NB)
    row <- data.table(
      frame = J$fr, cohort = J$co, horizon = as.integer(pz[2]),
      model = MLAB[[pz[1]]], n = J$n, events = sum(y), boots_ok = a$nok,
      mean_risk = mean(po[[k]]), mape_mean = mean(mape),
      mape_p95 = quantile(mape, 0.95, names = FALSE),
      mape_max = max(mape),
      inst10_mean = mean(inst[, 1]), inst10_share_over_20pct =
        mean(inst[, 1] > 0.2),
      inst20_mean = mean(inst[, 2]), inst20_share_over_20pct =
        mean(inst[, 2] > 0.2),
      ici_original = ici1(y, po[[k]]),
      ici_boot_median = median(a$ici, na.rm = TRUE),
      ici_boot_lo = quantile(a$ici, 0.025, na.rm = TRUE, names = FALSE),
      ici_boot_hi = quantile(a$ici, 0.975, na.rm = TRUE, names = FALSE))
    for (j in seq_along(TH_NB)) {
      tag <- sprintf("nb%02d", round(100 * TH_NB[j]))
      row[[paste0(tag, "_original")]] <- nbo[j]
      row[[paste0(tag, "_boot_lo")]] <-
        quantile(a$nb[, j], 0.025, names = FALSE)
      row[[paste0(tag, "_boot_hi")]] <-
        quantile(a$nb[, j], 0.975, names = FALSE)
    }
    SB[[length(SB) + 1]] <- row
    ps <- po[[k]][sub_idx]
    bin <- ceiling(20 * frank(ps, ties.method = "first") / length(ps))
    for (b in 1:20) {
      vals <- as.vector(a$sub[, bin == b, drop = FALSE])
      q <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
      BIN[[length(BIN) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, horizon = as.integer(pz[2]),
        model = MLAB[[pz[1]]], part = "instability plot (800 patients)",
        bin = b, mean_original = mean(ps[bin == b]), q025 = q[1],
        q25 = q[2], q50 = q[3], q75 = q[4], q975 = q[5],
        mape = NA_real_)
    }
    dec <- ceiling(10 * frank(po[[k]], ties.method = "first") / J$n)
    for (g in 1:10) {
      BIN[[length(BIN) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, horizon = as.integer(pz[2]),
        model = MLAB[[pz[1]]], part = "MAPE by risk decile (all patients)",
        bin = g, mean_original = mean(po[[k]][dec == g]), q025 = NA_real_,
        q25 = NA_real_, q50 = NA_real_, q75 = NA_real_, q975 = NA_real_,
        mape = mean(mape[dec == g]))
    }
    m <- a$sub[seq_len(min(NKEEP, nrow(a$sub))), , drop = FALSE]
    fwrite(data.table(
      model = k, boot = rep(seq_len(nrow(m)), times = ncol(m)),
      patient = rep(seq_len(ncol(m)), each = nrow(m)),
      p_original = rep(ps, each = nrow(m)), p_boot = as.vector(m)),
      file.path(SDIR, sprintf("stab_points_%s_%s_%s.csv", J$fr,
                              if (J$co == "EHR+ECG") "ehrecg" else "3mod",
                              gsub("|", "_", k, fixed = TRUE))))
  }
}
sb <- rbindlist(SB, fill = TRUE)
bn <- rbindlist(BIN)
fwrite(sb, file.path(RS, "stab_summary.csv"))
fwrite(bn, file.path(RS, "stab_bins.csv"))

# ---------------- Part C: subgroups ----------------
SC <- list()
for (S in SS) {
  cat(sprintf("=== PART C: Frame %s death %d d, %d bootstrap samples ===\n",
              S$fr, S$h, NBOOT_SUB))
  if (!length(S$G)) next
  B <- do.call(rbind, future_lapply(seq_len(NBOOT_SUB), sub_boot, S = S,
                                    future.seed = SEED,
                                    future.packages = PKG))
  pt <- sub_vec(S, seq_len(S$n))
  for (m in names(S$P)) for (gv in names(S$G)) for (lv in S$levels[[gv]]) {
    ii <- which(S$G[[gv]] == lv)
    row <- data.table(
      frame = S$fr, cohort = "EHR+ECG", horizon = S$h, model = m,
      variable = gv, level = lv, n = length(ii), events = sum(S$y[ii]),
      observed_risk = mean(S$y[ii]), mean_predicted = mean(S$P[[m]][ii]))
    for (ms in c("oe", "auc", "slope", "ici")) {
      k <- paste(m, gv, lv, ms, sep = "|")
      q <- quantile(B[, k], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      row[[ms]] <- pt[[k]]
      row[[paste0(ms, "_lo")]] <- q[1]
      row[[paste0(ms, "_hi")]] <- q[2]
    }
    SC[[length(SC) + 1]] <- row
  }
}
plan(sequential)
sc <- rbindlist(SC, fill = TRUE)
fwrite(sc, file.path(RS, "subgroups.csv"))

cat("\n=== BOOTSTRAP INSTABILITY (fusion stage) ===\n")
print(sb[, .(frame, cohort, horizon, model, mean_risk = round(mean_risk, 3),
             mape = round(mape_mean, 4), mape_p95 = round(mape_p95, 4),
             inst10 = round(inst10_mean, 3), inst20 = round(inst20_mean, 3),
             ici = round(ici_original, 4))], row.names = FALSE)
cat("\n=== SUBGROUPS ===\n")
if (nrow(sc)) {
  print(sc[, .(frame, model, variable, level, n, events,
               oe = round(oe, 2), auc = round(auc, 3),
               slope = round(slope, 2), ici = round(ici, 3))],
        row.names = FALSE)
}
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("stab_seed.csv", "stab_summary.csv", "stab_bins.csv",
            "subgroups.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", SDIR, " (patient-level points, local only)\n")
cat("  ", file.path(RS, "stab_log.txt"), " (terminal log)\n")
cat("\nOutputs: stab_seed.csv, stab_summary.csv, stab_bins.csv,",
    "subgroups.csv\n")