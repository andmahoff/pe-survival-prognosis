# 17a_gap_time.R
# Fusion gap mechanism along follow-up time (PREREGISTRATION_ADDENDUM_2.md,
# A9). Exploratory; intervals only.
#
# Examines how the gain from fusing two modalities depends on the gap
# between their single-modality AUROCs, over follow-up time. It follows the
# side project (pe-fusion-transfer), where gain fell with the gap across 12
# transfer and outcome cells (r = -0.902), and script 140 of the
# dissertation (pe-multimodal-prognosis-msc), where gain halved as the EHR
# modality strengthened. Theory: Su and Liu (1993), Pencina et al. (2012).
#
# Frame B (death after discharge) in the EHR+ECG cohort and the
# three-modality cohort, using the seed-averaged out-of-fold level-1 Cox
# scores and folds from 06a_level1_linear.py. A cell is a modality pair at
# one horizon: EHR+ECG (both cohorts), EHR+CTPA and ECG+CTPA (three-modality
# cohort), at 10, 14, 21, 30, 45, 60, 90, 120, 180, 270 and 365 days. In the
# training folds of every seed, the main model is a logistic stack per
# horizon on the scores (winsorised at the 0.5th and 99.5th percentiles,
# standardised) of each single modality and each pair. The sensitivity model
# is one Cox stack per score set, with Breslow risks at every horizon.
#
# The gap is the AUROC of the better single modality minus the other. The
# fusion gain is the pair minus the better single modality, measured by
# AUROC, IPA, a log-likelihood skill score (1 - LL / LL of the event-rate
# model) and the continuous NRI. The theoretical gain of the optimal linear
# combination under the binormal model is computed from the two AUROCs and
# the within-class correlation of the normal scores. Slope and correlation
# of gain on gap are estimated along time within each pair and pooled over
# all cells, with observed against theoretical gain and the same regression
# for the other gain measures. The incident/dynamic AUC (Heagerty and Zheng
# 2005) of the Cox-stack score is kernel-smoothed on log time (bandwidth
# 0.4) at the same 11 times and treated likewise.
#
# Intervals are 95% percentile intervals from 2,000 patient-level bootstrap
# resamples (seed 42), the three-modality cohort being a subset of each
# resample. The script stops if a check fails, runs in parallel with
# future.apply and begins with a timing probe that pauses 5 s.
#
# Inputs:  data/processed/frame_B.csv, level1/oof_linear_B_*.csv
# Outputs: results/gap_cells.csv, results/gap_regressions.csv,
#          results/gap_incident.csv, results/gap_checks.csv,
#          results/gap_log.txt
# Usage:   Rscript scripts/17a_gap_time.R [--probe-only]

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
NBOOT <- 2000L
SEED <- 42L
HZ <- c(10, 14, 21, 30, 45, 60, 90, 120, 180, 270, 365)
BW <- 0.4
MEAS <- c("auc", "ipa", "ll", "nri")
PKG <- c("survival", "data.table")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-60s %s\n", status, item, value))
}

# ---------------- stacking: one task = cohort x seed ----------------
stack_task <- function(tk) {
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
  n <- length(tk$t)
  lg <- list()
  cx <- list()
  lp <- list()
  for (sn in names(tk$sets)) {
    lp[[sn]] <- rep(NA_real_, n)
    for (h in tk$hz) {
      lg[[paste(sn, h, sep = "|")]] <- rep(NA_real_, n)
      cx[[paste(sn, h, sep = "|")]] <- rep(NA_real_, n)
    }
  }
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    for (sn in names(tk$sets)) {
      ms <- tk$sets[[sn]]
      Z <- matrix(sapply(ms, function(m) wz(tk$LP[[m]], tr)),
                  ncol = length(ms))
      X1 <- cbind(1, Z)
      for (h in tk$hz) {
        y <- as.numeric(tk$e == 1 & tk$t <= h + 0.5)
        f <- suppressWarnings(glm.fit(X1[tr, , drop = FALSE], y[tr],
                                      family = binomial()))
        b <- f$coefficients
        b[is.na(b)] <- 0
        lg[[paste(sn, h, sep = "|")]][te] <-
          plogis(as.vector(X1[te, , drop = FALSE] %*% b))
      }
      D <- data.frame(tt = tk$t, ee = tk$e, Z)
      zc <- setdiff(names(D), c("tt", "ee"))
      fc <- coxph(Surv(tt, ee) ~ ., data = D[tr, ], ties = "efron")
      bh <- basehaz(fc, centered = FALSE)
      l <- as.vector(as.matrix(D[te, zc, drop = FALSE]) %*% coef(fc)[zc])
      lp[[sn]][te] <- l
      for (h in tk$hz) {
        i <- findInterval(h + 0.5, bh$time)
        H0 <- if (i == 0) 0 else bh$hazard[i]
        cx[[paste(sn, h, sep = "|")]][te] <- 1 - exp(-H0 * exp(l))
      }
    }
  }
  list(co = tk$co, seed = tk$seed, lg = lg, cx = cx, lp = lp)
}

# ---------------- statistics on one resample ----------------
stat_all <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  ipa <- function(y, p) {
    pr <- mean(y)
    1 - mean((y - p)^2) / (pr * (1 - pr))
  }
  llsk <- function(y, p) {
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    pr <- mean(y)
    ll <- -mean(y * log(p) + (1 - y) * log(1 - p))
    l0 <- -(pr * log(pr) + (1 - pr) * log(1 - pr))
    1 - ll / l0
  }
  nri <- function(y, pn, po) {
    e <- y == 1
    (mean(pn[e] > po[e]) - mean(pn[e] < po[e])) +
      (mean(pn[!e] < po[!e]) - mean(pn[!e] > po[!e]))
  }
  wcor <- function(y, a, b) {
    na <- qnorm((data.table::frank(a) - 0.5) / length(a))
    nb <- qnorm((data.table::frank(b) - 0.5) / length(b))
    e <- y == 1
    ca <- c(na[e] - mean(na[e]), na[!e] - mean(na[!e]))
    cb <- c(nb[e] - mean(nb[e]), nb[!e] - mean(nb[!e]))
    sum(ca * cb) / sqrt(sum(ca^2) * sum(cb^2))
  }
  theory <- function(a1, a2, r) {
    d1 <- sqrt(2) * qnorm(a1)
    d2 <- sqrt(2) * qnorm(a2)
    dd <- (d1^2 + d2^2 - 2 * r * d1 * d2) / (1 - r^2)
    pnorm(sqrt(max(dd, 0)) / sqrt(2)) - max(a1, a2)
  }
  id_auc <- function(tt, ee, s, grid, bw) {
    dt <- sort(unique(tt[ee == 1]))
    ci <- numeric(0)
    ti <- numeric(0)
    for (d in dt) {
      sr <- sort(s[tt > d])
      if (!length(sr)) next
      sd_ <- s[tt == d & ee == 1]
      lo <- findInterval(sd_, sr, left.open = TRUE)
      le <- findInterval(sd_, sr)
      ci <- c(ci, (lo + 0.5 * (le - lo)) / length(sr))
      ti <- c(ti, rep(d, length(sd_)))
    }
    sapply(grid, function(g) {
      w <- dnorm((log(ti) - log(g)) / bw)
      sum(w * ci) / sum(w)
    })
  }
  slope_r <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0) return(c(NA_real_, NA_real_))
    c(unname(coef(lm(y[ok] ~ x[ok]))[2]), cor(x[ok], y[ok]))
  }
  v <- numeric(0)
  cells <- list()
  inc <- list()
  for (co in names(J$C)) {
    C <- J$C[[co]]
    ix <- if (co == "3-mod") idx[C$member[idx]] else idx
    tt <- J$t[ix]
    ee <- J$e[ix]
    for (sn in names(C$lp)) {
      a <- id_auc(tt, ee, C$lp[[sn]][ix], J$hz, J$bw)
      for (g in seq_along(J$hz)) {
        v[paste(co, "id", sn, J$hz[g], sep = "|")] <- a[g]
      }
      inc[[sn]] <- a
    }
    for (pp in C$pairs) {
      A <- pp[1]
      B <- pp[2]
      AB <- paste(A, B, sep = "+")
      pk <- paste(co, AB, sep = "|")
      ia <- inc[[A]]
      ib <- inc[[B]]
      bi <- ifelse(ia >= ib, ia, ib)
      igap <- abs(ia - ib)
      igain <- inc[[AB]] - bi
      for (g in seq_along(J$hz)) {
        v[paste(pk, "idgap", J$hz[g], sep = "|")] <- igap[g]
        v[paste(pk, "idgain", J$hz[g], sep = "|")] <- igain[g]
      }
      sr <- slope_r(igap, igain)
      v[paste(pk, "idreg", "slope", sep = "|")] <- sr[1]
      v[paste(pk, "idreg", "r", sep = "|")] <- sr[2]
      for (st in c("lg", "cx")) {
        P <- C[[st]]
        for (h in J$hz) {
          y <- as.numeric(ee == 1 & tt <= h + 0.5)
          pa <- P[[paste(A, h, sep = "|")]][ix]
          pb <- P[[paste(B, h, sep = "|")]][ix]
          pab <- P[[paste(AB, h, sep = "|")]][ix]
          aa <- auc(y, pa)
          ab <- auc(y, pb)
          best <- if (aa >= ab) pa else pb
          r <- wcor(y, C$l1[[A]][ix], C$l1[[B]][ix])
          k <- paste(pk, st, h, sep = "|")
          x <- c(auc_a = aa, auc_b = ab, auc_pair = auc(y, pab),
                 gap = abs(aa - ab), rho = r)
          x["gain_auc"] <- x[["auc_pair"]] - max(aa, ab)
          x["gain_ipa"] <- ipa(y, pab) - ipa(y, best)
          x["gain_ll"] <- llsk(y, pab) - llsk(y, best)
          x["gain_nri"] <- nri(y, pab, best)
          x["gain_theory"] <- theory(aa, ab, r)
          for (nm in names(x)) v[paste(k, nm, sep = "|")] <- x[[nm]]
          cells[[length(cells) + 1]] <- c(pair = pk, st = st, h = h, x)
        }
      }
    }
  }
  cm <- as.data.frame(do.call(rbind, cells), stringsAsFactors = FALSE)
  for (nm in setdiff(names(cm), c("pair", "st"))) {
    cm[[nm]] <- as.numeric(cm[[nm]])
  }
  for (st in c("lg", "cx")) {
    for (ms in c(J$meas, "theory")) {
      gcol <- if (ms == "theory") "gain_theory" else paste0("gain_", ms)
      for (pk in unique(cm$pair)) {
        s <- cm[cm$st == st & cm$pair == pk, ]
        sr <- slope_r(s$gap, s[[gcol]])
        v[paste(pk, st, "reg", ms, "slope", sep = "|")] <- sr[1]
        v[paste(pk, st, "reg", ms, "r", sep = "|")] <- sr[2]
      }
      s <- cm[cm$st == st, ]
      sr <- slope_r(s$gap, s[[gcol]])
      v[paste("pooled|all pairs", st, "reg", ms, "slope", sep = "|")] <- sr[1]
      v[paste("pooled|all pairs", st, "reg", ms, "r", sep = "|")] <- sr[2]
    }
    s <- cm[cm$st == st, ]
    sr <- slope_r(s$gain_theory, s$gain_auc)
    v[paste("pooled|all pairs", st, "obs_vs_theory", "slope", sep = "|")] <-
      sr[1]
    v[paste("pooled|all pairs", st, "obs_vs_theory", "r", sep = "|")] <-
      sr[2]
    v[paste("pooled|all pairs", st, "obs_vs_theory", "ratio", sep = "|")] <-
      sum(s$gain_auc) / sum(s$gain_theory)
  }
  v
}
boot_rep <- function(i, J) {
  stat_all(J, sample.int(J$n, J$n, replace = TRUE))
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fb <- fread(file.path(PR, "frame_B.csv"))
COH <- list(
  "EHR+ECG" = list(code = "ehrecg", mods = c("ehr", "ecg"),
                   pairs = list(c("EHR", "ECG"))),
  "3-mod" = list(code = "3mod", mods = c("ehr", "ecg", "ctpa"),
                 pairs = list(c("EHR", "ECG"), c("EHR", "CTPA"),
                              c("ECG", "CTPA"))))
tasks <- list()
L1 <- list()
for (co in names(COH)) {
  M <- COH[[co]]
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  lpc <- sprintf("%s.cox.death.lp", M$mods)
  oof <- fread(file.path(PR, "level1", sprintf("oof_linear_B_%s.csv",
                                               M$code)),
               select = c("seed", "hadm_id", "fold", lpc))
  seeds <- sort(unique(oof$seed))
  sets <- setNames(as.list(M$mods), toupper(M$mods))
  for (pp in M$pairs) {
    sets[[paste(pp, collapse = "+")]] <- tolower(pp)
  }
  L1[[co]] <- oof[, lapply(.SD, mean), by = hadm_id,
                  .SDcols = lpc][match(sub$hadm_id, hadm_id)]
  miss <- 0L
  for (s in seeds) {
    g <- oof[oof$seed == s][match(sub$hadm_id, hadm_id)]
    miss <- miss + sum(is.na(g$fold)) + sum(is.na(as.matrix(g[, ..lpc])))
    tasks[[length(tasks) + 1]] <- list(
      co = co, seed = s, fold = g$fold,
      LP = setNames(lapply(lpc, function(cc) g[[cc]]), M$mods),
      t = sub$time_death, e = as.integer(sub$event_death == 1),
      sets = sets, hz = HZ)
  }
  chk(sprintf("%s: missing level-1 scores or folds", co), miss,
      if (miss == 0) "PASS" else "FAIL")
  ev <- sapply(HZ, function(h) sum(sub$event_death == 1 &
                                     sub$time_death <= h + 0.5))
  chk(sprintf("%s: n / deaths at %s days", co,
              paste(range(HZ), collapse = " and ")),
      sprintf("%d / %d and %d", nrow(sub), ev[1], ev[length(ev)]),
      if (min(ev) >= 20) "PASS" else "WARN")
}
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "gap_checks.csv"))
  cat("\nSTOP: a check failed. see: gap_checks.csv\n")
  quit(status = 1)
}

cat(sprintf("\n=== STACKING (%d tasks, %d workers) ===\n", length(tasks),
            workers))
plan(multisession, workers = workers)
t0 <- Sys.time()
res <- future_lapply(tasks, stack_task, future.packages = PKG,
                     future.seed = TRUE)
cat(sprintf("  done in %.0f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

nb3 <- fb[, in_three_mod == 1]
J <- list(n = nrow(fb), t = fb$time_death,
          e = as.integer(fb$event_death == 1), hz = HZ, bw = BW,
          meas = MEAS, C = list())
nmax <- 0
rng <- c(Inf, -Inf)
for (co in names(COH)) {
  rr <- Filter(function(r) r$co == co, res)
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  pos <- match(sub$hadm_id, fb$hadm_id)
  full <- function(v) {
    out <- rep(NA_real_, nrow(fb))
    out[pos] <- v
    out
  }
  avgk <- function(part, key) {
    full(Reduce(`+`, lapply(rr, function(r) r[[part]][[key]])) / length(rr))
  }
  C <- list(member = if (co == "3-mod") nb3 else rep(TRUE, nrow(fb)),
            pairs = COH[[co]]$pairs, lg = list(), cx = list(), lp = list(),
            l1 = list())
  for (key in names(rr[[1]]$lg)) {
    C$lg[[key]] <- avgk("lg", key)
    C$cx[[key]] <- avgk("cx", key)
    nmax <- nmax + sum(is.na(C$lg[[key]][pos])) +
      sum(is.na(C$cx[[key]][pos]))
    rng <- c(min(rng[1], C$lg[[key]][pos], C$cx[[key]][pos]),
             max(rng[2], C$lg[[key]][pos], C$cx[[key]][pos]))
  }
  for (sn in names(rr[[1]]$lp)) C$lp[[sn]] <- avgk("lp", sn)
  for (m in COH[[co]]$mods) {
    C$l1[[toupper(m)]] <- full(L1[[co]][[sprintf("%s.cox.death.lp", m)]])
  }
  J$C[[co]] <- C
}
chk("missing stacked predictions", nmax, if (nmax == 0) "PASS" else "FAIL")
chk("range of stacked risks", sprintf("%.5f to %.5f", rng[1], rng[2]),
    if (rng[1] > 0 && rng[2] < 1) "PASS" else "FAIL")
if (any(rbindlist(CHK)$status == "FAIL")) {
  plan(sequential)
  fwrite(rbindlist(CHK), file.path(RS, "gap_checks.csv"))
  cat("\nSTOP: a check failed. see: gap_checks.csv\n")
  quit(status = 1)
}

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
nb <- workers * 2L
t0 <- Sys.time()
invisible(future_lapply(seq_len(nb), boot_rep, J = J, future.seed = 1L,
                        future.packages = PKG))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  %d replicates in %.1f s; %d: about %.1f min\n", nb, dt,
            NBOOT, dt * NBOOT / nb / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
pt <- stat_all(J, seq_len(J$n))
B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                  future.seed = SEED,
                                  future.packages = PKG))
plan(sequential)

# ---------------- assemble ----------------
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
row_ci <- function(k) {
  q <- q2(B[, k])
  data.table(estimate = pt[[k]], lo = q[1], hi = q[2])
}
STN <- c(lg = "logistic stack per horizon (main)",
         cx = "single Cox stack (sensitivity)")
CEL <- list()
REG <- list()
INC <- list()
for (k in names(pt)) {
  z <- strsplit(k, "|", fixed = TRUE)[[1]]
  is_cell <- length(z) == 5 && z[3] %in% c("lg", "cx") &&
    !is.na(suppressWarnings(as.numeric(z[4])))
  if (is_cell) {
    CEL[[length(CEL) + 1]] <- cbind(data.table(
      cohort = z[1], pair = z[2], stack = STN[[z[3]]],
      horizon = as.numeric(z[4]), measure = z[5]), row_ci(k))
  } else if (length(z) >= 5 && z[4] %in% c("reg", "obs_vs_theory")) {
    REG[[length(REG) + 1]] <- cbind(data.table(
      cohort = z[1], pair = z[2], stack = STN[[z[3]]],
      test = if (z[4] == "reg") sprintf("%s gain on AUROC gap", z[5]) else
        "observed on theoretical AUROC gain",
      statistic = z[length(z)]), row_ci(k))
  } else if (length(z) == 4 && z[3] == "idreg") {
    REG[[length(REG) + 1]] <- cbind(data.table(
      cohort = z[1], pair = z[2], stack = "Cox-stack score",
      test = "incident gain on incident gap", statistic = z[4]),
      row_ci(k))
  } else if (length(z) == 4 && z[2] == "id") {
    INC[[length(INC) + 1]] <- cbind(data.table(
      cohort = z[1], score = z[3], time = as.numeric(z[4]),
      measure = "incident/dynamic AUROC"), row_ci(k))
  } else if (length(z) == 4 && z[3] %in% c("idgap", "idgain")) {
    INC[[length(INC) + 1]] <- cbind(data.table(
      cohort = z[1], score = z[2], time = as.numeric(z[4]),
      measure = if (z[3] == "idgap") "incident gap" else
        "incident fusion gain"), row_ci(k))
  }
}
cel <- rbindlist(CEL)
reg <- rbindlist(REG)
inc <- rbindlist(INC)
reg[, excludes_zero := lo > 0 | hi < 0]
ea <- cel[grepl("^EHR", pair) & measure == "auc_a", estimate]
eb <- cel[grepl("^EHR", pair) & measure == "auc_b", estimate]
chk("EHR the better single modality in its pairs (share of cells)",
    sprintf("%.2f", mean(ea >= eb)), "INFO")
ck <- rbindlist(CHK)
fwrite(cel, file.path(RS, "gap_cells.csv"))
fwrite(reg, file.path(RS, "gap_regressions.csv"))
fwrite(inc, file.path(RS, "gap_incident.csv"))
fwrite(ck, file.path(RS, "gap_checks.csv"))

cat("\n=== GAP AND GAIN BY HORIZON, EHR+ECG cohort, main stack ===\n")
print(dcast(cel[cohort == "EHR+ECG" & grepl("logistic", stack) &
                  measure %in% c("gap", "gain_auc", "gain_theory",
                                 "gain_nri")],
            horizon ~ measure, value.var = "estimate")[, lapply(.SD, round,
                                                                4)],
      row.names = FALSE)
cat("\n=== REGRESSIONS (slopes and correlations) ===\n")
print(reg[, .(cohort, pair, stack = substr(stack, 1, 12), test,
              statistic, estimate = round(estimate, 3), lo = round(lo, 3),
              hi = round(hi, 3))], row.names = FALSE, nrows = 1000)
cat(sprintf("\n  elapsed %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("gap_cells.csv", "gap_regressions.csv", "gap_incident.csv",
            "gap_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "gap_log.txt"), " (terminal log)\n")
cat("\nOutputs: gap_cells.csv, gap_regressions.csv, gap_incident.csv,",
    "gap_checks.csv\n")