# 15c_ecg_features.R
# Which ECG features carry the early signal (PREREGISTRATION_ADDENDUM.md,
# A4). Exploratory; intervals only.
#
# Cohort: Frame B EHR+ECG, death within 365 days of discharge. Each of the
# 264 ECG features (71 SCP-statement logits and 193 measurements) is
# winsorised at its 0.5th and 99.5th percentiles and standardised, using
# complete cases per feature. The EHR score is the seed-averaged
# out-of-fold ridge-Cox death score, winsorised and standardised.
#
# Main model: Cox model with the feature's coefficient split at day 30 (cut
# at 30.5, the primary analysis's periods), adjusted for the EHR score,
# which is also split. Reported are the hazard ratio per SD early and late,
# their ratio with Wald 95% CI and p, and Benjamini-Hochberg q across
# features, together with the whole-year hazard ratio (constant effect,
# EHR-adjusted) with p and q. Each feature is oriented so that higher
# values mean higher one-year risk (the sign of its whole-year
# coefficient), so a ratio above 1 always means front-loaded.
#
# Sensitivity: the effect linear in log time (the best shape in A1.2),
# adjusted for the EHR score in the same way (slope per unit log(day),
# negative = front-loaded, with p and q), and the main model without EHR
# adjustment. Summaries: front- and late-loaded features (q < 0.05) by
# block (SCP statements, measurements and the right-ventricular strain
# statements named in A4), and the 20 most front-loaded features.
#
# Notes: the feature matrix is passed to the workers as Xm, not X, because
# X is the name of future_lapply's own first argument.
#
# Inputs:  data/processed/frame_B.csv, data/processed/feature_whitelist.csv,
#          data/processed/level1/oof_linear_B_ehrecg.csv,
#          data/inputs/ecg_opt.csv
# Outputs: results/ecg_features_tv.csv, results/ecg_features_blocks.csv,
#          results/ecg_features_log.txt
# Usage:   Rscript scripts/15c_ecg_features.R [--probe-only]

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
CUT <- 30.5
RV <- c("STACH", "CRBBB", "IRBBB", "RVH", "RAO", "INVT", "NT_", "TAB_")
PKG <- c("survival")
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()

wins_na <- function(v) {
  q <- quantile(v, c(0.005, 0.995), na.rm = TRUE, names = FALSE)
  w <- pmin(pmax(v, q[1]), q[2])
  (w - mean(w, na.rm = TRUE)) / sd(w, na.rm = TRUE)
}

feat_chunk <- function(js, d, Xm, cut) {
  d$rid <- seq_len(nrow(d))
  s2 <- survSplit(Surv(time_death, event_death) ~ rid + ehr, data = d,
                  cut = cut, episode = "per")
  s2$ehr_e <- s2$ehr * (s2$per == 1)
  s2$ehr_l <- s2$ehr * (s2$per == 2)
  ev <- sort(unique(d$time_death[d$event_death == 1]))
  sl <- survSplit(Surv(time_death, event_death) ~ rid + ehr, data = d,
                  cut = ev, episode = "iv")
  sl$lt <- log(sl$time_death)
  sl$ehr_lt <- sl$ehr * sl$lt
  ratio_of <- function(f) {
    b <- coef(f)
    V <- vcov(f)
    c(b[["x_e"]] - b[["x_l"]],
      sqrt(V["x_e", "x_e"] + V["x_l", "x_l"] - 2 * V["x_e", "x_l"]))
  }
  out <- list()
  for (j in js) {
    x <- Xm[, j]
    ok <- !is.na(x)
    base <- data.frame(feature = colnames(Xm)[j], n = sum(ok),
                       deaths = sum(d$event_death[ok]))
    r <- tryCatch({
      f0 <- coxph(Surv(time_death, event_death) ~ x + ehr,
                  data = data.frame(d, x = x)[ok, ], ties = "efron")
      s <- s2[ok[s2$rid], ]
      xs <- x[s$rid]
      s$x_e <- xs * (s$per == 1)
      s$x_l <- xs * (s$per == 2)
      fa <- coxph(Surv(tstart, time_death, event_death) ~
                    x_e + x_l + ehr_e + ehr_l, data = s, ties = "efron")
      fu <- coxph(Surv(tstart, time_death, event_death) ~ x_e + x_l,
                  data = s, ties = "efron")
      l <- sl[ok[sl$rid], ]
      l$x <- x[l$rid]
      l$x_lt <- l$x * l$lt
      fl <- coxph(Surv(tstart, time_death, event_death) ~
                    x + x_lt + ehr + ehr_lt, data = l, ties = "efron")
      b0 <- coef(f0)[["x"]]
      se0 <- sqrt(vcov(f0)["x", "x"])
      o <- if (b0 < 0) -1 else 1
      ra <- ratio_of(fa)
      ru <- ratio_of(fu)
      bl <- coef(fl)[["x_lt"]]
      sel <- sqrt(vcov(fl)["x_lt", "x_lt"])
      data.frame(
        orientation = o,
        hr_year = exp(o * b0), hr_year_lo = exp(o * b0 - 1.96 * se0),
        hr_year_hi = exp(o * b0 + 1.96 * se0),
        p_year = 2 * pnorm(-abs(b0 / se0)),
        hr_early = exp(o * coef(fa)[["x_e"]]),
        hr_late = exp(o * coef(fa)[["x_l"]]),
        ratio = exp(o * ra[1]), ratio_lo = exp(o * ra[1] - 1.96 * ra[2]),
        ratio_hi = exp(o * ra[1] + 1.96 * ra[2]),
        p_ratio = 2 * pnorm(-abs(ra[1] / ra[2])),
        ratio_unadj = exp(o * ru[1]),
        p_ratio_unadj = 2 * pnorm(-abs(ru[1] / ru[2])),
        slope_logt = o * bl, slope_lo = o * bl - 1.96 * sel,
        slope_hi = o * bl + 1.96 * sel,
        p_slope = 2 * pnorm(-abs(bl / sel)))
    }, error = function(e) NULL)
    if (!is.null(r)) out[[length(out) + 1]] <- cbind(base, r)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

cat("=== DATA ===\n")
fb <- fread(file.path(PR, "frame_B.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
fcols <- wl[modality == "ecg", column]
ecg <- unique(fread(file.path(IN, "ecg_opt.csv"),
                    select = c("hadm_id", fcols)), by = "hadm_id")
sc <- fread(file.path(PR, "level1", "oof_linear_B_ehrecg.csv"),
            select = c("hadm_id", "ehr.cox.death.lp"))
sc <- sc[, .(ehr = mean(`ehr.cox.death.lp`)), by = hadm_id]
dd <- merge(fb[, .(hadm_id, time_death, event_death)], sc, by = "hadm_id")
dd[, ehr := wins_na(ehr)]
Xall <- as.matrix(ecg[match(dd$hadm_id, hadm_id), ..fcols])
storage.mode(Xall) <- "double"
Xall <- apply(Xall, 2, wins_na)
colnames(Xall) <- fcols
d <- as.data.frame(dd[, .(time_death, event_death, ehr)])
cat(sprintf("  n=%d deaths=%d features=%d (SCP %d, measurements %d)\n",
            nrow(d), sum(d$event_death), ncol(Xall),
            sum(startsWith(fcols, "scp_")),
            sum(!startsWith(fcols, "scp_"))))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
t0 <- Sys.time()
invisible(feat_chunk(1:2, d, Xall, CUT))
t2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
proj <- t2 / 2 * ncol(Xall) / workers * 1.4
cat(sprintf("  two features (with setup) %.1f s; projected about %.1f min\n",
            t2, proj / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

plan(multisession, workers = workers)
chunks <- split(seq_len(ncol(Xall)),
                cut(seq_len(ncol(Xall)), workers * 3, labels = FALSE))
res <- future_lapply(chunks, feat_chunk, d = d, Xm = Xall, cut = CUT,
                     future.packages = PKG, future.seed = TRUE)
plan(sequential)
tv <- rbindlist(res)
tv[, block := ifelse(startsWith(feature, "scp_"), "SCP statement",
                     "measurement")]
short <- sub("^scp_", "", tv$feature)
tv[, rv_strain := block == "SCP statement" &
     Reduce(`|`, lapply(RV, function(p) grepl(p, short, fixed = TRUE)))]
tv[, q_year := p.adjust(p_year, "BH")]
tv[, q_ratio := p.adjust(p_ratio, "BH")]
tv[, q_ratio_unadj := p.adjust(p_ratio_unadj, "BH")]
tv[, q_slope := p.adjust(p_slope, "BH")]
setorder(tv, p_ratio)
fwrite(tv, file.path(RS, "ecg_features_tv.csv"))

BL <- list()
blocks <- list("SCP statement" = tv$block == "SCP statement",
               "measurement" = tv$block == "measurement",
               "right-ventricular-strain statements" = tv$rv_strain,
               "all features" = rep(TRUE, nrow(tv)))
for (bn in names(blocks)) {
  x <- tv[blocks[[bn]]]
  BL[[length(BL) + 1]] <- data.table(
    block = bn, features = nrow(x),
    whole_year_q05 = sum(x$q_year < 0.05),
    front_loaded_q05 = sum(x$q_ratio < 0.05 & x$ratio > 1),
    late_loaded_q05 = sum(x$q_ratio < 0.05 & x$ratio < 1),
    front_loaded_logt_q05 = sum(x$q_slope < 0.05 & x$slope_logt < 0),
    median_ratio = median(x$ratio),
    median_ratio_if_year_q05 = if (any(x$q_year < 0.05))
      median(x$ratio[x$q_year < 0.05]) else NA_real_)
}
bl <- rbindlist(BL)
fwrite(bl, file.path(RS, "ecg_features_blocks.csv"))

cat("\n=== BLOCKS ===\n")
print(bl, row.names = FALSE)
cat("\n=== 20 MOST FRONT-LOADED (ratio > 1, smallest p) ===\n")
print(head(tv[ratio > 1, .(feature, rv_strain, n, hr_year = round(hr_year, 2),
                           hr_early = round(hr_early, 2),
                           hr_late = round(hr_late, 2),
                           ratio = round(ratio, 2), lo = round(ratio_lo, 2),
                           hi = round(ratio_hi, 2), q = signif(q_ratio, 2),
                           q_logt = signif(q_slope, 2))], 20),
      row.names = FALSE)
cat(sprintf("\n  features analysed %d of %d; elapsed %.1f min\n", nrow(tv),
            ncol(Xall),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(RS, "ecg_features_tv.csv"), "\n")
cat("  ", file.path(RS, "ecg_features_blocks.csv"), "\n")
cat("  ", file.path(RS, "ecg_features_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: ecg_features_tv.csv, ecg_features_blocks.csv\n")