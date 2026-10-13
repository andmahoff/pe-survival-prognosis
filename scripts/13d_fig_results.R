# 13d_fig_results.R
# Core results figures (PREREGISTRATION_ADDENDUM_2.md A38).
#
# One panel per figure, with centred bold Title Case titles that state
# what the figure shows, direct labels instead of legends, Okabe-Ito
# colours with shapes or line types as a second channel, faint gridlines,
# Calibri where available, and no notes on the graphic (suggested captions
# go to a CSV). Titles that make a claim are worded from the numbers, so
# they stay true if the results change.
#
# Each figure is saved as a PNG file:
#   ecg_hr_early_late_B.png      the registered primary analysis (stage 7,
#                                Frame B, EHR+ECG cohort): hazard ratio per
#                                SD of the ECG and EHR scores for deaths in
#                                days 0-30 and 31-365 after discharge
#   ecg_ratio_sensitivity_B.png  the ECG's early / late ratio in the
#                                primary analysis and every sensitivity
#                                analysis (stages 7, 7b, 7c, 9a, 19a);
#                                missing files are skipped
#   cal_<P|B>_<outcome>_<cohort>.png
#                                calibration of each final model at 30 and
#                                365 days (stage 13b): loess curve with its
#                                bootstrap band and risk-group points; both
#                                cohorts (12 figures)
#   dca_<P|B>_<outcome>_<horizon>.png
#                                decision curves of each final model at 30
#                                and 365 days, EHR+ECG cohort (stage 13b):
#                                model, treat all, treat none and, at
#                                diagnosis for death, the sPESI and PESI-10
#                                rules (12 figures)
#
# Inputs:  results/primary.csv, primary_robust.csv, primary_hospice.csv,
#          s5_primary.csv, ab_ecg_ratio.csv, fin_metrics.csv,
#          fin_cal_bins.csv, fin_cal_curves.csv, fin_dca_curves.csv,
#          fin_dca_key.csv (aggregates only)
# Outputs: R Graphs/results_final/*.png (inside the project folder),
#          results/fin_fig_captions.csv, results/fin_fig_log.txt
#          (terminal log)
# Usage:   Rscript scripts/13d_fig_results.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
RS <- file.path(P2, "results")
FD <- file.path(P2, "R Graphs", "results_final")
dir.create(FD, showWarnings = FALSE, recursive = TRUE)
OI <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
        "#000000", "#999999")
fam <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  sf <- systemfonts::system_fonts()
  if (any(sf$family == "Calibri")) fam <- "Calibri"
}
dev_png <- function(path, w, h) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = w, height = h, units = "in", res = 300)
  } else {
    png(path, width = w, height = h, units = "in", res = 300)
  }
}
th <- theme_minimal(base_size = 12, base_family = fam) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92",
                                        linewidth = 0.3),
        legend.position = "none", plot.title.position = "plot",
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        plot.margin = margin(8, 14, 8, 8))
CAP <- list()
save_fig <- function(p, name, caption, w = 7, h = 4.5) {
  path <- file.path(FD, name)
  dev_png(path, w, h)
  print(p)
  invisible(dev.off())
  CAP[[length(CAP) + 1]] <<- data.table(file = name, caption = caption)
  cat("  saved", path, "\n")
}
# titles are written in Title Case here; this only wraps them
ttl <- function(s, w = 62) paste(strwrap(s, w), collapse = "\n")
hr_breaks <- function(lo, hi) {
  r <- range(c(lo, hi), na.rm = TRUE)
  ok <- function(b) b[b >= r[1] / 1.02 & b <= r[2] * 1.02]
  sets <- list(c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 7, 10, 20),
               c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4,
                 5),
               c(0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2))
  for (st in sets) if (length(ok(st)) >= 5) return(ok(st))
  ok(sets[[3]])
}
hr_lab <- function(x) formatC(x, format = "fg", digits = 3)
# a tidy 0-based axis whose upper limit is itself a labelled break
nice_axis <- function(v, cap = 1) {
  st <- diff(pretty(c(0, v), n = 5))[1]
  lim <- min(cap, st * ceiling(v / st - 1e-9))
  list(lim = lim, breaks = seq(0, lim + 1e-9, by = st))
}
spread <- function(y, gap) {
  o <- order(y)
  ys <- y[o]
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  out <- y
  out[o] <- ys
  out
}
f2 <- function(x) sprintf("%.2f", x)
f3 <- function(x) sprintf("%.3f", x)
ci2 <- function(e, l, h) sprintf("%.2f (%.2f to %.2f)", e, l, h)
pct <- function(x) paste0(trimws(formatC(100 * x, format = "fg",
                                         digits = 3)), "%")
FRT <- c(P = "At-Diagnosis", B = "At-Discharge")
FRS <- c(P = "at-diagnosis", B = "at-discharge")
CLOCK <- c(P = "admission", B = "discharge")
OCT <- c(death = "Death", comp = "Death or CV Readmission",
         cv = "CV Readmission")
OCS <- c(death = "death", comp = "death or CV readmission",
         cv = "CV readmission (death as a competing event)")
UNIT <- c(P = "admissions", B = "discharges")
cat(sprintf("font family: %s\n", fam))

# ================= ECG: the registered primary analysis =================
cat("\n=== ECG, PRIMARY ANALYSIS (FRAME B) ===\n")
fp <- file.path(RS, "primary.csv")
if (file.exists(fp)) {
  pr <- fread(fp)
  hr <- pr[grepl("^HR per SD", quantity) & method == "bootstrap percentile"]
  hr[, modality := ifelse(grepl("ECG", quantity), "ECG", "EHR")]
  hr[, per := ifelse(grepl("0-30", quantity), "Days 0 to 30",
                     "Days 31 to 365")]
  rat <- pr[grepl("^PRIMARY: ECG HR ratio", quantity)]
  npt <- pr[quantity == "patients"]$estimate
  lv <- c("EHR", "ECG")
  hr[, yy := match(modality, lv) + ifelse(per == "Days 0 to 30", 0.15,
                                          -0.15)]
  pc <- c("Days 0 to 30" = OI[5], "Days 31 to 365" = OI[2])
  ps <- c("Days 0 to 30" = 16, "Days 31 to 365" = 17)
  xr <- range(c(hr$lo, hr$hi))
  kx <- exp(log(xr[1]) + c(0.00, 0.30) * diff(log(xr)))
  ky <- length(lv) + 0.6
  t1 <- if (nrow(rat) && rat$lo > 1) {
    "The ECG Predicts Death More Strongly in the First 30 Days After Discharge"
  } else if (nrow(rat) && rat$hi < 1) {
    "The ECG Predicts Death Less Strongly in the First 30 Days After Discharge"
  } else {
    "The ECG's Hazard Ratio for Death in Early and Late Follow-Up After Discharge"
  }
  p <- ggplot(hr, aes(x = estimate, y = yy, colour = per, shape = per)) +
    annotate("segment", x = 1, xend = 1, y = 0.5, yend = length(lv) + 0.5,
             linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.1,
                  orientation = "y") +
    geom_point(size = 2.6) +
    annotate("point", x = kx, y = ky, colour = pc, shape = ps, size = 2.6) +
    annotate("text", x = kx * exp(0.03 * diff(log(xr))), y = ky,
             label = names(pc), hjust = 0, size = 3.5, family = fam,
             colour = pc) +
    scale_colour_manual(values = pc) + scale_shape_manual(values = ps) +
    scale_x_log10(breaks = hr_breaks(min(xr[1], 1), xr[2]), labels = hr_lab) +
    scale_y_continuous(breaks = seq_along(lv),
                       labels = c("EHR score", "ECG score"),
                       expand = expansion(add = c(0.4, 0.35))) +
    coord_cartesian(clip = "off") +
    labs(title = ttl(t1), x = "Hazard ratio per SD of the score (log scale)",
         y = NULL) + th
  save_fig(p, "ecg_hr_early_late_B.png", paste0(
    "Figure X: figure showing the hazard ratio of death per standard ",
    "deviation of the ECG and EHR scores in days 0 to 30 and days 31 to ",
    "365 after discharge, from the registered primary analysis (",
    format(npt, big.mark = ","), " discharges, EHR+ECG cohort), with 95% ",
    "patient-level bootstrap intervals from 2,000 refits. The ECG's early ",
    "to late ratio is ", ci2(rat$estimate, rat$lo, rat$hi), "."),
    w = 7.5, h = 3.6)

  # ---------- the ratio across sensitivity analyses ----------
  cat("\n=== ECG RATIO ACROSS SENSITIVITY ANALYSES ===\n")
  S <- list(data.table(analysis = "Registered primary analysis",
                       estimate = rat$estimate, lo = rat$lo, hi = rat$hi))
  add <- function(lab, d) {
    if (!is.null(d) && nrow(d) == 1) {
      S[[length(S) + 1]] <<- data.table(analysis = lab,
                                        estimate = d$estimate, lo = d$lo,
                                        hi = d$hi)
    } else {
      cat("   not found:", lab, "\n")
    }
  }
  rd <- function(f) {
    x <- file.path(RS, f)
    if (file.exists(x)) fread(x) else NULL
  }
  rb <- rd("primary_robust.csv")
  if (!is.null(rb)) {
    add("Scores winsorised", rb[variant == "winsorised" &
                                  grepl("^ECG HR ratio", quantity)])
    add("Scores rank-normalised", rb[variant == "rank_normal" &
                                       grepl("^ECG HR ratio", quantity)])
  }
  hp <- rd("primary_hospice.csv")
  if (!is.null(hp)) {
    add("Hospice discharges excluded",
        hp[grepl("^SENSITIVITY: ECG HR ratio", quantity)])
  }
  s5 <- rd("s5_primary.csv")
  if (!is.null(s5)) {
    add("Three-modality cohort", s5[startsWith(model, "(a)") &
                                      grepl("^ECG HR ratio", quantity)])
    add("Three-modality cohort, adjusted for CTPA",
        s5[startsWith(model, "(b)") & grepl("^ECG HR ratio", quantity)])
  }
  ab <- rd("ab_ecg_ratio.csv")
  if (!is.null(ab)) {
    add("Retuned EHR score, 28 features",
        ab[ehr_version == "EHR 28" & quantity == "ECG ratio"])
    add("Larger EHR score, 79 features",
        ab[ehr_version == "EHR 79" & quantity == "ECG ratio"])
    add("79 features plus missing-value indicators",
        ab[ehr_version == "EHR 79 + MI" & quantity == "ECG ratio"])
  }
  n3 <- ""
  if (!is.null(s5) && nrow(s5[quantity == "patients"])) {
    n3 <- sprintf(" The three-modality analyses use %s discharges.",
                  format(s5[quantity == "patients"]$estimate, big.mark = ","))
  }
  S <- rbindlist(S)
  print(S[, .(analysis, estimate = round(estimate, 3), lo = round(lo, 3),
              hi = round(hi, 3))], row.names = FALSE)
  S[, yy := rev(seq_len(.N))]
  S[, kind := ifelse(yy == max(yy), "primary", "sensitivity")]
  kc <- c(primary = OI[5], sensitivity = OI[4])
  ks <- c(primary = 16, sensitivity = 17)
  xr <- range(c(S$lo, S$hi, 1))
  t2 <- sprintf(paste("The ECG's Early to Late Hazard Ratio Stays Between",
                      "%.2f and %.2f Across Sensitivity Analyses"),
                min(S$estimate), max(S$estimate))
  if (min(S$estimate) <= 1) {
    t2 <- "The ECG's Early to Late Hazard Ratio Across Sensitivity Analyses"
  }
  p <- ggplot(S, aes(x = estimate, y = yy, colour = kind, shape = kind)) +
    annotate("segment", x = 1, xend = 1, y = 0.5, yend = nrow(S) + 0.5,
             linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25,
                  orientation = "y") +
    geom_point(size = 2.6) +
    scale_colour_manual(values = kc) + scale_shape_manual(values = ks) +
    scale_x_log10(breaks = hr_breaks(xr[1], xr[2]), labels = hr_lab) +
    scale_y_continuous(breaks = S$yy, labels = S$analysis,
                       expand = expansion(add = 0.5)) +
    labs(title = ttl(t2),
         x = "ECG hazard ratio, days 0 to 30 over days 31 to 365 (log scale)",
         y = NULL) + th
  save_fig(p, "ecg_ratio_sensitivity_B.png", paste0(
    "Figure X: figure showing the ratio of the ECG score's hazard ratio for ",
    "death in days 0 to 30 after discharge to that in days 31 to 365, in ",
    "the registered primary analysis and each sensitivity analysis, with ",
    "95% patient-level bootstrap intervals. Values above 1 mean the ECG's ",
    "signal is concentrated early.", n3, " The others use the EHR+ECG ",
    "cohort."),
    w = 7.5, h = 1.6 + 0.36 * nrow(S))
} else {
  cat("  ", fp, "not found: ECG figures skipped\n")
}

# ================= calibration (stage 13b) =================
cat("\n=== CALIBRATION ===\n")
need <- file.path(RS, c("fin_metrics.csv", "fin_cal_bins.csv",
                        "fin_cal_curves.csv", "fin_dca_curves.csv",
                        "fin_dca_key.csv"))
if (any(!file.exists(need))) {
  stop("run stage 13b first: missing ", basename(need[!file.exists(need)][1]))
}
MT <- fread(need[1])
BN <- fread(need[2])
CU <- fread(need[3])
DC <- fread(need[4])
KY <- fread(need[5])
H2 <- c(30, 365)
hc <- c("30" = OI[1], "365" = OI[4])
hs <- c("30" = 16, "365" = 17)
hl <- c("30" = "solid", "365" = "dashed")
for (fr in c("P", "B")) for (oc in c("death", "comp", "cv")) {
  for (co in c("EHR+ECG", "3-mod")) {
    cu <- CU[frame == fr & outcome == oc & cohort == co & horizon %in% H2]
    bn <- BN[frame == fr & outcome == oc & cohort == co & horizon %in% H2]
    mt <- MT[frame == fr & outcome == oc & cohort == co & horizon %in% H2]
    if (!nrow(cu)) next
    cu[, hz := as.character(horizon)]
    bn[, hz := as.character(horizon)]
    lim <- max(c(cu$risk, bn$mean_predicted, bn$observed), na.rm = TRUE)
    ax <- nice_axis(lim * 1.03)
    lim <- ax$lim
    cu[, `:=`(lo_c = pmin(pmax(lo, 0), lim), hi_c = pmin(pmax(hi, 0), lim),
              ob_c = pmin(pmax(observed, 0), lim))]
    # direct labels at the end of each curve
    lab <- cu[, .SD[which.max(risk)], by = hz]
    lab[, ylab := spread(ob_c, 0.05 * lim)]
    lab[, txt := sprintf("%s days", hz)]
    get <- function(h, m) mt[horizon == h & measure == m]$estimate
    close <- sapply(H2, function(h) {
      oe <- get(h, "oe")
      sl <- get(h, "slope")
      length(oe) && length(sl) && oe >= 0.9 && oe <= 1.1 && sl >= 0.9 &&
        sl <= 1.1
    })
    cpart <- if (co == "3-mod") " in the Three-Modality Cohort" else ""
    t3 <- if (all(close)) {
      sprintf("%s Predictions of %s Agree Closely With Observed Risk at 30 and 365 Days%s",
              FRT[[fr]], OCT[[oc]], cpart)
    } else if (close[1]) {
      sprintf("%s Predictions of %s Agree Closely With Observed Risk at 30 Days but Less at 365 Days%s",
              FRT[[fr]], OCT[[oc]], cpart)
    } else if (close[2]) {
      sprintf("%s Predictions of %s Agree Closely With Observed Risk at 365 Days but Less at 30 Days%s",
              FRT[[fr]], OCT[[oc]], cpart)
    } else {
      sprintf("%s Predictions of %s Depart From Observed Risk at 30 and 365 Days%s",
              FRT[[fr]], OCT[[oc]], cpart)
    }
    p <- ggplot() +
      annotate("segment", x = 0, y = 0, xend = lim, yend = lim,
               linetype = "dashed", colour = "grey55") +
      geom_ribbon(data = cu, aes(x = risk, ymin = lo_c, ymax = hi_c,
                                 fill = hz), alpha = 0.15) +
      geom_line(data = cu, aes(x = risk, y = ob_c, colour = hz,
                               linetype = hz), linewidth = 0.8) +
      geom_point(data = bn[!is.na(observed)],
                 aes(x = mean_predicted, y = pmin(observed, lim),
                     colour = hz, shape = hz), size = 2.2) +
      geom_text(data = lab, aes(x = risk + 0.02 * lim, y = ylab,
                                label = txt, colour = hz),
                hjust = 0, size = 3.6, family = fam) +
      scale_colour_manual(values = hc) + scale_fill_manual(values = hc) +
      scale_shape_manual(values = hs) + scale_linetype_manual(values = hl) +
      scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                         breaks = ax$breaks,
                         expand = expansion(mult = c(0.01, 0.16))) +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                         breaks = ax$breaks,
                         expand = expansion(mult = c(0.01, 0.03))) +
      coord_cartesian(xlim = c(0, lim), ylim = c(0, lim), clip = "off") +
      labs(title = ttl(t3, 56), x = "Predicted risk", y = "Observed risk") +
      th
    n <- mt$n[1]
    nm <- sprintf("cal_%s_%s_%s.png", fr, oc,
                  if (co == "EHR+ECG") "ehrecg" else "3mod")
    save_fig(p, nm, paste0(
      "Figure X: figure showing the calibration of the ", FRS[[fr]],
      " model for ", OCS[[oc]], " at 30 and 365 days after ", CLOCK[[fr]],
      " (", if (co == "EHR+ECG") "EHR+ECG cohort" else
        "three-modality cohort, with CTPA", ", ", format(n, big.mark = ","),
      " ", UNIT[[fr]], "). Lines are loess calibration curves with 95% ",
      "patient-clustered bootstrap bands; points are deciles of predicted ",
      "risk, merged from the lowest until each holds at least 10 events. ",
      "The dashed line marks perfect calibration. Calibration slope ",
      f2(get(30, "slope")), " at 30 days and ", f2(get(365, "slope")),
      " at 365 days; integrated calibration index ", f3(get(30, "ici")),
      " and ", f3(get(365, "ici")), "."), w = 6.2, h = 6.2)
  }
}

# ================= decision curves (stage 13b) =================
cat("\n=== DECISION CURVES (EHR+ECG cohort) ===\n")
SC <- c("Model" = OI[5], "Treat all" = OI[8], "Treat none" = OI[7],
        "sPESI 1 or more" = OI[4], "PESI-10 class III to V" = OI[3])
SL <- c("Model" = "solid", "Treat all" = "dashed", "Treat none" = "dotted",
        "sPESI 1 or more" = "dotdash", "PESI-10 class III to V" = "longdash")
SLAB <- c("Model" = "Model", "Treat all" = "Treat all",
          "Treat none" = "Treat none", "sPESI 1 or more" = "sPESI 1 or more",
          "PESI-10 class III to V" = "PESI-10 class III to V")
XEXP <- 0.36
# place each line's label where it is clearest: candidate spots above or
# below, right or left of, every point of the line are scored by their gap
# to every line and to the labels already placed (more is better, capped),
# with a small preference for the right; a label must fit inside the
# panel. cw and lh are a character's width and a label's height as shares
# of the panel.
place_labels <- function(d, labs, x0, x1, ylo, yhi, cw = 0.0115,
                         lh = 0.06) {
  R <- yhi - ylo
  W <- x1 - x0
  Y <- dcast(d, threshold ~ strategy, value.var = "net_benefit")
  xs <- Y$threshold
  curves <- setdiff(names(Y), "threshold")
  boxes <- list()
  out <- list()
  for (sn in names(labs)) {
    v <- Y[[sn]]
    if (is.null(v)) next
    wd <- nchar(labs[[sn]]) * cw * W
    ht <- lh * R
    best <- NULL
    bs <- -Inf
    for (i in seq_along(xs)) {
      if (!is.finite(v[i]) || v[i] < ylo || v[i] > yhi) next
      for (pos in 1:4) {
        # above right, below right, above left, below left of the point
        xa <- if (pos <= 2) xs[i] + 0.008 * W else xs[i] - 0.008 * W - wd
        yb <- if (pos %% 2 == 1) v[i] + 0.015 * R else v[i] - 0.015 * R - ht
        if (xa < x0 || xa + wd > x1 || yb < ylo + 0.005 * R ||
            yb + ht > yhi) next
        gx <- seq(xa, xa + wd, length.out = 15)
        gap <- Inf
        bad <- FALSE
        for (cn in curves) {
          cy <- suppressWarnings(approx(xs, Y[[cn]], gx)$y)
          cy <- cy[is.finite(cy) & cy >= ylo & cy <= yhi]
          if (!length(cy)) next
          if (any(cy >= yb & cy <= yb + ht)) {
            bad <- TRUE
            break
          }
          gap <- min(gap, abs(cy - yb), abs(cy - yb - ht))
        }
        if (bad) next
        for (bx in boxes) {
          if (xa < bx[2] && xa + wd > bx[1]) {
            if (yb < bx[4] && yb + ht > bx[3]) {
              bad <- TRUE
              break
            }
            gap <- min(gap, abs(yb - bx[4]), abs(bx[3] - yb - ht))
          }
        }
        if (bad) next
        sc <- min(gap / R, 0.08) + 0.03 * (xa - x0) / W
        if (sc > bs) {
          bs <- sc
          best <- c(xa, yb)
        }
      }
    }
    if (is.null(best)) {
      j <- max(which(is.finite(v) & v >= ylo & v <= yhi))
      best <- c(xs[j] + 0.008 * W, v[j])
    }
    boxes[[length(boxes) + 1]] <- c(best[1], best[1] + wd, best[2],
                                    best[2] + ht)
    out[[length(out) + 1]] <- data.table(strategy = sn, txt = labs[[sn]],
                                         x = best[1], y = best[2])
  }
  rbindlist(out)
}
for (fr in c("P", "B")) for (oc in c("death", "comp", "cv")) {
  for (h in H2) {
    d <- DC[frame == fr & outcome == oc & cohort == "EHR+ECG" & horizon == h]
    if (!nrow(d)) next
    prev <- d$prevalence[1]
    xmax <- min(0.5, max(0.10, ceiling(2 * prev / 0.05) * 0.05))
    d <- d[threshold <= xmax + 1e-9]
    d <- rbind(d, d[strategy == "Model"][, `:=`(strategy = "Treat none",
                                                net_benefit = 0)])
    ylo <- -0.15 * prev
    yhi <- 1.05 * prev
    d[, strategy := factor(strategy, levels = names(SC))]
    # where the model gives the most net benefit (thresholds of 1% or more)
    w <- dcast(d, threshold ~ strategy, value.var = "net_benefit")
    oth <- setdiff(names(w), c("threshold", "Model"))
    best <- w$Model >= do.call(pmax, w[, ..oth]) - 1e-6 & w$threshold >= 0.01
    r <- rle(best)
    ends <- cumsum(r$lengths)
    sts <- ends - r$lengths + 1
    k <- which(r$values)
    hzs <- sprintf("%d-Day %s", h, OCT[[oc]])
    if (length(k)) {
      k <- k[which.max(r$lengths[k])]
      a <- w$threshold[sts[k]]
      b <- w$threshold[ends[k]]
      t4 <- if (abs(b - xmax) < 1e-9 && a <= 0.0101) {
        sprintf("The %s Model Gives the Most Net Benefit for %s at Every Threshold From %s to %s",
                FRT[[fr]], hzs, pct(a), pct(b))
      } else {
        sprintf("The %s Model Gives the Most Net Benefit for %s at Thresholds From %s to %s",
                FRT[[fr]], hzs, pct(a), pct(b))
      }
    } else {
      t4 <- sprintf("The %s Model Does Not Give the Most Net Benefit for %s",
                    FRT[[fr]], hzs)
    }
    # direct labels where each is clearest (see place_labels)
    x1 <- xmax * (1 + XEXP)
    labs_txt <- SLAB[intersect(names(SLAB), as.character(unique(d$strategy)))]
    lab <- place_labels(d, labs_txt, -0.01 * xmax, x1, ylo, yhi)
    xbr <- pretty(c(0, xmax), n = 6)
    xbr <- xbr[xbr <= xmax + 1e-9]
    p <- ggplot(d, aes(x = threshold, y = net_benefit, colour = strategy,
                       linetype = strategy)) +
      geom_line(linewidth = 0.8) +
      geom_text(data = lab, aes(x = x, y = y, label = txt),
                hjust = 0, vjust = 0, size = 3.4, family = fam,
                inherit.aes = FALSE, colour = SC[lab$strategy]) +
      scale_colour_manual(values = SC) + scale_linetype_manual(values = SL) +
      scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                         breaks = xbr,
                         expand = expansion(mult = c(0.01, XEXP))) +
      coord_cartesian(xlim = c(0, xmax), ylim = c(ylo, yhi)) +
      labs(title = ttl(t4), x = "Risk threshold for treatment",
           y = "Net benefit") + th
    ky <- KY[frame == fr & outcome == oc & cohort == "EHR+ECG" &
               horizon == h & quantity == "model minus strategy" &
               strategy == "Treat all"]
    kt <- paste(sprintf("%s at %s", f3(ky$estimate), pct(ky$threshold)),
                collapse = ", ")
    who <- if (fr == "P" && oc == "death") {
      paste0(", of treating everyone, of treating no one, of the sPESI rule ",
             "(1 or more as high risk) and of the PESI-10 rule (class III ",
             "to V)")
    } else ", of treating everyone and of treating no one"
    nm <- sprintf("dca_%s_%s_%d.png", fr, oc, h)
    save_fig(p, nm, paste0(
      "Figure X: figure showing decision curves for ", OCS[[oc]],
      " within ", h, " days of ", CLOCK[[fr]], ": the net benefit of the ",
      FRS[[fr]], " model", who, ", across risk thresholds up to ",
      pct(xmax), " (EHR+ECG cohort; event in ", pct(prev), "). The ",
      "model's net benefit over treating everyone is ", kt, "."),
      w = 7.5, h = 4.6)
  }
}

caps <- rbindlist(CAP)
fwrite(caps, file.path(RS, "fin_fig_captions.csv"))
cat(sprintf("\n  %d figures written\n", nrow(caps)))
cat("\nOUTPUT FILES\n  ", FD, " (", nrow(caps), " PNGs)\n", sep = "")
cat("  ", file.path(RS, "fin_fig_captions.csv"), "\n")
cat("  ", file.path(RS, "fin_fig_log.txt"), " (terminal log)\n")
cat("\nOutputs: fin_fig_captions.csv and the PNGs in", FD, "\n")