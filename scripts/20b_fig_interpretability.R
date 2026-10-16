# 20b_fig_interpretability.R
# Figures for 20a (A23): partial effects of the spline Cox model and
# SurvSHAP(t) importance for the death model.
#
# Frame B, EHR+ECG cohort, death to 365 days after discharge. As in the
# other PE projects, each figure has one panel, a centred bold Title Case
# title that describes what it shows without stating the result, direct
# labels instead of legends, Okabe-Ito colours with shapes or line types
# as a second channel, faint gridlines and Calibri where available, with
# no notes on the graphic (suggested captions are written to a CSV).
# PNGs are written to R Graphs/interp/ in the project folder:
#   pe_B_<feature>.png         hazard ratio across the feature's 1st to
#                              99th percentile, relative to its median,
#                              from the spline Cox model; one figure for
#                              each of the 8 strongest terms
#   shap_share_B_<model>.png   each feature's share of SurvSHAP(t)
#                              importance at each horizon, 8 largest at
#                              365 days; one figure per forest
#   importance_agreement_B.png importance of every EHR feature in the
#                              spline model and the two forests, each
#                              scaled to its own maximum
#
# Inputs:  results/interp_spline_terms.csv, interp_spline_curves.csv,
#          interp_shap_importance.csv
# Outputs: R Graphs/interp/*.png, results/interp_fig_captions.csv,
#          results/interp_fig_log.txt
# Usage:   Rscript scripts/20b_fig_interpretability.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
RS <- file.path(P2, "results")
FD <- file.path(P2, "R Graphs", "interp")
dir.create(FD, showWarnings = FALSE, recursive = TRUE)
FR <- "B"
CO <- "EHR+ECG"
NTOP <- 8L
OI <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
        "#000000", "#999999")
SHAPES <- c(16, 17, 15, 18, 1, 2, 0, 5)
LTY <- c("solid", "dashed", "dotdash", "longdash", "solid", "dashed",
         "dotdash", "longdash")
MODEL_LONG <- c("RSF, Brier-tuned" = "random survival forest",
                "oblique RSF, Brier-tuned" =
                  "oblique random survival forest")
MODEL_SHORT <- c("RSF, Brier-tuned" = "RSF",
                 "oblique RSF, Brier-tuned" = "Oblique RSF")
# readable names and units for the 28 features of the dissertation's EHR,
# named as in stage 20q
LABEL <- c(
  age = "Age", aniongap = "Anion gap", bicarb = "Bicarbonate",
  bun = "Urea nitrogen", calcium = "Calcium", cancer = "Cancer",
  chloride = "Chloride", copd = "COPD", creatinine = "Creatinine",
  dbp = "Diastolic BP", glucose = "Glucose",
  hct = "Haematocrit", heart_failure = "Heart failure",
  afib = "Atrial fibrillation", hgb = "Haemoglobin", hr = "Heart rate",
  mch = "MCH", mchc = "MCHC", mcv = "MCV", plt = "Platelets",
  potassium = "Potassium", rbc = "Red cell count", rdw = "RDW",
  rr = "Respiratory rate", sbp = "Systolic BP",
  sodium = "Sodium", temp = "Temperature", wbc = "White cell count")
UNITS <- c(
  age = "years", aniongap = "mEq/L", bicarb = "mEq/L", bun = "mg/dL",
  calcium = "mg/dL", chloride = "mEq/L", creatinine = "mg/dL",
  dbp = "mmHg", glucose = "mg/dL", hct = "%", hgb = "g/dL",
  hr = "beats/min", mch = "pg", mchc = "g/dL", mcv = "fL", plt = "K/uL",
  potassium = "mEq/L", rbc = "m/uL", rdw = "%", rr = "breaths/min",
  sbp = "mmHg", sodium = "mEq/L", temp = "degrees C", wbc = "K/uL")
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
safe <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
# Title Case as in the other PE projects: every word capitalised except
# short articles, conjunctions and prepositions after the first; words that
# already hold a capital (EHR, pH) are left as they are
tc <- function(s) {
  small <- c("a", "an", "the", "and", "but", "or", "nor", "for", "of",
             "on", "in", "to", "by", "at", "per", "as", "via")
  w <- strsplit(s, " ", fixed = TRUE)[[1]]
  cap <- function(x) {
    if (grepl("[A-Z]", x) || !grepl("^[a-z]", x)) return(x)
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }
  for (i in seq_along(w)) {
    if (i > 1 && w[i] %in% small) next
    parts <- strsplit(w[i], "-", fixed = TRUE)[[1]]
    w[i] <- paste(vapply(parts, cap, ""), collapse = "-")
  }
  paste(w, collapse = " ")
}
# Title Case, wrapped into lines of similar length so no word is left on
# its own
ttl <- function(s, w = 72) {
  s <- tc(s)
  n <- length(strwrap(s, w))
  if (n > 1) {
    for (k in seq(ceiling(nchar(s) / n), w)) {
      if (length(strwrap(s, k)) == n) {
        w <- k
        break
      }
    }
  }
  paste(strwrap(s, w), collapse = "\n")
}
lab_of <- function(f) {
  out <- f
  hit <- f %in% names(LABEL)
  out[hit] <- unname(LABEL[f[hit]])
  out
}
# lower-case the first letter for use mid-sentence, unless it starts an
# abbreviation
lc1 <- function(s) {
  keep <- grepl("^[A-Z][A-Z0-9-]", s)
  ifelse(keep, s, paste0(tolower(substring(s, 1, 1)), substring(s, 2)))
}
# labelled ticks for a hazard-ratio axis on the log scale
hr_breaks <- function(lo, hi) {
  r <- range(c(lo, hi), na.rm = TRUE)
  ok <- function(b) b[b >= r[1] / 1.02 & b <= r[2] * 1.02]
  b <- ok(c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 7, 10, 20))
  if (length(b) < 3) {
    b <- ok(c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75, 2))
  }
  b
}
# spread direct labels apart vertically so none overlap (no leader lines),
# then centre each group of touching labels on its own lines, keeping the
# lowest label at or above lo
spread <- function(y, gap, lo = -Inf) {
  o <- order(y)
  yo <- y[o]
  ys <- yo
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  g <- cumsum(c(TRUE, diff(ys) > gap + 1e-9))
  for (k in unique(g)) {
    j <- which(g == k)
    ys[j] <- ys[j] - mean(ys[j] - yo[j])
  }
  ys[1] <- max(ys[1], lo)
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  out <- y
  out[o] <- ys
  out
}
tm <- fread(file.path(RS, "interp_spline_terms.csv"))
cv <- fread(file.path(RS, "interp_spline_curves.csv"))
im <- fread(file.path(RS, "interp_shap_importance.csv"))
cat(sprintf("font family: %s\n", fam))

# ---------- partial effects ----------
cat("\n=== PARTIAL EFFECTS ===\n")
tb <- tm[frame == FR & cohort == CO][order(-chi_sq)]
for (ft in tb$feature[seq_len(min(NTOP, nrow(tb)))]) {
  d <- cv[frame == FR & cohort == CO & feature == ft]
  lb <- lab_of(ft)
  un <- if (ft %in% names(UNITS)) sprintf(" (%s)", UNITS[[ft]]) else ""
  if (d$type[1] == "spline") {
    p <- ggplot(d, aes(x = x, y = exp(log_hr))) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), fill = OI[4],
                  alpha = 0.18) +
      geom_line(colour = OI[4], linewidth = 0.9) +
      geom_rug(data = data.frame(x = c(d$x_p10[1], d$x_median[1],
                                       d$x_p90[1])),
               aes(x = x), inherit.aes = FALSE, sides = "b",
               length = unit(0.04, "npc"), linewidth = 0.7)
    xl <- paste0(lb, un)
  } else {
    bin <- nrow(d) == 2 && all(sort(round(d$x, 8)) == c(0, 1))
    d[, xf := if (bin) factor(ifelse(round(x) == 1, "Yes", "No"),
                              levels = c("No", "Yes")) else
                                factor(signif(x, 3))]
    p <- ggplot(d, aes(x = xf, y = exp(log_hr))) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi)), width = 0.12,
                    colour = OI[4]) +
      geom_point(colour = OI[4], size = 2.6)
    xl <- if (bin) paste(lb, "recorded") else paste0(lb, un)
  }
  yb <- hr_breaks(exp(d$lo), exp(d$hi))
  p <- p + scale_y_log10(breaks = yb, labels = as.character(yb)) +
    labs(title = ttl(sprintf(paste("effect of %s on the hazard of death",
                                   "after discharge"), lc1(lb))),
         x = xl, y = "Hazard ratio relative to the median (log scale)") +
    th
  cap <- if (d$type[1] == "spline") {
    paste0(
      "Figure X: figure showing the hazard ratio for death within one ",
      "year of discharge across values of ", lc1(lb), ", relative to its ",
      "median. Estimates come from the penalised spline Cox model (Frame ",
      "B, EHR+ECG cohort) with a shaded 95% interval, and ticks mark the ",
      "10th percentile, median and 90th percentile.")
  } else {
    paste0(
      "Figure X: figure showing the hazard ratio for death within one ",
      "year of discharge by ", lc1(lb), ", relative to its median value. ",
      "Estimates come from the linear term of the spline Cox model ",
      "(Frame B, EHR+ECG cohort), with 95% intervals.")
  }
  save_fig(p, sprintf("pe_%s_%s.png", FR, safe(ft)), cap)
}

# ---------- SurvSHAP(t) share over time ----------
cat("\n=== SurvSHAP(t) SHARE OVER TIME ===\n")
for (mn in unique(im$model)) {
  x <- im[frame == FR & cohort == CO & model == mn]
  if (!nrow(x)) next
  ml <- if (mn %in% names(MODEL_LONG)) MODEL_LONG[[mn]] else mn
  hmax <- max(x$horizon)
  top <- x[horizon == hmax][order(rank)]$feature[seq_len(NTOP)]
  d <- x[feature %in% top]
  d[, grp := factor(feature, levels = top)]
  d[, label := lab_of(feature)]
  h <- 5.4
  lab <- d[horizon == hmax]
  gap <- 0.17 / (h - 1.1) * diff(range(c(0, 100 * d$share)))
  lab[, ylab := spread(100 * share, gap, lo = 0)]
  p <- ggplot(d, aes(x = horizon, y = 100 * share, colour = grp,
                     shape = grp, linetype = grp)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    geom_text(data = lab, aes(y = ylab, label = label), hjust = 0,
              nudge_x = 0.06, size = 3.5, family = fam,
              show.legend = FALSE) +
    scale_x_log10(breaks = sort(unique(d$horizon)),
                  expand = expansion(mult = c(0.03, 0.62))) +
    scale_colour_manual(values = OI) + scale_shape_manual(values = SHAPES) +
    scale_linetype_manual(values = LTY) +
    scale_y_continuous(breaks = function(l) pretty(l, n = 6)) +
    expand_limits(y = 0) +
    labs(title = ttl(sprintf(
           "share of EHR feature importance over time in the %s", ml)),
         x = "Days after discharge (log scale)",
         y = "Share of importance (%)") + th
  save_fig(p, sprintf("shap_share_%s_%s.png", FR, safe(mn)), paste0(
    "Figure X: figure showing each feature's share of mean absolute ",
    "SurvSHAP(t) value at 10, 30, 90, 180 and 365 days after discharge in ",
    "the ", ml, " (Frame B, EHR+ECG cohort), for the eight features with ",
    "the largest share at 365 days. SurvSHAP(t) values are on the risk ",
    "scale and were estimated by antithetic permutation sampling."),
    w = 7.5, h = h)
}

# ---------- agreement ----------
cat("\n=== IMPORTANCE AGREEMENT ===\n")
sp <- tm[frame == FR & cohort == CO, .(feature, v = term_sd / max(term_sd))]
sp[, model := "Spline Cox"]
hmax <- max(im[frame == FR]$horizon)
fo <- im[frame == FR & cohort == CO & horizon == hmax,
         .(v = mean_abs_shap / max(mean_abs_shap), feature), by = model]
fo[, model := ifelse(model %in% names(MODEL_SHORT),
                     MODEL_SHORT[model], model)]
d <- rbind(sp, fo[, .(feature, v, model)])
ordf <- lab_of(sp[order(v)]$feature)
d[, label := factor(lab_of(feature), levels = ordf)]
mods <- unique(d$model)
cols <- OI[c(4, 1, 3)][seq_along(mods)]
shp <- c(16, 17, 15)[seq_along(mods)]
ytop <- length(ordf) + 1.2
kx <- c(0.00, 0.30, 0.52)[seq_along(mods)]
p <- ggplot(d, aes(x = v, y = label, colour = model, shape = model)) +
  geom_point(size = 2.4) +
  annotate("point", x = kx, y = ytop, colour = cols, shape = shp,
           size = 2.4) +
  annotate("text", x = kx + 0.025, y = ytop, label = mods, hjust = 0,
           size = 3.2, family = fam, colour = cols) +
  scale_colour_manual(values = setNames(cols, mods)) +
  scale_shape_manual(values = setNames(shp, mods)) +
  scale_x_continuous(limits = c(0, 1.05)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.6))) +
  coord_cartesian(clip = "off") +
  labs(title = ttl("relative importance of the EHR features in three models"),
       x = "Relative importance (each model's largest = 1)", y = NULL) + th
save_fig(p, sprintf("importance_agreement_%s.png", FR), paste0(
  "Figure X: figure showing the importance of each of the 28 EHR features ",
  "for death after discharge (Frame B, EHR+ECG cohort), scaled to each ",
  "model's largest value. The spline Cox model is measured by the spread ",
  "of each term, and the random survival forest (RSF) and oblique RSF by ",
  sprintf("mean absolute SurvSHAP(t) value at %d days.", hmax)),
  w = 7.5, h = 2 + 0.24 * length(ordf))
caps <- rbindlist(CAP)
fwrite(caps, file.path(RS, "interp_fig_captions.csv"))
cat("\nOUTPUT FILES\n  ", FD, " (pe_*.png, shap_share_*.png,",
    "importance_agreement_*.png)\n")
cat("  ", file.path(RS, "interp_fig_captions.csv"), "\n")
cat("  ", file.path(RS, "interp_fig_log.txt"), " (terminal log)\n")
