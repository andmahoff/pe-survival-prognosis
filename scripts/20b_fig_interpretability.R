# 20b_fig_interpretability.R
# Figures for 20a (A23): partial effects of the spline Cox model and
# SurvSHAP(t) importance for the death model.
#
# Frame B, EHR+ECG cohort, death to 365 days after discharge. Each figure
# has one panel, direct labels instead of legends, and Okabe-Ito colours
# with a redundant channel (shape or line type), in Calibri where available.
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
# Outputs: R Graphs/interp/*.png, results/interp_fig_log.txt
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
fam <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  sf <- systemfonts::system_fonts()
  if (any(sf$family == "Calibri")) fam <- "Calibri"
}
dev_png <- function(path, w = 7, h = 4.5) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = w, height = h, units = "in", res = 300)
  } else {
    png(path, width = w, height = h, units = "in", res = 300)
  }
}
th <- theme_minimal(base_size = 12, base_family = fam) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10.5, colour = "grey25"),
        plot.caption = element_text(size = 9, colour = "grey35",
                                    hjust = 0),
        axis.title = element_text(size = 11))
save_fig <- function(p, name, w = 7, h = 4.5) {
  path <- file.path(FD, name)
  dev_png(path, w, h)
  print(p)
  invisible(dev.off())
  cat("  saved", path, "\n")
}
safe <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
wrapc <- function(s, w = 95) paste(strwrap(s, w), collapse = "\n")
# spread direct labels apart vertically so none overlap (no leader lines)
spread <- function(y, gap) {
  o <- order(y)
  ys <- y[o]
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
  r <- tb[feature == ft]
  sub <- sprintf("%s term, edf %.1f, chi-square %.1f, p %s", r$type, r$edf,
                 r$chi_sq, format.pval(r$p_value, digits = 2, eps = 1e-4))
  cap <- wrapc(paste("Penalised-spline Cox model, Frame B (death to 365",
                     "days after discharge), EHR+ECG cohort. Hazard ratio",
                     "relative to the median; shaded band and bars: 95%",
                     "interval. Ticks on the axis: 10th percentile, median",
                     "and 90th percentile."))
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
  } else {
    p <- ggplot(d, aes(x = factor(signif(x, 3)), y = exp(log_hr))) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi)), width = 0.12,
                    colour = OI[4]) +
      geom_point(colour = OI[4], size = 2.6)
  }
  p <- p + scale_y_log10() +
    labs(title = ft, subtitle = sub, caption = cap, x = ft,
         y = "Hazard ratio (log scale)") + th
  save_fig(p, sprintf("pe_%s_%s.png", FR, safe(ft)))
}

# ---------- SurvSHAP(t) share over time ----------
cat("\n=== SurvSHAP(t) SHARE OVER TIME ===\n")
for (mn in unique(im$model)) {
  x <- im[frame == FR & cohort == CO & model == mn]
  if (!nrow(x)) next
  hmax <- max(x$horizon)
  top <- x[horizon == hmax][order(rank)]$feature[seq_len(NTOP)]
  d <- x[feature %in% top]
  d[, feature := factor(feature, levels = top)]
  lab <- d[horizon == hmax]
  lab[, ylab := spread(100 * share, 0.045 * diff(range(100 * d$share)))]
  p <- ggplot(d, aes(x = horizon, y = 100 * share, colour = feature,
                     shape = feature, linetype = feature)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    geom_text(data = lab, aes(y = ylab, label = feature), hjust = 0,
              nudge_x = 0.06, size = 3.6, family = fam,
              show.legend = FALSE) +
    scale_x_log10(breaks = sort(unique(d$horizon)),
                  expand = expansion(mult = c(0.03, 0.30))) +
    scale_colour_manual(values = OI) + scale_shape_manual(values = SHAPES) +
    scale_linetype_manual(values = LTY) +
    labs(title = wrapc(sprintf("EHR feature importance over time: %s",
                               mn), 60),
         subtitle = paste("Share of mean |SurvSHAP(t)| at each horizon,",
                          "8 largest at 365 days"),
         caption = wrapc(paste("Frame B, EHR+ECG cohort, death after",
                               "discharge. SurvSHAP(t) on the risk scale,",
                               "estimated by antithetic permutation",
                               "sampling.")),
         x = "Days after discharge (log scale)",
         y = "Share of importance (%)") + th
  save_fig(p, sprintf("shap_share_%s_%s.png", FR, safe(mn)))
}

# ---------- agreement ----------
cat("\n=== IMPORTANCE AGREEMENT ===\n")
sp <- tm[frame == FR & cohort == CO, .(feature, v = term_sd / max(term_sd))]
sp[, model := "spline Cox (term SD)"]
hmax <- max(im[frame == FR]$horizon)
fo <- im[frame == FR & cohort == CO & horizon == hmax,
         .(v = mean_abs_shap / max(mean_abs_shap), feature), by = model]
d <- rbind(sp, fo[, .(feature, v, model)])
ordf <- sp[order(v)]$feature
d[, feature := factor(feature, levels = ordf)]
mods <- unique(d$model)
cols <- OI[c(4, 1, 3)][seq_along(mods)]
shp <- c(16, 17, 15)[seq_along(mods)]
ytop <- length(ordf) + 1.2
kx <- c(0.00, 0.36, 0.72)[seq_along(mods)]
p <- ggplot(d, aes(x = v, y = feature, colour = model, shape = model)) +
  geom_point(size = 2.4) +
  annotate("point", x = kx, y = ytop, colour = cols, shape = shp,
           size = 2.4) +
  annotate("text", x = kx + 0.025, y = ytop, label = mods, hjust = 0,
           size = 3.4, family = fam, colour = cols) +
  scale_colour_manual(values = setNames(cols, mods)) +
  scale_shape_manual(values = setNames(shp, mods)) +
  scale_x_continuous(limits = c(0, 1.05)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.6))) +
  coord_cartesian(clip = "off") +
  labs(title = "Do the spline model and the forests agree on what matters?",
       subtitle = sprintf(paste("Importance of each EHR feature, scaled to",
                                "each model's maximum (forests at %d days)"),
                          hmax),
       caption = "Frame B, EHR+ECG cohort, death after discharge.",
       x = "Relative importance", y = NULL) + th
save_fig(p, sprintf("importance_agreement_%s.png", FR), w = 7, h = 7.5)
cat("\nOUTPUT FILES\n  ", FD, " (pe_*.png, shap_share_*.png,",
    "importance_agreement_*.png)\n")
cat("  ", file.path(RS, "interp_fig_log.txt"), " (terminal log)\n")