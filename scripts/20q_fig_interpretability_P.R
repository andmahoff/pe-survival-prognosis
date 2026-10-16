# 20q_fig_interpretability_P.R
# Figures for 20p (A37): interpretability of the A36 at-diagnosis model, and
# for 16p, the CTPA block split at diagnosis.
#
# Frame P, EHR+ECG cohort, death within 365 days of admission. Each figure
# has one panel, a centred bold Title Case title that states what it shows,
# direct labels instead of legends, Okabe-Ito colours with shapes or line
# types as a second channel, faint gridlines and Calibri where available,
# with no notes on the graphic (suggested captions are written to a CSV).
# PNGs are written to R Graphs/interpP/ in the project folder:
#   pe_P_<feature>.png          hazard ratio relative to the median
#                               across the feature's 1st to 99th
#                               percentile, from the spline Cox model;
#                               one figure for each of the 8 largest
#                               terms (term SD)
#   shap_share_P_<model>.png    each feature's share of SurvSHAP(t)
#                               importance at each horizon, 8 largest
#                               at 365 days; one figure per forest
#   shap_domain_P_<model>.png   the same for the five EHR domains
#   importance_agreement_P.png  relative importance in the spline model
#                               and both forests, for the features in
#                               any model's top 10
#   ctpa_blocks_P.png           hazard ratio per SD of each CTPA block
#                               for deaths in days 0-30 and 31-365
#                               (16p; skipped if its results are
#                               missing)
#
# Charlson and Elixhauser indicators built from the same codes (for example
# metastatic cancer) are identical columns, so each model splits their
# credit in two. Such twins are found from the 20p spline terms (identical
# statistics in both cohorts) and, in the importance figures only, shown as
# one feature: the SurvSHAP(t) values of the pair are summed per admission
# before taking the mean absolute value, and the spline term SDs are summed.
# The twins found are listed in results/interpP_fig_twins.csv.
#
# Inputs:  results/interpP_{spline_terms,spline_curves,shap_importance,
#          shap_domains}.csv, results/ctpa_blocks_P_time.csv (optional),
#          data/processed/interpP/shap_P_EHRECG_*.csv (patient level,
#          not in the repo)
# Outputs: R Graphs/interpP/*.png, results/interpP_fig_captions.csv,
#          results/interpP_fig_twins.csv, results/interpP_fig_log.txt
# Usage:   Rscript scripts/20q_fig_interpretability_P.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
RS <- file.path(P2, "results")
OD <- file.path(P2, "data", "processed", "interpP")
FD <- file.path(P2, "R Graphs", "interpP")
dir.create(FD, showWarnings = FALSE, recursive = TRUE)
CO <- "EHR+ECG"
NTOP <- 8L
NAGR <- 10L
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
# clearer names for MIMIC items (stage 20p's labels otherwise)
RELAB <- c(
  lab_50820_ph_first = "Blood pH", lab_51491_ph_first = "Urine pH",
  lab_50861_alanine_aminotransferase_alt_first = "ALT",
  lab_50878_asparate_aminotransferase_as_first = "AST",
  lab_50911_creatine_kinase_mb_isoenzyme_first = "CK-MB",
  lab_50910_creatine_kinase_ck_first = "Creatine kinase",
  lab_50954_lactate_dehydrogenase_ld_first = "LDH",
  lab_51516_wbc_first = "Urine white cells",
  lab_51493_rbc_first = "Urine red cells",
  lab_51498_specific_gravity_first = "Urine specific gravity",
  lab_51476_epithelial_cells_first = "Urine epithelial cells",
  lab_51082_creatinine_urine_first = "Urine creatinine",
  lab_51492_protein_first = "Urine protein",
  lab_52172_rdw_sd_first = "RDW-SD",
  lab_51214_fibrinogen_functional_first = "Fibrinogen",
  lab_50804_calculated_total_co2_first = "Total CO2 (blood gas)",
  lab_50808_free_calcium_first = "Ionised calcium",
  lab_50953_iron_binding_capacity_total_first =
    "Total iron binding capacity")
# units for the partial-effect x-axes
UNITS <- c(
  age = "years", hx_charlson_score = "points",
  hx_elixhauser_vw = "points", creatinine_first = "mg/dL",
  bun_first = "mg/dL", sodium_first = "mEq/L", potassium_first = "mEq/L",
  bicarb_first = "mEq/L", glucose_first = "mg/dL", hgb_first = "g/dL",
  wbc_first = "K/uL", plt_first = "K/uL", hr_first = "beats/min",
  sbp_first = "mmHg", dbp_first = "mmHg", rr_first = "breaths/min",
  spo2_first = "%", temp_first = "degrees C",
  lab_50813_lactate_first = "mmol/L", lab_51003_troponin_t_first = "ng/mL",
  lab_50963_ntprobnp_first = "pg/mL",
  lab_50861_alanine_aminotransferase_alt_first = "IU/L",
  lab_50878_asparate_aminotransferase_as_first = "IU/L",
  lab_50863_alkaline_phosphatase_first = "IU/L",
  lab_50954_lactate_dehydrogenase_ld_first = "IU/L",
  lab_50911_creatine_kinase_mb_isoenzyme_first = "ng/mL",
  lab_50862_albumin_first = "g/dL",
  lab_51244_lymphocytes_first = "% of white cells",
  lab_50818_pco2_first = "mmHg", lab_50821_po2_first = "mmHg",
  lab_50802_base_excess_first = "mEq/L", rdw_first = "%")
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
# lower-case the first letter for use mid-sentence, unless it starts a
# proper noun or an abbreviation
lc1 <- function(s) {
  keep <- grepl("^(Charlson|Elixhauser|[A-Z][A-Z0-9-]|pH|pO2|pCO2)", s)
  ifelse(keep, s, paste0(tolower(substring(s, 1, 1)), substring(s, 2)))
}
relabel <- function(f, l) {
  out <- l
  hit <- f %in% names(RELAB)
  out[hit] <- unname(RELAB[f[hit]])
  out
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
# spread direct labels apart vertically so none overlap (no leader lines)
spread <- function(y, gap) {
  o <- order(y)
  ys <- y[o]
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  out <- y
  out[o] <- ys
  out
}
need <- file.path(RS, c("interpP_spline_terms.csv",
                        "interpP_spline_curves.csv",
                        "interpP_shap_importance.csv",
                        "interpP_shap_domains.csv"))
if (any(!file.exists(need))) {
  stop("run stage 20p first: missing ", basename(need[!file.exists(need)][1]))
}
tm <- fread(need[1])
cv <- fread(need[2])
im <- fread(need[3])
dm <- fread(need[4])
tm[, label := relabel(feature, label)]
cv[, label := relabel(feature, label)]
im[, label := relabel(feature, label)]
cat(sprintf("font family: %s\n", fam))

# ---------- twins: identical Charlson and Elixhauser indicators ----------
cat("\n=== TWIN INDICATORS ===\n")
hx <- tm[startsWith(feature, "hx_") & type != "spline"]
sig <- hx[order(cohort), .(s = paste(cohort, signif(term_sd, 10),
                                     signif(edf, 8), collapse = "|"),
                           label = label[1]), by = feature]
grp <- sig[, .(features = list(feature), labels = list(label), n = .N),
           by = s][n == 2]
TW <- rbindlist(lapply(seq_len(nrow(grp)), function(i) {
  f <- grp$features[[i]]
  b <- unique(sub(" \\((Charlson|Elixhauser)\\)$", "", grp$labels[[i]]))
  b <- c(b[1], lc1(b[-1]))
  data.table(feature_1 = f[1], feature_2 = f[2],
             label = sprintf("%s (Charlson, Elixhauser)",
                             paste(b, collapse = " or ")))
}))
if (nrow(TW)) {
  print(TW, row.names = FALSE)
} else {
  TW <- data.table(feature_1 = character(0), feature_2 = character(0),
                   label = character(0))
  cat("  none found\n")
}
fwrite(TW, file.path(RS, "interpP_fig_twins.csv"))
gkey <- function(f) {
  k <- f
  for (i in seq_len(nrow(TW))) {
    k[f %in% c(TW$feature_1[i], TW$feature_2[i])] <- TW$feature_1[i]
  }
  k
}
glab <- function(k, l) {
  out <- l
  m <- match(k, TW$feature_1)
  out[!is.na(m)] <- TW$label[m[!is.na(m)]]
  out
}
# SurvSHAP(t) importance with twins merged, from patient-level values
merged_imp <- function(mn) {
  f <- file.path(OD, sprintf("shap_P_EHRECG_%s.csv", safe(mn)))
  x <- im[cohort == CO & model == mn]
  if (!nrow(TW) || !file.exists(f)) {
    if (nrow(TW)) cat("  ", basename(f), "not found: twins not merged\n")
    return(x[, .(feature, label, domain, horizon, mean_abs_shap, share)])
  }
  L <- fread(f, select = c("hadm_id", "feature", "horizon", "shap_risk"))
  L[, key := gkey(feature)]
  S <- L[, .(v = sum(shap_risk)), by = .(hadm_id, key, horizon)]
  M <- S[, .(mean_abs_shap = mean(abs(v))), by = .(key, horizon)]
  M[, share := mean_abs_shap / sum(mean_abs_shap), by = horizon]
  lb <- unique(x[, .(feature, label, domain)])
  M[, label := glab(key, lb$label[match(key, lb$feature)])]
  M[, domain := lb$domain[match(key, lb$feature)]]
  setnames(M, "key", "feature")
  M[, .(feature, label, domain, horizon, mean_abs_shap, share)]
}
twin_note <- if (nrow(TW)) {
  paste(" Identical Charlson and Elixhauser indicators are shown as one",
        "feature.")
} else {
  ""
}

# ---------- partial effects ----------
cat("\n=== PARTIAL EFFECTS ===\n")
tb <- tm[cohort == CO][order(-term_sd)]
for (ft in tb$feature[seq_len(min(NTOP, nrow(tb)))]) {
  d <- cv[cohort == CO & feature == ft]
  r <- tb[feature == ft]
  lb <- r$label
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
    labs(title = ttl(sprintf("effect of %s on the hazard of death", lb)),
         x = xl, y = "Hazard ratio relative to the median (log scale)") +
    th
  cap <- if (d$type[1] == "spline") {
    paste0(
      "Figure X: figure showing the hazard ratio for death within one ",
      "year of admission across values of ", lc1(lb), " at presentation, ",
      "relative to its median. Estimates come from the penalised spline ",
      "Cox model with a shaded 95% interval, and ticks mark the 10th ",
      "percentile, median and 90th percentile.")
  } else {
    paste0(
      "Figure X: figure showing the hazard ratio for death within one ",
      "year of admission by ", lc1(lb), " at presentation, relative to ",
      "its median value. Estimates come from the ridge penalised linear ",
      "term of the spline Cox model, with 95% intervals.")
  }
  save_fig(p, sprintf("pe_P_%s.png", safe(ft)), cap)
}

# ---------- SurvSHAP(t) share over time, by feature ----------
cat("\n=== SurvSHAP(t) SHARE OVER TIME ===\n")
line_fig <- function(d, key, lab_col, file, title, cap, h = 4.5) {
  hmax <- max(d$horizon)
  lv <- d[horizon == hmax][order(-share)][[key]]
  d[, grp := factor(get(key), levels = lv)]
  lab <- d[horizon == hmax]
  gap <- 0.17 / (h - 1.1) * diff(range(c(0, 100 * d$share)))
  lab[, ylab := spread(100 * share, gap)]
  p <- ggplot(d, aes(x = horizon, y = 100 * share, colour = grp,
                     shape = grp, linetype = grp)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    geom_text(data = lab, aes(y = ylab, label = get(lab_col)), hjust = 0,
              nudge_x = 0.06, size = 3.5, family = fam,
              show.legend = FALSE) +
    scale_x_log10(breaks = sort(unique(d$horizon)),
                  expand = expansion(mult = c(0.03, 0.62))) +
    scale_colour_manual(values = OI) + scale_shape_manual(values = SHAPES) +
    scale_linetype_manual(values = LTY) +
    scale_y_continuous(breaks = function(l) pretty(l, n = 6)) +
    expand_limits(y = 0) +
    labs(title = title, x = "Days from admission (log scale)",
         y = "Share of importance (%)") + th
  save_fig(p, file, cap, w = 7.5, h = h)
}
MI <- list()
for (mn in unique(im$model)) {
  x <- merged_imp(mn)
  MI[[mn]] <- x
  if (!nrow(x)) next
  ml <- if (mn %in% names(MODEL_LONG)) MODEL_LONG[[mn]] else mn
  hm <- max(x$horizon)
  top <- x[horizon == hm][order(-share)]$feature[seq_len(NTOP)]
  line_fig(
    x[feature %in% top], "feature", "label",
    sprintf("shap_share_P_%s.png", safe(mn)),
    ttl(sprintf("share of EHR feature importance over time in the %s",
                ml)),
    paste0(
      "Figure X: figure showing each feature's share of mean absolute ",
      "SurvSHAP(t) value at 10, 30, 90, 180 and 365 days after admission ",
      "in the ", ml, ", for the eight features with the largest share at ",
      "365 days. Values come from 300 admissions explained against 100 ",
      "background admissions.", twin_note), h = 5.4)
}

# ---------- SurvSHAP(t) share over time, by domain ----------
cat("\n=== SurvSHAP(t) SHARE BY DOMAIN ===\n")
for (mn in unique(dm$model)) {
  x <- dm[cohort == CO & model == mn]
  if (!nrow(x)) next
  ml <- if (mn %in% names(MODEL_LONG)) MODEL_LONG[[mn]] else mn
  x[, dlab := sprintf("%s (%d)", domain, n_features)]
  line_fig(
    x, "domain", "dlab", sprintf("shap_domain_P_%s.png", safe(mn)),
    ttl(sprintf("share of importance by EHR domain over time in the %s",
                ml)),
    paste0(
      "Figure X: figure showing the share of SurvSHAP(t) importance held by ",
      "each EHR domain at each horizon in the ", ml, ". A domain's ",
      "importance is the mean absolute sum of its features' values, and ",
      "the number of features in each domain is shown in brackets."))
}

# ---------- agreement ----------
cat("\n=== IMPORTANCE AGREEMENT ===\n")
s0 <- tm[cohort == CO]
s0[, key := gkey(feature)]
sp <- s0[, .(v = sum(term_sd), label = label[1]), by = key]
sp[, label := glab(key, label)]
sp[, v := v / max(v)]
setnames(sp, "key", "feature")
sp[, model := "Spline Cox"]
fo <- rbindlist(lapply(names(MI), function(mn) {
  x <- MI[[mn]][horizon == max(horizon)]
  x[, .(feature, label, v = mean_abs_shap / max(mean_abs_shap),
        model = if (mn %in% names(MODEL_SHORT)) MODEL_SHORT[[mn]] else mn)]
}))
hmax <- max(im$horizon)
d <- rbind(sp[, .(feature, label, v, model)], fo)
keep <- unique(d[, .SD[order(-v)][seq_len(min(NAGR, .N))], by = model]$
                 feature)
d <- d[feature %in% keep]
ordf <- sp[feature %in% keep][order(v)]$label
d[, label := factor(label, levels = ordf)]
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
  labs(title = ttl(paste("relative importance of the leading EHR features",
                         "in three models")),
       x = "Relative importance (each model's largest = 1)", y = NULL) + th
save_fig(p, "importance_agreement_P.png", paste0(
  "Figure X: figure showing the importance of each EHR feature that ranks ",
  "in the top ten of any model, scaled to each model's largest value. ",
  "The spline Cox model is measured by the spread of each term, and the ",
  "random survival forest (RSF) and oblique RSF by mean absolute ",
  sprintf("SurvSHAP(t) value at %d days.", hmax), twin_note),
  w = 7.5, h = 2 + 0.24 * length(ordf))

# ---------- CTPA blocks at diagnosis (stage 16p) ----------
cat("\n=== CTPA BLOCKS ===\n")
fc <- file.path(RS, "ctpa_blocks_P_time.csv")
if (file.exists(fc)) {
  cb <- fread(fc)[model == "alone" & quantity %in%
                    c("HR days 0-30 (per SD)", "HR days 31-365 (per SD)")]
  ordb <- c("All CTPA", "PE descriptors", "PE: right-heart strain signs",
            "PE: clot location and extent", "Incidental findings",
            "Comorbidity", "Support devices", "Report metadata")
  ordb <- intersect(ordb, cb$block)
  nf <- cb[, .(n = n_features[1]), by = block]
  cb[, blab := sprintf("%s (%d)", block, nf$n[match(block, nf$block)])]
  lv <- rev(sprintf("%s (%d)", ordb, nf$n[match(ordb, nf$block)]))
  cb[, blab := factor(blab, levels = lv)]
  cb[, per := ifelse(grepl("0-30", quantity), "Days 0 to 30",
                     "Days 31 to 365")]
  cb[, yy := as.numeric(blab) + ifelse(per == "Days 0 to 30", 0.15, -0.15)]
  pc <- c("Days 0 to 30" = OI[5], "Days 31 to 365" = OI[2])
  ps <- c("Days 0 to 30" = 16, "Days 31 to 365" = 17)
  xr <- range(c(cb$lo, cb$hi), na.rm = TRUE)
  ky <- length(lv) + 1
  kx <- exp(log(xr[1]) + c(0.00, 0.30) * diff(log(xr)))
  p <- ggplot(cb, aes(x = estimate, y = yy, colour = per, shape = per)) +
    annotate("segment", x = 1, xend = 1, y = 0.5, yend = length(lv) + 0.5,
             linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.12,
                  orientation = "y") +
    geom_point(size = 2.4) +
    annotate("point", x = kx, y = ky, colour = pc, shape = ps, size = 2.4) +
    annotate("text", x = kx * exp(0.03 * diff(log(xr))), y = ky,
             label = names(pc), hjust = 0, size = 3.4, family = fam,
             colour = pc) +
    scale_colour_manual(values = pc) + scale_shape_manual(values = ps) +
    scale_x_log10(breaks = c(0.25, 0.33, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2,
                             3, 4, 6)) +
    scale_y_continuous(breaks = seq_along(lv), labels = lv,
                       expand = expansion(add = c(0.5, 0.9))) +
    coord_cartesian(clip = "off") +
    labs(title = ttl(paste("hazard ratio per SD of each CTPA block in early",
                           "and late follow-up after admission")),
         x = "Hazard ratio per SD of the block score (log scale)",
         y = NULL) + th
  save_fig(p, "ctpa_blocks_P.png", paste0(
    "Figure X: figure showing the hazard ratio per standard deviation of ",
    "each CTPA block's ridge Cox score for deaths in days 0 to 30 and 31 ",
    "to 365 after admission, with 95% patient-level bootstrap intervals. ",
    "The number of features in each block is shown in brackets."),
    w = 7.5, h = 4.8)
} else {
  cat("  ", fc, "not found: run stage 16p first; figure skipped\n")
}
caps <- rbindlist(CAP)
fwrite(caps, file.path(RS, "interpP_fig_captions.csv"))
cat("\nOUTPUT FILES\n  ", FD, " (pe_P_*.png, shap_share_P_*.png,",
    "shap_domain_P_*.png, importance_agreement_P.png, ctpa_blocks_P.png)\n")
cat("  ", file.path(RS, "interpP_fig_captions.csv"), "\n")
cat("  ", file.path(RS, "interpP_fig_twins.csv"), "\n")
cat("  ", file.path(RS, "interpP_fig_log.txt"), " (terminal log)\n")
cat("\nOutputs: interpP_fig_captions.csv, interpP_fig_twins.csv\n")