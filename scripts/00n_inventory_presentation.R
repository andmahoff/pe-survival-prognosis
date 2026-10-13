# 00n_inventory_presentation.R
# Inventory of the inputs for the at-diagnosis model (A31).
#
# Summarises input files (column names, row counts, coverage, missing
# values) and changes nothing. It answers four questions on which the
# design of the at-diagnosis model depends:
#   1. Can a 365-day death outcome be timed from admission for the Frame A
#      cohort (coverage of dod_days in mimic_labels_harmonised.csv)?
#   2. Which of the 28 EHR columns are condition flags, to be rebuilt from
#      previous admissions only?
#   3. Is the ECG modality one row per admission (an aggregate) or one row
#      per ECG, and does a per-ECG feature file with times exist anywhere?
#      A first-ECG-only modality would need one.
#   4. Does anything record which CTPA report window ctpa46.csv was built
#      from?
#
# Inputs:  data/inputs/{mimic_labels_harmonised,ehr28,ehr_extra,ecg_opt,
#          ecg_timing,ctpa46}.csv, data/processed/frame_A.csv,
#          data/processed/feature_whitelist.csv; file names (not
#          contents) under PROJECTS_DIR
# Outputs: results/inventory.csv, results/inventory_log.txt
# Usage:   Rscript scripts/00n_inventory_presentation.R

suppressPackageStartupMessages(library(data.table))
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
ROOT <- dirname(P2)
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OUT <- list()
put <- function(section, item, value) {
  OUT[[length(OUT) + 1]] <<- data.table(section = section, item = item,
                                        value = as.character(value))
  cat(sprintf("  [%s] %s: %s\n", section, item, value))
}
hdr <- function(f) {
  tryCatch(names(fread(f, nrows = 1, sep = ",")),
           error = function(e) character(0))
}

cat("=== 1. FRAME A AND DEATH TIMING ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
put("frame A", "rows / columns", sprintf("%d / %d", nrow(fa), ncol(fa)))
put("frame A", "column names", paste(names(fa), collapse = " "))
if ("time_death" %in% names(fa)) {
  put("frame A", "time_death max / deaths",
      sprintf("%.1f / %d", max(fa$time_death), sum(fa$event_death == 1)))
}
lab <- file.path(IN, "mimic_labels_harmonised.csv")
if (file.exists(lab)) {
  L <- fread(lab)
  put("labels", "column names", paste(names(L), collapse = " "))
  L <- L[hadm_id %in% fa$hadm_id]
  for (cc in intersect(c("dod_days", "inhosp_days"), names(L))) {
    v <- suppressWarnings(as.numeric(L[[cc]]))
    put("labels", sprintf("%s: non-missing / max / <= 365", cc),
        sprintf("%d / %.0f / %d", sum(!is.na(v)), max(v, na.rm = TRUE),
                sum(v <= 365, na.rm = TRUE)))
  }
}

cat("\n=== 2. THE 28 EHR COLUMNS ===\n")
e28 <- fread(file.path(IN, "ehr28.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
c28 <- intersect(wl[modality == "ehr", column], names(e28))
for (cc in c28) {
  v <- e28[[cc]]
  u <- length(unique(v[!is.na(v)]))
  put("ehr28", cc, sprintf("%s distinct values, %.1f%% missing",
                           if (u <= 2) "binary," else paste0(u),
                           100 * mean(is.na(v))))
}
ex <- file.path(IN, "ehr_extra.csv")
if (file.exists(ex)) {
  h <- hdr(ex)
  put("ehr_extra", "non-lab columns",
      paste(h[!startsWith(h, "lab_")], collapse = " "))
}

cat("\n=== 3. ECG MODALITY ===\n")
for (f in c("ecg_opt.csv", "ecg_timing.csv")) {
  p <- file.path(IN, f)
  if (!file.exists(p)) next
  E <- fread(p)
  put(f, "rows / distinct hadm_id",
      sprintf("%d / %d", nrow(E), uniqueN(E$hadm_id)))
  put(f, "first 15 column names", paste(head(names(E), 15),
                                        collapse = " "))
}

cat("\n=== 4. CTPA MODALITY ===\n")
ct <- file.path(IN, "ctpa46.csv")
if (file.exists(ct)) {
  C <- fread(ct)
  put("ctpa46", "rows / distinct hadm_id",
      sprintf("%d / %d", nrow(C), uniqueN(C$hadm_id)))
  put("ctpa46", "non-feature columns",
      paste(setdiff(names(C), wl[modality == "ctpa", column]),
            collapse = " "))
}

cat("\n=== 5. FILES ELSEWHERE THAT MAY HOLD PER-ECG OR WINDOW INFO ===\n")
pat <- paste0("(bench_feats|record_list|machine_measurements|ecg.*(feat|",
              "logit|time)|ctpa.*(window|present|note)|presentation)")
ff <- list.files(ROOT, pattern = "\\.(csv|parquet)$", recursive = TRUE,
                 full.names = TRUE, ignore.case = TRUE)
ff <- ff[grepl(pat, basename(ff), ignore.case = TRUE)]
ff <- ff[!grepl("/(venv|venv_surv|\\.git)/", ff)]
put("search", "matching files", length(ff))
for (f in head(ff, 60)) {
  put("search", sub(paste0("^", ROOT, "/"), "", f),
      sprintf("%.1f MB; columns: %s", file.size(f) / 1e6,
              paste(head(hdr(f), 8), collapse = " ")))
}
inv <- rbindlist(OUT)
fwrite(inv, file.path(RS, "inventory.csv"))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(RS, "inventory.csv"), "\n")
cat("  ", file.path(RS, "inventory_log.txt"), " (terminal log)\n")
cat("\nOutputs: inventory.csv\n")