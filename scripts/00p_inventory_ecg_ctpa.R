# 00p_inventory_ecg_ctpa.R
# Inventory for the first-ECG and presentation-CTPA modalities (A32).
#
# Reads only file headers, row counts and identifier columns (never report
# text) and changes nothing. It answers three questions:
#   1. Which per-ECG file holds each column of the ECG modality
#      (ecg_opt.csv)? Columns are matched exactly and after removing the
#      suffixes (_mean, _max, _min and others) added when the ECGs were
#      averaged per admission.
#   2. Do the per-ECG files share the record key with the ECG index that
#      carries each ECG's time?
#   3. Which columns of ctpa_notes_index.csv can time the CTPA report?
#
# Inputs:  data/inputs/ecg_opt.csv, data/inputs/ctpa46.csv;
#          side project: data/processed/exp0_logits_record.csv,
#          data/processed/ecg_derived_*record*.csv,
#          data/raw/ecg/ecg_record_index*.csv;
#          dissertation: figwork/data/ctpa_notes_index.csv
# Outputs: results/inventory_ecg_ctpa.csv,
#          results/inventory_ecg_ctpa_log.txt
# Usage:   Rscript scripts/00p_inventory_ecg_ctpa.R

suppressPackageStartupMessages(library(data.table))
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
ROOT <- dirname(P2)
SP <- file.path(ROOT, "SideProject")
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
nrow_ids <- function(f, col) {
  x <- tryCatch(fread(f, select = col), error = function(e) NULL)
  if (is.null(x)) return(c(NA, NA))
  c(nrow(x), uniqueN(x[[1]]))
}

cat("=== 1. THE ECG MODALITY'S COLUMNS ===\n")
eo <- file.path(P2, "data", "inputs", "ecg_opt.csv")
oc <- setdiff(hdr(eo), c("subject_id", "hadm_id"))
put("ecg_opt", "feature columns", length(oc))
put("ecg_opt", "all column names", paste(oc, collapse = " "))
base <- unique(sub("_(mean|max|min|sd|first|last)$", "", oc))
put("ecg_opt", "names after removing _mean/_max/_min/_sd suffixes",
    length(base))

cat("\n=== 2. PER-ECG FILES ===\n")
rf <- c(file.path(SP, "data", "processed", "exp0_logits_record.csv"),
        Sys.glob(file.path(SP, "data", "processed",
                           "ecg_derived_*record*.csv")))
rf <- rf[file.exists(rf)]
put("per-ECG", "files found", length(rf))
for (f in rf) {
  h <- hdr(f)
  key <- intersect(c("rec", "path", "study_id", "ecg_study_id"), h)
  ni <- if (length(key)) nrow_ids(f, key[1]) else c(NA, NA)
  put("per-ECG", basename(f), sprintf(paste(
    "%s rows, %s distinct %s; %d columns; ecg_opt columns matched",
    "exactly %d, after suffix removal %d"),
    ni[1], ni[2], if (length(key)) key[1] else "(no key)", length(h),
    sum(oc %in% h), sum(base %in% h)))
}
cat("\n=== 3. ECG TIMES ===\n")
ix <- Sys.glob(file.path(SP, "data", "raw", "ecg",
                         "ecg_record_index*.csv"))
for (f in ix) {
  h <- hdr(f)
  ni <- if ("rec" %in% h) nrow_ids(f, "rec") else c(NA, NA)
  put("ECG index", basename(f), sprintf("%s rows, %s distinct rec; columns: %s",
                                        ni[1], ni[2], paste(h, collapse = " ")))
}
lg <- file.path(SP, "data", "processed", "exp0_logits_record.csv")
if (file.exists(lg) && length(ix)) {
  a <- fread(lg, select = "rec")$rec
  b <- fread(ix[1], select = "rec")$rec
  put("ECG index", "per-ECG logits whose rec is in the index",
      sprintf("%.1f%%", 100 * mean(a %in% b)))
}

cat("\n=== 4. CTPA TIMING ===\n")
cn <- file.path(ROOT, "Dissertation", "figwork", "data",
                "ctpa_notes_index.csv")
put("CTPA", "ctpa_notes_index.csv found", file.exists(cn))
if (file.exists(cn)) {
  h <- hdr(cn)
  put("CTPA", "columns (no text read)", paste(h, collapse = " "))
  keep <- intersect(c("hadm_id", "note_id", "charttime", "admittime"), h)
  x <- fread(cn, select = keep)
  put("CTPA", "rows / distinct hadm_id", sprintf("%d / %d", nrow(x),
                                                  uniqueN(x$hadm_id)))
  ct <- unique(fread(file.path(P2, "data", "inputs", "ctpa46.csv"),
                     select = "hadm_id"))
  put("CTPA", "ctpa46 admissions with a note in the index",
      sprintf("%.1f%%", 100 * mean(ct$hadm_id %in% x$hadm_id)))
  if (all(c("charttime", "admittime") %in% names(x))) {
    hh <- as.numeric(difftime(as.POSIXct(x$charttime, tz = "UTC"),
                              as.POSIXct(x$admittime, tz = "UTC"),
                              units = "hours"))
    inw <- x$hadm_id %in% ct$hadm_id
    put("CTPA", "ctpa46 notes: hours from admission (min / median / max)",
        sprintf("%.1f / %.1f / %.1f", min(hh[inw], na.rm = TRUE),
                median(hh[inw], na.rm = TRUE), max(hh[inw], na.rm = TRUE)))
    put("CTPA", "ctpa46 notes inside -48 h to +24 h",
        sprintf("%.1f%%", 100 * mean(hh[inw] >= -48 & hh[inw] <= 24,
                                     na.rm = TRUE)))
  }
}
inv <- rbindlist(OUT)
fwrite(inv, file.path(RS, "inventory_ecg_ctpa.csv"))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(RS, "inventory_ecg_ctpa.csv"), "\n")
cat("  ", file.path(RS, "inventory_ecg_ctpa_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: inventory_ecg_ctpa.csv\n")