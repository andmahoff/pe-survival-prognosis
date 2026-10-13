# 00q_ecg_time_check.R
# Two checks on the first-ECG modality before any model is fitted (A32).
#
# Prints aggregates only (no dates, no patient rows) and changes nothing.
#   1. Time resolution. Do ECG times and admission times carry clock
#      times, or only dates? 04p_build_presentation.R found the first ECG
#      at exactly 0.0 hours from admission for at least three quarters of
#      admissions, which is only possible if one or both are dates.
#   2. Unmatched columns. 04p_build_presentation.R found 166 of the 264
#      columns of ecg_opt.csv in a per-ECG file. Which suffixes do the
#      other 98 carry, and does each have a matched sibling with the same
#      base name?
#
# Inputs:  data/inputs/ecg_opt.csv, data/inputs/index_admission_times.csv,
#          results/presentation_ecg_sources.csv; side project:
#          data/raw/ecg/ecg_record_index.csv
# Outputs: results/ecg_time_check.csv, results/ecg_time_check_log.txt
# Usage:   Rscript scripts/00q_ecg_time_check.R

suppressPackageStartupMessages(library(data.table))
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
SP <- file.path(dirname(P2), "SideProject")
IN <- file.path(P2, "data", "inputs")
RS <- file.path(P2, "results")
OUT <- list()
put <- function(item, value) {
  OUT[[length(OUT) + 1]] <<- data.table(item = item,
                                        value = as.character(value))
  cat(sprintf("  %-60s %s\n", item, value))
}
tparse <- function(x) {
  x <- sub("T", " ", as.character(x), fixed = TRUE)
  as.POSIXct(x, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%OS",
                                           "%Y-%m-%d %H:%M", "%Y-%m-%d"),
             optional = TRUE)
}
shape <- function(x) {
  x <- as.character(x)
  sprintf("%.1f%% contain a clock time; lengths %s",
          100 * mean(grepl(":", x)),
          paste(head(names(sort(table(nchar(x)), decreasing = TRUE)), 2),
                collapse = " / "))
}
cat("=== 1. TIME RESOLUTION ===\n")
ix <- fread(file.path(SP, "data", "raw", "ecg", "ecg_record_index.csv"))
adm <- unique(fread(file.path(IN, "index_admission_times.csv"),
                    select = c("hadm_id", "admittime")), by = "hadm_id")
put("ECG times (ecg_charttime)", shape(ix$ecg_charttime))
put("admission times (admittime)", shape(adm$admittime))
tE <- tparse(ix$ecg_charttime)
put("ECG times at exactly midnight", sprintf(
  "%.1f%%", 100 * mean(format(tE, "%H:%M:%S") == "00:00:00",
                       na.rm = TRUE)))
tA <- tparse(adm$admittime[match(ix$hadm_id, adm$hadm_id)])
h <- as.numeric(difftime(tE, tA, units = "hours"))
put("ECG minus admission (h): share exactly 0",
    sprintf("%.1f%%", 100 * mean(h == 0, na.rm = TRUE)))
put("ECG minus admission (h): share a whole multiple of 24",
    sprintf("%.1f%%", 100 * mean(abs(h / 24 - round(h / 24)) < 1e-9,
                                 na.rm = TRUE)))
put("ECG minus admission (h): 5%, 25%, 50%, 75%, 95%",
    paste(round(quantile(h, c(0.05, 0.25, 0.5, 0.75, 0.95),
                         na.rm = TRUE), 2), collapse = ", "))
put("ECG times identical to the admission time",
    sprintf("%.1f%%", 100 * mean(tE == tA, na.rm = TRUE)))

cat("\n=== 2. UNMATCHED ECG COLUMNS ===\n")
oc <- setdiff(names(fread(file.path(IN, "ecg_opt.csv"), nrows = 1)),
              c("subject_id", "hadm_id"))
mt <- fread(file.path(RS, "presentation_ecg_sources.csv"))$column
um <- setdiff(oc, mt)
suf <- ifelse(grepl("_[a-z]+$", um), sub(".*_([a-z]+)$", "\\1", um), "")
tb <- sort(table(suf), decreasing = TRUE)
put("unmatched columns", length(um))
put("their suffixes (count)", paste(sprintf("%s %d", names(tb), tb),
                                    collapse = "; "))
base_um <- sub("_[a-z]+$", "", um)
base_mt <- sub("_[a-z]+$", "", mt)
put("unmatched columns whose base name is matched under another suffix",
    sprintf("%d of %d", sum(base_um %in% base_mt), length(um)))
put("unmatched columns with no matched sibling (first 10)",
    paste(head(um[!(base_um %in% base_mt)], 10), collapse = " "))
fwrite(rbindlist(OUT), file.path(RS, "ecg_time_check.csv"))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(RS, "ecg_time_check.csv"), "\n")
cat("  ", file.path(RS, "ecg_time_check_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: ecg_time_check.csv\n")