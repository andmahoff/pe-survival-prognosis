# 13r_compare_retuned.R
# Comparison of the results of the registered scripts in a working copy
# with a base set of results (A17, A20).
#
# The working copy is built by 00m_make_retuned_mirror.R. Every results CSV
# written in the copy after it was built (newer than its .mirror_stamp
# file) and present in the base is compared cell by cell. The registered
# scripts are unchanged, so each file has the same rows in the same order.
# Numbers are compared as base, copy and difference. TRUE/FALSE and text
# columns that are identical in both files label the rows; the others
# (decisions, "excludes zero" flags, interval labels) are values, and every
# cell that changed is listed. Files with different columns or row counts
# are listed, not compared, and bootstrap-replicate files (*_boot.csv) are
# skipped.
#
# --variant=retuned (default) compares the copy SurvivalProject_retuned
# (A17) and --variant=retuned2 the copy SurvivalProject_retuned2 (A20).
# --base=registered (default) uses the main project's registered results;
# --base=retuned uses the A17 copy's results, which with
# --variant=retuned2 isolates what retuning the binary twins did. Outputs
# are named {variant}_vs_{base}, for example retuned2_vs_registered. The
# script stops if the variant or base is invalid (they must differ), the
# base results folder is missing or the working copy has not been built.
# Version 4: adds --variant and --base, and names the value columns base
# and copy.
#
# Inputs:  results/*.csv in the main project and in the working copy
# Outputs: results/{variant}_vs_{base}.csv, {variant}_vs_{base}_text.csv,
#          {variant}_vs_{base}_summary.csv, {variant}_vs_{base}_log.txt
#          (terminal log), all in the main project's results folder
# Usage:   Rscript scripts/13r_compare_retuned.R [--variant=retuned2]
#          [--base=retuned]

suppressPackageStartupMessages(library(data.table))
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
args <- commandArgs(trailingOnly = TRUE)
arg <- function(nm, def) {
  a <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", nm, "="), "", a) else def
}
VAR <- arg("variant", "retuned")
BASE <- arg("base", "registered")
if (!(VAR %in% c("retuned", "retuned2")) ||
    !(BASE %in% c("registered", "retuned")) || VAR == BASE) {
  cat("STOP: --variant must be retuned or retuned2; --base registered or",
      "retuned, and different from --variant\n")
  quit(status = 1)
}
DST <- paste0(P2, "_", VAR)
RS0 <- if (BASE == "registered") file.path(P2, "results") else
  file.path(paste0(P2, "_", BASE), "results")
RS1 <- file.path(DST, "results")
OUT <- file.path(P2, "results")
TAG <- sprintf("%s_vs_%s", VAR, BASE)
if (!dir.exists(RS0)) {
  cat("STOP: no results folder for --base:", RS0, "\n")
  quit(status = 1)
}
options(width = 220)
st <- file.path(RS1, ".mirror_stamp")
if (!file.exists(st)) {
  cat("STOP: no working copy (run 00m_make_retuned_mirror.R first)\n")
  quit(status = 1)
}
stamp <- file.info(st)$mtime
f1 <- list.files(RS1, pattern = "\\.csv$")
f1 <- f1[file.info(file.path(RS1, f1))$mtime > stamp]
f1 <- f1[!grepl("_boot\\.csv$", f1) & file.exists(file.path(RS0, f1))]
cat(sprintf("=== %s: %d results files rewritten in the copy and present",
            TAG, length(f1)), "in the base ===\n")
LONG <- list()
TXT <- list()
SUM <- list()
for (f in f1) {
  a <- fread(file.path(RS0, f))
  b <- fread(file.path(RS1, f))
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) {
    SUM[[length(SUM) + 1]] <- data.table(
      file = f, rows = nrow(b), numeric_columns = NA_integer_,
      cells_changed = NA_integer_, max_abs_diff = NA_real_,
      text_cells_changed = NA_integer_, compared = FALSE)
    next
  }
  num <- names(a)[sapply(a, is.numeric)]
  oth <- setdiff(names(a), num)
  same_col <- vapply(oth, function(k) {
    identical(as.character(a[[k]]), as.character(b[[k]]))
  }, logical(1))
  key <- oth[same_col]
  tval <- oth[!same_col]
  kh <- intersect(c("horizon", "cut", "h"), num)
  lab <- if (length(c(key, kh))) {
    do.call(paste, c(lapply(c(key, kh), function(k) a[[k]]), sep = " | "))
  } else {
    as.character(seq_len(nrow(a)))
  }
  vals <- setdiff(num, kh)
  for (v in vals) {
    ra <- as.numeric(a[[v]])
    rb <- as.numeric(b[[v]])
    LONG[[length(LONG) + 1]] <- data.table(
      file = f, row = seq_len(nrow(a)), row_key = lab, column = v,
      base = ra, copy = rb, difference = rb - ra)
  }
  nt <- 0L
  for (v in tval) {
    ch <- which(as.character(a[[v]]) != as.character(b[[v]]) |
                  xor(is.na(a[[v]]), is.na(b[[v]])))
    nt <- nt + length(ch)
    if (length(ch)) {
      TXT[[length(TXT) + 1]] <- data.table(
        file = f, row = ch, row_key = lab[ch], column = v,
        base = as.character(a[[v]][ch]),
        copy = as.character(b[[v]][ch]))
    }
  }
  dd <- as.numeric(unlist(lapply(vals, function(v) {
    as.numeric(b[[v]]) - as.numeric(a[[v]])
  })))
  SUM[[length(SUM) + 1]] <- data.table(
    file = f, rows = nrow(b), numeric_columns = length(vals),
    cells_changed = sum(abs(dd) > 1e-12, na.rm = TRUE),
    max_abs_diff = if (any(is.finite(dd))) max(abs(dd), na.rm = TRUE)
    else NA_real_,
    text_cells_changed = nt, compared = TRUE)
}
long <- rbindlist(LONG)
txt <- rbindlist(TXT)
if (!nrow(txt)) {
  txt <- data.table(file = character(0), row = integer(0),
                    row_key = character(0), column = character(0),
                    base = character(0), copy = character(0))
}
sm <- rbindlist(SUM)
fwrite(long, file.path(OUT, paste0(TAG, ".csv")))
fwrite(txt, file.path(OUT, paste0(TAG, "_text.csv")))
fwrite(sm, file.path(OUT, paste0(TAG, "_summary.csv")))
print(sm, row.names = FALSE)
for (f in intersect(c("primary.csv", "s5_primary.csv"), sm$file)) {
  cat(sprintf("\n=== %s: %s (base) vs %s (copy) ===\n", f, BASE, VAR))
  print(long[file == f & column %in% c("estimate", "lo", "hi"),
             .(row_key, column, base = round(base, 4),
               copy = round(copy, 4))], row.names = FALSE)
}
cat("\n=== TEXT AND TRUE / FALSE CELLS THAT CHANGED ===\n")
print(txt, row.names = FALSE, nrows = 500)
cat("\nOUTPUT FILES\n")
for (f in paste0(TAG, c(".csv", "_text.csv", "_summary.csv"))) {
  cat("  ", file.path(OUT, f), "\n")
}
cat("  ", file.path(OUT, paste0(TAG, "_log.txt")),
    " (terminal log)\n")
cat(sprintf("\nOutputs: %s.csv, %s_text.csv, %s_summary.csv\n", TAG, TAG,
            TAG))