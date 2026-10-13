# 18z_diag_preds.R
# Diagnostic for the fusion architecture predictions written from 18a.
#
# For every file, model, set and horizon in data/processed/fusion_arch it
# counts missing risks, risks at exactly 0 or 1 and risks within 1e-8 of 0
# or 1, and reports the smallest and largest risk and the number of patients
# affected. Output is aggregate only. It runs in seconds, so it has no
# timing probe and no parallelism.
#
# Inputs:  data/processed/fusion_arch/arch_preds_*.csv (patient level,
#          not in the repo)
# Outputs: results/arch_pred_diag.csv, results/arch_pred_diag_log.txt
# Usage:   Rscript scripts/18z_diag_preds.R

suppressPackageStartupMessages(library(data.table))
options(width = 220)
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
RS <- file.path(P2, "results")
OD <- file.path(PR, "fusion_arch")
files <- list.files(OD, pattern = "^arch_preds_.*\\.csv$", full.names = TRUE)
cat(sprintf("files found: %d\n", length(files)))
D <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  s <- x[, .(
    n = .N,
    missing = sum(is.na(risk)),
    patients_missing = uniqueN(hadm_id[is.na(risk)]),
    exact_0 = sum(risk == 0, na.rm = TRUE),
    exact_1 = sum(risk == 1, na.rm = TRUE),
    below_1e8 = sum(risk < 1e-8, na.rm = TRUE),
    above_1m1e8 = sum(risk > 1 - 1e-8, na.rm = TRUE),
    min_risk = suppressWarnings(min(risk, na.rm = TRUE)),
    max_risk = suppressWarnings(max(risk, na.rm = TRUE))),
    by = .(model, set, horizon)]
  s[, file := basename(f)]
  s
}))
setcolorder(D, "file")
D[, problem := missing > 0 | below_1e8 > 0 | above_1m1e8 > 0]
fwrite(D, file.path(RS, "arch_pred_diag.csv"))
cat("\n=== MODELS WITH MISSING OR SATURATED RISKS ===\n")
print(D[problem == TRUE, .(file, model, set, horizon, missing,
                           patients_missing, exact_0, exact_1, below_1e8,
                           above_1m1e8)], row.names = FALSE, nrows = 500)
cat(sprintf("\n%d of %d model x set x horizon cells affected\n",
            sum(D$problem), nrow(D)))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(RS, "arch_pred_diag.csv"), "\n")
cat("  ", file.path(RS, "arch_pred_diag_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: arch_pred_diag.csv\n")