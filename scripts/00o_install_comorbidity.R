# 00o_install_comorbidity.R
# Installs the comorbidity package, which 06j_level1_ehrplus.R needs (A30).
#
# Nothing already installed is upgraded: install.packages() adds only the
# missing dependencies (checkmate, data.table, stringi). The script tries
# three CRAN mirrors with R's default downloader (libcurl), then with the
# Windows downloader (wininet), and stops at the first combination that
# works. On success it runs a self-test that scores a three-row table of
# ICD-10 codes with the Charlson index; the expected scores are 8 and 2.
#
# Inputs:  none
# Outputs: none (console output only)
# Usage:   Rscript scripts/00o_install_comorbidity.R

deps <- c("checkmate", "data.table", "stringi")
have <- rownames(installed.packages())
cat("dependencies already installed:",
    paste(intersect(deps, have), collapse = ", "), "\n")
cat("dependencies missing (will be added):",
    if (length(setdiff(deps, have))) paste(setdiff(deps, have),
                                           collapse = ", ") else "none",
    "\n")
mirrors <- c("https://cran.ma.imperial.ac.uk",
             "https://www.stats.bris.ac.uk/R",
             "https://cloud.r-project.org")
ok <- requireNamespace("comorbidity", quietly = TRUE)
if (ok) cat("\ncomorbidity is already installed.\n")
for (m in c("libcurl", "wininet")) {
  for (r in mirrors) {
    if (ok) break
    cat(sprintf("\n--- trying %s with the %s downloader ---\n", r, m))
    try(install.packages("comorbidity", repos = r, method = m,
                         dependencies = NA))
    ok <- requireNamespace("comorbidity", quietly = TRUE)
  }
}
if (ok) {
  cat(sprintf("\nSUCCESS: comorbidity %s installed and loads.\n",
              as.character(packageVersion("comorbidity"))))
  suppressPackageStartupMessages(library(comorbidity))
  x <- data.frame(id = c(1, 1, 2), code = c("C787", "I509", "E119"))
  s <- score(comorbidity(x, id = "id", code = "code",
                         map = "charlson_icd10_quan", assign0 = FALSE),
             weights = "charlson", assign0 = TRUE)
  cat("self-test Charlson scores (expect 8 and 2):",
      paste(as.numeric(s), collapse = " "), "\n")
} else {
  cat("\nFAILED on every mirror and downloader.",
      "Use the manual route in the instructions.\n")
}