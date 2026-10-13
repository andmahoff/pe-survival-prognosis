# 00b_check_packages.R
# Installs missing dependencies of rstpm2, randomForestSRC and glmnet.
#
# rstpm2 is used by 14a_rp_curve.R, randomForestSRC by
# 06b_level1_forests.R and glmnet by 16b_transfer_s6.R. The script lists
# the packages and dependencies that are not installed, installs them as
# Windows binaries from CRAN, and never upgrades a package that is already
# installed. It then checks that each of the three packages loads and
# prints its version.
#
# Inputs:  none
# Outputs: none (console output only)
# Usage:   Rscript scripts/00b_check_packages.R

repos <- "https://cloud.r-project.org"
want <- c("rstpm2", "randomForestSRC", "glmnet")
ap <- available.packages(repos = repos)
deps <- unique(unlist(tools::package_dependencies(
  want, db = ap, which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE)))
have <- rownames(installed.packages())
miss <- setdiff(c(want, deps), have)
cat("R", as.character(getRversion()), "\n")
cat("missing packages:",
    if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")
if (length(miss)) {
  install.packages(miss, repos = repos, type = "binary")
}
still <- setdiff(c(want, deps), rownames(installed.packages()))
if (length(still)) {
  cat("STILL MISSING after install:", paste(still, collapse = ", "), "\n")
}
for (p in want) {
  ok <- tryCatch({
    suppressWarnings(suppressPackageStartupMessages(
      library(p, character.only = TRUE)))
    TRUE
  }, error = function(e) {
    cat("  ", p, "FAILED TO LOAD:", conditionMessage(e), "\n")
    FALSE
  })
  if (ok) {
    cat(sprintf("  %-16s %-10s loads OK\n", p,
                as.character(packageVersion(p))))
  }
}