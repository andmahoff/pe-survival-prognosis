# 00_r_packages.R
# Installs the R packages the survival analysis needs, if missing.
#
# Installs only packages (and dependencies) that are not already
# installed, as pre-built Windows binaries, so nothing is compiled and no
# installed package is upgraded or reinstalled. Required: survival,
# riskRegression, prodlim, cmprsk, future.apply and data.table. The script
# sets R_LIBCURL_SSL_REVOKE_BEST_EFFORT, which relaxes Windows
# certificate-revocation checking (a common cause of "SSL connect error"),
# and tries four CRAN mirrors in turn until nothing is missing. It then
# loads each package and prints its version. It exits with status 1 if a
# package still fails to load after every mirror has been tried.
#
# Inputs:  none
# Outputs: none (console output only)
# Usage:   Rscript scripts/00_r_packages.R

Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT = "TRUE")

need <- c("survival", "riskRegression", "prodlim", "cmprsk",
          "future.apply", "data.table")
mirrors <- c("https://cloud.r-project.org",
             "https://www.stats.bris.ac.uk/R",
             "https://cran.ma.imperial.ac.uk",
             "https://packagemanager.posit.co/cran/latest")

cat(R.version.string, "\n")
cat("library folders:\n")
cat(paste0("  ", .libPaths()), sep = "\n")

missing_from <- function(db) {
  deps <- unique(unlist(tools::package_dependencies(
    need, db = db, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE)))
  todo <- setdiff(c(need, deps), rownames(installed.packages()))
  intersect(todo, rownames(db))
}

for (repo in mirrors) {
  cat("\nmirror:", repo, "\n")
  db <- tryCatch(suppressWarnings(
    available.packages(repos = repo, type = "binary")),
    error = function(e) NULL)
  if (is.null(db) || nrow(db) == 0) {
    cat("  could not read the package index, trying next\n")
    next
  }
  todo <- missing_from(db)
  if (length(todo) == 0) {
    cat("  nothing missing\n")
    break
  }
  cat("  missing, to install:", paste(todo, collapse = ", "), "\n")
  tryCatch(
    install.packages(todo, repos = repo, type = "binary",
                     dependencies = FALSE),
    error = function(e) cat("  install error:",
                            conditionMessage(e), "\n"))
  still <- setdiff(todo, rownames(installed.packages()))
  if (length(still) == 0) {
    cat("  all installed from this mirror\n")
    break
  }
  cat("  still missing:", paste(still, collapse = ", "), "\n")
}

cat("\nload check:\n")
ok <- sapply(need, function(p) {
  r <- tryCatch({
    suppressPackageStartupMessages(
      library(p, character.only = TRUE))
    "loaded"
  }, error = function(e) conditionMessage(e))
  cat(sprintf("  %-15s %s\n", p, r))
  identical(r, "loaded")
})
if (!all(ok)) {
  cat("\nStill failing on every mirror. The usual cause is an\n")
  cat("antivirus that scans secure web traffic, or a VPN.\n")
  cat("Send me this log.\n")
  quit(status = 1)
}
cat("\nversions:\n")
for (p in need) {
  cat(sprintf("  %-15s %s\n", p,
              as.character(packageVersion(p))))
}
cat("\nALL LOADED. Now re-run 05b_gates.R\n")