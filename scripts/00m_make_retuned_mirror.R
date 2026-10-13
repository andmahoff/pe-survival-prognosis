# 00m_make_retuned_mirror.R
# Working copy of the project with retuned level-1 scores (A17, A20;
# PREREGISTRATION_ADDENDUM_2.md).
#
# The level-1 ridge-Cox file of each frame and cohort is replaced by a
# retuned one in a copy of the project, so the registered scripts can be
# rerun unchanged on retuned scores without touching any registered
# output. The script:
#   1. copies the project to SurvivalProject_<variant>, leaving out .git,
#      venv_surv, checkpoint folders and the large output folders that the
#      registered scripts do not read (fusion_arch, super_learner,
#      transfer_ext, ctpa_time, spline, analysis_b);
#   2. replaces each level1/oof_linear_<code>.csv (codes A_ehrecg, A_3mod,
#      B_ehrecg, B_3mod) with oof_linear_<variant>_<code>.csv;
#   3. in the copy's .R and .py scripts, replaces the project path
#      SurvivalProject by SurvivalProject_<variant>, byte for byte, so
#      nothing else in a script changes;
#   4. writes MIRROR_README.txt and results/.mirror_stamp. The comparison
#      script 13r_compare_retuned.R uses only results written after the
#      stamp.
# Variants: retuned (default, A17) takes the scores from
# 06e_level1_linear_retuned.R (ridge-Cox retuned); retuned2 (A20) takes
# them from 06f_level1_binary_retuned.R (ridge-Cox and binary twins
# retuned). The script stops if a retuned file is missing, or if the copy
# already exists and --rebuild is not given.
#
# Inputs:  the whole project; data/processed/level1/
#          oof_linear_<variant>_<code>.csv
# Outputs: SurvivalProject_<variant>/ (sibling folder, the working copy)
# Usage:   Rscript scripts/00m_make_retuned_mirror.R
#              [--variant=retuned|retuned2] [--rebuild]

args <- commandArgs(trailingOnly = TRUE)
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
va <- grep("^--variant=", args, value = TRUE)
VAR <- if (length(va)) sub("^--variant=", "", va) else "retuned"
if (!(VAR %in% c("retuned", "retuned2"))) {
  cat("STOP: --variant must be retuned or retuned2\n")
  quit(status = 1)
}
DST <- paste0(P2, "_", VAR)
EXCL <- c(".git", "venv_surv", "fusion_arch", "super_learner",
          "transfer_ext", "ctpa_time", "spline", "analysis_b")
CODES <- c("A_ehrecg", "A_3mod", "B_ehrecg", "B_3mod")
t0 <- Sys.time()
fail <- function(msg) {
  cat("\nSTOP:", msg, "\n")
  quit(status = 1)
}
cat("=== CHECK ===\n")
cat(sprintf("  variant %s: working copy %s\n", VAR, DST))
rt <- file.path(P2, "data", "processed", "level1",
                sprintf("oof_linear_%s_%s.csv", VAR, CODES))
if (!all(file.exists(rt))) {
  fail(paste("retuned level-1 file missing:",
             paste(basename(rt[!file.exists(rt)]), collapse = ", ")))
}
cat("  retuned files found:", paste(basename(rt), collapse = ", "), "\n")
if (dir.exists(DST)) {
  if (!("--rebuild" %in% args)) {
    fail(paste(DST, "already exists; rerun with --rebuild to replace it"))
  }
  cat("  removing the existing copy (--rebuild)\n")
  unlink(DST, recursive = TRUE, force = TRUE)
}

cat("\n=== COPY ===\n")
rel <- list.files(P2, recursive = TRUE, all.files = TRUE, no.. = TRUE)
parts <- strsplit(rel, "/", fixed = TRUE)
skip <- sapply(parts, function(p) {
  d <- p[-length(p)]
  any(d %in% EXCL) || any(grepl("_ckpt$", d))
})
rel <- rel[!skip]
sz <- sum(file.info(file.path(P2, rel))$size, na.rm = TRUE)
cat(sprintf("  %d files, %.0f MB to copy\n", length(rel), sz / 2^20))
for (d in unique(dirname(rel))) {
  dir.create(file.path(DST, d), recursive = TRUE, showWarnings = FALSE)
}
ok <- file.copy(file.path(P2, rel), file.path(DST, rel),
                overwrite = TRUE, copy.date = TRUE)
if (!all(ok)) fail(sprintf("%d files could not be copied", sum(!ok)))
cat("  copied\n")

cat("\n=== REPLACE THE LEVEL-1 RIDGE-COX FILES ===\n")
for (cd in CODES) {
  src <- file.path(P2, "data", "processed", "level1",
                   sprintf("oof_linear_%s_%s.csv", VAR, cd))
  dst <- file.path(DST, "data", "processed", "level1",
                   sprintf("oof_linear_%s.csv", cd))
  if (!file.copy(src, dst, overwrite = TRUE)) {
    fail(paste("could not write", dst))
  }
  cat(sprintf("  %s  <-  %s\n", basename(dst), basename(src)))
}

cat("\n=== POINT THE COPY'S SCRIPTS AT THE COPY ===\n")
sc <- list.files(file.path(DST, "scripts"), pattern = "\\.(R|py)$",
                 full.names = TRUE)
nch <- 0L
for (f in sc) {
  raw <- readBin(f, "raw", file.info(f)$size)
  txt <- rawToChar(raw)
  new <- gsub("SurvivalProject(?!_retuned)", paste0("SurvivalProject_", VAR),
              txt, perl = TRUE, useBytes = TRUE)
  if (!identical(new, txt)) {
    writeBin(charToRaw(new), f)
    nch <- nch + 1L
  }
}
cat(sprintf("  %d of %d scripts rewritten\n", nch, length(sc)))
dir.create(file.path(DST, "results"), showWarnings = FALSE)
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           file.path(DST, "results", ".mirror_stamp"))
writeLines(c(
  sprintf("WORKING COPY (%s) for PREREGISTRATION_ADDENDUM_2.md %s.", VAR,
          if (VAR == "retuned") "A17" else "A20"),
  "Do not edit.",
  sprintf("Built %s from %s.", format(Sys.time(), "%Y-%m-%d %H:%M"), P2),
  sprintf("level1/oof_linear_*.csv here are copies of oof_linear_%s_*.csv.",
          VAR),
  "Scripts here point at this folder; registered outputs are untouched."),
  file.path(DST, "MIRROR_README.txt"))
cat(sprintf("\n  done in %.1f min. Working copy: %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), DST))