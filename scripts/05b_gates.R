# 05b_gates.R
# Pre-registered gates G3 to G5 with the R evaluator (stage 5b).
#
# Repeats G3 (the time-dependent AUC at 30.5 days equals the binary AUROC),
# G4 (Harrell's C is above 0.5 and matches Python) and G5 (Uno's C equals
# Harrell's C) with riskRegression::Score and survival::concordance, the
# functions every later evaluation uses. Each R value is also checked
# against the Python value from 05_gates.py, within 5e-4 for G3 and 1e-4 for
# G4 and G5. The script stops if survival, riskRegression or prodlim is
# missing.
#
# Inputs:  data/processed/gate_scores.csv and results/gates.csv (both from
#          05_gates.py)
# Outputs: results/gates_r.csv, results/gates_r_log.txt
# Usage:   Rscript scripts/05b_gates.R

need <- c("survival", "riskRegression", "prodlim")
have <- sapply(need, requireNamespace, quietly = TRUE)
if (!all(have)) {
  cat("Missing R packages:", need[!have], "\n")
  cat("Install them with the R block from the setup step.\n")
  quit(status = 1)
}
suppressPackageStartupMessages({
  library(survival)
  library(riskRegression)
  library(prodlim)
})
cat(sprintf("R %s.%s  survival %s  riskRegression %s\n",
            R.version$major, R.version$minor,
            as.character(packageVersion("survival")),
            as.character(packageVersion("riskRegression"))))

P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
sc <- read.csv(file.path(P2, "data", "processed",
                         "gate_scores.csv"))
py <- read.csv(file.path(P2, "results", "gates.csv"))
CUT <- 30.5
TAU_C <- 31

pyv <- function(item) {
  v <- py$value[py$item == item]
  if (length(v) == 0) NA_real_ else as.numeric(v[1])
}

rows <- list()
add <- function(gate, item, r_value, py_value, ok) {
  st <- if (isTRUE(ok)) "PASS" else "FAIL"
  rows[[length(rows) + 1]] <<- data.frame(
    gate = gate, item = item, r_value = r_value,
    py_value = py_value, status = st)
  cat(sprintf("  %-4s %-3s %-40s R %.6f   Python %s\n",
              st, gate, item, r_value,
              ifelse(is.na(py_value), "NA",
                     sprintf("%.6f", py_value))))
}

cat(sprintf("  rows %d, deaths by day 30 %d\n", nrow(sc),
            sum(sc$event_death == 1 & sc$time_death <= CUT)))

# G3: binary AUROC by rank formula, then riskRegression
yb <- as.integer(sc$event_death == 1 & sc$time_death <= CUT)
r_ <- rank(sc$score_ecg_death)
n1 <- sum(yb)
n0 <- length(yb) - n1
auc_bin <- (sum(r_[yb == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
p_bin <- pyv("binary AUROC death by 30.5 d")
add("G3", "binary AUROC death by 30.5 d", auc_bin, p_bin,
    abs(auc_bin - p_bin) < 1e-6)

sco <- Score(list(ecg = matrix(sc$score_ecg_death, ncol = 1)),
             formula = Hist(time_death, event_death) ~ 1,
             data = sc, times = CUT, metrics = "auc",
             null.model = FALSE, conf.int = FALSE,
             contrasts = FALSE)
tauc <- as.numeric(sco$AUC$score$AUC[1])
p_t <- pyv("time-dependent AUC at 30.5 d")
add("G3", "time-dependent AUC at 30.5 d (Score)", tauc, p_t,
    abs(tauc - auc_bin) < 5e-4 && abs(tauc - p_t) < 5e-4)

# G4: Harrell's C; reverse = TRUE because higher score = higher risk
h <- as.numeric(concordance(
  Surv(time_death, event_death) ~ score_ecg_death,
  data = sc, reverse = TRUE)$concordance)
p_h <- pyv("Harrell C")
add("G4", "Harrell C (concordance, reverse)", h, p_h,
    h > 0.5 && abs(h - p_h) < 1e-4)

# G5: Uno-type weights, truncated at day 31
u <- as.numeric(concordance(
  Surv(time_death, event_death) ~ score_ecg_death,
  data = sc, reverse = TRUE, timewt = "n/G2",
  ymax = TAU_C)$concordance)
p_u <- pyv("Uno C tau 31")
add("G5", "Uno C (n/G2 weights, ymax 31)", u, p_u,
    abs(u - h) < 1e-4 && abs(u - p_u) < 1e-4)

out <- do.call(rbind, rows)
write.csv(out, file.path(P2, "results", "gates_r.csv"),
          row.names = FALSE)
nf <- sum(out$status == "FAIL")
cat(sprintf("\n  PASS %d   FAIL %d\n",
            sum(out$status == "PASS"), nf))
cat("\nOUTPUT FILES\n")
cat("  ", file.path(P2, "results", "gates_r.csv"), "\n")
cat("  ", file.path(P2, "results", "gates_r_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: gates.csv and gates_r.csv\n")
if (nf > 0) quit(status = 1)