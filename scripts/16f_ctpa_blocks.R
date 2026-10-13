# 16f_ctpa_blocks.R
# Which content of the CTPA report carries its short- and long-term signal
# for death (PREREGISTRATION_ADDENDUM_2.md, A28, item 2).
#
# Blocks: the 46 CTPA features in the five content blocks of the extractor
# in the dissertation (pe-multimodal-prognosis-msc): PE descriptors (13),
# incidental findings (7), comorbidity (cm_*, 16), support devices (dv_*, 8)
# and report metadata (txt_len, n_sent). Columns are matched on the name
# after any "ctpa_" or "ctpa." prefix, and every CTPA column must fall in a
# block.
#
# Scores: Frame B, three-modality cohort (death within 365 days of
# discharge). For each block, and for all 46 features, a ridge-Cox death
# score fitted as in 16e (18a's preprocessing, glmnet alpha 0, D11 path,
# inner 3-fold patient-grouped deviance), out of fold over the folds of 04
# and six seeds and seed-averaged. The EHR score is the retuned ridge-Cox
# score of 06d, seed-averaged.
#
# Time pattern: per standardised score, as in the primary analysis and 16e,
# the hazard ratios for death in days 0-30 and 31-365 and their ratio (early
# / late; below 1 means late-loaded), alone and with the EHR score (both
# with early and late terms). Landmark discrimination: AUROC for death by
# day 30, and for death in days 31-365 among patients alive at day 30.
# Contrasts: the ratio of early/late ratios of each block over all 46
# features and over the PE descriptors, and the landmark AUROC (days 31-365)
# of each block minus the PE descriptors, with a 2,000-sample bootstrap
# (seed 42) and the scores held fixed.
#
# Checks:  stops if inputs or scores are missing, a CTPA column is in no
#          block or the D2 per-fold rule fails; warns on block sizes other
#          than 13/7/16/8/2 and on over 20% of penalty choices at the
#          least-penalised end.
# Notes:   a timing probe runs first; --workers=N sets the worker count.
#
# Inputs:  data/processed/frame_B.csv, folds.csv, feature_whitelist.csv,
#          level1/oof_extra_B_3mod.csv; data/inputs/ctpa46.csv
# Outputs: results/ctpa_blocks_time.csv, ctpa_blocks_contrasts.csv,
#          ctpa_blocks_checks.csv, ctpa_blocks_log.txt
# Usage:   Rscript scripts/16f_ctpa_blocks.R [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CUT <- 30.5
TAU <- 365.5
NLAM <- 80L
LMR <- 1e-7
NBOOT <- 2000L
SEED <- 42L
PE_DESC <- c("pe_pos", "pe_neg", "saddle", "central", "lobar", "segmental",
             "subseg", "bilateral", "rv_strain", "septal_bow", "reflux",
             "mpa_enlarge", "infarct", "chronic_pe")
INCID <- c("effusion", "malignancy", "consolid", "atelect", "edema",
           "cardiomeg", "adenopathy")
META <- c("txt_len", "n_sent")
EXPECT <- c("PE descriptors" = 13, "Incidental findings" = 7,
            "Comorbidity" = 16, "Support devices" = 8,
            "Report metadata" = 2)
PKG <- c("survival", "glmnet", "data.table")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 200)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}
stop_if_fail <- function() {
  ck <- rbindlist(CHK)
  if (any(ck$status == "FAIL")) {
    fwrite(ck, file.path(RS, "ctpa_blocks_checks.csv"))
    cat("\nSTOP: a check failed. see: ctpa_blocks_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- one scoring task = block x seed ----------------
score_task <- function(tk) {
  prep <- function(X, tr) {
    for (j in seq_len(ncol(X))) {
      m <- median(X[tr, j], na.rm = TRUE)
      if (!is.finite(m)) m <- 0
      X[is.na(X[, j]), j] <- m
      if (length(unique(X[tr, j])) > 2) {
        q <- quantile(X[tr, j], c(0.005, 0.995), names = FALSE)
        X[, j] <- pmin(pmax(X[, j], q[1]), q[2])
      }
    }
    mu <- colMeans(X[tr, , drop = FALSE])
    s <- apply(X[tr, , drop = FALSE], 2, sd)
    keep <- is.finite(s) & s > 0
    sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
  }
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  lp <- rep(NA_real_, length(tk$t))
  least <- logical(0)
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- prep(tk$X, tr)
    if (ncol(X) == 1) X <- cbind(X, 0 * X)
    cv <- glmnet::cv.glmnet(X[tr, , drop = FALSE],
                            survival::Surv(tk$t[tr], tk$e[tr]),
                            family = "cox", alpha = 0,
                            foldid = grouped(tr, tk$seed * 100 + k, 3),
                            nlambda = tk$nlam, lambda.min.ratio = tk$lmr,
                            type.measure = "deviance",
                            penalty.factor = rep(1, ncol(X)),
                            standardize = FALSE)
    least <- c(least, which(cv$lambda == cv$lambda.min) ==
                 length(cv$lambda))
    lp[te] <- as.vector(predict(cv, X[te, , drop = FALSE],
                                s = "lambda.min"))
  }
  list(block = tk$block, seed = tk$seed, lp = lp, least = least)
}

# ---------------- statistics for one resample ----------------
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  tt <- J$t[idx]
  ee <- J$e[idx]
  z <- function(v) (v - mean(v)) / sd(v)
  ep <- survival::survSplit(Surv(tt, ee) ~ ., cut = J$cut,
                            episode = "per",
                            data = data.frame(tt = tt, ee = ee,
                                              id = seq_along(tt)))
  early <- as.numeric(ep$per == 1)
  ze <- z(J$ehr[idx])[ep$id]
  v <- numeric(0)
  y30 <- ee == 1 & tt <= J$cut
  alive <- tt > J$cut
  ylate <- ee[alive] == 1 & tt[alive] <= J$tau
  for (b in names(J$S)) {
    s <- J$S[[b]][idx]
    zs <- z(s)[ep$id]
    fit1 <- tryCatch(coef(survival::coxph(
      survival::Surv(ep$tstart, ep$tt, ep$ee) ~
        I(zs * early) + I(zs * (1 - early)), ties = "efron")),
      error = function(e) c(NA, NA))
    fit2 <- tryCatch(coef(survival::coxph(
      survival::Surv(ep$tstart, ep$tt, ep$ee) ~
        I(zs * early) + I(zs * (1 - early)) + I(ze * early) +
        I(ze * (1 - early)), ties = "efron")),
      error = function(e) c(NA, NA, NA, NA))
    v[paste(b, "alone", "early_logHR", sep = "|")] <- fit1[[1]]
    v[paste(b, "alone", "late_logHR", sep = "|")] <- fit1[[2]]
    v[paste(b, "alone", "log_ratio", sep = "|")] <- fit1[[1]] - fit1[[2]]
    v[paste(b, "with EHR", "early_logHR", sep = "|")] <- fit2[[1]]
    v[paste(b, "with EHR", "late_logHR", sep = "|")] <- fit2[[2]]
    v[paste(b, "with EHR", "log_ratio", sep = "|")] <- fit2[[1]] -
      fit2[[2]]
    v[paste(b, "alone", "auc_day30", sep = "|")] <- auc(y30, s)
    v[paste(b, "alone", "auc_days31_365", sep = "|")] <- auc(ylate,
                                                             s[alive])
  }
  v
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PR, c("frame_B.csv", "folds.csv",
                          "feature_whitelist.csv",
                          "level1/oof_extra_B_3mod.csv")),
          file.path(IN, "ctpa46.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fb <- fread(file.path(PR, "frame_B.csv"))
sub <- fb[in_three_mod == 1]
folds <- fread(file.path(PR, "folds.csv"))[frame == "B" & cohort == "3-mod"]
wl <- fread(file.path(PR, "feature_whitelist.csv"))
ct <- unique(fread(file.path(IN, "ctpa46.csv")), by = "hadm_id")
cc <- intersect(wl[modality == "ctpa", column], names(ct))
base <- sub("^ctpa[._]", "", cc)
blk <- ifelse(base %in% PE_DESC, "PE descriptors",
              ifelse(base %in% INCID, "Incidental findings",
                     ifelse(base %in% META, "Report metadata",
                            ifelse(startsWith(base, "cm_"), "Comorbidity",
                                   ifelse(startsWith(base, "dv_"),
                                          "Support devices", NA)))))
un <- cc[is.na(blk)]
chk("CTPA columns assigned to a block (unmatched)",
    sprintf("%d of %d%s", sum(!is.na(blk)), length(cc),
            if (length(un)) paste0(" (", paste(head(un, 6),
                                                collapse = ", "), ")") else
              ""), if (length(un) == 0) "PASS" else "FAIL")
stop_if_fail()
sz <- table(factor(blk, levels = names(EXPECT)))
chk("block sizes (PE / incidental / comorbidity / devices / metadata)",
    paste(as.integer(sz), collapse = " / "),
    if (all(as.integer(sz) == EXPECT)) "PASS" else "WARN")
XA <- as.matrix(ct[match(sub$hadm_id, hadm_id), ..cc])
storage.mode(XA) <- "double"
BLOCKS <- c(list("All CTPA" = cc), split(cc, blk)[names(EXPECT)])
XT <- fread(file.path(L1, "oof_extra_B_3mod.csv"),
            select = c("seed", "hadm_id", "ehr.cox2.death.lp"))
ehr <- XT[, .(v = mean(ehr.cox2.death.lp)), by = hadm_id][
  match(sub$hadm_id, hadm_id), v]
seeds <- sort(unique(folds$seed))
tasks <- list()
for (b in names(BLOCKS)) for (s in seeds) {
  f <- folds[folds$seed == s][match(sub$hadm_id, hadm_id), fold]
  tasks[[length(tasks) + 1]] <- list(
    block = b, seed = s, X = XA[, BLOCKS[[b]], drop = FALSE],
    t = sub$time_death, e = as.integer(sub$event_death == 1),
    subj = sub$subject_id, fold = as.integer(f), nlam = NLAM, lmr = LMR)
}
nm <- sum(!is.finite(ehr)) + sum(sapply(tasks, function(t) anyNA(t$fold)))
chk("EHR scores and folds complete (missing)", nm,
    if (nm == 0) "PASS" else "FAIL")
stop_if_fail()
cat(sprintf("  three-modality cohort n=%d, deaths by 365 %d; %d blocks x %d",
            nrow(sub), sum(sub$event_death == 1 & sub$time_death <= TAU),
            length(BLOCKS), length(seeds)), "seeds\n")

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
t0 <- Sys.time()
invisible(score_task(tasks[[1]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  one scoring task %.0f s; %d tasks on %d workers: about",
            tp, length(tasks), workers),
    sprintf("%.1f min\n", tp * ceiling(length(tasks) / workers) * 1.2 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, score_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1L)

cat("\n=== CHECK SCORES ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
S <- list()
nbad <- 0L
nmiss <- 0L
lst <- c()
for (b in names(BLOCKS)) {
  rr <- Filter(function(r) r$block == b, res)
  for (r in rr) {
    f <- tasks[[which(sapply(tasks, function(t) t$block == b &&
                               t$seed == r$seed))]]$fold
    v <- sapply(sort(unique(f)), function(k) {
      ii <- which(f == k)
      cidx(sub$time_death[ii], sub$event_death[ii] == 1, r$lp[ii])
    })
    if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
    lst <- c(lst, r$least)
  }
  S[[b]] <- Reduce(`+`, lapply(rr, function(r) r$lp)) / length(rr)
  nmiss <- nmiss + sum(!is.finite(S[[b]]))
}
chk("block scores missing", nmiss, if (nmiss == 0) "PASS" else "FAIL")
chk("block scores failing the per-fold rule (D2)", nbad,
    if (nbad == 0) "PASS" else "FAIL")
chk("penalty at the least-penalised end", sprintf("%.1f%%",
                                                  100 * mean(lst)),
    if (mean(lst) <= 0.20) "PASS" else "WARN")
stop_if_fail()
J <- list(S = S, ehr = ehr, t = sub$time_death,
          e = as.integer(sub$event_death == 1), cut = CUT, tau = TAU)

cat(sprintf("\n=== BOOTSTRAP (%d resamples) ===\n", NBOOT))
pt <- stat_vec(J, seq_len(nrow(sub)))
B <- do.call(rbind, future_lapply(seq_len(NBOOT), function(i) {
  stat_vec(J, sample.int(nrow(sub), nrow(sub), replace = TRUE))
}, future.seed = SEED, future.packages = PKG))
plan(sequential)
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
OUT <- list()
for (k in names(pt)) {
  z <- strsplit(k, "|", fixed = TRUE)[[1]]
  q <- q2(B[, k])
  isl <- grepl("logHR|log_ratio", z[3])
  OUT[[length(OUT) + 1]] <- data.table(
    block = z[1], model = z[2],
    quantity = sub("_logHR$|^log_", "", z[3]),
    estimate = if (isl) exp(pt[[k]]) else pt[[k]],
    lo = if (isl) exp(q[1]) else q[1],
    hi = if (isl) exp(q[2]) else q[2],
    n_features = if (z[1] == "All CTPA") length(cc) else
      length(BLOCKS[[z[1]]]))
}
tm <- rbindlist(OUT)
tm[quantity == "ratio", quantity := "early / late ratio"]
tm[quantity == "early", quantity := "HR days 0-30 (per SD)"]
tm[quantity == "late", quantity := "HR days 31-365 (per SD)"]
CON <- list()
for (b in setdiff(names(BLOCKS), "All CTPA")) {
  for (ref in c("All CTPA", "PE descriptors")) {
    if (b == ref) next
    for (md in c("alone", "with EHR")) {
      ka <- paste(b, md, "log_ratio", sep = "|")
      kb <- paste(ref, md, "log_ratio", sep = "|")
      q <- q2(B[, ka] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        contrast = sprintf("%s over %s", b, ref), model = md,
        quantity = "ratio of early/late ratios",
        estimate = exp(pt[[ka]] - pt[[kb]]), lo = exp(q[1]),
        hi = exp(q[2]), excludes_one = q[1] > 0 | q[2] < 0)
    }
  }
  if (b != "PE descriptors") {
    ka <- paste(b, "alone", "auc_days31_365", sep = "|")
    kb <- paste("PE descriptors", "alone", "auc_days31_365", sep = "|")
    q <- q2(B[, ka] - B[, kb])
    CON[[length(CON) + 1]] <- data.table(
      contrast = sprintf("%s minus PE descriptors", b), model = "alone",
      quantity = "AUROC, deaths in days 31-365 among day-30 survivors",
      estimate = pt[[ka]] - pt[[kb]], lo = q[1], hi = q[2],
      excludes_one = q[1] > 0 | q[2] < 0)
  }
}
cn <- rbindlist(CON)
setnames(cn, "excludes_one", "excludes_null")
fwrite(tm, file.path(RS, "ctpa_blocks_time.csv"))
fwrite(cn, file.path(RS, "ctpa_blocks_contrasts.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "ctpa_blocks_checks.csv"))
cat("\n=== TIME PATTERN BY BLOCK (CTPA alone) ===\n")
print(dcast(tm[model == "alone"], block + n_features ~ quantity,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE)
cat("\n=== CONTRASTS ===\n")
print(cn[, .(contrast, model, quantity, est = round(estimate, 3),
             lo = round(lo, 3), hi = round(hi, 3), excludes_null)],
      row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("ctpa_blocks_time.csv", "ctpa_blocks_contrasts.csv",
            "ctpa_blocks_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "ctpa_blocks_log.txt"), " (terminal log)\n")
cat("\nOutputs: ctpa_blocks_time.csv, ctpa_blocks_contrasts.csv,",
    "ctpa_blocks_checks.csv\n")