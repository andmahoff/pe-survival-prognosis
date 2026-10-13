# 16p_ctpa_blocks_P.R
# The CTPA block split at diagnosis (PREREGISTRATION_ADDENDUM_2.md, A37,
# part 2): script 16f (A28, item 4) repeated on the A36 at-diagnosis model,
# for death from admission.
#
# Blocks: the five content blocks of 16f (PE descriptors 13, incidental
# findings 7, comorbidity 16, support devices 8, report metadata 2), plus
# two sub-blocks of the PE descriptors (A37): right-heart strain signs
# (rv_strain, septal_bow, reflux, mpa_enlarge) and clot location and extent
# (the other nine). Every CTPA column must fall in a block.
#
# Scores: Frame P, three-modality cohort (admissions with a presentation
# CTPA; death within 365 days of admission). For each block, sub-block and
# all 46 features, a ridge-Cox death score fitted exactly as the CTPA
# ridge-Cox of 06p: 18a's preprocessing inside each training fold, glmnet
# alpha 0, inner 3-fold patient-grouped deviance, and the path extended
# 100-fold if its smallest penalty is chosen (at most twice). Scores are out
# of fold over Frame A's three-modality folds and six seeds, then
# seed-averaged. The EHR score is the retuned ridge-Cox EHR score of 06p
# (A36 history), seed-averaged.
#
# Time pattern and contrasts, as in 16f but from admission: per standardised
# score, hazard ratios for death in days 0-30 and 31-365 and their ratio
# (early / late; below 1 means late-loaded), alone and with the EHR score
# (both with early and late terms), and AUROC for death by day 30 and for
# death in days 31-365 among admissions alive at day 30. Added contrasts:
# the AUROC for death by day 30 minus the PE descriptors, and right-heart
# strain signs against clot location and extent. The bootstrap has 2,000
# resamples of patients (seed 42), keeping repeat admissions together, with
# the scores held fixed.
#
# Checks:  stops if inputs or scores are missing, a CTPA column is in no
#          block, the D2 per-fold rule fails or "All CTPA" does not
#          reproduce 06p to 1e-6; warns on block sizes other than
#          13/7/16/8/2 and 4/9.
# Notes:   a timing probe runs first; --workers=N sets the worker count.
#
# Inputs:  data/processed/presentation/frame_P.csv, ctpaP.csv;
#          data/processed/folds.csv, feature_whitelist.csv,
#          level1/oof_P_3mod.csv; results/ctpa_blocks_time.csv (optional)
# Outputs: results/ctpa_blocks_P_time.csv, ctpa_blocks_P_contrasts.csv,
#          ctpa_blocks_P_checks.csv, ctpa_blocks_P_log.txt,
#          ctpa_blocks_P_vs_B.csv (if the Frame B results exist)
# Usage:   Rscript scripts/16p_ctpa_blocks_P.R [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CUT <- 30.5
TAU <- 365.5
NLAM <- 80L
LMR <- 1e-7
NBOOT <- 2000L
SEED <- 42L
TOLREP <- 1e-6
PE_DESC <- c("pe_pos", "pe_neg", "saddle", "central", "lobar", "segmental",
             "subseg", "bilateral", "rv_strain", "septal_bow", "reflux",
             "mpa_enlarge", "infarct", "chronic_pe")
RV <- c("rv_strain", "septal_bow", "reflux", "mpa_enlarge")
INCID <- c("effusion", "malignancy", "consolid", "atelect", "edema",
           "cardiomeg", "adenopathy")
META <- c("txt_len", "n_sent")
EXPECT <- c("PE descriptors" = 13, "Incidental findings" = 7,
            "Comorbidity" = 16, "Support devices" = 8,
            "Report metadata" = 2)
SUBS <- c("PE: right-heart strain signs" = 4,
          "PE: clot location and extent" = 9)
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
    fwrite(ck, file.path(RS, "ctpa_blocks_P_checks.csv"))
    cat("\nSTOP: a check failed. see: ctpa_blocks_P_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- one scoring task = block x seed ----------------
# stage 6p's CTPA ridge-Cox, line for line
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
  nconst <- 0L
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- prep(tk$X, tr)
    if (ncol(X) == 0) {
      lp[te] <- 0
      nconst <- nconst + 1L
      next
    }
    if (ncol(X) == 1) X <- cbind(X, 0 * X)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    lmr <- tk$lmr
    for (pass in 1:3) {
      cv <- glmnet::cv.glmnet(X[tr, , drop = FALSE],
                              survival::Surv(tk$t[tr], tk$e[tr]),
                              family = "cox", alpha = 0, foldid = inner,
                              nlambda = tk$nlam, lambda.min.ratio = lmr,
                              type.measure = "deviance",
                              penalty.factor = rep(1, ncol(X)),
                              standardize = FALSE)
      at_end <- which(cv$lambda == cv$lambda.min) == length(cv$lambda)
      if (!at_end || pass == 3) break
      lmr <- lmr / 100
    }
    least <- c(least, at_end)
    lp[te] <- as.vector(predict(cv, X[te, , drop = FALSE],
                                s = "lambda.min"))
  }
  list(block = tk$block, seed = tk$seed, lp = lp, least = least,
       nconst = nconst)
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
# resample patients, keeping each patient's admissions together
boot_idx <- function(subj) {
  us <- unique(subj)
  rows <- split(seq_along(subj), match(subj, us))
  pick <- sample.int(length(us), length(us), replace = TRUE)
  unlist(rows[pick], use.names = FALSE)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PP, c("frame_P.csv", "ctpaP.csv")),
          file.path(PR, c("folds.csv", "feature_whitelist.csv")),
          file.path(L1, "oof_P_3mod.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- fread(file.path(PP, "frame_P.csv"))
sub <- fp[in_three_mod_p == 1]
folds <- fread(file.path(PR, "folds.csv"))[frame == "A" &
                                             cohort == "3-mod"]
wl <- fread(file.path(PR, "feature_whitelist.csv"))
ct <- unique(fread(file.path(PP, "ctpaP.csv")), by = "hadm_id")
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
pe <- cc[blk %in% "PE descriptors"]
rv <- pe[sub("^ctpa[._]", "", pe) %in% RV]
clot <- setdiff(pe, rv)
chk("PE sub-blocks (right-heart strain / clot location and extent)",
    sprintf("%d / %d", length(rv), length(clot)),
    if (length(rv) == SUBS[[1]] && length(clot) == SUBS[[2]]) "PASS" else
      "WARN")
XA <- as.matrix(ct[match(sub$hadm_id, hadm_id), ..cc])
storage.mode(XA) <- "double"
chk("three-modality admissions with a CTPA row (missing)",
    sum(!sub$hadm_id %in% ct$hadm_id),
    if (all(sub$hadm_id %in% ct$hadm_id)) "PASS" else "FAIL")
stop_if_fail()
BLOCKS <- c(list("All CTPA" = cc), split(cc, blk)[names(EXPECT)],
            list("PE: right-heart strain signs" = rv,
                 "PE: clot location and extent" = clot))
BLOCKS <- BLOCKS[sapply(BLOCKS, length) > 0]
XT <- fread(file.path(L1, "oof_P_3mod.csv"),
            select = c("seed", "hadm_id", "fold", "ehr.cox2.death.lp",
                       "ctpa.cox2.death.lp"))
ehr <- XT[, .(v = mean(ehr.cox2.death.lp)), by = hadm_id][
  match(sub$hadm_id, hadm_id), v]
seeds <- sort(unique(folds$seed))
tasks <- list()
nmf <- 0L
for (b in names(BLOCKS)) for (s in seeds) {
  f <- folds[folds$seed == s][match(sub$hadm_id, hadm_id), fold]
  nmf <- nmf + sum(is.na(f))
  tasks[[length(tasks) + 1]] <- list(
    block = b, seed = s, X = XA[, BLOCKS[[b]], drop = FALSE],
    t = sub$time_death, e = as.integer(sub$event_death == 1),
    subj = sub$subject_id, fold = as.integer(f), nlam = NLAM, lmr = LMR)
}
nm <- sum(!is.finite(ehr)) + nmf
chk("EHR scores and folds complete (missing)", nm,
    if (nm == 0) "PASS" else "FAIL")
stop_if_fail()
d30 <- sum(sub$event_death == 1 & sub$time_death <= CUT)
d365 <- sum(sub$event_death == 1 & sub$time_death <= TAU)
cat(sprintf(paste("  three-modality cohort: %d admissions (%d patients),",
                  "deaths by day 30 %d, by day 365 %d\n"), nrow(sub),
            length(unique(sub$subject_id)), d30, d365))
cat(sprintf("  %d blocks x %d seeds = %d scoring tasks\n", length(BLOCKS),
            length(seeds), length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
t0 <- Sys.time()
invisible(score_task(tasks[[1]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
J0 <- list(S = setNames(rep(list(ehr), length(BLOCKS)), names(BLOCKS)),
           ehr = ehr, t = sub$time_death,
           e = as.integer(sub$event_death == 1), cut = CUT, tau = TAU)
set.seed(1)
t0 <- Sys.time()
for (i in 1:10) invisible(stat_vec(J0, boot_idx(sub$subject_id)))
tb <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10
proj <- (tp * ceiling(length(tasks) / workers) +
           tb * NBOOT / workers) * 1.2 / 60
cat(sprintf(paste("  one scoring task %.0f s; one bootstrap resample %.2f",
                  "s; about %.1f min on %d workers\n"), tp, tb, proj,
            workers))
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
ncon <- 0L
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
    ncon <- ncon + r$nconst
  }
  S[[b]] <- Reduce(`+`, lapply(rr, function(r) r$lp)) / length(rr)
  nmiss <- nmiss + sum(!is.finite(S[[b]]))
}
chk("block scores missing", nmiss, if (nmiss == 0) "PASS" else "FAIL")
chk("block scores failing the per-fold rule (D2)", nbad,
    if (nbad == 0) "PASS" else "FAIL")
chk("folds where a block had no varying feature", ncon,
    if (ncon == 0) "PASS" else "WARN")
chk("penalty at the least-penalised end, after extension",
    sprintf("%.1f%%", 100 * mean(lst)),
    if (mean(lst) <= 0.20) "PASS" else "WARN")
# "All CTPA" must reproduce stage 6p's CTPA ridge-Cox, seed by seed
dmax <- 0
for (r in Filter(function(r) r$block == "All CTPA", res)) {
  y <- XT[XT$seed == r$seed][match(sub$hadm_id, hadm_id)]
  dmax <- max(dmax, abs(y$ctpa.cox2.death.lp - r$lp), na.rm = FALSE)
}
chk("all 46 features reproduce stage 6p's CTPA ridge-Cox (max |diff|)",
    sprintf("%.1e", dmax), if (is.finite(dmax) && dmax <= TOLREP) "PASS"
    else "FAIL")
stop_if_fail()
nc <- sapply(S, function(s) sd(s) > 0)
if (!all(nc)) {
  chk("blocks with a constant score (left out)",
      paste(names(S)[!nc], collapse = ", "), "WARN")
  S <- S[nc]
}
J <- list(S = S, ehr = ehr, t = sub$time_death,
          e = as.integer(sub$event_death == 1), cut = CUT, tau = TAU)

cat(sprintf("\n=== BOOTSTRAP (%d resamples of patients) ===\n", NBOOT))
pt <- stat_vec(J, seq_len(nrow(sub)))
subj <- sub$subject_id
B <- do.call(rbind, future_lapply(seq_len(NBOOT), function(i) {
  stat_vec(J, boot_idx(subj))
}, future.seed = SEED, future.packages = PKG,
future.globals = c("J", "subj", "stat_vec", "boot_idx")))
plan(sequential)
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
nfeat <- function(b) length(BLOCKS[[b]])
OUT <- list()
for (k in names(pt)) {
  z <- strsplit(k, "|", fixed = TRUE)[[1]]
  q <- q2(B[, k])
  isl <- grepl("logHR|log_ratio", z[3])
  OUT[[length(OUT) + 1]] <- data.table(
    frame = "P", block = z[1], model = z[2],
    quantity = sub("_logHR$|^log_", "", z[3]),
    estimate = if (isl) exp(pt[[k]]) else pt[[k]],
    lo = if (isl) exp(q[1]) else q[1],
    hi = if (isl) exp(q[2]) else q[2],
    n_features = nfeat(z[1]))
}
tm <- rbindlist(OUT)
tm[quantity == "ratio", quantity := "early / late ratio"]
tm[quantity == "early", quantity := "HR days 0-30 (per SD)"]
tm[quantity == "late", quantity := "HR days 31-365 (per SD)"]
tm[quantity == "auc_day30", quantity := "AUROC, death by day 30"]
tm[quantity == "auc_days31_365",
   quantity := "AUROC, deaths in days 31-365 among day-30 survivors"]
CON <- list()
add_ratio <- function(b, ref) {
  for (md in c("alone", "with EHR")) {
    ka <- paste(b, md, "log_ratio", sep = "|")
    kb <- paste(ref, md, "log_ratio", sep = "|")
    if (!(ka %in% names(pt) && kb %in% names(pt))) next
    q <- q2(B[, ka] - B[, kb])
    CON[[length(CON) + 1]] <<- data.table(
      contrast = sprintf("%s over %s", b, ref), model = md,
      quantity = "ratio of early/late ratios",
      estimate = exp(pt[[ka]] - pt[[kb]]), lo = exp(q[1]),
      hi = exp(q[2]), excludes_null = q[1] > 0 | q[2] < 0)
  }
}
add_auc <- function(b, ref, what, label) {
  ka <- paste(b, "alone", what, sep = "|")
  kb <- paste(ref, "alone", what, sep = "|")
  if (!(ka %in% names(pt) && kb %in% names(pt))) return(invisible())
  q <- q2(B[, ka] - B[, kb])
  CON[[length(CON) + 1]] <<- data.table(
    contrast = sprintf("%s minus %s", b, ref), model = "alone",
    quantity = label, estimate = pt[[ka]] - pt[[kb]], lo = q[1],
    hi = q[2], excludes_null = q[1] > 0 | q[2] < 0)
}
LA <- "AUROC, deaths in days 31-365 among day-30 survivors"
L30 <- "AUROC, death by day 30"
main5 <- intersect(names(EXPECT), names(S))
for (b in main5) {
  for (ref in c("All CTPA", "PE descriptors")) {
    if (b != ref) add_ratio(b, ref)
  }
  if (b != "PE descriptors") {
    add_auc(b, "PE descriptors", "auc_days31_365", LA)
    add_auc(b, "PE descriptors", "auc_day30", L30)
  }
}
rvb <- "PE: right-heart strain signs"
clb <- "PE: clot location and extent"
add_ratio(rvb, clb)
add_auc(rvb, clb, "auc_day30", L30)
add_auc(rvb, clb, "auc_days31_365", LA)
cn <- rbindlist(CON)
cn[, frame := "P"]
setcolorder(cn, "frame")
fwrite(tm, file.path(RS, "ctpa_blocks_P_time.csv"))
fwrite(cn, file.path(RS, "ctpa_blocks_P_contrasts.csv"))

# ---------------- beside Frame B (stage 16f) ----------------
fb16 <- file.path(RS, "ctpa_blocks_time.csv")
VS <- NULL
if (file.exists(fb16)) {
  tb16 <- fread(fb16)
  tb16[quantity == "auc_day30", quantity := "AUROC, death by day 30"]
  tb16[quantity == "auc_days31_365",
       quantity := "AUROC, deaths in days 31-365 among day-30 survivors"]
  qs <- c("early / late ratio", "AUROC, death by day 30",
          "AUROC, deaths in days 31-365 among day-30 survivors")
  a <- tm[model == "alone" & quantity %in% qs,
          .(block, quantity, P = estimate, P_lo = lo, P_hi = hi)]
  b <- tb16[model == "alone" & quantity %in% qs,
            .(block, quantity, B = estimate, B_lo = lo, B_hi = hi)]
  VS <- merge(a, b, by = c("block", "quantity"), all.x = TRUE)
  fwrite(VS, file.path(RS, "ctpa_blocks_P_vs_B.csv"))
  chk("Frame B's stage 16f results found (written side by side)", "yes",
      "INFO")
} else {
  chk("Frame B's stage 16f results (results/ctpa_blocks_time.csv)",
      "not found, comparison skipped", "INFO")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "ctpa_blocks_P_checks.csv"))
cat("\n=== TIME PATTERN BY BLOCK, AT DIAGNOSIS (CTPA alone) ===\n")
print(dcast(tm[model == "alone"], block + n_features ~ quantity,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 3) else v
            })], row.names = FALSE)
cat("\n=== CONTRASTS ===\n")
print(cn[, .(contrast, model, quantity, est = round(estimate, 3),
             lo = round(lo, 3), hi = round(hi, 3), excludes_null)],
      row.names = FALSE, nrows = 200)
if (!is.null(VS)) {
  cat("\n=== AT DIAGNOSIS (P) BESIDE AFTER DISCHARGE (B), CTPA alone ===\n")
  print(VS[, .(block, quantity, P = round(P, 3), B = round(B, 3))],
        row.names = FALSE, nrows = 200)
}
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("ctpa_blocks_P_time.csv", "ctpa_blocks_P_contrasts.csv",
            if (!is.null(VS)) "ctpa_blocks_P_vs_B.csv",
            "ctpa_blocks_P_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "ctpa_blocks_P_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: ctpa_blocks_P_time.csv, ctpa_blocks_P_contrasts.csv,",
    if (!is.null(VS)) "ctpa_blocks_P_vs_B.csv,", "ctpa_blocks_P_checks.csv\n")