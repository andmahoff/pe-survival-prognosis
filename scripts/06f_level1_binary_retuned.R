# 06f_level1_binary_retuned.R
# Level-1 binary learners refitted like the retuned ridge-Cox
# (PREREGISTRATION_ADDENDUM_2.md A20).
#
# Stage 6a's binary twins (L2 logistic regression, one per label and
# horizon) are refitted as in stage 6e, so that the registered
# Cox-versus-binary comparisons (S1, S3, S2, calibration) are like for like:
# the same preprocessing, then glmnet binomial, alpha 0, 80 penalties down
# to 1e-7 of the largest, lambda by inner 3-fold cross-validated binomial
# deviance (log-loss), standardize = FALSE, predicted probabilities. Inner
# folds group each patient's admissions and are stratified by the label, so
# rare labels (CV readmission by day 10) have events in every inner fold.
# Folds are the stage 4 folds, for every seed. A label with fewer than 3
# events (or non-events) in a training fold, or a fit glmnet cannot
# complete, gets the training fold's event rate as its prediction; these
# cases are counted in the checks.
#
# Labels are every {modality}.lr.{label}_h{h} column of the stage 6a file
# (death_h, composite_h, cvu_h, cvc_h). oof_linear_retuned2_*.csv is stage
# 6e's file with every binary column replaced, same rows and columns in the
# same order, so the registered scripts can read it in place of 6a's. No
# performance numbers are written. The checks FAIL if a binary column is not
# replaced or has missing values, if the ridge-Cox columns differ from stage
# 6e, if the per-fold AUROC rule fails (mean below 0.40, as D2) or if over
# 5% of fits fall back to the event rate (WARN if any). They WARN when over
# 20% of penalty choices are at the least-penalised end.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_linear_retuned_*.csv (stage
#          6e); data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: results/level1_binary_retuned_checks.csv,
#          results/level1_binary_retuned_log.txt;
#          data/processed/level1/oof_linear_retuned2_*.csv (patient level,
#          not in the repo)
# Usage:   Rscript scripts/06f_level1_binary_retuned.R [--probe-only]
#          [--workers=N]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. --workers=N sets the worker
#          count (default: physical cores - 1).

suppressPackageStartupMessages({
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
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
NLAM <- 80L
LMR <- 1e-7
PKG <- c("glmnet", "data.table")
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
    fwrite(ck, file.path(RS, "level1_binary_retuned_checks.csv"))
    cat("\nSTOP: a check failed. see: level1_binary_retuned_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- one task = frame, cohort, modality, seed ----------------
run_task <- function(tk) {
  n <- length(tk$subj)
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
    X <- sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
    X
  }
  strat_grouped <- function(rows, y, seed, k) {
    sj <- tk$subj[rows]
    us <- unique(sj)
    pos <- as.vector(tapply(y[rows], sj, max)[as.character(us)])
    set.seed(seed)
    f <- integer(length(us))
    for (g in c(0, 1)) {
      ii <- which(pos == g)
      if (length(ii)) f[ii] <- sample(rep_len(seq_len(k), length(ii)))
    }
    f[match(sj, us)]
  }
  cols <- list()
  info <- list()
  for (lb in names(tk$L)) {
    cols[[sprintf("%s.lr.%s", tk$mod, lb)]] <- rep(NA_real_, n)
  }
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- prep(tk$X, tr)
    for (lb in names(tk$L)) {
      y <- tk$L[[lb]]$y
      trl <- tr[tk$L[[lb]]$keep[tr]]
      cn <- sprintf("%s.lr.%s", tk$mod, lb)
      ne <- sum(y[trl])
      fit <- if (ne < 3 || length(trl) - ne < 3) {
        "fewer than 3 events or non-events"
      } else {
        inner <- strat_grouped(trl, y, tk$seed * 100 + k, 3)
        tryCatch(suppressWarnings(cv.glmnet(
          X[trl, , drop = FALSE], y[trl], family = "binomial", alpha = 0,
          foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
          type.measure = "deviance", standardize = FALSE)),
          error = function(e) conditionMessage(e))
      }
      if (is.character(fit)) {
        cols[[cn]][te] <- mean(y[trl])
        info[[length(info) + 1]] <- data.table(
          label = lb, fold = k, ok = FALSE, least = NA,
          msg = substr(fit, 1, 100))
        next
      }
      i <- which(fit$lambda == fit$lambda.min)[1]
      info[[length(info) + 1]] <- data.table(
        label = lb, fold = k, ok = TRUE, least = i == length(fit$lambda),
        msg = "")
      cols[[cn]][te] <- as.vector(predict(fit, X[te, , drop = FALSE],
                                          s = "lambda.min",
                                          type = "response"))
    }
  }
  list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
       cols = as.data.table(cols), info = rbindlist(info))
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- lapply(MATS, function(f) unique(fread(file.path(IN, f)),
                                        by = "hadm_id"))
ldef <- function(sub, lb) {
  base <- sub("_h[0-9]+$", "", lb)
  cut <- as.numeric(sub("^.*_h", "", lb)) + 0.5
  allr <- rep(TRUE, nrow(sub))
  cvb <- sub$event_cr == 1 & sub$time_first <= cut
  dfb <- sub$event_cr == 2 & sub$time_first <= cut
  switch(base,
         death = list(y = as.integer(sub$event_death == 1 &
                                       sub$time_death <= cut), keep = allr),
         composite = list(y = as.integer(sub$event_comp == 1 &
                                           sub$time_first <= cut),
                          keep = allr),
         cvu = list(y = as.integer(cvb), keep = allr),
         cvc = list(y = as.integer(cvb), keep = !dfb),
         NULL)
}
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  Fd <- if (fr == "A") fa else fb
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  key <- paste(fr, co)
  f6 <- file.path(L1, sprintf("oof_linear_retuned_%s_%s.csv", fr, code))
  if (!file.exists(f6)) {
    chk(sprintf("%s: stage 6e file", key), basename(f6), "FAIL")
    stop_if_fail()
  }
  hdr <- names(fread(f6, nrows = 0))
  bc <- grep("\\.lr\\.", hdr, value = TRUE)
  mods <- unique(sub("\\..*$", "", bc))
  lbs <- unique(sub("^[a-z]+\\.lr\\.", "", bc))
  unk <- lbs[sapply(lbs, function(l) is.null(ldef(sub, l)))]
  chk(sprintf("%s: binary labels / columns in the file", key),
      sprintf("%d / %d%s", length(lbs), length(bc),
              if (length(unk)) paste0(" (unknown: ", unk[1], ")") else ""),
      if (length(unk) == 0) "PASS" else "FAIL")
  stop_if_fail()
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  for (m in mods) {
    LB <- list()
    for (l in lbs) {
      if (sprintf("%s.lr.%s", m, l) %in% bc) LB[[l]] <- ldef(sub, l)
    }
    vc <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..vc])
    storage.mode(Xm) <- "double"
    for (s in seeds) {
      fs <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
      if (anyNA(fs)) {
        chk(sprintf("%s seed %s: folds complete", key, s), sum(is.na(fs)),
            "FAIL")
        stop_if_fail()
      }
      tasks[[length(tasks) + 1]] <- list(
        fr = fr, co = co, mod = m, seed = s, X = Xm, L = LB,
        subj = sub$subject_id, fold = as.integer(fs))
    }
  }
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub, f6 = f6,
                      bc = bc, seeds = seeds, mods = mods, lbs = lbs)
}
cat(sprintf("  %d tasks (frame x cohort x modality x seed)\n",
            length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) nrow(t$X) * ncol(t$X) *
                          length(t$L)))
t0 <- Sys.time()
invisible(run_task(tasks[[big]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(paste("  largest task %.0f s; %d tasks on %d workers: at most",
                  "about %.1f min\n"), tp, length(tasks), workers,
            tp * ceiling(length(tasks) / workers) * 1.2 / 60))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, run_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1L)
plan(sequential)

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE AND CHECK ===\n")
auc <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- frank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
INFO <- rbindlist(lapply(res, function(r) r$info[, mod := r$mod]))
nf <- sum(!INFO$ok)
chk("fits replaced by the training event rate (first reason)",
    sprintf("%d of %d%s", nf, nrow(INFO),
            if (nf) paste0(" (", INFO$msg[!INFO$ok][1], ")") else ""),
    if (nf == 0) "PASS" else if (nf / nrow(INFO) <= 0.05) "WARN" else
      "FAIL")
if (nf) {
  ft <- INFO[ok == FALSE, .N, by = .(mod, label)]
  chk("  where", paste(sprintf("%s %s x%d", ft$mod, ft$label, ft$N),
                       collapse = "; "), "INFO")
}
for (key in names(META)) {
  M <- META[[key]]
  L6 <- fread(M$f6)
  out <- copy(L6)
  for (s in M$seeds) {
    rr <- Filter(function(r) r$fr == M$fr && r$co == M$co &&
                   r$seed == s, res)
    rows <- which(out$seed == s)
    pos <- match(out$hadm_id[rows], M$sub$hadm_id)
    for (r in rr) {
      for (cn in names(r$cols)) {
        if (cn %in% M$bc) set(out, rows, cn, r$cols[[cn]][pos])
      }
    }
  }
  nmiss <- sum(!is.finite(as.matrix(out[, M$bc, with = FALSE])))
  nsame <- sum(sapply(M$bc, function(cn) {
    isTRUE(all.equal(out[[cn]], L6[[cn]], tolerance = 0))
  }))
  chk(sprintf("%s: binary columns replaced / unchanged / missing", key),
      sprintf("%d / %d / %d", length(M$bc), nsame, nmiss),
      if (nmiss == 0 && nsame == 0) "PASS" else "FAIL")
  oth <- setdiff(names(L6), M$bc)
  same_oth <- all(sapply(oth, function(cn) identical(out[[cn]], L6[[cn]])))
  chk(sprintf("%s: ridge-Cox and other columns unchanged (%d)", key,
              length(oth)), same_oth, if (same_oth) "PASS" else "FAIL")
  nbad <- 0L
  for (cn in M$bc) {
    lb <- sub("^[a-z]+\\.lr\\.", "", cn)
    ld <- ldef(M$sub, lb)
    v <- c()
    for (s in M$seeds) {
      os <- out[out$seed == s]
      pos <- match(os$hadm_id, M$sub$hadm_id)
      for (k in sort(unique(os$fold))) {
        ii <- which(os$fold == k)
        v <- c(v, auc(ld$y[pos[ii]], os[[cn]][ii]))
      }
    }
    if (!is.finite(mean(v, na.rm = TRUE)) || mean(v, na.rm = TRUE) < 0.40) {
      nbad <- nbad + 1L
    }
  }
  chk(sprintf("%s: binary scores failing the per-fold AUROC rule", key),
      nbad, if (nbad == 0) "PASS" else "FAIL")
  fwrite(out, file.path(L1, sprintf("oof_linear_retuned2_%s_%s.csv",
                                    M$fr, M$code)))
}
INFO[, base := sub("_h[0-9]+$", "", label)]
es <- INFO[ok == TRUE, .(least = mean(least)), by = .(mod, base)]
for (i in seq_len(nrow(es))) {
  chk(sprintf("%s %s: penalty at the least-penalised end",
              toupper(es$mod[i]), es$base[i]),
      sprintf("%.1f%%", 100 * es$least[i]),
      if (es$least[i] <= 0.20) "PASS" else "WARN")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_binary_retuned_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(L1, sprintf("oof_linear_retuned2_%s_%s.csv",
                                  META[[key]]$fr, META[[key]]$code)),
      " (local only)\n")
}
cat("  ", file.path(RS, "level1_binary_retuned_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_binary_retuned_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: level1_binary_retuned_checks.csv\n")