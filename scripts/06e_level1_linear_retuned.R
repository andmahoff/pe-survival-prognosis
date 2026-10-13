# 06e_level1_linear_retuned.R
# Level-1 ridge-Cox refitted with the D11 settings
# (PREREGISTRATION_ADDENDUM_2.md A17).
#
# Stage 6a's ridge-Cox is refitted for every target, modality, frame, cohort
# and seed, so that the registered analyses can be rerun on retuned scores.
# The method is that of stage 6d's ridge-Cox, which is stage 18a's EHR
# alone: median imputation, winsorising of non-binary features at the
# training fold's 0.5th and 99.5th percentiles, standardisation; glmnet,
# alpha 0, 80 penalties down to 1e-7 of the largest (D11); lambda by inner
# 3-fold subject-grouped cross-validated partial-likelihood deviance (18a's
# inner seeds), standardize = FALSE; risks from the centred Breslow baseline
# of the same target. Folds are the stage 4 folds, for every seed.
#
# oof_linear_retuned_*.csv copies the stage 6a file (same rows and columns,
# same order) with every ridge-Cox column replaced
# ({modality}.cox.{target}.lp and .r{h}; targets death, composite, cs_cv,
# cs_death). The binary columns are copied unchanged, so the registered
# scripts can read it in place of 6a's file. No performance numbers are
# written. The checks FAIL if a ridge-Cox column is not replaced or has
# missing values, if the death columns differ from stage 6d's ridge-Cox
# (beyond 1e-9) or if the D2 per-fold rule fails. They WARN when over 20% of
# penalty choices are at the least-penalised end.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, level1/oof_linear_*.csv,
#          level1/oof_extra_*.csv; data/inputs/ehr28.csv, ecg_opt.csv,
#          ctpa46.csv
# Outputs: results/level1_retuned_checks.csv,
#          results/level1_retuned_log.txt;
#          data/processed/level1/oof_linear_retuned_*.csv (patient level,
#          not in the repo)
# Usage:   Rscript scripts/06e_level1_linear_retuned.R [--probe-only]
#          [--workers=N]
# Notes:   --probe-only prints a run-time projection and stops; otherwise
#          the run starts after a 5 s pause. --workers=N sets the worker
#          count (default: physical cores - 1). Tasks run in the workers, as
#          in stage 18a, so the random inner folds match it.

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
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
MATS <- c(ehr = "ehr28.csv", ecg = "ecg_opt.csv", ctpa = "ctpa46.csv")
NLAM <- 80L
LMR <- 1e-7
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
    fwrite(ck, file.path(RS, "level1_retuned_checks.csv"))
    cat("\nSTOP: a check failed. see: level1_retuned_checks.csv\n")
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
  grouped <- function(rows, seed, k) {
    us <- unique(tk$subj[rows])
    set.seed(seed)
    f <- sample(rep_len(seq_len(k), length(us)))
    f[match(tk$subj[rows], us)]
  }
  breslow <- function(tt, ee, lp, cut) {
    o <- order(tt)
    tt <- tt[o]
    ee <- ee[o]
    w <- exp(lp[o] - max(lp))
    rs <- rev(cumsum(rev(w)))
    ut <- unique(tt[ee == 1 & tt <= cut])
    if (!length(ut)) return(0)
    d <- sapply(ut, function(u) sum(ee[tt == u]))
    sum(d / rs[match(ut, tt)])
  }
  cols <- list()
  edges <- list()
  for (tg in names(tk$T)) {
    cols[[sprintf("%s.cox.%s.lp", tk$mod, tg)]] <- rep(NA_real_, n)
    for (h in tk$T[[tg]]$hz) {
      cols[[sprintf("%s.cox.%s.r%d", tk$mod, tg, h)]] <- rep(NA_real_, n)
    }
  }
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    inner <- grouped(tr, tk$seed * 100 + k, 3)
    X <- prep(tk$X, tr)
    for (tg in names(tk$T)) {
      tt <- tk$T[[tg]]$t
      ee <- tk$T[[tg]]$e
      cv <- cv.glmnet(X[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                      family = "cox", alpha = 0, foldid = inner,
                      nlambda = NLAM, lambda.min.ratio = LMR,
                      type.measure = "deviance",
                      penalty.factor = rep(1, ncol(X)),
                      standardize = FALSE)
      i <- which(cv$lambda == cv$lambda.min)[1]
      edges[[length(edges) + 1]] <- data.table(
        target = tg, least = i == length(cv$lambda), most = i == 1)
      l_tr <- as.vector(predict(cv, X[tr, , drop = FALSE],
                                s = "lambda.min"))
      l_te <- as.vector(predict(cv, X[te, , drop = FALSE],
                                s = "lambda.min"))
      cols[[sprintf("%s.cox.%s.lp", tk$mod, tg)]][te] <- l_te
      c0 <- max(l_tr)
      for (h in tk$T[[tg]]$hz) {
        hs <- breslow(tt[tr], ee[tr], l_tr, h + 0.5)
        cols[[sprintf("%s.cox.%s.r%d", tk$mod, tg, h)]][te] <-
          if (hs == 0) 0 else -expm1(-hs * exp(l_te - c0))
      }
    }
  }
  list(fr = tk$fr, co = tk$co, mod = tk$mod, seed = tk$seed,
       cols = as.data.table(cols), edges = rbindlist(edges))
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- lapply(MATS, function(f) unique(fread(file.path(IN, f)),
                                        by = "hadm_id"))
tdef <- function(sub, tg) {
  switch(tg,
         death = list(t = sub$time_death,
                      e = as.integer(sub$event_death == 1)),
         composite = list(t = sub$time_first,
                          e = as.integer(sub$event_comp == 1)),
         cs_cv = list(t = sub$time_first,
                      e = as.integer(sub$event_cr == 1)),
         cs_death = list(t = sub$time_first,
                         e = as.integer(sub$event_cr == 2)),
         NULL)
}
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  Fd <- if (fr == "A") fa else fb
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  key <- paste(fr, co)
  f6 <- file.path(L1, sprintf("oof_linear_%s_%s.csv", fr, code))
  if (!file.exists(f6)) {
    chk(sprintf("%s: stage 6a file", key), basename(f6), "FAIL")
    stop_if_fail()
  }
  hdr <- names(fread(f6, nrows = 0))
  cc <- grep("\\.cox\\.", hdr, value = TRUE)
  mods <- unique(sub("\\..*$", "", cc))
  tgs <- unique(sub("^[a-z]+\\.cox\\.([a-z_]+)\\..*$", "\\1", cc))
  unk <- tgs[sapply(tgs, function(g) is.null(tdef(sub, g)))]
  chk(sprintf("%s: ridge-Cox targets / columns in stage 6a's file", key),
      sprintf("%s / %d", paste(tgs, collapse = ", "), length(cc)),
      if (length(unk) == 0) "PASS" else "FAIL")
  stop_if_fail()
  fo <- folds[frame == fr & cohort == co]
  seeds <- sort(unique(fo$seed))
  for (m in mods) {
    TG <- list()
    for (g in tgs) {
      if (!(sprintf("%s.cox.%s.lp", m, g) %in% cc)) next
      rh <- grep(sprintf("^%s\\.cox\\.%s\\.r[0-9]+$", m, g), cc,
                 value = TRUE)
      TG[[g]] <- c(tdef(sub, g),
                   list(hz = sort(as.integer(sub(".*\\.r", "", rh)))))
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
        fr = fr, co = co, mod = m, seed = s, X = Xm, T = TG,
        subj = sub$subject_id, fold = as.integer(fs))
    }
  }
  META[[key]] <- list(fr = fr, co = co, code = code, sub = sub, f6 = f6,
                      cc = cc, seeds = seeds, mods = mods)
}
cat(sprintf("  %d tasks (frame x cohort x modality x seed)\n",
            length(tasks)))

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
big <- which.max(sapply(tasks, function(t) nrow(t$X) * ncol(t$X) *
                          length(t$T)))
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
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
ED <- rbindlist(lapply(res, function(r) r$edges[, mod := r$mod]))
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
        if (cn %in% M$cc) set(out, rows, cn, r$cols[[cn]][pos])
      }
    }
  }
  nmiss <- sum(is.na(as.matrix(out[, M$cc, with = FALSE])))
  nsame <- sum(sapply(M$cc, function(cn) {
    isTRUE(all.equal(out[[cn]], L6[[cn]], tolerance = 0))
  }))
  chk(sprintf("%s: ridge-Cox columns replaced / unchanged / missing", key),
      sprintf("%d / %d / %d", length(M$cc), nsame, nmiss),
      if (nmiss == 0 && nsame == 0) "PASS" else "FAIL")
  oth <- setdiff(names(L6), M$cc)
  same_oth <- all(sapply(oth, function(cn) identical(out[[cn]], L6[[cn]])))
  chk(sprintf("%s: other columns copied unchanged (%d)", key, length(oth)),
      same_oth, if (same_oth) "PASS" else "FAIL")
  fx <- file.path(L1, sprintf("oof_extra_%s_%s.csv", M$fr, M$code))
  if (file.exists(fx)) {
    X6 <- fread(fx)
    dm <- 0
    for (m in M$mods) {
      for (cn in grep(paste0("^", m, "\\.cox\\.death\\."), M$cc,
                      value = TRUE)) {
        xn <- sub("\\.cox\\.", ".cox2.", cn)
        if (!(xn %in% names(X6))) next
        v <- X6[[xn]][match(paste(out$seed, out$hadm_id),
                            paste(X6$seed, X6$hadm_id))]
        dm <- max(dm, max(abs(v - out[[cn]])))
      }
    }
    chk(sprintf("%s: death columns equal stage 6d (max |diff|)", key),
        sprintf("%.1e", dm), if (is.finite(dm) && dm <= 1e-9) "PASS" else
          "FAIL")
  } else {
    chk(sprintf("%s: stage 6d file for the reproduction check", key),
        "missing", "WARN")
  }
  nbad <- 0L
  for (cn in grep("\\.lp$", M$cc, value = TRUE)) {
    g <- sub("^[a-z]+\\.cox\\.([a-z_]+)\\.lp$", "\\1", cn)
    td <- tdef(M$sub, g)
    v <- c()
    for (s in M$seeds) {
      os <- out[out$seed == s]
      pos <- match(os$hadm_id, M$sub$hadm_id)
      for (k in sort(unique(os$fold))) {
        ii <- which(os$fold == k)
        v <- c(v, cidx(td$t[pos[ii]], td$e[pos[ii]], os[[cn]][ii]))
      }
    }
    if (!is.finite(mean(v, na.rm = TRUE)) || mean(v, na.rm = TRUE) < 0.40) {
      nbad <- nbad + 1L
    }
  }
  chk(sprintf("%s: scores failing the per-fold rule (D2)", key), nbad,
      if (nbad == 0) "PASS" else "FAIL")
  fwrite(out, file.path(L1, sprintf("oof_linear_retuned_%s_%s.csv", M$fr,
                                    M$code)))
}
es <- ED[, .(least = mean(least), most = mean(most)), by = .(mod, target)]
for (i in seq_len(nrow(es))) {
  chk(sprintf("%s %s: penalty at the least / most penalised end",
              toupper(es$mod[i]), es$target[i]),
      sprintf("%.1f%% / %.1f%%", 100 * es$least[i], 100 * es$most[i]),
      if (es$least[i] <= 0.20) "PASS" else "WARN")
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "level1_retuned_checks.csv"))
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (key in names(META)) {
  cat("  ", file.path(L1, sprintf("oof_linear_retuned_%s_%s.csv",
                                  META[[key]]$fr, META[[key]]$code)),
      " (local only)\n")
}
cat("  ", file.path(RS, "level1_retuned_checks.csv"), "\n")
cat("  ", file.path(RS, "level1_retuned_log.txt"),
    " (terminal log)\n")
cat("\nOutputs: level1_retuned_checks.csv\n")