# 06j_level1_ehrplus.R
# Level-1 death scores from "EHR+", the 79-feature EHR plus information
# known at discharge (PREREGISTRATION_ADDENDUM_2.md, A30 and A31).
#
# Frame B only, since every addition is known by discharge (in Frame A the
# diagnosis codes would leak the outcome). EHR+ has four blocks:
#   1. EHR 79, built as in stage 19a and 06i
#   2. measurements pulled by stage 3b but not used before: spo2, nlr,
#      hr_first, sbp_first, spo2_first, sex_f, n_prior, n_codes
#   3. discharge destination (followup_365.csv) as indicators, home being
#      the reference: home health care, skilled nursing, rehabilitation,
#      long-term acute care, hospice, other
#   4. comorbidity from the index admission's diagnosis codes (stage 3d):
#      the 17 Charlson and 31 Elixhauser categories (R package comorbidity,
#      Quan's ICD-9-CM and ICD-10 maps, hierarchy applied when scoring),
#      each kept if its prevalence in Frame B is between 1% and 99%, plus
#      the Charlson score (Charlson weights) and the Elixhauser score (van
#      Walraven weights)
# Families, tuning, folds, checks, timing probe and checkpointing are as in
# 06i: ridge-Cox retuned, spline Cox, and Brier-tuned RSF and oblique RSF.
#
# Spline Cox, version 2 (D16): terms entered linearly (binary indicators
# and features with fewer than 10 distinct values) carry a ridge penalty
# (mgcv paraPen, strength estimated by REML with the smooths). A rare
# indicator with no deaths, or only deaths, in a training fold (such as
# hospice discharge) has no finite unpenalised coefficient, and the penalty
# keeps every coefficient finite. The other families are unchanged.
#
# --refit-spline refits only the spline family, reading the other scores
# from the first run's checkpoints (ehrplus_ckpt/) and saving to
# ehrplus_ckpt_v2/. --no-hospice (A31) omits the hospice destination
# indicator, so hospice discharges look like home, and fits every family
# afresh, with "ehrplus_nh" in place of "ehrplus" in every file name.
#
# Inputs:  data/processed/frame_B.csv, folds.csv, feature_whitelist.csv
#          data/inputs/ehr28.csv, ehr_extra.csv, followup_365.csv,
#          dx_index.csv, scripts/tune_expand.R
# Outputs: results/ehrplus_checks.csv, results/ehrplus_features.csv,
#          results/ehrplus_tuning_summary.csv, results/ehrplus_log.txt
#          data/processed/level1/oof_ehrplus_B_*.csv (patient level, not in
#          the repo; read by stage 18o)
# Usage:   Rscript scripts/06j_level1_ehrplus.R [--probe-only]
#          [--workers=N] [--refit-spline | --no-hospice]

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(mgcv)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf not installed")
if (!requireNamespace("comorbidity", quietly = TRUE)) {
  stop("package comorbidity not installed: install.packages(\"comorbidity\")")
}
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
refit_spline <- "--refit-spline" %in% args
no_hospice <- "--no-hospice" %in% args
TAG <- if (no_hospice) "ehrplus_nh" else "ehrplus"
if (no_hospice && refit_spline) stop("--refit-spline is for the main run")
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
IN <- file.path(P2, "data", "inputs")
PR <- file.path(P2, "data", "processed")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CKPT_OLD <- file.path(L1, paste0(TAG, "_ckpt"))
CKPT <- file.path(L1, paste0(TAG, "_ckpt_v2"))
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
source(file.path(P2, "scripts", "tune_expand.R"))
HZ <- c(10, 30, 90, 180, 365)
EXCL <- c("vancomycin", "_50825_")
EXTRA <- c("spo2", "nlr", "hr_first", "sbp_first", "spo2_first", "sex_f",
           "n_prior", "n_codes")
DEST <- c("HOME HEALTH CARE" = "dest_home_health",
          "SKILLED NURSING FACILITY" = "dest_skilled_nursing",
          "REHAB" = "dest_rehab",
          "CHRONIC/LONG TERM ACUTE CARE" = "dest_long_term_care",
          "HOSPICE" = "dest_hospice")
PREV <- c(0.01, 0.99)
NLAM <- 80L
LMR <- 1e-7
KNOTS <- 5L
MINU <- 10L
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
NODES <- c(5, 15, 40, 80, 160)
MTRYF <- c(0.5, 1, 3)
LEAF <- c(5, 20, 60)
OMTRYF <- c(0.5, 1, 2, 4)
BATCH <- 50L
TOL <- 0.001
MAXEXT <- 3L
NFOLD <- 5L
PKG <- c("survival", "glmnet", "mgcv", "data.table", "randomForestSRC",
         "aorsf")
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
    fwrite(ck, file.path(RS, paste0(TAG, "_checks.csv")))
    cat(sprintf("\nSTOP: a check failed. see: %s_checks.csv\n", TAG))
    quit(status = 1)
  }
}

# ---------------- one task = cohort x seed ----------------
run_task <- function(tk, folds_to_run = sort(unique(tk$fold)),
                     probe = FALSE) {
  op <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(op), add = TRUE)
  t_task <- proc.time()[["elapsed"]]
  n <- length(tk$t)
  H <- tk$hz
  tt <- tk$t
  ee <- tk$e
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
  old <- if (tk$mode == "spline") readRDS(tk$ckpt_old) else NULL
  out <- if (is.null(old)) {
    list(cox2 = rep(NA_real_, n), spline = rep(NA_real_, n),
         rsfb = rep(NA_real_, n), orsfb = rep(NA_real_, n))
  } else {
    old$out
  }
  prm <- list()
  exts <- list()
  fails <- list()
  for (k in folds_to_run) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    Xs <- prep(tk$X, tr)
    colnames(Xs) <- sprintf("f%d", seq_len(ncol(Xs)))
    if (is.null(old)) {
      # ridge-Cox retuned, path extended if the smallest penalty is chosen
      inner <- grouped(tr, tk$seed * 100 + k, 3)
      lmr <- tk$lmr
      for (pass in 1:3) {
        cv <- cv.glmnet(Xs[tr, , drop = FALSE], Surv(tt[tr], ee[tr]),
                        family = "cox", alpha = 0, foldid = inner,
                        nlambda = tk$nlam, lambda.min.ratio = lmr,
                        type.measure = "deviance",
                        penalty.factor = rep(1, ncol(Xs)),
                        standardize = FALSE)
        at_end <- which(cv$lambda == cv$lambda.min) == length(cv$lambda)
        if (!at_end || pass == 3) break
        lmr <- lmr / 100
      }
      out$cox2[te] <- as.vector(predict(cv, Xs[te, , drop = FALSE],
                                        s = "lambda.min"))
    }
    # spline Cox (D14's fallback rule); linear terms ridge-penalised (D16)
    nu <- apply(Xs[tr, , drop = FALSE], 2, function(v) length(unique(v)))
    sm <- colnames(Xs)[nu >= tk$minu]
    li <- setdiff(colnames(Xs), sm)
    dsp <- data.frame(t = tt[tr], e = ee[tr], Xs[tr, sm, drop = FALSE])
    dnw <- data.frame(Xs[te, sm, drop = FALSE])
    pp <- NULL
    if (length(li)) {
      dsp$L <- Xs[tr, li, drop = FALSE]
      dnw$L <- Xs[te, li, drop = FALSE]
      pp <- list(L = list(diag(length(li))))
    }
    sk <- 0L
    for (kk in unique(c(tk$knots, 4L, 3L))) {
      fml <- as.formula(paste("t ~", paste(c(sprintf(
        "s(%s, k = %d, bs = \"cr\")", sm, kk), if (length(li)) "L"),
        collapse = " + ")))
      g <- tryCatch(mgcv::gam(fml, family = mgcv::cox.ph(), weights = e,
                              data = dsp, method = "REML", select = TRUE,
                              optimizer = "efs", paraPen = pp),
                    error = function(err) conditionMessage(err))
      if (is.character(g)) {
        fails[[length(fails) + 1]] <- data.table(
          fold = k, family = "spline Cox", setting = sprintf("k = %d", kk),
          message = substr(g, 1, 150))
        next
      }
      out$spline[te] <- as.vector(predict(g, newdata = dnw,
                                          type = "link"))
      sk <- kk
      break
    }
    if (sk == 0L) out$spline[te] <- out$cox2[te]
    if (!is.null(old)) {
      pr <- old$prm[old$prm$fold == k]
      pr$spline_k <- sk
      prm[[length(prm) + 1]] <- pr
      next
    }
    Ytr <- matrix(sapply(H, function(h) as.numeric(ee[tr] == 1 &
                                                      tt[tr] <= h + 0.5)),
                  ncol = length(H))
    # RSF, Brier-tuned, grid extended past its edges
    X <- tk$X
    med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    na <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na) > 0) X[na] <- med[na[, 2]]
    colnames(X) <- sprintf("f%d", seq_len(ncol(X)))
    ev <- sort(unique(tt[tr][ee[tr] == 1]))
    at_h <- sapply(H, function(h) {
      v <- ev[ev <= h + 0.5]
      if (length(v)) max(v) else NA_real_
    })
    pick <- ev[unique(round(seq(1, length(ev),
                                length.out = min(tk$ntime, length(ev)))))]
    nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
    dtr <- data.frame(time = tt[tr], status = ee[tr], X[tr, , drop = FALSE])
    dte <- as.data.frame(X[te, , drop = FALSE])
    p <- ncol(X)
    gi <- 0L
    rsf_eval <- function(par) {
      gi <<- gi + 1L
      f <- rfsrc(Surv(time, status) ~ ., data = dtr, ntree = tk$ntree,
                 nodesize = par$node, mtry = par$mtry, nsplit = tk$nsplit,
                 ntime = nt, splitrule = "logrank", save.memory = TRUE,
                 seed = -as.integer(tk$seed * 1000 + k * 100 + gi))
      ro <- sapply(H, function(h) {
        i <- findInterval(h + 0.5, f$time.interest)
        if (i == 0) rep(0, length(tr)) else 1 - f$survival.oob[, i]
      })
      br <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2))
      sc <- predict(f, newdata = dte)$predicted
      rm(f)
      invisible(gc(verbose = FALSE))
      list(score = br, pred = sc)
    }
    rr <- tune_expand(
      list(node = tk$nodes,
           mtry = unique(pmin(p, ceiling(tk$mtryf * sqrt(p))))),
      rsf_eval,
      list(node = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) v * 2, lo = 1,
                       hi = floor(length(tr) / 4)),
           mtry = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) ceiling(v * 1.75), lo = 1, hi = p)),
      tol = tk$tol, max_ext = tk$maxext)
    if (!is.null(rr$fit)) out$rsfb[te] <- rr$fit$pred
    # oblique RSF, Brier-tuned, grid extended past its edges
    ps <- ncol(Xs)
    otr <- data.frame(time = tt[tr], status = ee[tr], Xs[tr, , drop = FALSE])
    ote <- as.data.frame(Xs[te, , drop = FALSE])
    nb <- max(1L, tk$ntree %/% tk$batch)
    oi <- 0L
    orsf_eval <- function(par) {
      oi <<- oi + 1L
      ro <- list()
      ms <- list()
      for (bb in seq_len(nb)) {
        f <- aorsf::orsf(otr, Surv(time, status) ~ ., n_tree = tk$batch,
                         mtry = par$mtry, leaf_min_obs = par$leaf,
                         split_min_obs = 2L * par$leaf,
                         oobag_pred_type = "mort",
                         tree_seeds = as.integer(tk$seed * 10000 +
                                                   k * 1000 + oi * 20 + bb),
                         n_thread = 1)
        ro[[bb]] <- matrix(predict(f, oobag = TRUE, pred_type = "risk",
                                   pred_horizon = H + 0.5,
                                   boundary_checks = FALSE),
                           ncol = length(H))
        ms[[bb]] <- as.vector(predict(f, new_data = ote,
                                      pred_type = "mort"))
        rm(f)
      }
      invisible(gc(verbose = FALSE))
      br <- sum(colMeans((Ytr - Reduce(`+`, ro) / nb)^2, na.rm = TRUE))
      list(score = br, pred = rowMeans(do.call(cbind, ms)))
    }
    ro <- tune_expand(
      list(leaf = tk$leaf,
           mtry = unique(pmax(1, pmin(ps, ceiling(tk$omtryf * sqrt(ps)))))),
      orsf_eval,
      list(leaf = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) v * 2, lo = 1,
                       hi = floor(length(tr) / 8)),
           mtry = list(down = function(v) max(1, floor(v / 2)),
                       up = function(v) ceiling(v * 1.75), lo = 1,
                       hi = ps)),
      tol = tk$tol, max_ext = tk$maxext)
    if (!is.null(ro$fit)) out$orsfb[te] <- ro$fit$pred
    for (fp in list(list("RSF", rr$path, "node"),
                    list("oblique RSF", ro$path, "leaf"))) {
      pth <- fp[[2]]
      if (is.null(pth) || !nrow(pth)) next
      bad <- pth[pth$error != ""]
      if (!nrow(bad)) next
      fails[[length(fails) + 1]] <- data.table(
        fold = k, family = fp[[1]],
        setting = sprintf("%s %g, mtry %g", fp[[3]], bad[[fp[[3]]]],
                          bad$mtry), message = bad$error)
    }
    prm[[length(prm) + 1]] <- data.table(
      fold = k, ridge_path_passes = pass, spline_k = sk,
      rsf_node = if (is.null(rr$best)) NA else rr$best$node,
      rsf_mtry = if (is.null(rr$best)) NA else rr$best$mtry,
      rsf_edge_node = if (is.null(rr$edge)) NA else rr$edge[["node"]],
      rsf_edge_mtry = if (is.null(rr$edge)) NA else rr$edge[["mtry"]],
      rsf_fails = rr$nfail,
      orsf_leaf = if (is.null(ro$best)) NA else ro$best$leaf,
      orsf_mtry = if (is.null(ro$best)) NA else ro$best$mtry,
      orsf_edge_leaf = if (is.null(ro$edge)) NA else ro$edge[["leaf"]],
      orsf_edge_mtry = if (is.null(ro$edge)) NA else ro$edge[["mtry"]],
      orsf_fails = ro$nfail, p = p)
    if (nrow(rr$ext)) {
      exts[[length(exts) + 1]] <- cbind(data.table(fold = k,
                                                   family = "RSF"), rr$ext)
    }
    if (nrow(ro$ext)) {
      exts[[length(exts) + 1]] <- cbind(data.table(fold = k,
                                                   family = "oblique RSF"),
                                        ro$ext)
    }
  }
  if (!is.null(old)) {
    keep <- if (nrow(old$fails)) old$fails[old$fails$family !=
                                              "spline Cox"] else old$fails
    res <- list(co = tk$co, seed = tk$seed, out = out,
                prm = rbindlist(prm), ext = old$ext,
                fails = rbindlist(list(keep, rbindlist(fails, fill = TRUE)),
                                  fill = TRUE),
                sec = proc.time()[["elapsed"]] - t_task)
  } else {
    res <- list(co = tk$co, seed = tk$seed, out = out,
                prm = rbindlist(prm), ext = rbindlist(exts, fill = TRUE),
                fails = rbindlist(fails, fill = TRUE),
                sec = proc.time()[["elapsed"]] - t_task)
  }
  if (!probe) {
    tmp <- paste0(tk$ckpt, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, tk$ckpt)
  }
  res
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PR, c("frame_B.csv", "folds.csv",
                          "feature_whitelist.csv")),
          file.path(IN, c("ehr28.csv", "ehr_extra.csv", "followup_365.csv",
                          "dx_index.csv")))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fb <- fread(file.path(PR, "frame_B.csv"))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
e28 <- unique(fread(file.path(IN, "ehr28.csv")), by = "hadm_id")
ex <- unique(fread(file.path(IN, "ehr_extra.csv")), by = "hadm_id")
labs <- grep("^lab_", names(ex), value = TRUE)
drop <- labs[grepl(paste(EXCL, collapse = "|"), labs)]
labs <- setdiff(labs, drop)
c28 <- intersect(wl[modality == "ehr", column], names(e28))
xtra <- intersect(EXTRA, names(ex))
chk("EHR 79: registered / extra labs (excluded)",
    sprintf("%d / %d (%s)", length(c28), length(labs),
            paste(drop, collapse = ", ")),
    if (length(labs) >= 40) "PASS" else "FAIL")
chk("block 2, unused measurements found",
    sprintf("%d of %d (%s)", length(xtra), length(EXTRA),
            paste(xtra, collapse = ", ")),
    if (length(xtra) >= 6) "PASS" else "FAIL")
stop_if_fail()
# block 3: discharge destination
fu <- unique(fread(file.path(IN, "followup_365.csv")), by = "hadm_id")
loc <- toupper(trimws(as.character(
  fu$discharge_location[match(fb$hadm_id, fu$hadm_id)])))
loc[is.na(loc) | loc == "NA"] <- ""
DS <- sapply(names(DEST), function(k) as.numeric(loc == k))
colnames(DS) <- DEST
DS <- cbind(DS, dest_other = as.numeric(!(loc %in% c("HOME", names(DEST)))))
if (no_hospice) {
  DS <- DS[, colnames(DS) != "dest_hospice", drop = FALSE]
  chk("A31: hospice destination indicator removed (hospice discharges)",
      sum(loc == "HOSPICE"), "INFO")
}
chk("block 3, discharge destination (home / other share)",
    sprintf("%.1f%% / %.1f%%", 100 * mean(loc == "HOME"),
            100 * mean(DS[, "dest_other"])), "INFO")
# block 4: comorbidity from the index admission's codes
dx <- fread(file.path(IN, "dx_index.csv"), colClasses = list(
  character = "icd_code"))
dx <- dx[hadm_id %in% fb$hadm_id]
dx[, icd_code := gsub("[^A-Za-z0-9]", "", toupper(icd_code))]
cov_dx <- mean(fb$hadm_id %in% dx$hadm_id)
chk("block 4, admissions with diagnosis codes", sprintf("%.1f%%",
                                                       100 * cov_dx),
    if (cov_dx >= 0.99) "PASS" else "FAIL")
stop_if_fail()
cmb <- function(map9, map10) {
  out <- NULL
  for (v in c(9, 10)) {
    d <- as.data.frame(dx[icd_version == v, .(hadm_id, icd_code)])
    if (!nrow(d)) next
    r <- comorbidity::comorbidity(d, id = "hadm_id", code = "icd_code",
                                  map = if (v == 9) map9 else map10,
                                  assign0 = FALSE)
    out <- rbind(out, as.data.frame(r))
  }
  out <- as.data.table(out)[, lapply(.SD, max), by = hadm_id]
  out <- out[match(fb$hadm_id, hadm_id)]
  out[is.na(hadm_id), hadm_id := fb$hadm_id[is.na(out$hadm_id)]]
  for (k in setdiff(names(out), "hadm_id")) {
    set(out, which(is.na(out[[k]])), k, 0)
  }
  out
}
CH <- cmb("charlson_icd9_quan", "charlson_icd10_quan")
EL <- cmb("elixhauser_icd9_quan", "elixhauser_icd10_quan")
sc_ch <- as.numeric(comorbidity::score(
  structure(as.data.frame(CH), map = "charlson_icd10_quan",
            class = c("comorbidity", "data.frame")),
  weights = "charlson", assign0 = TRUE))
sc_el <- as.numeric(comorbidity::score(
  structure(as.data.frame(EL), map = "elixhauser_icd10_quan",
            class = c("comorbidity", "data.frame")),
  weights = "vw", assign0 = TRUE))
ind <- cbind(setNames(as.data.frame(CH)[, -1],
                      paste0("ch_", names(CH)[-1])),
             setNames(as.data.frame(EL)[, -1],
                      paste0("el_", names(EL)[-1])))
pv <- colMeans(ind)
keep <- pv >= PREV[1] & pv <= PREV[2]
CM <- cbind(as.matrix(ind[, keep, drop = FALSE]), charlson_score = sc_ch,
            elixhauser_vw_score = sc_el)
chk("block 4, Charlson / Elixhauser categories kept (1%-99%)",
    sprintf("%d of %d / %d of %d", sum(keep[startsWith(names(keep), "ch_")]),
            ncol(CH) - 1, sum(keep[startsWith(names(keep), "el_")]),
            ncol(EL) - 1), "INFO")
chk("block 4, median Charlson / Elixhauser score",
    sprintf("%g / %g", median(sc_ch), median(sc_el)), "INFO")
FT <- rbind(
  data.table(block = "1 EHR 79", feature = c(c28, labs), prevalence = NA),
  data.table(block = "2 unused measurements", feature = xtra,
             prevalence = NA),
  data.table(block = "3 discharge destination", feature = colnames(DS),
             prevalence = colMeans(DS)),
  data.table(block = "4 comorbidity", feature = names(pv),
             prevalence = pv)[, kept := keep][],
  data.table(block = "4 comorbidity", feature = c("charlson_score",
                                                   "elixhauser_vw_score"),
             prevalence = NA), fill = TRUE)
FT[is.na(kept), kept := TRUE]
fwrite(FT, file.path(RS, paste0(TAG, "_features.csv")))
XALL <- cbind(as.matrix(e28[match(fb$hadm_id, hadm_id), ..c28]),
              as.matrix(ex[match(fb$hadm_id, hadm_id), ..labs]),
              as.matrix(ex[match(fb$hadm_id, hadm_id), ..xtra]),
              DS, CM)
storage.mode(XALL) <- "double"
chk("EHR+ features in total", ncol(XALL), "INFO")
tasks <- list()
META <- list()
for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  X <- XALL[match(sub$hadm_id, fb$hadm_id), , drop = FALSE]
  fo <- folds[frame == "B" & cohort == co]
  seeds <- sort(unique(fo$seed))
  fmap <- list()
  for (s in seeds) {
    f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(f)) {
      chk(sprintf("%s seed %s: folds complete", co, s), sum(is.na(f)),
          "FAIL")
      stop_if_fail()
    }
    fmap[[as.character(s)]] <- as.integer(f)
    tasks[[length(tasks) + 1]] <- list(
      co = co, seed = s, X = X, t = sub$time_death,
      e = as.integer(sub$event_death == 1), subj = sub$subject_id,
      fold = as.integer(f), hz = HZ, nlam = NLAM, lmr = LMR,
      knots = KNOTS, minu = MINU, ntree = NTREE, nsplit = NSPLIT,
      ntime = NTIME, nodes = NODES, mtryf = MTRYF, leaf = LEAF,
      omtryf = OMTRYF, batch = BATCH, tol = TOL, maxext = MAXEXT,
      mode = if (refit_spline) "spline" else "full",
      ckpt_old = file.path(CKPT_OLD, sprintf("B_%s_%d.rds", code, s)),
      ckpt = file.path(CKPT, sprintf("B_%s_%d.rds", code, s)))
  }
  META[[co]] <- list(co = co, code = code, sub = sub, seeds = seeds,
                     fmap = fmap)
  cat(sprintf("  B %-7s n=%d  features %d  seeds %d\n", co, nrow(sub),
              ncol(X), length(seeds)))
}
kind <- sapply(tasks, function(t) t$co)
if (refit_spline) {
  nold <- sum(sapply(tasks, function(t) file.exists(t$ckpt_old)))
  chk("first-run checkpoints found for the spline refit",
      sprintf("%d of %d", nold, length(tasks)),
      if (nold == length(tasks)) "PASS" else "FAIL")
  stop_if_fail()
  cat("\n=== TIMING PROBE (spline refit, fold 0, largest task) ===\n")
  big <- which.max(sapply(tasks, function(t) length(t$t)))
  pr0 <- run_task(tasks[[big]], folds_to_run = 0L, probe = TRUE)
  costs <- rep(pr0$sec * NFOLD, length(tasks))
  cat(sprintf("  %.1f min per fold\n", pr0$sec / 60))
  if (probe_only) {
    cat(sprintf("  %d tasks on %d workers: about %.1f hours\n",
                length(tasks), workers,
                max(sum(costs) / workers, max(costs)) * 1.1 / 3600))
    cat("--probe-only: stopping here (nothing saved).\n")
    quit(status = 0)
  }
} else if (file.exists(PLANF)) {
  PL <- readRDS(PLANF)
  costs <- PL$costs
  cat(sprintf("\n=== RESUMING (plan of %s); probe skipped ===\n", PL$when))
} else {
  cat("\n=== TIMING PROBE (real data, fold 0, in parallel) ===\n")
  reps <- tasks[!duplicated(kind)]
  names(reps) <- unique(kind)
  plan(multisession, workers = min(workers, length(reps)))
  pres <- future_lapply(reps, run_task, folds_to_run = 0L, probe = TRUE,
                        future.packages = PKG, future.seed = TRUE,
                        future.chunk.size = 1L)
  plan(sequential)
  for (kd in names(pres)) {
    cat(sprintf("  %-8s %6.1f min per fold\n", kd, pres[[kd]]$sec / 60))
  }
  costs <- sapply(pres, function(r) r$sec)[kind] * NFOLD
  if (probe_only) {
    cat(sprintf("  %d tasks on %d workers: about %.1f hours\n",
                length(tasks), workers,
                max(sum(costs) / workers, max(costs)) * 1.1 / 3600))
    cat("--probe-only: stopping here (nothing saved).\n")
    quit(status = 0)
  }
  saveRDS(list(costs = costs, when = format(Sys.time(), "%Y-%m-%d %H:%M")),
          PLANF)
}
ok_ck <- function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(r) && all(sapply(r$out, function(v) all(is.finite(v))))
}
done <- which(sapply(tasks, function(t) ok_ck(t$ckpt)))
todo <- setdiff(seq_along(tasks), done)
proj <- max(sum(costs[todo]) / workers, max(c(0, costs[todo]))) * 1.1
cat(sprintf(paste("  %d of %d tasks already saved; %d to fit, about %.1f",
                  "hours on %d workers\n"), length(done), length(tasks),
            length(todo), proj / 3600, workers))
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
if (length(todo)) {
  plan(multisession, workers = workers)
  ord <- todo[order(-costs[todo])]
  invisible(future_lapply(tasks[ord], run_task, future.packages = PKG,
                          future.seed = TRUE, future.chunk.size = 1L))
  plan(sequential)
}
res <- lapply(tasks, function(t) readRDS(t$ckpt))

# ---------------- assemble and check ----------------
cat("\n=== ASSEMBLE AND CHECK ===\n")
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
PRM <- list()
EXT <- list()
FL <- list()
for (co in names(META)) {
  M <- META[[co]]
  rr <- Filter(function(r) r$co == co, res)
  parts <- list()
  nbad <- 0L
  nmiss <- 0L
  for (r in rr) {
    f <- M$fmap[[as.character(r$seed)]]
    d <- data.table(seed = r$seed, hadm_id = M$sub$hadm_id, fold = f,
                    ehrp.cox2.death.lp = r$out$cox2,
                    ehrp.spline.death.lp = r$out$spline,
                    ehrp.rsfb.death.score = r$out$rsfb,
                    ehrp.orsfb.death.score = r$out$orsfb)
    parts[[length(parts) + 1]] <- d
    for (cn in names(d)[-(1:3)]) {
      nmiss <- nmiss + sum(!is.finite(d[[cn]]))
      v <- sapply(sort(unique(f)), function(k) {
        ii <- which(f == k)
        cidx(M$sub$time_death[ii], M$sub$event_death[ii] == 1,
             d[[cn]][ii])
      })
      if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
    }
    PRM[[length(PRM) + 1]] <- cbind(data.table(cohort = co,
                                               seed = r$seed), r$prm)
    if (nrow(r$ext)) {
      EXT[[length(EXT) + 1]] <- cbind(data.table(cohort = co,
                                                 seed = r$seed), r$ext)
    }
    if (nrow(r$fails)) {
      FL[[length(FL) + 1]] <- cbind(data.table(cohort = co,
                                               seed = r$seed), r$fails)
    }
  }
  fwrite(rbindlist(parts), file.path(L1, sprintf("oof_%s_B_%s.csv", TAG,
                                                 M$code)))
  chk(sprintf("%s: missing scores", co), nmiss,
      if (nmiss == 0) "PASS" else "FAIL")
  chk(sprintf("%s: scores failing the per-fold rule (D2)", co), nbad,
      if (nbad == 0) "PASS" else "FAIL")
}
P <- rbindlist(PRM)
XT <- rbindlist(EXT, fill = TRUE)
FLS <- rbindlist(FL, fill = TRUE)
fwrite(P, file.path(L1, paste0("params_", TAG, ".csv")))
nfo <- sum(P$rsf_fails) + sum(P$orsf_fails)
chk("forest fits that failed (best remaining one used)", nfo,
    if (nfo == 0) "PASS" else "WARN")
chk("spline Cox fits: k = 5 / smaller basis / replaced by ridge-Cox",
    sprintf("%d / %d / %d", sum(P$spline_k == KNOTS),
            sum(P$spline_k %in% c(3L, 4L)), sum(P$spline_k == 0L)),
    if (all(P$spline_k == KNOTS)) "PASS" else "WARN")
chk("ridge-Cox: folds whose path was extended",
    sum(P$ridge_path_passes > 1), "INFO")
for (pr in c("rsf_edge_node", "rsf_edge_mtry", "orsf_edge_leaf",
             "orsf_edge_mtry")) {
  sh <- mean(P[[pr]], na.rm = TRUE)
  chk(sprintf("share of folds still at an edge after extension: %s",
              sub("_edge_", " ", pr)), sprintf("%.0f%%", 100 * sh),
      if (sh <= 0.20) "PASS" else "WARN")
}
SUMM <- list()
if (nrow(XT)) {
  SUMM[[1]] <- XT[, .(tried = .N, kept = sum(kept)),
                  by = .(family, param, direction)]
  for (i in seq_len(nrow(SUMM[[1]]))) {
    z <- SUMM[[1]][i]
    chk(sprintf("extensions %s %s %s: tried / kept", z$family, z$param,
                z$direction), sprintf("%d / %d", z$tried, z$kept), "INFO")
  }
}
cs <- rbind(
  P[, .(family = "RSF", setting = sprintf("node %g, mtry %g", rsf_node,
                                          rsf_mtry))],
  P[, .(family = "oblique RSF", setting = sprintf("leaf %g, mtry %g",
                                                  orsf_leaf, orsf_mtry))])
cs <- cs[, .N, by = .(family, setting)][order(family, -N)]
cs[, share := round(N / sum(N), 3), by = family]
fwrite(rbind(cbind(type = "chosen setting", cs[, .(family, setting, N,
                                                    share)]),
             if (length(SUMM)) cbind(type = "extensions",
                                     SUMM[[1]][, .(family,
                                                   setting = paste(param,
                                                                   direction),
                                                   N = tried,
                                                   share = round(kept /
                                                                   tried,
                                                                 3))])),
       file.path(RS, paste0(TAG, "_tuning_summary.csv")))
if (nrow(FLS)) {
  cat("\n  fit failures:\n")
  print(FLS, row.names = FALSE)
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, paste0(TAG, "_checks.csv")))
cat("\n=== CHOSEN FOREST SETTINGS (share of folds) ===\n")
print(cs, row.names = FALSE, nrows = 100)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (co in names(META)) {
  cat("  ", file.path(L1, sprintf("oof_%s_B_%s.csv", TAG, META[[co]]$code)),
      " (local only)\n")
}
cat("  ", file.path(L1, paste0("params_", TAG, ".csv")), " (local only)\n")
for (f in c("_checks.csv", "_tuning_summary.csv", "_features.csv")) {
  cat("  ", file.path(RS, paste0(TAG, f)), "\n")
}
cat("  ", file.path(RS, paste0(TAG, "_log.txt")),
    " (terminal log)\n")
cat(sprintf("\nOutputs: %s_checks.csv, %s_features.csv,", TAG, TAG),
    sprintf("%s_tuning_summary.csv\n", TAG))