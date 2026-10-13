# 06q_level1_optionB.R
# Level-1 scores for the composite and CV readmission outcomes, out of
# fold, in both final pipelines (PREREGISTRATION_ADDENDUM_2.md, A33 Option
# B; A34 for --no-hospice).
#
# Frames:
#   --frame=P  at diagnosis: stage 4r's frame_P_outcomes.csv (clock at
#     admission); EHR at presentation, first ECG and presentation CTPA
#     (stages 4p, 6p); Frame A's folds.
#   --frame=B  at discharge: frame_B.csv (clock at discharge); EHR+, ECG
#     and CTPA as stage 4r wrote them; Frame B's folds.
# Outcomes:
#   --outcome=comp  the composite: time to the first of CV readmission or
#     death (time_first, event_comp), a single event.
#   --outcome=cv  CV readmission first, death competing (event_cr 1 = CV
#     first, 2 = death first). Ridge-Cox, spline Cox and oblique RSF are
#     cause-specific (deaths censored at death). The RSF is a
#     competing-risks forest (randomForestSRC, status 0/1/2), Brier-tuned
#     on its out-of-bag cumulative incidence of CV readmission at each
#     horizon; its test score is its CV-readmission prediction.
# Cohorts are EHR+ECG and three-modality, with horizons of 10, 30, 90, 180
# and 365 days. Families, tuning, checks, timing probe and checkpointing
# are as in 06p: ridge-Cox retuned; spline Cox (EHR only, D16); RSF and
# oblique RSF, Brier-tuned with tune_expand.R against the outcome by each
# horizon (either event for the composite, CV readmission first for CV).
#
# --no-hospice (A34, --frame=B only) refits the EHR modality alone on EHR+
# without the hospice destination indicator. ECG and CTPA are unaffected
# and reused from the main runs. Files are tagged Bcomp_nh or Bcv_nh.
#
# Inputs:  data/processed/presentation/frame_P_outcomes.csv, ehrP.csv,
#          ecgP.csv, ctpaP.csv (frame P); data/processed/frame_B.csv and
#          discharge/ehrplusB.csv, ecgB.csv, ctpaB.csv (frame B);
#          data/processed/folds.csv, feature_whitelist.csv;
#          scripts/tune_expand.R
# Outputs: results/<TAG>_level1_checks.csv,
#          results/<TAG>_level1_tuning_summary.csv,
#          results/<TAG>_level1_log.txt
#          data/processed/level1/oof_<TAG>_<cohort>.csv (patient level, not
#          in the repo), with TAG = Pcomp, Pcv, Bcomp or Bcv (plus _nh)
# Usage:   Rscript scripts/06q_level1_optionB.R --frame=P|B
#          --outcome=comp|cv [--no-hospice] [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(mgcv)
  library(data.table)
  library(future.apply)
  library(randomForestSRC)
})
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf not installed")
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
geta <- function(k, ok) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args,
                                           value = TRUE))
  if (length(v) != 1 || !(v %in% ok)) {
    stop(sprintf("give --%s= one of %s", k, paste(ok, collapse = ", ")))
  }
  v
}
FR <- geta("frame", c("P", "B"))
OC <- geta("outcome", c("comp", "cv"))
no_hosp <- "--no-hospice" %in% args
if (no_hosp && FR != "B") stop("--no-hospice applies to --frame=B only")
TAG <- paste0(FR, OC, if (no_hosp) "_nh" else "")
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
CKPT <- file.path(L1, paste0(TAG, "_ckpt"))
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
PLANF <- file.path(CKPT, "plan.rds")
source(file.path(P2, "scripts", "tune_expand.R"))
HZ <- c(10, 30, 90, 180, 365)
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
    fwrite(ck, file.path(RS, paste0(TAG, "_level1_checks.csv")))
    cat(sprintf("\nSTOP: a check failed. see: %s_level1_checks.csv\n",
                TAG))
    quit(status = 1)
  }
}

# ---------------- one task = cohort x seed x modality ----------------
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
    # spline Cox (D14's fallback rule); linear terms ridge-penalised (D16);
    # fitted for the EHR only (ECG and CTPA use ridge-Cox, as in v4)
    nu <- apply(Xs[tr, , drop = FALSE], 2, function(v) length(unique(v)))
    if (!tk$spline) nu[] <- 0L
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
    for (kk in if (tk$spline) unique(c(tk$knots, 4L, 3L)) else integer(0)) {
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
    if (!tk$spline) sk <- NA_integer_
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
    st_tr <- if (is.null(tk$cr)) ee[tr] else tk$cr[tr]
    dtr <- data.frame(time = tt[tr], status = st_tr, X[tr, , drop = FALSE])
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
        if (i == 0) return(rep(0, length(tr)))
        # competing risks (CV): out-of-bag CIF of CV readmission (cause 1)
        if (is.null(tk$cr)) 1 - f$survival.oob[, i] else f$cif.oob[, i, 1]
      })
      br <- sum(colMeans((Ytr - matrix(ro, ncol = length(H)))^2))
      pp <- predict(f, newdata = dte)$predicted
      sc <- if (is.null(tk$cr)) pp else pp[, 1]
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
    res <- list(co = tk$co, seed = tk$seed, mod = tk$mod, out = out,
                prm = rbindlist(prm), ext = old$ext,
                fails = rbindlist(list(keep, rbindlist(fails, fill = TRUE)),
                                  fill = TRUE),
                sec = proc.time()[["elapsed"]] - t_task)
  } else {
    res <- list(co = tk$co, seed = tk$seed, mod = tk$mod, out = out,
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
cat(sprintf("=== BUILD (frame %s, outcome %s) ===\n", FR, OC))
DIR <- if (FR == "P") PP else file.path(PR, "discharge")
fframe <- if (FR == "P") file.path(PP, "frame_P_outcomes.csv") else
  file.path(PR, "frame_B.csv")
fm <- if (FR == "P") {
  c(ehr = "ehrP.csv", ecg = "ecgP.csv", ctpa = "ctpaP.csv")
} else {
  c(ehr = "ehrplusB.csv", ecg = "ecgB.csv", ctpa = "ctpaB.csv")
}
need <- c(fframe, file.path(DIR, fm),
          file.path(PR, c("folds.csv", "feature_whitelist.csv")))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- fread(fframe)
if (FR == "B") fp[, in_three_mod_p := in_three_mod]
folds <- fread(file.path(PR, "folds.csv"))[frame == (if (FR == "P") "A"
                                                     else "B")]
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mat <- function(f, rows, keep = NULL) {
  D <- unique(fread(file.path(DIR, f)), by = "hadm_id")
  cc <- setdiff(names(D), c("hadm_id", "subject_id"))
  if (!is.null(keep)) cc <- intersect(keep, cc)
  if (no_hosp) cc <- setdiff(cc, "dest_hospice")
  X <- as.matrix(D[match(fp$hadm_id[rows], hadm_id), ..cc])
  storage.mode(X) <- "double"
  list(X = X, miss = sum(is.na(match(fp$hadm_id[rows], D$hadm_id))))
}
ok_ev <- all(c("time_first", "event_cr", "event_comp") %in% names(fp))
chk("first-event columns present (time_first, event_cr, event_comp)",
    if (ok_ev) "yes" else "no", if (ok_ev) "PASS" else "FAIL")
stop_if_fail()
fp[, e_out := if (OC == "comp") as.integer(event_comp == 1) else
  as.integer(event_cr == 1)]
COH <- list("EHR+ECG" = list(code = "ehrecg", mods = c("ehr", "ecg"),
                             rows = rep(TRUE, nrow(fp))),
            "3-mod" = list(code = "3mod", mods = c("ehr", "ecg", "ctpa"),
                           rows = fp$in_three_mod_p == 1))
if (no_hosp) for (co in names(COH)) COH[[co]]$mods <- "ehr"
tasks <- list()
META <- list()
for (co in names(COH)) {
  C <- COH[[co]]
  sub <- fp[C$rows]
  MX <- list()
  for (m in C$mods) {
    mm <- mat(fm[[m]], C$rows, if (m == "ctpa") wl[modality == "ctpa",
                                                   column] else NULL)
    if (no_hosp && m == "ehr") {
      chk(sprintf("%s EHR+: hospice indicator removed", co),
          if ("dest_hospice" %in% colnames(mm$X)) "no" else "yes",
          if ("dest_hospice" %in% colnames(mm$X)) "FAIL" else "PASS")
    }
    chk(sprintf("%s %s: features (admissions missing)", co, toupper(m)),
        sprintf("%d (%d)", ncol(mm$X), mm$miss),
        if (ncol(mm$X) > 0 && mm$miss == 0) "PASS" else "FAIL")
    MX[[m]] <- mm$X
  }
  stop_if_fail()
  fo <- folds[cohort == co]
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
    for (m in C$mods) {
      tasks[[length(tasks) + 1]] <- list(
        co = co, seed = s, mod = m, X = MX[[m]], t = sub$time_first,
        e = sub$e_out, cr = if (OC == "cv") as.integer(sub$event_cr) else
          NULL, subj = sub$subject_id,
        fold = as.integer(f), hz = HZ, nlam = NLAM, lmr = LMR,
        knots = KNOTS, minu = MINU, ntree = NTREE, nsplit = NSPLIT,
        ntime = NTIME, nodes = NODES, mtryf = MTRYF, leaf = LEAF,
        omtryf = OMTRYF, batch = BATCH, tol = TOL, maxext = MAXEXT,
        spline = (m == "ehr"), mode = "full", ckpt_old = NA_character_,
        ckpt = file.path(CKPT, sprintf("%s_%s_%d.rds", C$code, m, s)))
    }
  }
  META[[co]] <- list(co = co, code = C$code, sub = sub, seeds = seeds,
                     fmap = fmap, mods = C$mods)
  ev <- sapply(HZ, function(h) sum(sub$e_out == 1 &
                                     sub$time_first <= h + 0.5))
  cat(sprintf("  %-7s n=%d  events by %s days: %s\n", co, nrow(sub),
              paste(HZ, collapse = "/"), paste(ev, collapse = "/")))
  chk(sprintf("%s: %s events by 10 / 365 days", co, OC),
      sprintf("%d / %d", ev[1], ev[length(ev)]), "INFO")
}
kind <- sapply(tasks, function(t) paste(t$co, t$mod))
if (file.exists(PLANF)) {
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
    cat(sprintf("  %-14s %6.1f min per fold\n", kd, pres[[kd]]$sec / 60))
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
  parts <- list()
  nbad <- 0L
  nmiss <- 0L
  for (s in M$seeds) {
    f <- M$fmap[[as.character(s)]]
    d <- data.table(seed = s, hadm_id = M$sub$hadm_id, fold = f)
    for (r in Filter(function(r) r$co == co && r$seed == s, res)) {
      m <- r$mod
      d[[sprintf("%s.cox2.%s.lp", m, OC)]] <- r$out$cox2
      if (m == "ehr") {
        d[[sprintf("%s.spline.%s.lp", m, OC)]] <- r$out$spline
      }
      d[[sprintf("%s.rsfb.%s.score", m, OC)]] <- r$out$rsfb
      d[[sprintf("%s.orsfb.%s.score", m, OC)]] <- r$out$orsfb
      PRM[[length(PRM) + 1]] <- cbind(data.table(cohort = co, seed = s,
                                                 modality = m), r$prm)
      if (nrow(r$ext)) {
        EXT[[length(EXT) + 1]] <- cbind(data.table(cohort = co, seed = s,
                                                   modality = m), r$ext)
      }
      if (nrow(r$fails)) {
        FL[[length(FL) + 1]] <- cbind(data.table(cohort = co, seed = s,
                                                 modality = m), r$fails)
      }
    }
    for (cn in names(d)[-(1:3)]) {
      nmiss <- nmiss + sum(!is.finite(d[[cn]]))
      v <- sapply(sort(unique(f)), function(k) {
        ii <- which(f == k)
        cidx(M$sub$time_first[ii], M$sub$e_out[ii] == 1, d[[cn]][ii])
      })
      if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
    }
    parts[[length(parts) + 1]] <- d
  }
  fwrite(rbindlist(parts), file.path(L1, sprintf("oof_%s_%s.csv", TAG,
                                                 M$code)))
  chk(sprintf("%s: missing scores", co), nmiss,
      if (nmiss == 0) "PASS" else "FAIL")
  chk(sprintf("%s: scores failing the per-fold rule (D2)", co), nbad,
      if (nbad == 0) "PASS" else "FAIL")
}
P <- rbindlist(PRM, fill = TRUE)
XT <- rbindlist(EXT, fill = TRUE)
FLS <- rbindlist(FL, fill = TRUE)
fwrite(P, file.path(L1, paste0("params_", TAG, ".csv")))
nfo <- sum(P$rsf_fails) + sum(P$orsf_fails)
chk("forest fits that failed (best remaining one used)", nfo,
    if (nfo == 0) "PASS" else "WARN")
pe <- P[modality == "ehr"]
chk("EHR spline Cox fits: k = 5 / smaller basis / replaced by ridge-Cox",
    sprintf("%d / %d / %d", sum(pe$spline_k == KNOTS),
            sum(pe$spline_k %in% c(3L, 4L)), sum(pe$spline_k == 0L)),
    if (all(pe$spline_k == KNOTS)) "PASS" else "WARN")
chk("ridge-Cox: folds whose path was extended",
    sum(P$ridge_path_passes > 1), "INFO")
for (pr in c("rsf_edge_node", "rsf_edge_mtry", "orsf_edge_leaf",
             "orsf_edge_mtry")) {
  sh <- mean(P[[pr]], na.rm = TRUE)
  chk(sprintf("share of folds still at an edge after extension: %s",
              sub("_edge_", " ", pr)), sprintf("%.0f%%", 100 * sh),
      if (sh <= 0.20) "PASS" else "WARN")
}
cs <- rbind(
  P[, .(family = "RSF", modality, setting = sprintf("node %g, mtry %g",
                                                    rsf_node, rsf_mtry))],
  P[, .(family = "oblique RSF", modality,
        setting = sprintf("leaf %g, mtry %g", orsf_leaf, orsf_mtry))])
cs <- cs[, .N, by = .(family, modality, setting)][order(family, modality,
                                                        -N)]
cs[, share := round(N / sum(N), 3), by = .(family, modality)]
ex2 <- if (nrow(XT)) {
  XT[, .(N = .N, share = round(mean(kept), 3)),
     by = .(family, modality, setting = paste(param, direction))]
} else NULL
fwrite(rbind(cbind(type = "chosen setting", cs),
             if (!is.null(ex2)) cbind(type = "extensions", ex2),
             fill = TRUE),
       file.path(RS, paste0(TAG, "_level1_tuning_summary.csv")))
if (nrow(FLS)) {
  cat("\n  fit failures:\n")
  print(FLS, row.names = FALSE)
}
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, paste0(TAG, "_level1_checks.csv")))
cat("\n=== CHOSEN FOREST SETTINGS (share of folds) ===\n")
print(cs[share >= 0.2], row.names = FALSE, nrows = 100)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f hours\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
cat("\nOUTPUT FILES\n")
for (co in names(META)) {
  cat("  ", file.path(L1, sprintf("oof_%s_%s.csv", TAG, META[[co]]$code)),
      " (local only)\n")
}
cat("  ", file.path(L1, paste0("params_", TAG, ".csv")), " (local only)\n")
cat("  ", file.path(RS, paste0(TAG, "_level1_checks.csv")), "\n")
cat("  ", file.path(RS, paste0(TAG, "_level1_tuning_summary.csv")), "\n")
cat("  ", file.path(RS, paste0(TAG, "_level1_log.txt")),
    " (terminal log)\n")
cat(sprintf("\nOutputs: %s_level1_checks.csv, %s_level1_tuning_summary.csv\n",
            TAG, TAG))