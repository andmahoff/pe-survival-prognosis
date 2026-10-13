# 12_sensitivity.R
# Registered sensitivity analyses not run elsewhere (PREREGISTRATION.md
# section 8; PREREGISTRATION_ADDENDUM_2.md A25).
#
# Part 1, final data year (Frame B, CV readmission and composite).
# Patients whose index admission may fall in MIMIC's final data year
# (final_year_possible, stage 4) are excluded, since a later CV
# readmission would be missed. The registered out-of-fold predictions of
# S2's four CV formulations (stage 9) and the composite models of stage
# 11b (Cox stack, binary stack, recalibrated late fusion where present)
# are re-scored on the rest with the models unchanged (AUROC, Brier score,
# observed/expected, contrasts), beside the full cohort.
#
# Part 2, all admissions (Frame B, primary analysis). Frame B is rebuilt
# with every admission of the EHR+ECG cohort rather than each patient's
# first, with stage 4's exclusions and time coding. EHR and ECG death
# scores are ridge-Cox retuned (D11; stage 18a preprocessing and inner
# folds), out of fold over 5 folds grouped by patient and stratified by
# death, for the six stage 4 seeds. The primary estimand (ECG hazard ratio
# in days 0-30 over days 31-365) is shown beside the same estimand in the
# first-admission frame on the stage 6d scores, so only the cohort differs.
#
# Part 3, 48-hour landmark (Frame A). Patients with no death and no CV
# readmission on calendar days 0 or 1 (event times above 2) are kept, with
# outcomes by day 30 as registered. The S1 scores (stage 8 fused, stage 6a
# per modality) are re-scored on them as AUROC at 30.5 days and Uno's C to
# 31 days, with the S1 contrasts, beside the full Frame A cohort.
#
# Bootstrap: 2,000 samples, seed 42, clustered by patient in Parts 2 and 3.
# Checks (the script stops on failure): inputs present, rebuilt first
# admissions equal frame_B.csv, no missing scores, D2 per-fold rule on the
# all-admission scores. Version 2: Part 1 uses the composite models in the
# stage 11b file, which must include the Cox and binary stacks.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, s2/s2_preds_B_*.csv,
#          calibration/cal_preds_B_*.csv, fusion/fused_A_*.csv,
#          level1/oof_linear_A_*.csv, level1/oof_extra_B_ehrecg.csv,
#          data/inputs/cohorts.csv, mimic_labels_harmonised.csv,
#          index_admission_times.csv, ecg_timing.csv, followup_365.csv,
#          ehr28.csv, ecg_opt.csv
# Outputs: results/sens_final_year.csv, sens_all_admissions.csv,
#          sens_landmark.csv, sens_checks.csv, sens_log.txt (terminal log)
# Usage:   Rscript scripts/12_sensitivity.R [--probe-only] [--workers=N]

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
RS <- file.path(P2, "results")
HB <- c(10, 30, 90, 180, 365)
TAU_B <- 365
CUT <- 30.5
LMK <- 2
NLAM <- 80L
LMR <- 1e-7
NBOOT <- 2000L
SEED <- 42L
PKG <- c("survival", "data.table")
wk_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(wk_arg)) {
  as.integer(sub("^--workers=", "", wk_arg))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
options(width = 220)
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
    fwrite(ck, file.path(RS, "sens_checks.csv"))
    cat("\nSTOP: a check failed. see: sens_checks.csv\n")
    quit(status = 1)
  }
}
need <- c(file.path(PR, c("frame_A.csv", "frame_B.csv", "folds.csv",
                          "feature_whitelist.csv",
                          "level1/oof_extra_B_ehrecg.csv")),
          file.path(IN, c("cohorts.csv", "mimic_labels_harmonised.csv",
                          "index_admission_times.csv", "ecg_timing.csv",
                          "followup_365.csv", "ehr28.csv", "ecg_opt.csv")))
for (co in c("ehrecg", "3mod")) {
  need <- c(need, file.path(PR, c(sprintf("s2/s2_preds_B_%s.csv", co),
                                  sprintf("calibration/cal_preds_B_%s.csv",
                                          co),
                                  sprintf("fusion/fused_A_%s.csv", co),
                                  sprintf("level1/oof_linear_A_%s.csv",
                                          co))))
}
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
auc <- function(y, s) {
  n1 <- sum(y)
  n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- data.table::frank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
unoc <- function(tt, ee, s, ymax) {
  as.numeric(tryCatch(survival::concordance(
    survival::Surv(tt, ee) ~ s, reverse = TRUE, timewt = "n/G2",
    ymax = ymax)$concordance, error = function(e) NA_real_))
}
avg_cols <- function(path, cols, keys) {
  x <- fread(path, select = c("hadm_id", cols))
  x[, lapply(.SD, mean), by = hadm_id, .SDcols = cols][
    match(keys, hadm_id)]
}
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
fa <- fread(file.path(PR, "frame_A.csv"))
fb <- fread(file.path(PR, "frame_B.csv"))
JOBS <- list()

# ================= PART 1: final data year =================
cat("\n=== PART 1: FINAL DATA YEAR (Frame B) ===\n")
if (!("final_year_possible" %in% names(fb))) {
  chk("frame_B.csv has final_year_possible", "no", "FAIL")
  stop_if_fail()
}
for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fb else fb[in_three_mod == 1]
  keep <- sub$final_year_possible == 0
  chk(sprintf("%s: patients kept / removed (CV events by 365 removed)", co),
      sprintf("%d / %d (%d)", sum(keep), sum(!keep),
              sum(!keep & sub$event_cr == 1 & sub$time_first <= 365.5)),
      "INFO")
  s2c <- as.vector(outer(c("F1", "F2", "F3", "F4"), HB, function(f, h)
    sprintf("%s_h%d", f, h)))
  S2 <- avg_cols(file.path(PR, "s2", sprintf("s2_preds_B_%s.csv", code)),
                 s2c, sub$hadm_id)
  cpath <- file.path(PR, "calibration", sprintf("cal_preds_B_%s.csv", code))
  ccl0 <- as.vector(outer(c("coxstack", "binstack", "latefusion"),
                          c(30, 365), function(m, h)
                            sprintf("composite|%s|%d", m, h)))
  ccl <- intersect(ccl0, names(fread(cpath, nrows = 1, sep = ",")))
  absent <- unique(sub("^composite\\|([a-z]+)\\|.*$", "\\1",
                       setdiff(ccl0, ccl)))
  need_c <- sprintf("composite|%s|%d", rep(c("coxstack", "binstack"), 2),
                    rep(c(30, 365), each = 2))
  chk(sprintf("%s: composite models in stage 11b's file", co),
      if (length(absent)) paste("absent:", paste(absent, collapse = ", "))
      else "all", if (all(need_c %in% ccl)) "PASS" else "FAIL")
  stop_if_fail()
  CA <- avg_cols(cpath, ccl, sub$hadm_id)
  nb <- sum(is.na(as.matrix(S2[, ..s2c]))) +
    sum(is.na(as.matrix(CA[, ..ccl])))
  chk(sprintf("%s: missing S2 or composite predictions", co), nb,
      if (nb == 0) "PASS" else "FAIL")
  stop_if_fail()
  ycv <- lapply(HB, function(h) sub$event_cr == 1 & sub$time_first <=
                  h + 0.5)
  dfi <- lapply(HB, function(h) sub$event_cr == 2 & sub$time_first <=
                  h + 0.5)
  yco <- lapply(c(30, 365), function(h) sub$event_comp == 1 &
                  sub$time_first <= h + 0.5)
  names(ycv) <- names(dfi) <- HB
  names(yco) <- c(30, 365)
  S2l <- as.list(S2[, ..s2c])
  CAl <- as.list(CA[, ..ccl])
  met <- local({
    # local copies travel with the closure to the parallel workers
    HB <- HB
    auc <- auc
    ycv <- ycv
    dfi <- dfi
    yco <- yco
    S2l <- S2l
    CAl <- CAl
    function(idx) {
      v <- numeric(0)
      for (h in HB) {
        hk <- as.character(h)
        y <- ycv[[hk]][idx]
        kp <- !dfi[[hk]][idx]
        for (f in c("F1", "F2", "F3", "F4")) {
          p <- S2l[[sprintf("%s_h%d", f, h)]][idx]
          v[paste("S2", h, f, "auc_full", sep = "|")] <- auc(y, p)
          v[paste("S2", h, f, "auc_cond", sep = "|")] <- auc(y[kp], p[kp])
          v[paste("S2", h, f, "brier", sep = "|")] <- mean((y - p)^2)
          v[paste("S2", h, f, "oe", sep = "|")] <- mean(y) / mean(p)
        }
      }
      for (h in c(30, 365)) {
        y <- yco[[as.character(h)]][idx]
        pr <- mean(y)
        for (m in c("coxstack", "binstack", "latefusion")) {
          if (is.null(CAl[[sprintf("composite|%s|%d", m, h)]])) next
          p <- CAl[[sprintf("composite|%s|%d", m, h)]][idx]
          v[paste("composite", h, m, "auc", sep = "|")] <- auc(y, p)
          v[paste("composite", h, m, "ipa", sep = "|")] <-
            1 - mean((y - p)^2) / (pr * (1 - pr))
          v[paste("composite", h, m, "oe", sep = "|")] <- pr / mean(p)
        }
      }
      v
    }
  })
  prs <- list()
  for (h in HB) {
    for (pp in list(c("F2", "F1"), c("F3", "F1"), c("F4", "F1"),
                    c("F4", "F3"))) {
      for (ms in c("auc_full", "auc_cond", "brier")) {
        prs[[length(prs) + 1]] <- c(paste("S2", h, pp[1], ms, sep = "|"),
                                    paste("S2", h, pp[2], ms, sep = "|"))
      }
    }
  }
  for (h in c(30, 365)) {
    for (pp in list(c("coxstack", "binstack"),
                    c("coxstack", "latefusion"))) {
      for (ms in c("auc", "ipa")) {
        prs[[length(prs) + 1]] <- c(
          paste("composite", h, pp[1], ms, sep = "|"),
          paste("composite", h, pp[2], ms, sep = "|"))
      }
    }
  }
  for (ss in c("full cohort", "final year excluded")) {
    rows <- if (ss == "full cohort") seq_len(nrow(sub)) else which(keep)
    JOBS[[length(JOBS) + 1]] <- list(
      part = "final_year", cohort = co, subset = ss, met = met,
      draw = local({
        r <- rows
        function() sample(r, length(r), replace = TRUE)
      }), all = rows, pairs = prs)
  }
}

# ================= PART 2: all admissions =================
cat("\n=== PART 2: ALL ADMISSIONS (Frame B, the primary) ===\n")
flt <- function(x) suppressWarnings(as.numeric(x))
coh <- fread(file.path(IN, "cohorts.csv"))
lab <- fread(file.path(IN, "mimic_labels_harmonised.csv"))
for (cc in c("dod_days", "inhosp_days", "cv_days")) lab[[cc]] <- flt(
  lab[[cc]])
adm <- unique(fread(file.path(IN, "index_admission_times.csv")),
              by = "hadm_id")
adm[, at := as.POSIXct(admittime, tz = "UTC")]
adm[, dt := as.POSIXct(dischtime, tz = "UTC")]
adm[, los_d := as.numeric(difftime(as.Date(dt), as.Date(at),
                                   units = "days"))]
adm[, los_h := as.numeric(difftime(dt, at, units = "hours"))]
ecgt <- fread(file.path(IN, "ecg_timing.csv"))
fu <- fread(file.path(IN, "followup_365.csv"))
fu[, cv_days_all := flt(cv_days_all)]
d <- merge(coh, lab, by = c("subject_id", "hadm_id"), all.x = TRUE)
d <- merge(d, adm[, .(hadm_id, at, los_d, los_h)], by = "hadm_id",
           all.x = TRUE)
d <- merge(d, ecgt[, .(hadm_id, ecg_h_last)], by = "hadm_id", all.x = TRUE)
d <- merge(d, fu[, .(hadm_id, discharge_location, cv_days_all)],
           by = "hadm_id", all.x = TRUE)
A <- d[in_ehr_ecg == 1][order(subject_id, at, hadm_id)]
A[, first_adm := as.integer(!duplicated(subject_id))]
r_ih <- !is.na(A$inhosp_days)
r_early <- !is.na(A$dod_days) & A$dod_days < A$los_d
pre <- A[!r_ih & !r_early]
r_died <- toupper(as.character(pre$discharge_location)) %in% "DIED"
r_ecg <- !is.na(pre$ecg_h_last) & pre$ecg_h_last > pre$los_h
BA <- pre[!r_died & !r_ecg]
code_t <- function(t_d, t_cv, tau) {
  de <- !is.na(t_d) & t_d >= 0 & t_d <= tau
  cv <- !is.na(t_cv) & t_cv >= 0 & t_cv <= tau
  cvf <- cv & (!de | (t_cv <= t_d))
  dthf <- de & !cvf
  list(time_death = ifelse(de, t_d + 0.5, tau + 1),
       event_death = as.integer(de),
       time_first = ifelse(cvf, t_cv + 0.5, ifelse(dthf, t_d + 0.5,
                                                   tau + 1)))
}
ct <- code_t(BA$dod_days - BA$los_d, BA$cv_days_all - BA$los_d, TAU_B)
BA[, `:=`(time_death = ct$time_death, event_death = ct$event_death)]
chk("all admissions: no follow-up row (discharge location unknown)",
    sum(is.na(BA$discharge_location)), "INFO")
f1 <- BA[first_adm == 1]
same <- nrow(f1) == nrow(fb) && setequal(f1$hadm_id, fb$hadm_id)
if (same) {
  m <- match(fb$hadm_id, f1$hadm_id)
  same <- isTRUE(all.equal(f1$time_death[m], fb$time_death)) &&
    identical(as.integer(f1$event_death[m]), as.integer(fb$event_death))
}
chk("rebuilt first admissions reproduce frame_B.csv (rows, times)",
    sprintf("%d vs %d, %s", nrow(f1), nrow(fb), same),
    if (same) "PASS" else "FAIL")
stop_if_fail()
chk("all-admission Frame B: admissions / patients / deaths by 365",
    sprintf("%d / %d / %d", nrow(BA), uniqueN(BA$subject_id),
            sum(BA$event_death)), "INFO")
wl <- fread(file.path(PR, "feature_whitelist.csv"))
folds <- fread(file.path(PR, "folds.csv"))
SEEDS <- sort(unique(folds$seed))
XM <- lapply(c(ehr = "ehr28.csv", ecg = "ecg_opt.csv"), function(f) {
  m <- if (f == "ehr28.csv") "ehr" else "ecg"
  M <- unique(fread(file.path(IN, f)), by = "hadm_id")
  cc <- intersect(wl[modality == m, column], names(M))
  X <- as.matrix(M[match(BA$hadm_id, hadm_id), ..cc])
  storage.mode(X) <- "double"
  X
})
gfold <- function(subj, ev, seed, k = 5) {
  us <- unique(subj)
  pos <- tapply(ev, subj, max)[as.character(us)]
  set.seed(seed)
  f <- integer(length(us))
  for (g in c(0, 1)) {
    ii <- which(pos == g)
    if (length(ii)) f[ii] <- sample(rep_len(seq_len(k), length(ii)))
  }
  f[match(subj, us)]
}
ridge_task <- function(tk) {
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
  for (k in sort(unique(tk$fold))) {
    tr <- which(tk$fold != k)
    te <- which(tk$fold == k)
    X <- prep(tk$X, tr)
    cv <- glmnet::cv.glmnet(X[tr, , drop = FALSE],
                            survival::Surv(tk$t[tr], tk$e[tr]),
                            family = "cox", alpha = 0,
                            foldid = grouped(tr, tk$seed * 100 + k, 3),
                            nlambda = tk$nlam, lambda.min.ratio = tk$lmr,
                            type.measure = "deviance",
                            penalty.factor = rep(1, ncol(X)),
                            standardize = FALSE)
    lp[te] <- as.vector(predict(cv, X[te, , drop = FALSE],
                                s = "lambda.min"))
  }
  list(mod = tk$mod, seed = tk$seed, lp = lp)
}
tasks <- list()
FM <- list()
for (s in SEEDS) {
  FM[[as.character(s)]] <- gfold(BA$subject_id, BA$event_death, s)
  for (m in names(XM)) {
    tasks[[length(tasks) + 1]] <- list(
      mod = m, seed = s, X = XM[[m]], t = BA$time_death,
      e = BA$event_death, subj = BA$subject_id,
      fold = FM[[as.character(s)]], nlam = NLAM, lmr = LMR)
  }
}
t0 <- Sys.time()
invisible(ridge_task(tasks[[2]]))
tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  level-1 probe: %.0f s per task; %d tasks\n", tp,
            length(tasks)))
ratio_vec <- function(R, w) {
  k <- w[R$id] > 0
  X <- cbind(R$ze * R$early, R$ze * (1 - R$early), R$zc * R$early,
             R$zc * (1 - R$early))[k, , drop = FALSE]
  b <- tryCatch(coef(coxph(Surv(R$tstart[k], R$tstop[k], R$ev[k]) ~ X,
                           weights = w[R$id][k], ties = "efron")),
                error = function(e) rep(NA_real_, 4))
  c(ecg_early = b[[3]], ecg_late = b[[4]], ecg_ratio = b[[3]] - b[[4]],
    ehr_ratio = b[[1]] - b[[2]])
}
make_R <- function(t, e, ze, zc) {
  z <- function(v) (v - mean(v)) / sd(v)
  ep <- survSplit(Surv(t, e) ~ ., cut = CUT, episode = "per",
                  data = data.frame(t = t, e = e, id = seq_along(t)))
  list(ze = z(ze)[ep$id], zc = z(zc)[ep$id], early = as.numeric(ep$per ==
                                                                   1),
       tstart = ep$tstart, tstop = ep$t, ev = ep$e, id = ep$id)
}

# ================= PART 3: 48-hour landmark =================
cat("\n=== PART 3: 48-HOUR LANDMARK (Frame A) ===\n")
for (co in c("EHR+ECG", "3-mod")) {
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") fa else fa[in_three_mod == 1]
  E <- list(
    death = list(cox = "cox.death.lp", lr = "lr.death_h30",
                 t = sub$time_death, e = sub$event_death == 1),
    composite = list(cox = "cox.composite.lp", lr = "lr.composite_h30",
                     t = sub$time_first, e = sub$event_comp == 1),
    cv = list(cox = "cox.cs_cv.lp", lr = "lr.cvu_h30", t = sub$time_first,
              e = sub$event_cr == 1))
  fcols <- as.vector(outer(c("lr", "cox"), names(E), function(l1, en)
    sprintf("fuse.%s.auc.%s", l1, en)))
  ucols <- unlist(lapply(names(E), function(en)
    sprintf("ehr.%s", c(E[[en]]$cox, E[[en]]$lr))))
  FU <- avg_cols(file.path(PR, "fusion", sprintf("fused_A_%s.csv", code)),
                 fcols, sub$hadm_id)
  UL <- avg_cols(file.path(PR, "level1", sprintf("oof_linear_A_%s.csv",
                                                 code)), ucols, sub$hadm_id)
  nb <- sum(is.na(as.matrix(FU[, ..fcols]))) +
    sum(is.na(as.matrix(UL[, ..ucols])))
  chk(sprintf("%s: missing S1 scores", co), nb,
      if (nb == 0) "PASS" else "FAIL")
  stop_if_fail()
  lmk <- sub$time_death > LMK & sub$time_first > LMK
  chk(sprintf("%s: admissions kept / removed by the landmark", co),
      sprintf("%d / %d (deaths %d, CV %d removed)", sum(lmk), sum(!lmk),
              sum(!lmk & sub$event_death == 1),
              sum(!lmk & sub$event_cr == 1)), "INFO")
  SC <- c(as.list(FU[, ..fcols]), as.list(UL[, ..ucols]))
  met <- local({
    auc <- auc
    unoc <- unoc
    E <- E
    SC <- SC
    function(idx) {
      v <- numeric(0)
      for (en in names(E)) {
        y <- (E[[en]]$e & E[[en]]$t <= 30.5)[idx]
        tt <- E[[en]]$t[idx]
        ee <- E[[en]]$e[idx]
        sc <- c("fused lr" = sprintf("fuse.lr.auc.%s", en),
                "fused cox" = sprintf("fuse.cox.auc.%s", en),
                "EHR lr" = sprintf("ehr.%s", E[[en]]$lr),
                "EHR cox" = sprintf("ehr.%s", E[[en]]$cox))
        for (nm in names(sc)) {
          s <- SC[[sc[[nm]]]][idx]
          v[paste("S1", en, nm, "auc30", sep = "|")] <- auc(y, s)
          v[paste("S1", en, nm, "unoC", sep = "|")] <- unoc(tt, ee, s, 31)
        }
      }
      v
    }
  })
  us <- unique(sub$subject_id)
  for (ss in c("full cohort", "48-hour landmark")) {
    rows <- if (ss == "full cohort") seq_len(nrow(sub)) else which(lmk)
    rb <- split(rows, match(sub$subject_id[rows],
                            unique(sub$subject_id[rows])))
    prs <- list()
    for (en in names(E)) for (ms in c("auc30", "unoC")) {
      K <- function(nm) paste("S1", en, nm, ms, sep = "|")
      prs <- c(prs, list(c(K("fused lr"), K("fused cox")),
                         c(K("fused cox"), K("EHR cox")),
                         c(K("fused lr"), K("EHR lr")),
                         c(K("EHR lr"), K("EHR cox"))))
    }
    JOBS[[length(JOBS) + 1]] <- list(
      part = "landmark", cohort = co, subset = ss, met = met,
      draw = local({
        b <- rb
        function() unlist(b[sample.int(length(b), length(b),
                                       replace = TRUE)], use.names = FALSE)
      }), all = rows, pairs = prs)
  }
}

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
plan(multisession, workers = workers)
proj <- tp * length(tasks) / workers
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(workers), function(i) J$met(J$draw()),
                          future.seed = 1L, future.packages = PKG))
  proj <- proj + as.numeric(difftime(Sys.time(), t0, units = "secs")) *
    NBOOT / workers
}
cat(sprintf("  projected total: about %.1f min\n", proj / 60 * 1.3))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---- part 2: fit and bootstrap ----
res <- future_lapply(tasks, ridge_task, future.packages = c(PKG, "glmnet"),
                     future.seed = TRUE, future.chunk.size = 1L)
cidx <- function(t, e, s) {
  as.numeric(tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE)$
                        concordance, error = function(err) NA_real_))
}
nbad <- 0L
LP <- list()
for (m in names(XM)) {
  rr <- Filter(function(r) r$mod == m, res)
  for (r in rr) {
    f <- FM[[as.character(r$seed)]]
    v <- sapply(sort(unique(f)), function(k) {
      ii <- which(f == k)
      cidx(BA$time_death[ii], BA$event_death[ii], r$lp[ii])
    })
    if (mean(v, na.rm = TRUE) < 0.40) nbad <- nbad + 1L
  }
  LP[[m]] <- Reduce(`+`, lapply(rr, function(r) r$lp)) / length(rr)
}
chk("all-admission scores failing the per-fold rule (D2)", nbad,
    if (nbad == 0) "PASS" else "FAIL")
nm <- sum(sapply(LP, function(v) sum(!is.finite(v))))
chk("all-admission scores missing", nm, if (nm == 0) "PASS" else "FAIL")
stop_if_fail()
XT <- fread(file.path(PR, "level1", "oof_extra_B_ehrecg.csv"))
x6 <- XT[, lapply(.SD, mean), by = hadm_id,
         .SDcols = c("ehr.cox2.death.lp", "ecg.cox2.death.lp")][
           match(fb$hadm_id, hadm_id)]
VERS <- list(
  "first admissions (registered frame), ridge-Cox retuned (stage 6d)" =
    list(R = make_R(fb$time_death, fb$event_death, x6$ehr.cox2.death.lp,
                    x6$ecg.cox2.death.lp), subj = fb$subject_id,
         n = nrow(fb), d = sum(fb$event_death)),
  "all admissions, ridge-Cox retuned refitted" =
    list(R = make_R(BA$time_death, BA$event_death, LP$ehr, LP$ecg),
         subj = BA$subject_id, n = nrow(BA), d = sum(BA$event_death)))
AA <- list()
for (vn in names(VERS)) {
  V <- VERS[[vn]]
  n <- V$n
  pt <- ratio_vec(V$R, rep(1, n))
  rb <- split(seq_len(n), match(V$subj, unique(V$subj)))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), function(i) {
    idx <- unlist(rb[sample.int(length(rb), length(rb), replace = TRUE)],
                  use.names = FALSE)
    ratio_vec(V$R, tabulate(idx, n))
  }, future.seed = SEED, future.packages = PKG))
  for (q in names(pt)) {
    ci <- q2(B[, q])
    AA[[length(AA) + 1]] <- data.table(
      version = vn, admissions = n, patients = length(rb), deaths = V$d,
      quantity = q, estimate = exp(pt[[q]]), lo = exp(ci[1]),
      hi = exp(ci[2]))
  }
}
aa <- rbindlist(AA)
fwrite(aa, file.path(RS, "sens_all_admissions.csv"))
cat("\n  primary estimand, first vs all admissions:\n")
print(aa[quantity == "ecg_ratio", .(version, admissions, patients,
                                    estimate = round(estimate, 3),
                                    lo = round(lo, 3), hi = round(hi, 3))],
      row.names = FALSE)

# ---- parts 1 and 3: bootstrap ----
OUT <- list()
for (J in JOBS) {
  pt <- J$met(J$all)
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), function(i) {
    J$met(J$draw())
  }, future.seed = SEED, future.packages = PKG))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    ci <- q2(B[, k])
    OUT[[length(OUT) + 1]] <- data.table(
      part = J$part, cohort = J$cohort, subset = J$subset,
      analysis = z[1], horizon_or_endpoint = z[2], model = z[3],
      measure = z[4], type = "score", estimate = pt[[k]], lo = ci[1],
      hi = ci[2])
  }
  for (pp in J$pairs) {
    if (!all(pp %in% names(pt))) next
    za <- strsplit(pp[1], "|", fixed = TRUE)[[1]]
    zb <- strsplit(pp[2], "|", fixed = TRUE)[[1]]
    ci <- q2(B[, pp[1]] - B[, pp[2]])
    OUT[[length(OUT) + 1]] <- data.table(
      part = J$part, cohort = J$cohort, subset = J$subset,
      analysis = za[1], horizon_or_endpoint = za[2],
      model = sprintf("%s - %s", za[3], zb[3]), measure = za[4],
      type = "contrast", estimate = pt[[pp[1]]] - pt[[pp[2]]],
      lo = ci[1], hi = ci[2])
  }
  cat(sprintf("  %s %s %s done\n", J$part, J$cohort, J$subset))
}
plan(sequential)
o <- rbindlist(OUT)
o[, excludes_zero := type == "contrast" & (lo > 0 | hi < 0)]
fwrite(o[part == "final_year"], file.path(RS, "sens_final_year.csv"))
fwrite(o[part == "landmark"], file.path(RS, "sens_landmark.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "sens_checks.csv"))
cat("\n=== S1 CONTRASTS, AUROC at 30 days: full cohort vs landmark ===\n")
print(dcast(o[part == "landmark" & type == "contrast" &
                measure == "auc30"],
            cohort + horizon_or_endpoint + model ~ subset,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 4) else v
            })], row.names = FALSE)
cat("\n=== S2 AUROC (full cohort), 365 days: full vs final year out ===\n")
print(dcast(o[part == "final_year" & analysis == "S2" & type == "score" &
                measure == "auc_full" & horizon_or_endpoint == "365"],
            cohort + model ~ subset, value.var = "estimate")[
              , lapply(.SD, function(v) if (is.numeric(v)) round(v, 4)
                       else v)], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("sens_final_year.csv", "sens_all_admissions.csv",
            "sens_landmark.csv", "sens_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "sens_log.txt"), " (terminal log)\n")
cat("\nOutputs: sens_final_year.csv, sens_all_admissions.csv,",
    "sens_landmark.csv, sens_checks.csv\n")