# 20p_interpretability_P.R
# Interpretability of the A36 at-diagnosis model (A37, part 1), as 20a (A23)
# did for Frames A and B.
#
# The outcome is death within 365 days of admission. Part 1 refits the EHR
# spline Cox model of 06p on the presentation EHR (ehrP.csv) once on every
# admission of each cohort (EHR+ECG: all; three-modality: those with a
# presentation CTPA), with the 18a preprocessing fitted on all of them, for
# interpretation only. The model is mgcv cox.ph with cubic regression
# splines (k = 5; features with at least 10 distinct values), other features
# linear with a ridge penalty (D16), REML (optimizer "efs", select = TRUE)
# and D14's fallback to k = 4, then 3. Per feature it reports effective
# degrees of freedom, Wald chi-square, p-value, the log hazard ratio against
# the median from the 1st to the 99th percentile with 95% intervals, and the
# SD of the term across admissions. The SD ranks the features, because the
# chi-square of a heavily shrunk term can be large while its p-value is not.
#
# Part 2 computes SurvSHAP(t) for the Brier-tuned RSF and oblique RSF of 06p
# on the presentation EHR (EHR+ECG cohort), each refitted once on every
# admission at the settings 06p chose most often. As in 20a, interventional
# Shapley values of the risk at 10, 30, 90, 180 and 365 days use 16
# antithetic pairs of permutations against 100 background admissions, for
# 300 random admissions (seed 2026). A positive value raises risk, and each
# admission's values must add up to its risk minus the background mean (the
# script stops beyond 1e-8). Per feature and horizon the summaries are the
# mean |SurvSHAP(t)|, its share of the total, rank and mean signed value;
# the same is done per domain (demographics, vital signs, laboratory tests,
# history, earlier care). Agreement is the Spearman correlation of each
# forest's importance at 365 days with the spline term SD, and between the
# forests.
#
# SHAP batches in data/processed/interpP/ckpt are reused on a rerun if the
# refitted forest reproduces the saved background predictions (to 1e-10);
# delete the folder to start afresh. The script also stops on missing inputs
# or modal settings, or if no spline fit succeeds for a cohort.
#
# Inputs:  data/processed/presentation/frame_P.csv, ehrP.csv;
#          data/processed/level1/params_P.csv
# Outputs: results/interpP_{spline_terms,spline_curves,shap_importance,
#          shap_domains,agreement,checks}.csv, results/interpP_log.txt
# Usage:   Rscript scripts/20p_interpretability_P.R [--probe-only]
#          [--workers=N]
# Notes:   each worker holds one forest (about 1 to 1.5 GB).

suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(data.table)
  library(parallel)
  library(randomForestSRC)
})
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf not installed")
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
PP <- file.path(PR, "presentation")
RS <- file.path(P2, "results")
OD <- file.path(PR, "interpP")
CK <- file.path(OD, "ckpt")
dir.create(CK, showWarnings = FALSE, recursive = TRUE)
HZ <- c(10, 30, 90, 180, 365)
KNOTS <- 5L
MINU <- 10L
NGRID <- 50L
NEXPL <- 300L
NBG <- 100L
NPERM <- 16L
ROWCAP <- 40000L
BATCHK <- 4L
NTREE <- 500L
NSPLIT <- 10L
NTIME <- 60L
BATCH <- 50L
SEED <- 2026L
TOLFP <- 1e-10
TOLADD <- 1e-8
PKG <- c("survival", "data.table", "randomForestSRC", "aorsf")
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
    fwrite(ck, file.path(RS, "interpP_checks.csv"))
    cat("\nSTOP: a check failed. see: interpP_checks.csv\n")
    quit(status = 1)
  }
}

# ---------------- readable names ----------------
LABMAP <- c(
  age = "Age", sex_f = "Female sex", no_pres_vitals =
    "Missing vital signs at presentation", hr = "Heart rate",
  heart_rate = "Heart rate", sbp = "Systolic BP", dbp = "Diastolic BP",
  mbp = "Mean arterial pressure", map = "Mean arterial pressure",
  rr = "Respiratory rate", resp_rate = "Respiratory rate",
  spo2 = "Oxygen saturation", temp = "Temperature",
  temperature = "Temperature", gcs = "Glasgow Coma Scale",
  creatinine = "Creatinine", bun = "Urea nitrogen", urea = "Urea nitrogen",
  sodium = "Sodium", potassium = "Potassium", bicarb = "Bicarbonate",
  bicarbonate = "Bicarbonate", glucose = "Glucose", hgb = "Haemoglobin",
  hemoglobin = "Haemoglobin", wbc = "White cell count",
  plt = "Platelets", platelets = "Platelets", lactate = "Lactate",
  troponin_t = "Troponin T", troponin_i = "Troponin I",
  ntprobnp = "NT-proBNP", aniongap = "Anion gap", anion_gap = "Anion gap",
  chloride = "Chloride", calcium = "Calcium",
  calcium_total = "Calcium", magnesium = "Magnesium",
  phosphate = "Phosphate", hct = "Haematocrit",
  hematocrit = "Haematocrit", rdw = "RDW", mcv = "MCV", mch = "MCH",
  mchc = "MCHC", rbc = "Red cell count", inr = "INR",
  inr_pt = "INR", pt = "Prothrombin time", ptt = "PTT",
  albumin = "Albumin", alt = "ALT", ast = "AST",
  alp = "Alkaline phosphatase", alkaline_phosphatase =
    "Alkaline phosphatase", bilirubin = "Bilirubin",
  bilirubin_total = "Bilirubin", ldh = "LDH", ld_ldh = "LDH",
  ck = "Creatine kinase", ck_cpk = "Creatine kinase", ckmb = "CK-MB",
  ck_mb = "CK-MB", crp = "CRP", c_reactive_protein = "CRP", ph = "pH",
  po2 = "pO2", pco2 = "pCO2", d_dimer = "D-dimer", ddimer = "D-dimer",
  fibrinogen = "Fibrinogen", lipase = "Lipase",
  base_excess = "Base excess")
HXMAP <- c(
  mi = "Myocardial infarction", chf = "Heart failure",
  pvd = "Peripheral vascular disease", cevd = "Cerebrovascular disease",
  dementia = "Dementia", cpd = "Chronic lung disease",
  rheumd = "Rheumatic disease", pud = "Peptic ulcer",
  mld = "Mild liver disease", diab = "Diabetes",
  diabwc = "Diabetes with complications", hp = "Hemiplegia",
  rend = "Renal disease", canc = "Cancer", msld = "Severe liver disease",
  metacanc = "Metastatic cancer", aids = "HIV/AIDS",
  carit = "Arrhythmia", valv = "Valvular disease",
  pcd = "Pulmonary circulation disorder", hypunc = "Hypertension",
  hypc = "Hypertension with complications", para = "Paralysis",
  ond = "Other neurological disorder", diabunc = "Diabetes",
  diabc = "Diabetes with complications", hypothy = "Hypothyroidism",
  rf = "Renal failure", ld = "Liver disease", lymph = "Lymphoma",
  solidtum = "Solid tumour", coag = "Coagulopathy", obes = "Obesity",
  wloss = "Weight loss", fed = "Fluid or electrolyte disorder",
  blane = "Blood loss anaemia", dane = "Deficiency anaemia",
  alcohol = "Alcohol misuse", drug = "Drug misuse",
  psycho = "Psychosis", depre = "Depression")
label_of <- function(nm) {
  sapply(nm, function(x) {
    if (x == "hx_charlson_score") return("Charlson score")
    if (x == "hx_elixhauser_vw") return("Elixhauser score")
    if (x == "n_prior") return("Earlier admissions")
    if (x == "n_codes") return("Earlier diagnosis codes")
    if (grepl("^hx_(ch|el)_", x)) {
      k <- sub("^hx_(ch|el)_", "", x)
      ix <- if (startsWith(x, "hx_ch_")) "Charlson" else "Elixhauser"
      nmk <- if (k %in% names(HXMAP)) HXMAP[[k]] else k
      return(sprintf("%s (%s)", nmk, ix))
    }
    b <- sub("_first$", "", x)
    b <- sub("^lab_[0-9]+_", "", b)
    k <- tolower(b)
    if (k %in% names(LABMAP)) return(LABMAP[[k]])
    s <- gsub("_", " ", b)
    paste0(toupper(substring(s, 1, 1)), substring(s, 2))
  }, USE.NAMES = FALSE)
}
domain_of <- function(nm) {
  iv <- which(nm == "no_pres_vitals")
  is <- which(nm == "sex_f")
  vit <- if (length(iv) && length(is) && iv > is) {
    nm[(is + 1):iv]
  } else {
    c("no_pres_vitals", grep("^(hr|sbp|dbp|rr|spo2|temp)_first$", nm,
                             value = TRUE))
  }
  ifelse(nm %in% c("age", "sex_f"), "Demographics",
         ifelse(nm %in% vit, "Vital signs",
                ifelse(startsWith(nm, "hx_"), "History",
                       ifelse(nm %in% c("n_prior", "n_codes"),
                              "Earlier care", "Laboratory tests"))))
}

# ---------------- shared helpers ----------------
# stage 18a preprocessing, fitted on the rows `tr`
prep_fit <- function(X, tr = seq_len(nrow(X))) {
  med <- q1 <- q2 <- mu <- s <- numeric(ncol(X))
  bin <- logical(ncol(X))
  Z <- X
  for (j in seq_len(ncol(X))) {
    med[j] <- median(X[tr, j], na.rm = TRUE)
    if (!is.finite(med[j])) med[j] <- 0
    Z[is.na(Z[, j]), j] <- med[j]
    bin[j] <- length(unique(Z[tr, j])) <= 2
    q1[j] <- if (bin[j]) -Inf else quantile(Z[tr, j], 0.005, names = FALSE)
    q2[j] <- if (bin[j]) Inf else quantile(Z[tr, j], 0.995, names = FALSE)
    Z[, j] <- pmin(pmax(Z[, j], q1[j]), q2[j])
    mu[j] <- mean(Z[tr, j])
    s[j] <- sd(Z[tr, j])
  }
  keep <- is.finite(s) & s > 0
  list(med = med, q1 = q1, q2 = q2, mu = mu, s = s, keep = keep,
       bin = bin, names = colnames(X))
}
prep_apply <- function(P, X) {
  for (j in seq_len(ncol(X))) {
    X[is.na(X[, j]), j] <- P$med[j]
    X[, j] <- pmin(pmax(X[, j], P$q1[j]), P$q2[j])
  }
  sweep(sweep(X[, P$keep, drop = FALSE], 2, P$mu[P$keep]), 2,
        P$s[P$keep], "/")
}
modal <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  as.numeric(names(sort(table(v), decreasing = TRUE))[1])
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
need <- c(file.path(PP, c("frame_P.csv", "ehrP.csv")),
          file.path(PR, "level1", "params_P.csv"))
miss <- need[!file.exists(need)]
chk("input files found", if (length(miss)) basename(miss[1]) else "all",
    if (length(miss) == 0) "PASS" else "FAIL")
stop_if_fail()
fp <- fread(file.path(PP, "frame_P.csv"))
E <- unique(fread(file.path(PP, "ehrP.csv")), by = "hadm_id")
cE <- setdiff(names(E), c("hadm_id", "subject_id"))
XALL <- as.matrix(E[match(fp$hadm_id, hadm_id), ..cE])
storage.mode(XALL) <- "double"
nmiss_row <- sum(!fp$hadm_id %in% E$hadm_id)
chk("EHR at presentation: features (admissions without a row)",
    sprintf("%d (%d)", length(cE), nmiss_row),
    if (nmiss_row == 0) "PASS" else "FAIL")
stop_if_fail()
LAB <- label_of(cE)
DOM <- domain_of(cE)
chk("features by domain (demog / vitals / labs / history / care)",
    paste(sapply(c("Demographics", "Vital signs", "Laboratory tests",
                   "History", "Earlier care"),
                 function(d) sum(DOM == d)), collapse = " / "), "INFO")
COH <- list("EHR+ECG" = rep(TRUE, nrow(fp)),
            "3-mod" = fp$in_three_mod_p == 1)

# ---------------- part 1: spline Cox partial effects ----------------
cat("\n=== PART 1: SPLINE COX PARTIAL EFFECTS ===\n")
spline_fit <- function(Z, tt, ee) {
  nu <- apply(Z, 2, function(v) length(unique(v)))
  sm <- colnames(Z)[nu >= MINU]
  li <- setdiff(colnames(Z), sm)
  D <- data.frame(t = tt, e = ee, Z[, sm, drop = FALSE])
  pp <- NULL
  if (length(li)) {
    D$L <- Z[, li, drop = FALSE]
    pp <- list(L = list(diag(length(li))))
  }
  msg <- character(0)
  for (kk in unique(c(KNOTS, 4L, 3L))) {
    fml <- as.formula(paste("t ~", paste(c(sprintf(
      "s(%s, k = %d, bs = \"cr\")", sm, kk), if (length(li)) "L"),
      collapse = " + ")))
    g <- tryCatch(gam(fml, family = cox.ph(), weights = e, data = D,
                      method = "REML", select = TRUE, optimizer = "efs",
                      paraPen = pp),
                  error = function(err) conditionMessage(err))
    if (!is.character(g)) {
      return(list(fit = g, k = kk, sm = sm, li = li, msg = msg, D = D))
    }
    msg <- c(msg, sprintf("k = %d: %s", kk, substr(g, 1, 100)))
  }
  list(fit = NULL, k = 0L, sm = sm, li = li, msg = msg, D = D)
}
TERMS <- list()
CURV <- list()
for (co in names(COH)) {
  rows <- which(COH[[co]])
  sub <- fp[rows]
  X <- XALL[rows, , drop = FALSE]
  P <- prep_fit(X)
  Z <- prep_apply(P, X)
  kj <- which(P$keep)
  colnames(Z) <- sprintf("f%d", kj)
  if (any(!P$keep)) {
    chk(sprintf("%s: constant features left out", co),
        paste(cE[!P$keep], collapse = ", "), "INFO")
  }
  tt <- sub$time_death
  ee <- as.integer(sub$event_death == 1)
  t0 <- Sys.time()
  S1 <- spline_fit(Z, tt, ee)
  if (is.null(S1$fit)) {
    chk(sprintf("%s: spline Cox fit", co), S1$msg[1], "FAIL")
    stop_if_fail()
  }
  fit <- S1$fit
  cvg <- fit$outer.info$conv
  chk(sprintf("%s: spline Cox fit, k = %d (%.1f min, total edf %.1f)",
              co, S1$k, as.numeric(difftime(Sys.time(), t0,
                                            units = "mins")),
              sum(fit$edf)), if (is.null(cvg)) "converged" else cvg,
      if (S1$k == KNOTS && (is.null(cvg) ||
                            grepl("full convergence", cvg))) "PASS" else
        "WARN")
  b <- coef(fit)
  Vp <- fit$Vp
  st <- summary(fit)$s.table
  sid <- setNames(seq_along(fit$smooth),
                  sapply(fit$smooth, function(s) s$term))
  for (fj in colnames(Z)) {
    j <- as.integer(sub("^f", "", fj))
    smooth <- fj %in% S1$sm
    z <- Z[, fj]
    zmed <- median(z)
    if (smooth) {
      smo <- fit$smooth[[sid[[fj]]]]
      cols <- smo$first.para:smo$last.para
      rn <- sprintf("s(%s)", fj)
      ed <- st[rn, "edf"]
      chi <- st[rn, "Chi.sq"]
      pv <- st[rn, "p-value"]
      tsd <- sd(as.vector(PredictMat(smo, S1$D) %*% b[cols]))
      zg <- seq(quantile(z, 0.01), quantile(z, 0.99), length.out = NGRID)
      Bg <- PredictMat(smo, setNames(data.frame(c(zg, zmed)), fj))
      Dm <- Bg[seq_along(zg), , drop = FALSE] -
        matrix(Bg[length(zg) + 1, ], length(zg), ncol(Bg), byrow = TRUE)
      est <- as.vector(Dm %*% b[cols])
      se <- sqrt(pmax(0, rowSums((Dm %*% Vp[cols, cols, drop = FALSE]) *
                                   Dm)))
    } else {
      jj <- which(names(b) == paste0("L", fj))
      ed <- fit$edf[jj]
      chi <- b[jj]^2 / Vp[jj, jj]
      pv <- pchisq(chi, 1, lower.tail = FALSE)
      tsd <- sd(z * b[jj])
      zg <- sort(unique(z))
      est <- b[jj] * (zg - zmed)
      se <- abs(zg - zmed) * sqrt(Vp[jj, jj])
    }
    TERMS[[length(TERMS) + 1]] <- data.table(
      frame = "P", cohort = co, feature = cE[j], label = LAB[j],
      domain = DOM[j], type = if (smooth) "spline" else "linear (ridge)",
      edf = ed, chi_sq = chi, p_value = pv, term_sd = tsd, k = S1$k,
      n = nrow(sub), patients = length(unique(sub$subject_id)),
      deaths = sum(ee))
    CURV[[length(CURV) + 1]] <- data.table(
      frame = "P", cohort = co, feature = cE[j], label = LAB[j],
      type = if (smooth) "spline" else "linear (ridge)",
      x = P$mu[j] + P$s[j] * zg, x_median = P$mu[j] + P$s[j] * zmed,
      x_p10 = quantile(X[, j], 0.10, na.rm = TRUE, names = FALSE),
      x_p90 = quantile(X[, j], 0.90, na.rm = TRUE, names = FALSE),
      log_hr = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }
  rm(fit, S1)
  invisible(gc(verbose = FALSE))
}
terms <- rbindlist(TERMS)
setorder(terms, frame, cohort, -term_sd)
fwrite(terms, file.path(RS, "interpP_spline_terms.csv"))
fwrite(rbindlist(CURV), file.path(RS, "interpP_spline_curves.csv"))
cat("\n  largest terms, EHR+ECG cohort (term SD):\n")
print(terms[cohort == "EHR+ECG"][1:12,
            .(label, domain, type, edf = round(edf, 1),
              chi_sq = round(chi_sq, 1), p_value = signif(p_value, 2),
              term_sd = round(term_sd, 3))], row.names = FALSE)

# ---------------- part 2: SurvSHAP(t) ----------------
cat("\n=== PART 2: SurvSHAP(t) FOR THE FORESTS ===\n")
fit_rsf <- function(X, tt, ee, node, mtryv) {
  P <- prep_fit(X)
  Xi <- X
  for (j in seq_len(ncol(X))) Xi[is.na(Xi[, j]), j] <- P$med[j]
  colnames(Xi) <- sprintf("f%d", seq_len(ncol(Xi)))
  ev <- sort(unique(tt[ee == 1]))
  at_h <- sapply(HZ, function(h) {
    v <- ev[ev <= h + 0.5]
    if (length(v)) max(v) else NA_real_
  })
  pick <- ev[unique(round(seq(1, length(ev),
                              length.out = min(NTIME, length(ev)))))]
  nt <- sort(unique(c(pick, at_h[is.finite(at_h)])))
  f <- rfsrc(Surv(time, status) ~ ., ntree = NTREE, nodesize = node,
             mtry = mtryv, nsplit = NSPLIT, ntime = nt,
             splitrule = "logrank", save.memory = TRUE, seed = -SEED,
             data = data.frame(time = tt, status = ee, Xi))
  list(kind = "rsf", f = f, P = P, H = HZ)
}
fit_orsf <- function(X, tt, ee, leaf, mtry) {
  P <- prep_fit(X)
  Z <- prep_apply(P, X)
  colnames(Z) <- sprintf("f%d", seq_len(ncol(Z)))
  mtry <- max(1L, min(ncol(Z), as.integer(mtry)))
  d <- data.frame(time = tt, status = ee, Z)
  fs <- lapply(seq_len(NTREE %/% BATCH), function(bb) {
    aorsf::orsf(d, Surv(time, status) ~ ., n_tree = BATCH, mtry = mtry,
                leaf_min_obs = leaf, split_min_obs = 2L * leaf,
                oobag_pred_type = "none",
                tree_seeds = as.integer(SEED * 10 + bb), n_thread = 1)
  })
  list(kind = "orsf", f = fs, P = P, H = HZ, mtry = mtry)
}
predict_risk <- function(M, X) {
  if (M$kind == "rsf") {
    Xi <- X
    for (j in seq_len(ncol(X))) Xi[is.na(Xi[, j]), j] <- M$P$med[j]
    colnames(Xi) <- sprintf("f%d", seq_len(ncol(Xi)))
    pr <- predict(M$f, newdata = as.data.frame(Xi))
    R <- sapply(M$H, function(h) {
      i <- findInterval(h + 0.5, pr$time.interest)
      if (i == 0) rep(0, nrow(X)) else 1 - pr$survival[, i]
    })
    matrix(R, nrow = nrow(X))
  } else {
    Z <- prep_apply(M$P, X)
    colnames(Z) <- sprintf("f%d", seq_len(ncol(Z)))
    Zd <- as.data.frame(Z)
    Reduce(`+`, lapply(M$f, function(f) {
      matrix(predict(f, new_data = Zd, pred_type = "risk",
                     pred_horizon = M$H + 0.5, boundary_checks = FALSE),
             ncol = length(M$H))
    })) / length(M$f)
  }
}
# antithetic permutation SHAP of the risk vector, for one admission
shap_one <- function(M, x, BG, nperm, seed, chunk) {
  p <- length(x)
  K <- length(M$H)
  set.seed(seed)
  perms <- lapply(seq_len(nperm), function(i) sample.int(p))
  perms <- c(perms, lapply(perms, rev))
  phi <- matrix(0, p, K)
  nb <- nrow(BG)
  for (ch in split(seq_along(perms), ceiling(seq_along(perms) / chunk))) {
    Xall <- do.call(rbind, lapply(perms[ch], function(pm) {
      do.call(rbind, lapply(0:p, function(s) {
        B <- BG
        if (s > 0) {
          cs <- pm[seq_len(s)]
          B[, cs] <- matrix(x[cs], nb, s, byrow = TRUE)
        }
        B
      }))
    }))
    R <- predict_risk(M, Xall)
    blk <- rep(seq_len(length(ch) * (p + 1)), each = nb)
    V <- rowsum(R, blk) / nb
    for (i in seq_along(ch)) {
      pm <- perms[[ch[i]]]
      v <- V[(i - 1) * (p + 1) + seq_len(p + 1), , drop = FALSE]
      phi[pm, ] <- phi[pm, ] + (v[-1, , drop = FALSE] -
                                  v[-(p + 1), , drop = FALSE])
    }
  }
  phi / length(perms)
}
# run on a worker, which holds the forest and data in JOB
shap_worker <- function(i) {
  shap_one(JOB$M, JOB$XE[i, ], JOB$BG, JOB$nperm, 1000L + i, JOB$chunk)
}
set_job <- function(obj) {
  assign("JOB", obj, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
  NULL
}
pp <- fread(file.path(PR, "level1", "params_P.csv"))
q <- pp[cohort == "EHR+ECG" & modality == "ehr"]
node <- modal(q$rsf_node)
mtry_r <- modal(q$rsf_mtry)
leaf <- modal(q$orsf_leaf)
mtry_o <- modal(q$orsf_mtry)
okm <- nrow(q) > 0 && all(is.finite(c(node, mtry_r, leaf, mtry_o)))
chk("stage 6p modal settings: RSF node/mtry, oblique leaf/mtry",
    sprintf("%g / %g, %g / %g", node, mtry_r, leaf, mtry_o),
    if (okm) "PASS" else "FAIL")
stop_if_fail()
X <- XALL
tt <- fp$time_death
ee <- as.integer(fp$event_death == 1)
set.seed(SEED)
bg <- sample.int(nrow(X), min(NBG, nrow(X)))
ex <- sample.int(nrow(X), min(NEXPL, nrow(X)))
BG <- X[bg, , drop = FALSE]
XE <- X[ex, , drop = FALSE]
p <- ncol(X)
CHUNK <- max(1L, floor(ROWCAP / ((p + 1) * nrow(BG))))
cat(sprintf(paste("  %d admissions; %d features; %d explained, %d",
                  "background; %d permutation orders, %d per",
                  "prediction call\n"), nrow(X), p, length(ex), nrow(BG),
            2L * NPERM, CHUNK))
MODELS <- list()
t0 <- Sys.time()
MODELS[["RSF, Brier-tuned"]] <- fit_rsf(X, tt, ee, node, mtry_r)
cat(sprintf("  RSF fitted in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
t0 <- Sys.time()
MODELS[["oblique RSF, Brier-tuned"]] <- fit_orsf(X, tt, ee, leaf, mtry_o)
cat(sprintf("  oblique RSF fitted in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
key_of <- function(mn) paste0("P_EHRECG_", gsub("[^A-Za-z]+", "_", mn))
done_ids <- function(key) {
  fs <- list.files(CK, pattern = paste0("^", key, "_i[0-9]+\\.rds$"),
                   full.names = TRUE)
  unlist(lapply(fs, function(f) readRDS(f)$ids))
}
# fingerprint: refits must reproduce the saved background predictions
FP <- list()
for (mn in names(MODELS)) {
  M <- MODELS[[mn]]
  key <- key_of(mn)
  f0 <- colMeans(predict_risk(M, BG))
  fx5 <- predict_risk(M, XE[1:5, , drop = FALSE])
  fpf <- file.path(CK, paste0(key, "_fp.rds"))
  nd <- length(done_ids(key))
  if (file.exists(fpf) && nd > 0) {
    old <- readRDS(fpf)
    d <- max(abs(c(old$f0 - f0, old$fx5 - fx5)))
    if (is.finite(d) && d <= TOLFP) {
      chk(sprintf("%s: saved SHAP batches reused (refit reproduces)", mn),
          sprintf("%d admissions, max |diff| %.1e", nd, d), "PASS")
    } else {
      chk(sprintf("%s: refit differs from the saved one, restarted", mn),
          sprintf("max |diff| %.1e", d), "WARN")
      unlink(list.files(CK, pattern = paste0("^", key, "_i"),
                        full.names = TRUE))
    }
  }
  saveRDS(list(f0 = f0, fx5 = fx5), fpf)
  FP[[mn]] <- f0
}
cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
op <- options(rf.cores = 1L, mc.cores = 1L)
tot <- 0
for (mn in names(MODELS)) {
  todo <- setdiff(seq_along(ex), done_ids(key_of(mn)))
  t0 <- Sys.time()
  invisible(shap_one(MODELS[[mn]], XE[1, ], BG, NPERM, 1001L, CHUNK))
  tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tot <- tot + tp * length(todo)
  cat(sprintf("  %-26s %6.1f s per admission, %d admissions to explain\n",
              mn, tp, length(todo)))
}
options(op)
cat(sprintf("  about %.1f hours on %d workers\n",
            tot * 1.15 / workers / 3600, workers))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
cl <- NULL
for (mn in names(MODELS)) {
  key <- key_of(mn)
  todo <- setdiff(seq_along(ex), done_ids(key))
  if (!length(todo)) next
  if (is.null(cl)) {
    cl <- makePSOCKcluster(workers)
    clusterCall(cl, function(pk) {
      for (k in pk) suppressPackageStartupMessages(
        library(k, character.only = TRUE))
      options(rf.cores = 1L, mc.cores = 1L)
      NULL
    }, PKG)
    clusterExport(cl, c("shap_one", "predict_risk", "prep_apply",
                        "shap_worker"))
  }
  clusterCall(cl, set_job, list(M = MODELS[[mn]], XE = XE, BG = BG,
                                nperm = NPERM, chunk = CHUNK))
  t0 <- Sys.time()
  bts <- split(todo, ceiling(seq_along(todo) / (workers * BATCHK)))
  nd <- length(ex) - length(todo)
  for (bt in bts) {
    ph <- parLapplyLB(cl, bt, shap_worker)
    f <- file.path(CK, sprintf("%s_i%04d.rds", key, min(bt)))
    saveRDS(list(ids = bt, phi = ph), paste0(f, ".tmp"))
    file.rename(paste0(f, ".tmp"), f)
    nd <- nd + length(bt)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %s: %d of %d explained (%.0f min; about %.0f min left)\n",
                mn, nd, length(ex), el,
                el / (nd - (length(ex) - length(todo))) *
                  (length(ex) - nd)))
  }
  clusterCall(cl, set_job, NULL)
}
if (!is.null(cl)) stopCluster(cl)

cat("\n=== SUMMARIES ===\n")
IMP <- list()
DMS <- list()
AGR <- list()
MA <- list()
doms <- c("Demographics", "Vital signs", "Laboratory tests", "History",
          "Earlier care")
for (mn in names(MODELS)) {
  M <- MODELS[[mn]]
  key <- key_of(mn)
  PH <- vector("list", length(ex))
  for (f in list.files(CK, pattern = paste0("^", key, "_i[0-9]+\\.rds$"),
                       full.names = TRUE)) {
    r <- readRDS(f)
    PH[r$ids] <- r$phi
  }
  nmis <- sum(sapply(PH, is.null))
  chk(sprintf("%s: admissions explained", mn), sprintf(
    "%d of %d", length(ex) - nmis, length(ex)),
    if (nmis == 0) "PASS" else "FAIL")
  stop_if_fail()
  fx <- predict_risk(M, XE)
  f0 <- FP[[mn]]
  add <- max(sapply(seq_along(PH), function(i) {
    max(abs(colSums(PH[[i]]) - (fx[i, ] - f0)))
  }))
  chk(sprintf("%s: SHAP additivity (max |diff|)", mn),
      sprintf("%.1e", add), if (add <= TOLADD) "PASS" else "FAIL")
  A <- array(unlist(PH), c(p, length(HZ), length(PH)))
  long <- rbindlist(lapply(seq_along(PH), function(i) {
    data.table(hadm_id = fp$hadm_id[ex[i]], feature = rep(cE, length(HZ)),
               horizon = rep(HZ, each = p),
               value = rep(XE[i, ], length(HZ)),
               shap_risk = as.vector(PH[[i]]))
  }))
  fwrite(long, file.path(OD, sprintf("shap_P_EHRECG_%s.csv",
                                     gsub("[^A-Za-z]+", "_", mn))))
  ma <- apply(abs(A), c(1, 2), mean)
  MA[[mn]] <- ma
  sh <- sweep(ma, 2, colSums(ma), "/")
  for (hi in seq_along(HZ)) {
    IMP[[length(IMP) + 1]] <- data.table(
      frame = "P", cohort = "EHR+ECG", model = mn, feature = cE,
      label = LAB, domain = DOM, horizon = HZ[hi],
      mean_abs_shap = ma[, hi], share = sh[, hi],
      rank = rank(-ma[, hi], ties.method = "first"),
      mean_shap = apply(A[, hi, , drop = FALSE], 1, mean),
      n_explained = length(PH), n_background = nrow(BG))
  }
  # domain attribution: sum of the domain's values per admission
  IND <- sapply(doms, function(d) as.numeric(DOM == d))
  G <- sapply(seq_along(PH), function(i) {
    abs(crossprod(IND, PH[[i]]))
  }, simplify = "array")
  gm <- apply(G, c(1, 2), mean)
  gs <- sweep(gm, 2, colSums(gm), "/")
  for (hi in seq_along(HZ)) {
    DMS[[length(DMS) + 1]] <- data.table(
      frame = "P", cohort = "EHR+ECG", model = mn, domain = doms,
      n_features = as.integer(table(factor(DOM, levels = doms))),
      horizon = HZ[hi], mean_abs_shap = gm[, hi], share = gs[, hi])
  }
  sp <- terms[cohort == "EHR+ECG"]
  K <- length(HZ)
  r <- suppressWarnings(cor(ma[, K], sp$term_sd[match(cE, sp$feature)],
                            method = "spearman", use = "complete.obs"))
  AGR[[length(AGR) + 1]] <- data.table(
    frame = "P", cohort = "EHR+ECG", model = mn,
    comparison = "importance at 365 days vs spline term SD",
    spearman = r, forest_top5 = paste(LAB[order(-ma[, K])][1:5],
                                      collapse = "; "),
    spline_top5 = paste(sp[order(-term_sd)]$label[1:5], collapse = "; "))
}
K <- length(HZ)
r2 <- cor(MA[[1]][, K], MA[[2]][, K], method = "spearman")
AGR[[length(AGR) + 1]] <- data.table(
  frame = "P", cohort = "EHR+ECG", model = paste(names(MA),
                                                 collapse = " vs "),
  comparison = "importance at 365 days, forest vs forest", spearman = r2,
  forest_top5 = paste(LAB[order(-MA[[1]][, K])][1:5], collapse = "; "),
  spline_top5 = paste(LAB[order(-MA[[2]][, K])][1:5], collapse = "; "))
imp <- rbindlist(IMP)
dms <- rbindlist(DMS)
agr <- rbindlist(AGR)
fwrite(imp, file.path(RS, "interpP_shap_importance.csv"))
fwrite(dms, file.path(RS, "interpP_shap_domains.csv"))
fwrite(agr, file.path(RS, "interpP_agreement.csv"))
ck <- rbindlist(CHK)
fwrite(ck, file.path(RS, "interpP_checks.csv"))
cat("\n=== AGREEMENT ===\n")
print(agr[, .(model, comparison, spearman = round(spearman, 2))],
      row.names = FALSE)
cat("\n=== SHARE OF IMPORTANCE BY DOMAIN AND HORIZON ===\n")
print(dcast(dms, model + domain ~ horizon, value.var = "share")[
  , lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)],
  row.names = FALSE)
cat("\n=== SHARE BY HORIZON, top 10 features at 365 days ===\n")
for (mn in unique(imp$model)) {
  x <- imp[model == mn]
  top <- x[horizon == max(horizon)][order(rank)]$feature[1:10]
  cat("  ", mn, "\n")
  y <- dcast(x[feature %in% top], label ~ horizon, value.var = "share")
  print(y[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)],
        row.names = FALSE)
}
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("interpP_spline_terms.csv", "interpP_spline_curves.csv",
            "interpP_shap_importance.csv", "interpP_shap_domains.csv",
            "interpP_agreement.csv", "interpP_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", OD, " (shap_*.csv, patient-level, local only)\n")
cat("  ", file.path(RS, "interpP_log.txt"), " (terminal log)\n")
cat("\nOutputs: interpP_spline_terms.csv, interpP_shap_importance.csv,",
    "interpP_shap_domains.csv, interpP_agreement.csv,",
    "interpP_checks.csv\n")