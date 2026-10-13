# 11e_eval_s4_neural.R
# S4, neural learners against the linear learners (PREREGISTRATION.md
# section 6; PREREGISTRATION_ADDENDUM_2.md A26).
#
# Compares the registered neural learners (stage 10) with the linear
# learners, per modality and fused, using seed-averaged out-of-fold scores
# for every frame, cohort and modality. DeepSurv is compared with ridge-Cox
# (stage 6a) for death, composite and CV readmission (cause-specific) and,
# for death, also with ridge-Cox retuned (stage 6d), by Harrell's C over
# the frame's follow-up, AUROC at every horizon, and IPA and
# log-likelihood skill of the horizon risks (death, composite). The MLP
# twin is compared with L2 logistic (stage 6a) by AUROC, IPA and
# log-likelihood skill at every horizon, for death, composite and
# unconditional CV readmission. DeepHit's cumulative incidence of CV
# readmission is compared with the cause-specific ridge-Cox score by AUROC
# at every horizon (deaths first count as no event), with DeepHit's IPA and
# log-likelihood skill as point estimates.
#
# Fused models (death, each cohort's full modality set) are, per seed and
# outer fold, a Cox stack on rank-based normal scores (as in S4 and A14) of
# each family's death scores (ridge-Cox, ridge-Cox retuned, DeepSurv), with
# Breslow risks averaged over seeds.
#
# Contrasts are DeepSurv minus ridge-Cox (and minus ridge-Cox retuned for
# death), MLP twin minus L2 logistic, DeepHit minus cause-specific
# ridge-Cox, and fused DeepSurv minus each fused ridge-Cox. Intervals come
# from 2,000 bootstrap samples (seed 42), clustered by patient in Frame A.
#
# Checks: inputs present and aligned, no missing scores (FAIL).
# Version 2 corrects the sign of the log-likelihood skill score (AUROC, C
# and IPA were unaffected).
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv,
#          data/processed/level1/oof_linear_*.csv, oof_extra_*.csv,
#          oof_deep_*.csv
# Outputs: results/s4n_metrics.csv, s4n_contrasts.csv, s4n_checks.csv,
#          s4n_log.txt (terminal log)
# Usage:   Rscript scripts/11e_eval_s4_neural.R [--probe-only] [--workers=N]

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(future.apply)
})
args <- commandArgs(trailingOnly = TRUE)
probe_only <- "--probe-only" %in% args
P2 <- file.path(Sys.getenv("PROJECTS_DIR", unset = dirname(getwd())), "SurvivalProject")
PR <- file.path(P2, "data", "processed")
L1 <- file.path(PR, "level1")
RS <- file.path(P2, "results")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
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
    fwrite(ck, file.path(RS, "s4n_checks.csv"))
    cat("\nSTOP: a check failed. see: s4n_checks.csv\n")
    quit(status = 1)
  }
}
nscore <- function(x, ref) {
  qnorm((findInterval(x, sort(ref)) + 0.5) / (length(ref) + 1))
}
bres_cuts <- function(tt, ee, lp, cuts) {
  o <- order(tt)
  tt <- tt[o]
  ee <- ee[o]
  w <- exp(lp[o] - max(lp))
  rs <- rev(cumsum(rev(w)))
  ut <- unique(tt[ee == 1])
  if (!length(ut)) return(rep(0, length(cuts)))
  d <- tabulate(match(tt[ee == 1], ut), length(ut))
  ch <- cumsum(d / rs[match(ut, tt)])
  sapply(cuts, function(cu) {
    j <- sum(ut <= cu)
    if (j == 0) 0 else ch[j]
  })
}
stack_ph <- function(S, tr, te, tt, ee, cuts) {
  Ztr <- sapply(S, function(v) nscore(v[tr], v[tr]))
  Zte <- sapply(S, function(v) nscore(v[te], v[tr]))
  D <- data.frame(t = tt[tr], e = ee[tr], Ztr)
  b <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
  b[is.na(b)] <- 0
  lt <- as.vector(as.matrix(Ztr) %*% b)
  lv <- as.vector(as.matrix(Zte) %*% b)
  H <- bres_cuts(tt[tr], ee[tr], lt, cuts)
  sapply(seq_along(cuts), function(i) -expm1(-H[i] * exp(lv - max(lt))))
}
# one statistic per item: auc (s, label), cidx (s, time, event), ipa and
# ll (risk, label)
stat_vec <- function(J, idx) {
  auc <- function(y, s) {
    n1 <- sum(y)
    n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- data.table::frank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  v <- numeric(length(J$items))
  names(v) <- names(J$items)
  for (k in names(J$items)) {
    it <- J$items[[k]]
    s <- J$S[[it$s]][idx]
    if (it$type == "cidx") {
      tm <- J$T[[it$en]]$t[idx]
      ev <- J$T[[it$en]]$e[idx]
      v[[k]] <- as.numeric(tryCatch(survival::concordance(
        survival::Surv(tm, ev) ~ s, reverse = TRUE)$concordance,
        error = function(e) NA_real_))
    } else {
      y <- J$Y[[it$lab]][idx]
      if (it$type == "auc") {
        v[[k]] <- auc(y, s)
      } else {
        p <- pmin(pmax(s, 1e-6), 1 - 1e-6)
        pr <- mean(y)
        v[[k]] <- if (it$type == "ipa") {
          1 - mean((y - p)^2) / (pr * (1 - pr))
        } else {
          1 - mean(y * log(p) + (1 - y) * log(1 - p)) /
            (pr * log(pr) + (1 - pr) * log(1 - pr))
        }
      }
    }
  }
  v
}
boot_rep <- function(i, J) {
  idx <- unlist(J$rows_by[sample.int(length(J$rows_by),
                                     length(J$rows_by), replace = TRUE)],
                use.names = FALSE)
  stat_vec(J, idx)
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
JOBS <- list()
for (fr in c("A", "B")) for (co in c("EHR+ECG", "3-mod")) {
  key <- paste(fr, co)
  code <- if (co == "EHR+ECG") "ehrecg" else "3mod"
  sub <- if (co == "EHR+ECG") frames[[fr]] else
    frames[[fr]][in_three_mod == 1]
  n <- nrow(sub)
  hz <- HZ[[fr]]
  mods <- c("ehr", "ecg", if (co == "3-mod") "ctpa")
  ff <- file.path(L1, sprintf("oof_%s_%s_%s.csv",
                              c("linear", "extra", "deep"), fr, code))
  ok <- all(file.exists(ff))
  chk(sprintf("%s: stage 6a, 6d and 10 files found", key),
      paste(basename(ff), collapse = ", "), if (ok) "PASS" else "FAIL")
  stop_if_fail()
  LT <- fread(ff[1])
  XT <- fread(ff[2])
  DT <- fread(ff[3])
  avg <- function(D, cols) {
    miss <- setdiff(cols, names(D))
    if (length(miss)) {
      chk(sprintf("%s: columns present", key), miss[1], "FAIL")
      stop_if_fail()
    }
    D[, lapply(.SD, mean), by = hadm_id, .SDcols = cols][
      match(sub$hadm_id, hadm_id)]
  }
  TG <- c(death = "death", composite = "composite", cv = "cs_cv")
  LB <- c(death = "death", composite = "composite", cv = "cvu")
  lcols <- c(
    unlist(lapply(mods, function(m) sprintf("%s.cox.%s.lp", m, TG))),
    unlist(lapply(mods, function(m) unlist(lapply(c("death",
      "composite"), function(g) sprintf("%s.cox.%s.r%d", m, g, hz))))),
    unlist(lapply(mods, function(m) unlist(lapply(LB, function(g)
      sprintf("%s.lr.%s_h%d", m, g, hz))))))
  dcols <- c(unlist(lapply(mods, function(m) sprintf("%s.dsurv.%s.lp", m,
                                                     TG))),
    unlist(lapply(mods, function(m) unlist(lapply(c("death",
      "composite"), function(g) sprintf("%s.dsurv.%s.r%d", m, g, hz))))),
    unlist(lapply(mods, function(m) sprintf("%s.dhit.cv.r%d", m, hz))),
    unlist(lapply(mods, function(m) unlist(lapply(LB, function(g)
      sprintf("%s.mlp.%s_h%d", m, g, hz))))))
  xcols <- c(sprintf("%s.cox2.death.lp", mods),
             unlist(lapply(mods, function(m)
               sprintf("%s.cox2.death.r%d", m, hz))))
  A1 <- avg(LT, lcols)
  A2 <- avg(DT, dcols)
  A3 <- avg(XT, xcols)
  S <- c(as.list(A1[, ..lcols]), as.list(A2[, ..dcols]),
         as.list(A3[, ..xcols]))
  nmiss <- sum(sapply(S, function(v) sum(!is.finite(v))))
  chk(sprintf("%s: patients / missing scores", key),
      sprintf("%d / %d", n, nmiss), if (nmiss == 0) "PASS" else "FAIL")
  stop_if_fail()
  TV <- list(death = list(t = sub$time_death, e = sub$event_death == 1),
             composite = list(t = sub$time_first, e = sub$event_comp == 1),
             cv = list(t = sub$time_first, e = sub$event_cr == 1))
  Y <- list()
  for (en in names(TV)) for (h in hz) {
    Y[[paste(en, h)]] <- TV[[en]]$e & TV[[en]]$t <= h + 0.5
  }
  items <- list()
  add <- function(k, type, s, lab = NULL, en = NULL) {
    items[[k]] <<- list(type = type, s = s, lab = lab, en = en)
  }
  for (m in mods) {
    u <- toupper(m)
    for (en in names(TG)) {
      g <- TG[[en]]
      fam <- list("ridge-Cox" = sprintf("%s.cox.%s.lp", m, g),
                  "DeepSurv" = sprintf("%s.dsurv.%s.lp", m, g))
      if (en == "death") {
        fam[["ridge-Cox retuned"]] <- sprintf("%s.cox2.death.lp", m)
      }
      for (lr in names(fam)) {
        add(paste(u, en, "survival", lr, "all", "C", sep = "|"), "cidx",
            fam[[lr]], en = en)
        for (h in hz) {
          add(paste(u, en, "survival", lr, h, "AUROC", sep = "|"), "auc",
              fam[[lr]], lab = paste(en, h))
        }
      }
      for (h in hz) {
        if (en != "cv") {
          rk <- list("ridge-Cox" = sprintf("%s.cox.%s.r%d", m, g, h),
                     "DeepSurv" = sprintf("%s.dsurv.%s.r%d", m, g, h))
          if (en == "death") {
            rk[["ridge-Cox retuned"]] <- sprintf("%s.cox2.death.r%d", m, h)
          }
          for (lr in names(rk)) for (ms in c("IPA", "LL")) {
            add(paste(u, en, "survival", lr, h, ms, sep = "|"),
                tolower(ifelse(ms == "IPA", "ipa", "ll")), rk[[lr]],
                lab = paste(en, h))
          }
        } else {
          add(paste(u, en, "competing", "DeepHit", h, "AUROC", sep = "|"),
              "auc", sprintf("%s.dhit.cv.r%d", m, h), lab = paste(en, h))
          for (ms in c("IPA", "LL")) {
            add(paste(u, en, "competing", "DeepHit", h, ms, sep = "|"),
                tolower(ifelse(ms == "IPA", "ipa", "ll")),
                sprintf("%s.dhit.cv.r%d", m, h), lab = paste(en, h))
          }
        }
        bl <- list("L2 logistic" = sprintf("%s.lr.%s_h%d", m, LB[[en]], h),
                   "MLP twin" = sprintf("%s.mlp.%s_h%d", m, LB[[en]], h))
        for (lr in names(bl)) for (ms in c("AUROC", "IPA", "LL")) {
          add(paste(u, en, "binary", lr, h, ms, sep = "|"),
              c(AUROC = "auc", IPA = "ipa", LL = "ll")[[ms]], bl[[lr]],
              lab = paste(en, h))
        }
      }
    }
  }
  # fused death: Cox stack on each family's death scores
  fo <- LT[, .(seed, hadm_id, fold)]
  seeds <- sort(unique(fo$seed))
  famcols <- list("ridge-Cox" = "%s.cox.death.lp",
                  "ridge-Cox retuned" = "%s.cox2.death.lp",
                  "DeepSurv" = "%s.dsurv.death.lp")
  srcs <- list("ridge-Cox" = LT, "ridge-Cox retuned" = XT, "DeepSurv" = DT)
  for (fm in names(famcols)) {
    P <- matrix(0, n, length(hz))
    for (s in seeds) {
      x <- srcs[[fm]][srcs[[fm]]$seed == s][match(sub$hadm_id, hadm_id)]
      f <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
      SS <- lapply(mods, function(m) x[[sprintf(famcols[[fm]], m)]])
      for (k in sort(unique(f))) {
        tr <- which(f != k)
        te <- which(f == k)
        P[te, ] <- P[te, ] + matrix(stack_ph(SS, tr, te, sub$time_death,
                                             as.integer(sub$event_death ==
                                                          1), hz + 0.5),
                                    ncol = length(hz)) / length(seeds)
      }
    }
    for (h in seq_along(hz)) {
      nm <- sprintf("fused|%s|%d", fm, hz[h])
      S[[nm]] <- P[, h]
      for (ms in c("AUROC", "IPA", "LL")) {
        add(paste("FUSED", "death", "survival", fm, hz[h], ms, sep = "|"),
            c(AUROC = "auc", IPA = "ipa", LL = "ll")[[ms]], nm,
            lab = paste("death", hz[h]))
      }
    }
  }
  JOBS[[key]] <- list(fr = fr, co = co, n = n, hz = hz, S = S, T = TV,
                      Y = Y, items = items, mods = mods,
                      rows_by = split(seq_len(n), match(
                        sub$subject_id, unique(sub$subject_id))))
  cat(sprintf("  %s: %d patients, %d statistics\n", key, n, length(items)))
}

cat(sprintf("\n=== TIMING PROBE (%d workers) ===\n", workers))
plan(multisession, workers = workers)
proj <- 0
for (J in JOBS) {
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(workers * 2), boot_rep, J = J,
                          future.seed = 1L, future.packages = PKG))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  proj <- proj + dt * NBOOT / (workers * 2)
}
cat(sprintf("  projected bootstrap: about %.1f min\n", proj / 60))
if (probe_only) {
  plan(sequential)
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)

# ---------------- bootstrap ----------------
q2 <- function(v) unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
MET <- list()
CON <- list()
for (J in JOBS) {
  pt <- stat_vec(J, seq_len(J$n))
  B <- do.call(rbind, future_lapply(seq_len(NBOOT), boot_rep, J = J,
                                    future.seed = SEED,
                                    future.packages = PKG))
  cat(sprintf("  Frame %s %s done\n", J$fr, J$co))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    q <- q2(B[, k])
    MET[[length(MET) + 1]] <- data.table(
      frame = J$fr, cohort = J$co, modality = z[1], endpoint = z[2],
      family = z[3], learner = z[4], horizon = z[5], measure = z[6],
      estimate = pt[[k]], lo = q[1], hi = q[2])
  }
  pairs <- list(c("DeepSurv", "ridge-Cox", "survival"),
                c("DeepSurv", "ridge-Cox retuned", "survival"),
                c("MLP twin", "L2 logistic", "binary"))
  for (k in names(pt)) {
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    for (pp in pairs) {
      if (z[3] != pp[3] || z[4] != pp[1]) next
      kb <- paste(z[1], z[2], z[3], pp[2], z[5], z[6], sep = "|")
      if (!(kb %in% names(pt))) next
      q <- q2(B[, k] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, modality = z[1], endpoint = z[2],
        contrast = sprintf("%s - %s", pp[1], pp[2]), horizon = z[5],
        measure = z[6], estimate = pt[[k]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
    if (z[3] == "competing" && z[6] == "AUROC") {
      kb <- paste(z[1], "cv", "survival", "ridge-Cox", z[5], "AUROC",
                  sep = "|")
      q <- q2(B[, k] - B[, kb])
      CON[[length(CON) + 1]] <- data.table(
        frame = J$fr, cohort = J$co, modality = z[1], endpoint = "cv",
        contrast = "DeepHit - cause-specific ridge-Cox", horizon = z[5],
        measure = "AUROC", estimate = pt[[k]] - pt[[kb]], lo = q[1],
        hi = q[2], excludes_zero = q[1] > 0 | q[2] < 0)
    }
  }
}
plan(sequential)
mt <- rbindlist(MET)
cn <- rbindlist(CON)
ck <- rbindlist(CHK)
fwrite(mt, file.path(RS, "s4n_metrics.csv"))
fwrite(cn, file.path(RS, "s4n_contrasts.csv"))
fwrite(ck, file.path(RS, "s4n_checks.csv"))
cat("\n=== DEATH, AUROC: DeepSurv minus ridge-Cox (6a and retuned) ===\n")
x <- cn[endpoint == "death" & measure == "AUROC" &
          grepl("^DeepSurv", contrast)]
x[, horizon := factor(horizon, levels = as.character(
  sort(unique(as.numeric(horizon)))))]
print(dcast(x, frame + cohort + modality + contrast ~ horizon,
            value.var = "estimate")[, lapply(.SD, function(v) {
              if (is.numeric(v)) round(v, 4) else v
            })], row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
for (f in c("s4n_metrics.csv", "s4n_contrasts.csv", "s4n_checks.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "s4n_log.txt"), " (terminal log)\n")
cat("\nOutputs: s4n_metrics.csv, s4n_contrasts.csv, s4n_checks.csv\n")