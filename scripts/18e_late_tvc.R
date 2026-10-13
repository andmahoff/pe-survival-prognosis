# 18e_late_tvc.R
# Late-fusion Cox stack whose modality weights change with follow-up time
# (PREREGISTRATION_ADDENDUM_2.md, A13). Writes out-of-fold predictions only;
# 18c_arch_eval.R evaluates them with every other architecture.
#
# Frame A (death at 10 and 30 days) and Frame B (death at 10, 30, 90, 180
# and 365 days), both cohorts, A10 modality sets. The nested modality
# scores are computed as in 18a_early_linear.R (D9, D10, D11), so the script
# first rebuilds 18a's Late PH and Late TV and checks them.
#
# The stack is a Cox model on the winsorised modality scores Z, with hazard
# h0(t) exp(Z b + Z c g(t)): every score has a constant coefficient b and a
# coefficient c on a shared function of time g(t). It is fitted with coxph
# on data split at the interval boundaries (Efron ties). The candidates for
# g are none (the proportional-hazards stack, c = 0), step7 (0 before day 7,
# 1 from day 7), step30 (the same at day 30; Frame B only) and log (log
# time, constant within intervals bounded at 1, 2, 4, 7, 14, 30, 60, 120 and
# 240 days, Frame A up to 14, valued at each interval's geometric midpoint).
# Late TVC chooses g inside each training fold by the cross-validated
# partial log-likelihood (Verweij and van Houwelingen 1993) over 18a's
# inner folds. Late TVC log-time fixes g at log time (sensitivity analysis).
# EHR alone TVC applies the chosen-g model to the EHR score only, for the
# fusion gain. Risks at each horizon use the Breslow baseline with the
# time-varying linear predictor, centred so that it cannot overflow.
#
# Checks:  18a's Late PH and Late TV rebuilt, largest absolute risk
#          difference below 1e-6 (FAIL, and nothing is written); no missing
#          risks (FAIL); risks of exactly 0 or 1 (WARN).
# Notes:   stops if 18a predictions are missing or 18e predictions from an
#          earlier run remain. Tasks (frame, cohort, seed) run in parallel
#          with future.apply, heaviest first so that the inner folds match
#          18a, after a timing probe that pauses 5 s.
#
# Inputs:  data/processed/frame_A.csv, frame_B.csv, folds.csv,
#          feature_whitelist.csv, fusion_arch/arch_preds_linear_*.csv;
#          data/inputs/ehr28.csv, ecg_opt.csv, ctpa46.csv
# Outputs: data/processed/fusion_arch/arch_preds_tvc_*.csv (patient level,
#          not in the repo); results/tvc_forms.csv (share of folds
#          choosing each g), results/tvc_coefs.csv (b and c of each
#          modality in the log-time model), results/tvc_checks.csv,
#          results/tvc_warnings.csv, results/tvc_log.txt
# Usage:   Rscript scripts/18e_late_tvc.R [--probe-only]

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
OD <- file.path(PR, "fusion_arch")
HZ <- list(A = c(10, 30), B = c(10, 30, 90, 180, 365))
TCUTS <- list(A = c(1, 2, 4, 7, 14), B = c(1, 2, 4, 7, 14, 30, 60, 120, 240))
TMAX <- c(A = 31, B = 366)
FORMS <- list(A = c("none", "step7", "log"),
              B = c("none", "step7", "step30", "log"))
NLAM <- 80L
LMR <- 1e-7
PKG <- c("survival", "glmnet", "data.table")
options(width = 220)
workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t_start <- Sys.time()
CHK <- list()
chk <- function(item, value, status) {
  CHK[[length(CHK) + 1]] <<- data.table(item = item,
                                        value = as.character(value),
                                        status = status)
  cat(sprintf("  %-5s %-62s %s\n", status, item, value))
}

# ---------------- time-varying Cox stack (sent to workers) ----------------
gvals <- function(form, lb, mid) {
  switch(form, none = rep(0, length(lb)), step7 = as.numeric(lb >= 7),
         step30 = as.numeric(lb >= 30), log = log(mid))
}
tvc_design <- function(Z, ep, g) {
  Zr <- Z[ep$id, , drop = FALSE]
  if (all(g == 0)) Zr else cbind(Zr, Zr * g[ep$ep])
}
tvc_fit <- function(Z, ep, g, rows = NULL) {
  X <- tvc_design(Z, ep, g)
  k <- if (is.null(rows)) rep(TRUE, length(ep$id)) else rows
  f <- coxph(Surv(ep$tstart[k], ep$tstop[k], ep$ev[k]) ~
               X[k, , drop = FALSE], ties = "efron")
  b <- coef(f)
  b[is.na(b)] <- 0
  list(b = b, loglik = f$loglik[2], X = X)
}
tvc_loglik <- function(X, ep, b) {
  f <- suppressWarnings(coxph(Surv(ep$tstart, ep$tstop, ep$ev) ~ X,
                              ties = "efron", init = b,
                              control = coxph.control(iter.max = 0)))
  f$loglik[2]
}
tvc_cvpl <- function(Z, ep, g, fold_of_row) {
  tot <- 0
  for (j in sort(unique(fold_of_row))) {
    keep <- fold_of_row[ep$id] != j
    ft <- tvc_fit(Z, ep, g, keep)
    tot <- tot + tvc_loglik(ft$X, ep, ft$b) - ft$loglik
  }
  tot
}
tvc_risk <- function(Ztr, Zte, tt, ee, b, g, lb, cuts) {
  q <- ncol(Ztr)
  b0 <- b[seq_len(q)]
  b1 <- if (length(b) > q) b[q + seq_len(q)] else rep(0, q)
  l0 <- as.vector(Ztr %*% b0)
  l1 <- as.vector(Ztr %*% b1)
  n0 <- as.vector(Zte %*% b0)
  n1 <- as.vector(Zte %*% b1)
  o <- order(tt)
  ts <- tt[o]
  es <- ee[o]
  ut <- unique(ts[es == 1])
  d <- tabulate(match(ts[es == 1], ut), length(ut))
  pos <- match(ut, ts)
  pint <- findInterval(ut, lb, left.open = TRUE)
  inc <- numeric(length(ut))
  cp <- numeric(length(lb))
  for (p in unique(pint)) {
    lp <- l0[o] + l1[o] * g[p]
    cp[p] <- max(lp)
    rs <- rev(cumsum(rev(exp(lp - cp[p]))))
    sel <- which(pint == p)
    inc[sel] <- d[sel] / rs[pos[sel]]
  }
  sapply(cuts, function(cu) {
    Hx <- 0
    for (p in unique(pint)) {
      sel <- which(pint == p & ut <= cu)
      if (length(sel)) {
        Hx <- Hx + sum(inc[sel]) * exp(n0 + n1 * g[p] - cp[p])
      }
    }
    -expm1(-Hx)
  })
}

# ---------------- one task = frame x cohort x seed (all 5 folds) --------
tvc_task <- function(tk) {
  n <- length(tk$t)
  H <- tk$hz
  out <- list()
  forms <- list()
  coefs <- list()
  put <- function(model, set, h, idx, val) {
    k <- paste(model, set, h, sep = "|")
    if (is.null(out[[k]])) out[[k]] <<- rep(NA_real_, n)
    out[[k]][idx] <<- val
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
  wz <- function(v, rows) {
    q <- quantile(v[rows], c(0.005, 0.995), names = FALSE)
    w <- pmin(pmax(v, q[1]), q[2])
    (w - mean(w[rows])) / sd(w[rows])
  }
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
  lb <- c(0, tk$tcuts)
  ub <- c(tk$tcuts, tk$tmax)
  mid <- c(0.5, sqrt(lb[-1] * ub[-1]))
  warn <- character(0)
  for (k in sort(unique(tk$fold))) {
    t_fold <- Sys.time()
    nw0 <- length(warn)
    withCallingHandlers({
      tr <- which(tk$fold != k)
      te <- which(tk$fold == k)
      inner <- grouped(tr, tk$seed * 100 + k, 3)
      B <- lapply(tk$X, prep, tr = tr)
      tt <- tk$t
      ee <- tk$e
      cox_cv <- function(Xm, rows, fid) {
        cv.glmnet(Xm[rows, , drop = FALSE], Surv(tt[rows], ee[rows]),
                  family = "cox", alpha = 0, foldid = fid, nlambda = NLAM,
                  lambda.min.ratio = LMR, type.measure = "deviance",
                  penalty.factor = rep(1, ncol(Xm)), standardize = FALSE)
      }
      lpred <- function(cv, Xm, rows) {
        as.vector(predict(cv, Xm[rows, , drop = FALSE], s = "lambda.min"))
      }
      ph_risk <- function(lp_tr, lp_te) {
        c0 <- max(lp_tr)
        sapply(H, function(h) {
          hs <- breslow(tt[tr], ee[tr], lp_tr, h + 0.5)
          if (hs == 0) return(rep(0, length(lp_te)))
          -expm1(-hs * exp(lp_te - c0))
        })
      }
      # modality scores exactly as stage 18a's late arms
      mods <- unique(unlist(tk$sets))
      oof <- list()
      lte <- list()
      for (m in mods) {
        v <- rep(NA_real_, n)
        for (j in 1:3) {
          a <- tr[inner != j]
          b <- tr[inner == j]
          cv <- cox_cv(B[[m]], a, grouped(a, tk$seed * 1000 + k * 10 + j, 3))
          v[b] <- lpred(cv, B[[m]], b)
        }
        oof[[m]] <- v
        lte[[m]] <- lpred(cox_cv(B[[m]], tr, inner), B[[m]], te)
      }
      # training data split at the time boundaries
      ep <- survSplit(Surv(t, e) ~ ., cut = tk$tcuts, episode = "ep",
                      data = data.frame(t = tt[tr], e = ee[tr],
                                        id = seq_along(tr)))
      ep <- list(tstart = ep$tstart, tstop = ep$t, ev = ep$e, id = ep$id,
                 ep = ep$ep)
      fold_row <- inner
      sets <- c(tk$sets, list("EHR" = "ehr"))
      for (sn in names(sets)) {
        ms <- sets[[sn]]
        Ztr <- matrix(sapply(ms, function(m) wz(oof[[m]], tr)[tr]),
                      ncol = length(ms))
        Zte <- matrix(sapply(ms, function(m) {
          q <- quantile(oof[[m]][tr], c(0.005, 0.995), names = FALSE)
          w <- pmin(pmax(oof[[m]][tr], q[1]), q[2])
          (pmin(pmax(lte[[m]], q[1]), q[2]) - mean(w)) / sd(w)
        }), ncol = length(ms))
        late <- sn != "EHR"
        if (late) {
          # 18a's Late PH and Late TV, rebuilt as a check
          D <- data.frame(t = tt[tr], e = ee[tr], Ztr)
          bc <- coef(coxph(Surv(t, e) ~ ., data = D, ties = "efron"))
          rk <- matrix(ph_risk(as.vector(Ztr %*% bc),
                               as.vector(Zte %*% bc)), ncol = length(H))
          for (i in seq_along(H)) {
            y <- as.numeric(ee[tr] == 1 & tt[tr] <= H[i] + 0.5)
            gf <- suppressWarnings(glm.fit(cbind(1, Ztr), y,
                                           family = binomial()))
            cf <- gf$coefficients
            cf[is.na(cf)] <- 0
            put("check 18a Late PH", sn, H[i], te, rk[, i])
            put("check 18a Late TV", sn, H[i], te,
                plogis(as.vector(cbind(1, Zte) %*% cf)))
          }
        }
        # g chosen by cross-validated partial log-likelihood
        cvp <- sapply(tk$forms, function(fm) {
          tvc_cvpl(Ztr, ep, gvals(fm, lb, mid), fold_row)
        })
        best <- tk$forms[which.max(cvp)]
        forms[[length(forms) + 1]] <- data.table(
          set = sn, fold = k, form = best,
          cvpl_gain_over_none = max(cvp) - cvp[1])
        gb <- gvals(best, lb, mid)
        fb <- tvc_fit(Ztr, ep, gb)
        rb <- tvc_risk(Ztr, Zte, tt[tr], ee[tr], fb$b, gb, lb, H + 0.5)
        rb <- matrix(rb, ncol = length(H))
        lab <- if (late) "Late TVC" else "EHR alone TVC"
        for (i in seq_along(H)) put(lab, sn, H[i], te, rb[, i])
        if (late) {
          gl <- gvals("log", lb, mid)
          fl <- tvc_fit(Ztr, ep, gl)
          rl <- matrix(tvc_risk(Ztr, Zte, tt[tr], ee[tr], fl$b, gl, lb,
                                H + 0.5), ncol = length(H))
          for (i in seq_along(H)) {
            put("Late TVC log-time", sn, H[i], te, rl[, i])
          }
          q <- length(ms)
          coefs[[length(coefs) + 1]] <- data.table(
            set = sn, fold = k, modality = toupper(ms),
            b_constant = fl$b[seq_len(q)], c_log_time = fl$b[q + seq_len(q)])
        }
      }
    }, warning = function(w) {
      msg <- conditionMessage(w)
      if (!grepl("tie-handling", msg, fixed = TRUE)) warn <<- c(warn, msg)
      invokeRestart("muffleWarning")
    })
    if (!is.null(tk$progress)) {
      cat(sprintf("%s  %s %-7s seed %s fold %d done in %.0f s, %d warnings\n",
                  format(Sys.time(), "%H:%M:%S"), tk$fr, tk$co, tk$seed, k,
                  as.numeric(difftime(Sys.time(), t_fold, units = "secs")),
                  length(warn) - nw0),
          file = tk$progress, append = TRUE)
    }
  }
  list(fr = tk$fr, co = tk$co, seed = tk$seed, out = out, warn = warn,
       forms = rbindlist(forms)[, `:=`(frame = tk$fr, cohort = tk$co,
                                       seed = tk$seed)],
       coefs = rbindlist(coefs)[, `:=`(frame = tk$fr, cohort = tk$co,
                                       seed = tk$seed)])
}

# ---------------- build ----------------
cat("=== BUILD ===\n")
frames <- list(A = fread(file.path(PR, "frame_A.csv")),
               B = fread(file.path(PR, "frame_B.csv")))
folds <- fread(file.path(PR, "folds.csv"))
wl <- fread(file.path(PR, "feature_whitelist.csv"))
mats <- list(ehr = fread(file.path(IN, "ehr28.csv")),
             ecg = fread(file.path(IN, "ecg_opt.csv")),
             ctpa = fread(file.path(IN, "ctpa46.csv")))
mats <- lapply(mats, function(x) unique(x, by = "hadm_id"))
COH <- list("EHR+ECG" = list(code = "ehrecg",
                             sets = list("EHR+ECG" = c("ehr", "ecg"))),
            "3-mod" = list(code = "3mod",
                           sets = list("EHR+ECG" = c("ehr", "ecg"),
                                       "EHR+CTPA" = c("ehr", "ctpa"),
                                       "EHR+ECG+CTPA" =
                                         c("ehr", "ecg", "ctpa"))))
old <- list.files(OD, pattern = "^arch_preds_tvc_.*\\.csv$")
chk("previous stage 18e prediction files (archive them)", length(old),
    if (length(old) == 0) "PASS" else "FAIL")
tasks <- list()
META <- list()
for (fr in c("A", "B")) for (co in names(COH)) {
  Fd <- frames[[fr]]
  sub <- if (co == "EHR+ECG") Fd else Fd[in_three_mod == 1]
  lf <- file.path(OD, sprintf("arch_preds_linear_%s_%s.csv", fr,
                              COH[[co]]$code))
  chk(sprintf("%s %s: stage 18a predictions found", fr, co), basename(lf),
      if (file.exists(lf)) "PASS" else "FAIL")
  mods <- unique(unlist(COH[[co]]$sets))
  X <- list()
  for (m in mods) {
    cols <- intersect(wl[modality == m, column], names(mats[[m]]))
    Xm <- as.matrix(mats[[m]][match(sub$hadm_id, hadm_id), ..cols])
    storage.mode(Xm) <- "double"
    X[[m]] <- Xm
  }
  fo <- folds[frame == fr & cohort == co]
  META[[paste(fr, co)]] <- list(fr = fr, co = co, sub = sub, lf = lf,
                                code = COH[[co]]$code)
  for (s in sort(unique(fo$seed))) {
    fs <- fo[fo$seed == s][match(sub$hadm_id, hadm_id), fold]
    if (anyNA(fs)) stop("fold missing: ", fr, " ", co, " seed ", s)
    tasks[[length(tasks) + 1]] <- list(
      fr = fr, co = co, seed = s, fold = fs, subj = sub$subject_id,
      t = sub$time_death, e = as.integer(sub$event_death == 1),
      X = X, sets = COH[[co]]$sets, hz = HZ[[fr]], tcuts = TCUTS[[fr]],
      tmax = TMAX[[fr]], forms = FORMS[[fr]])
  }
}
cat(sprintf("  %d tasks (frame x cohort x seed), 5 folds each\n",
            length(tasks)))
if (any(rbindlist(CHK)$status == "FAIL")) {
  fwrite(rbindlist(CHK), file.path(RS, "tvc_checks.csv"))
  cat("\nSTOP: a check failed. see: tvc_checks.csv\n")
  quit(status = 1)
}

cat("\n=== TIMING PROBE: the slowest kinds of fit, per frame x cohort ===\n")
prep_fold <- function(X, tr) {
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
probe_kind <- function(tk) {
  k <- sort(unique(tk$fold))[1]
  tr <- which(tk$fold != k)
  us <- unique(tk$subj[tr])
  set.seed(tk$seed * 100 + k)
  f <- sample(rep_len(1:3, length(us)))
  inner <- f[match(tk$subj[tr], us)]
  B <- lapply(tk$X, prep_fold, tr = tr)
  tm <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    proc.time()[["elapsed"]] - t0
  }
  cx <- function(X) {
    suppressWarnings(cv.glmnet(
      X[tr, , drop = FALSE], Surv(tk$t[tr], tk$e[tr]), family = "cox",
      alpha = 0, foldid = inner, nlambda = NLAM, lambda.min.ratio = LMR,
      type.measure = "deviance", standardize = FALSE))
  }
  lb <- c(0, tk$tcuts)
  mid <- c(0.5, sqrt(lb[-1] * c(tk$tcuts[-1], tk$tmax)))
  ep <- survSplit(Surv(t, e) ~ ., cut = tk$tcuts, episode = "ep",
                  data = data.frame(t = tk$t[tr], e = tk$e[tr],
                                    id = seq_along(tr)))
  ep <- list(tstart = ep$tstart, tstop = ep$t, ev = ep$e, id = ep$id,
             ep = ep$ep)
  Z <- matrix(rnorm(3 * length(tr)), ncol = 3)
  tv <- tm(for (fm in tk$forms) tvc_cvpl(Z, ep, gvals(fm, lb, mid), inner))
  c(big = if (!is.null(B$ecg)) tm(cx(B$ecg)) else 0,
    small = tm(cx(B$ehr)), tvc = tv)
}
fits_per_fold <- function(tk) {
  mods <- unique(unlist(tk$sets))
  c(big = 4 * ("ecg" %in% mods),
    small = 4 + 4 * ("ctpa" %in% mods),
    tvc = 2 * length(tk$sets) + 1)
}
kind <- sapply(tasks, function(t) paste(t$fr, t$co))
reps <- tasks[!duplicated(kind)]
names(reps) <- unique(kind)
plan(multisession, workers = min(workers, length(reps)))
t0 <- Sys.time()
pt <- future_lapply(reps, probe_kind, future.packages = PKG,
                    future.seed = TRUE, future.chunk.size = 1)
plan(sequential)
fold_s <- c()
for (kd in names(reps)) {
  u <- pt[[kd]]
  nf <- fits_per_fold(reps[[kd]])
  fold_s[kd] <- sum(u * nf[names(u)])
  cat(sprintf(paste("  %-10s ECG Cox fit %5.1f s (x%d), EHR %4.1f s (x%d),",
                    "time-varying stack %4.1f s (x%d): a fold %4.1f min\n"),
              kd, u[["big"]], nf[["big"]], u[["small"]], nf[["small"]],
              u[["tvc"]], nf[["tvc"]], fold_s[kd] / 60))
}
cost <- fold_s[kind] * 5
tasks <- tasks[order(-cost)]
cost <- sort(cost, decreasing = TRUE)
cat(sprintf(paste("  probe took %.1f min; %d tasks on %d workers: about",
                  "%.1f hours (longest task %.1f hours)\n"),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            length(tasks), workers,
            max(sum(cost) / workers, max(cost)) * 1.2 / 3600,
            max(cost) / 3600))
if (probe_only) {
  cat("--probe-only: stopping here.\n")
  quit(status = 0)
}
PROG <- file.path(OD, "progress_18e.txt")
cat(sprintf("%s  started, %d tasks x 5 folds\n",
            format(Sys.time(), "%H:%M:%S"), length(tasks)), file = PROG)
for (i in seq_along(tasks)) tasks[[i]]$progress <- PROG
cat("  progress: one line per finished fold in", PROG, "\n")
cat("  starting in 5 s (Ctrl+C to stop)...\n")
Sys.sleep(5)
plan(multisession, workers = workers)
res <- future_lapply(tasks, tvc_task, future.packages = PKG,
                     future.seed = TRUE, future.chunk.size = 1)
plan(sequential)

# ---------------- assemble and check ----------------
cat("\n=== CHECKS ===\n")
NEW <- list()
for (key in names(META)) {
  M <- META[[key]]
  rr <- Filter(function(r) r$fr == M$fr && r$co == M$co, res)
  rr <- rr[order(sapply(rr, function(r) r$seed))]
  keys <- names(rr[[1]]$out)
  L <- rbindlist(lapply(keys, function(k) {
    v <- Reduce(`+`, lapply(rr, function(r) r$out[[k]])) / length(rr)
    z <- strsplit(k, "|", fixed = TRUE)[[1]]
    data.table(hadm_id = M$sub$hadm_id, model = z[1], set = z[2],
               horizon = as.numeric(z[3]), risk = v)
  }))
  old <- fread(M$lf)
  for (mn in c("Late PH", "Late TV")) {
    cmp <- merge(L[model == paste("check 18a", mn),
                   .(hadm_id, set, horizon, new = risk)],
                 old[model == mn, .(hadm_id, set, horizon, old = risk)],
                 by = c("hadm_id", "set", "horizon"))
    dmax <- max(abs(cmp$new - cmp$old))
    chk(sprintf("%s: stage 18a %s rebuilt (rows, max |diff|)", key, mn),
        sprintf("%d, %.1e", nrow(cmp), dmax),
        if (nrow(cmp) == nrow(old[model == mn]) && dmax < 1e-6) "PASS"
        else "FAIL")
  }
  L <- L[!grepl("^check ", model)]
  nna <- sum(is.na(L$risk))
  nsat <- sum(L$risk <= 0 | L$risk >= 1, na.rm = TRUE)
  rg <- range(L$risk, na.rm = TRUE)
  chk(sprintf("%s: missing predictions", key), nna,
      if (nna == 0) "PASS" else "FAIL")
  chk(sprintf("%s: risks of exactly 0 or 1 (range)", key),
      sprintf("%d (%.2e to %.6f)", nsat, rg[1], rg[2]),
      if (nsat == 0) "PASS" else "WARN")
  NEW[[key]] <- L
}
WRN <- rbindlist(lapply(res, function(r) {
  if (!length(r$warn)) return(NULL)
  data.table(frame = r$fr, cohort = r$co, message = substr(r$warn, 1, 150))
}))
if (nrow(WRN)) {
  wsum <- WRN[, .(count = .N), by = .(frame, cohort, message)]
  setorder(wsum, frame, cohort, -count)
} else {
  wsum <- data.table(frame = character(0), cohort = character(0),
                     message = character(0), count = integer(0))
}
fwrite(wsum, file.path(RS, "tvc_warnings.csv"))
chk("warnings from model fitting (see tvc_warnings.csv)",
    sum(wsum$count), "INFO")
ck <- rbindlist(CHK)
if (any(ck$status == "FAIL")) {
  fwrite(ck, file.path(RS, "tvc_checks.csv"))
  cat("\nSTOP: a check failed; no predictions written.",
      "Outputs: tvc_checks.csv\n")
  quit(status = 1)
}
for (key in names(META)) {
  M <- META[[key]]
  fwrite(NEW[[key]], file.path(OD, sprintf("arch_preds_tvc_%s_%s.csv",
                                           M$fr, M$code)))
}
FRM <- rbindlist(lapply(res, function(r) r$forms))
fsum <- FRM[, .(folds = .N, median_cvpl_gain = median(cvpl_gain_over_none)),
            by = .(frame, cohort, set, form)]
fsum[, share := folds / sum(folds), by = .(frame, cohort, set)]
setorder(fsum, frame, cohort, set, -folds)
fwrite(fsum, file.path(RS, "tvc_forms.csv"))
CF <- rbindlist(lapply(res, function(r) r$coefs))
csum <- CF[, .(folds = .N, b_constant = median(b_constant),
               c_log_time = median(c_log_time),
               c_q25 = quantile(c_log_time, 0.25),
               c_q75 = quantile(c_log_time, 0.75)),
           by = .(frame, cohort, set, modality)]
fwrite(csum, file.path(RS, "tvc_coefs.csv"))
fwrite(ck, file.path(RS, "tvc_checks.csv"))
cat("\n=== TIME FUNCTION CHOSEN (share of folds) ===\n")
print(fsum[, .(frame, cohort, set, form, share = round(share, 2))],
      row.names = FALSE)
cat("\n=== LOG-TIME MODEL: median coefficients per modality ===\n")
print(csum[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)],
      row.names = FALSE)
cat(sprintf("\n  checks: %d PASS, %d WARN, %d FAIL; elapsed %.1f min\n",
            sum(ck$status == "PASS"), sum(ck$status == "WARN"),
            sum(ck$status == "FAIL"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nOUTPUT FILES\n")
cat("  ", OD, " (arch_preds_tvc_*.csv, patient-level, local only)\n")
for (f in c("tvc_forms.csv", "tvc_coefs.csv", "tvc_checks.csv",
            "tvc_warnings.csv")) {
  cat("  ", file.path(RS, f), "\n")
}
cat("  ", file.path(RS, "tvc_log.txt"), " (terminal log)\n")
cat("\nOutputs: tvc_forms.csv, tvc_coefs.csv, tvc_checks.csv,",
    "tvc_warnings.csv\n")