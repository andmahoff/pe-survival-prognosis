# tune_expand.R
# Grid search that extends itself past an edge (A28).
#
# Helper sourced by the level-1 tuning scripts 06i, 06j, 06p, 06q and 06r
# and by 21b. It is used inside each training fold only, so test folds never
# influence it.
#
# tune_expand(axes, eval_fun, rules, tol = 0.001, max_ext = 3)
#   axes      named list: each hyperparameter's starting grid (numbers)
#   eval_fun  function(par) -> list(score = <lower is better>, ...); par
#             is a named list with one value per hyperparameter. An error,
#             NULL or a non-finite score counts as a failed fit.
#   rules     named list, one entry per hyperparameter:
#             list(down = function(v), up = function(v), lo = , hi = )
#             giving the next value beyond each edge and hard bounds
#   tol       an extension is kept only if it improves the best score by
#             more than tol x |best score| (so flat surfaces are not
#             chased)
#   max_ext   at most this many kept extensions per hyperparameter and
#             direction
#
# The full starting grid is evaluated first. Then, while the best value of a
# hyperparameter is the smallest (largest) value tried for it, the next
# value beyond it is tried with the other hyperparameters held at their
# best. It is kept if it beats the best by more than the tolerance, and the
# search continues in that direction; otherwise that direction stops. The
# search ends when no hyperparameter is at an edge, or a hard bound or
# max_ext is reached.
#
# The return value is list(best = par, score, fit = eval_fun's result for
# the best, path = every evaluation (with any error message), ext = every
# extension tried, edge = named logical, TRUE where the final best is still
# at an edge of the values tried (a hard bound or max_ext was reached while
# still improving), nfail = failed evaluations).
#
# Usage:   source(file.path(P2, "scripts", "tune_expand.R"))

tune_expand <- function(axes, eval_fun, rules, tol = 0.001, max_ext = 3L) {
  nm <- names(axes)
  path <- list()
  ext <- list()
  nfail <- 0L
  best <- NULL
  bscore <- Inf
  bfit <- NULL
  tried <- lapply(axes, function(v) sort(unique(v)))
  evaluate <- function(par, step) {
    msg <- ""
    r <- tryCatch(eval_fun(par), error = function(e) {
      msg <<- substr(conditionMessage(e), 1, 120)
      NULL
    })
    sc <- if (is.null(r) || is.null(r$score) || !is.finite(r$score)) {
      NA_real_
    } else {
      r$score
    }
    if (is.na(sc) && msg == "") msg <- "non-finite score"
    path[[length(path) + 1]] <<- data.table::as.data.table(
      c(par, list(score = sc, step = step, error = msg)))
    if (is.na(sc)) {
      nfail <<- nfail + 1L
      return(NULL)
    }
    list(score = sc, fit = r)
  }
  grid <- expand.grid(axes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    par <- as.list(grid[i, , drop = FALSE])
    names(par) <- nm
    e <- evaluate(par, "grid")
    if (!is.null(e) && e$score < bscore) {
      bscore <- e$score
      best <- par
      bfit <- e$fit
    }
  }
  if (is.null(best)) {
    return(list(best = NULL, score = NA_real_, fit = NULL,
                path = data.table::rbindlist(path, fill = TRUE),
                ext = data.table::data.table(), edge = NULL,
                nfail = nfail))
  }
  count <- lapply(nm, function(p) c(down = 0L, up = 0L))
  done <- lapply(nm, function(p) c(down = FALSE, up = FALSE))
  names(count) <- names(done) <- nm
  repeat {
    moved <- FALSE
    for (p in nm) for (dir in c("down", "up")) {
      if (done[[p]][[dir]]) next
      v <- best[[p]]
      at_edge <- if (dir == "down") v <= min(tried[[p]]) else
        v >= max(tried[[p]])
      if (!at_edge) next
      if (count[[p]][[dir]] >= max_ext) {
        done[[p]][[dir]] <- TRUE
        next
      }
      nv <- rules[[p]][[dir]](v)
      nv <- min(max(nv, rules[[p]]$lo), rules[[p]]$hi)
      if (!is.finite(nv) || nv %in% tried[[p]]) {
        done[[p]][[dir]] <- TRUE
        next
      }
      par <- best
      par[[p]] <- nv
      tried[[p]] <- sort(unique(c(tried[[p]], nv)))
      e <- evaluate(par, sprintf("extend %s %s", p, dir))
      keep <- !is.null(e) && e$score < bscore - tol * abs(bscore)
      ext[[length(ext) + 1]] <- data.table::data.table(
        param = p, direction = dir, from = v, to = nv,
        score = if (is.null(e)) NA_real_ else e$score,
        best_before = bscore, kept = keep)
      if (keep) {
        count[[p]][[dir]] <- count[[p]][[dir]] + 1L
        bscore <- e$score
        best <- par
        bfit <- e$fit
        moved <- TRUE
      } else {
        done[[p]][[dir]] <- TRUE
      }
    }
    if (!moved) break
  }
  edge <- sapply(nm, function(p) {
    best[[p]] <= min(tried[[p]]) || best[[p]] >= max(tried[[p]])
  })
  list(best = best, score = bscore, fit = bfit,
       path = data.table::rbindlist(path, fill = TRUE),
       ext = data.table::rbindlist(ext), edge = edge, nfail = nfail)
}