#' Compare second-order (mstate2) and first-order (mstate) n-step predictions
#'
#' Quantifies how much the second-order, history-conditional n-step prediction
#' \eqn{P_{hj\ell}(1,n)} diverges from the first-order Aalen-Johansen
#' prediction \eqn{P(X_{n+1} = \ell \mid X_1 = j)}, which ignores the preceding
#' state \eqn{h} entirely. This is the comparison methodology motivated by the
#' fact that a second-order model needs \eqn{M} matrices of dimension
#' \eqn{M \times M} instead of one first-order matrix: before paying that cost,
#' it is useful to see, transition by transition, how much (and for how long)
#' knowing \eqn{h} actually changes the forecast relative to the first-order
#' baseline.
#'
#' \code{pt} must be an \code{mstate::probtrans} object for the first-order
#' model, fit with landmark time \code{predt = 1}, i.e.
#' \code{mstate::probtrans(msfit_object, predt = 1)}, so that its curves start
#' at the same "time 1" as \code{ckequations()}'s \code{n = 1}, and its list
#' components are indexed by starting state in the same order as
#' \code{x$states} (the usual \code{mstate} convention when the same
#' \code{transMat} / state ordering is used throughout; see
#' \code{\link{as_tmat}}). Only the S3 class and structure of \code{pt} are
#' used here (via list-indexing and its \code{pstate*}/\code{se*} columns), so
#' this function does not itself require the \pkg{mstate} package to be
#' installed -- only that \code{pt} was produced by it.
#'
#' @param x A "P2est" object (second-order fit).
#' @param pt An \code{mstate::probtrans} object for the first-order model (see
#'   Details).
#' @param h Vector of preceding states to compare against the first-order
#'   baseline (labels or indices).
#' @param j Current state shared across the compared paths.
#' @param l Target state.
#' @param nsteps Number of steps (default 9).
#' @param conf.level Confidence level for both curves. Defaults to
#'   \code{x$conf.level}.
#' @return An object of class "msm2pred" (as returned by \code{\link{compare2}}):
#'   one group per preceding state \code{h}, plus a \code{"1st order (mstate)"}
#'   baseline group, so that \code{plot()} and, for a single \code{h},
#'   \code{\link{overlap_step}} work directly. See also \code{\link{divergence}}.
#' @export
compare_order <- function(x, pt, h, j, l, nsteps = 9L, conf.level = NULL) {
  stopifnot(inherits(x, "P2est"))
  if (!inherits(pt, "probtrans"))
    stop("`pt` must be an mstate::probtrans object (see ?compare_order).", call. = FALSE)
  if (is.null(conf.level)) conf.level <- x$conf.level
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  states <- x$states
  ji <- .resolve(j, states); li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)
  if (ji > length(pt))
    stop("`pt` has ", length(pt), " starting-state component(s), but state '",
         j, "' resolves to index ", ji, ". Check `pt` was built over the same ",
         "state space/order as `x$states`.", call. = FALSE)

  pt_j <- pt[[ji]]
  pcol <- paste0("pstate", li); scol <- paste0("se", li)
  miss <- setdiff(c("time", pcol, scol), names(pt_j))
  if (length(miss))
    stop("`pt[[j]]` is missing expected mstate::probtrans column(s): ",
         paste(miss, collapse = ", "),
         ". Was `pt` computed with predt = 1, i.e. mstate::probtrans(msf, predt = 1)?",
         call. = FALSE)

  predt <- pt_j$time[1L]
  target_t <- predt + seq_len(nsteps)
  idx <- vapply(target_t, function(tt) {
    cand <- which(pt_j$time <= tt)
    if (length(cand)) max(cand) else NA_integer_
  }, integer(1))
  if (anyNA(idx))
    stop("`pt` does not extend far enough in time to cover nsteps = ", nsteps,
         " step(s) from predt = ", predt, ".", call. = FALSE)

  first <- data.frame(h = "1st order (mstate)", n = seq_len(nsteps),
                      estimate = pt_j[[pcol]][idx], se = pt_j[[scol]][idx])
  first$lower <- pmax(0, first$estimate - z * first$se)
  first$upper <- pmin(1, first$estimate + z * first$se)
  first$se <- NULL

  second <- do.call(rbind, lapply(h, function(hh) {
    out <- ckequations(x, hh, j, l, nsteps = nsteps, bounds = TRUE)
    out$h <- if (is.numeric(hh)) states[hh] else as.character(hh)
    out
  }))

  df <- rbind(first[, c("h", "n", "estimate", "lower", "upper")],
             second[, c("h", "n", "estimate", "lower", "upper")])
  rownames(df) <- NULL

  jlab <- if (is.numeric(j)) states[j] else as.character(j)
  llab <- if (is.numeric(l)) states[l] else as.character(l)
  structure(df, class = c("msm2pred", "data.frame"),
            j = jlab, l = llab, bounds = TRUE, estimator = x$estimator,
            conf.level = conf.level, baseline = "1st order (mstate)")
}

#' Memory horizon: how long second-order predictions diverge from mstate
#'
#' For each preceding state in an object returned by \code{\link{compare_order}},
#' reports how many leading steps its second-order evolution interval and the
#' first-order (mstate) confidence interval fail to overlap -- the "memory
#' horizon" over which knowing the preceding state changes the forecast -- by
#' calling \code{\link{overlap_step}} on each (baseline, h) pair in turn.
#'
#' @param x An object returned by \code{\link{compare_order}}.
#' @return A data frame with one row per preceding state \code{h} (excluding
#'   the baseline), sorted by descending \code{separated_steps}, with columns
#'   \code{h, n, s, separated_steps, max_abs_diff} (the largest absolute
#'   difference in point estimates between that history and the first-order
#'   baseline, over the whole horizon).
#' @export
divergence <- function(x) {
  stopifnot(inherits(x, "msm2pred"))
  base_lab <- attr(x, "baseline")
  if (is.null(base_lab))
    stop("`x` has no first-order baseline group; build it with compare_order().",
         call. = FALSE)
  groups <- setdiff(unique(x$h), base_lab)
  if (!length(groups))
    stop("No second-order preceding-state group(s) found besides the baseline.",
         call. = FALSE)

  rows <- lapply(groups, function(g) {
    sub <- x[x$h %in% c(base_lab, g), ]
    sub <- structure(sub, class = c("msm2pred", "data.frame"),
                     j = attr(x, "j"), l = attr(x, "l"), bounds = TRUE,
                     estimator = attr(x, "estimator"), conf.level = attr(x, "conf.level"))
    os <- overlap_step(sub)
    data.frame(h = g, n = os$n, s = os$s, separated_steps = os$separated_steps,
              max_abs_diff = max(abs(sub$estimate[sub$h == g] -
                                     sub$estimate[sub$h == base_lab])))
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$separated_steps), , drop = FALSE]
  rownames(out) <- NULL
  out
}
