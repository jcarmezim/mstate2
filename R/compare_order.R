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
#' @param pt The first-order baseline: either an \code{mstate::probtrans}
#'   object (see Details), or a \code{\link{P1est}} fit with bootstrap
#'   replicates, the discrete-time first-order chain estimated from the same
#'   data (no \pkg{mstate} fit needed; the baseline group is then labelled
#'   \code{"1st order (discrete)"}).
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
#' @examples
#' st   <- c("A", "B", "C")                                 # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' fit <- P2est(prep2(panel))
#'
#' # A minimal stand-in for an mstate::probtrans object, built without the
#' # mstate package: it pools the RPE over the preceding state h (i.e.
#' # deliberately discards the history mstate2 otherwise keeps). With real,
#' # progressive multistate data, build `pt` from an actual first-order
#' # mstate fit instead (see Details):
#' #   tmat  <- as_tmat(d)
#' #   msdat <- mstate::msprep(time = ..., status = ..., data = ...,
#' #                          trans = tmat, keep = "age")
#' #   cx    <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~
#' #                            strata(trans), data = msdat, method = "breslow")
#' #   msf   <- mstate::msfit(cx, trans = tmat)
#' #   pt    <- mstate::probtrans(msf, predt = 1)
#' est   <- fit$estimate
#' denom <- sum(unique(est[est$j == "B", c("h", "at.risk")])$at.risk)
#' p_BB  <- sum(est$n.trans[est$j == "B" & est$l == "B"]) / denom
#' p_BC  <- sum(est$n.trans[est$j == "B" & est$l == "C"]) / denom
#' se_BB <- sqrt(p_BB * (1 - p_BB) / denom)
#' se_BC <- sqrt(p_BC * (1 - p_BC) / denom)
#'
#' P1 <- matrix(0, 3, 3, dimnames = list(st, st))
#' P1["A", "A"] <- 1; P1["B", "B"] <- p_BB; P1["B", "C"] <- p_BC; P1["C", "C"] <- 1
#' P1_tensor <- array(rep(as.vector(P1), 3), c(3, 3, 3), dimnames = list(st, st, st))
#'
#' nsteps <- 6
#' curve <- ckequations(P1_tensor, h = "A", j = "B", l = NULL, nsteps = nsteps)
#' pt_B  <- data.frame(time = c(1, 1 + seq_len(nsteps)))
#' mat   <- rbind(c(A = 0, B = 1, C = 0), curve[, st])
#' pt_B$pstate1 <- mat[, "A"]; pt_B$se1 <- 0
#' pt_B$pstate2 <- mat[, "B"]; pt_B$se2 <- se_BB
#' pt_B$pstate3 <- mat[, "C"]; pt_B$se3 <- se_BC
#' pt <- vector("list", 3); pt[[2]] <- pt_B; class(pt) <- "probtrans"
#'
#' cmp2 <- compare_order(fit, pt, h = c("A", "B"), j = "B", l = "B", nsteps = nsteps)
#' cmp2
#' @export
compare_order <- function(x, pt, h, j, l, nsteps = 9L, conf.level = NULL) {
  stopifnot(inherits(x, "P2est"))
  ## A P1est fit is turned into the equivalent probtrans object (first-order
  ## n-step predictions from j, with bootstrap SEs) on the same time grid.
  base_lab <- "1st order (mstate)"
  p1fit <- if (inherits(pt, "P1est")) pt else NULL
  if (inherits(pt, "P1est")) {
    if (is.null(pt$boot))
      stop("`pt` (P1est) has no bootstrap replicates, so the first-order ",
           "predictions have no standard errors; use P1est(..., B > 0).", call. = FALSE)
    pt <- probtrans2(pt, j = j, nsteps = nsteps)
    base_lab <- "1st order (discrete)"
  }
  if (!inherits(pt, "probtrans"))
    stop("`pt` must be an mstate::probtrans object or a P1est fit (see ?compare_order).",
         call. = FALSE)
  if (is.null(conf.level)) conf.level <- x$conf.level
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  states <- x$states
  ji <- .resolve(j, states); li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)
  ## `pt` is a list indexed by *starting* state, mstate's usual convention
  ## for a probtrans object -- pt[[ji]] should be the first-order model's
  ## predicted trajectory starting from state j. If pt has fewer components
  ## than x$states does, j's index into x$states can't possibly line up with
  ## pt's own indexing, which almost always means the two were built over
  ## different (or differently-ordered) state spaces.
  if (ji > length(pt))
    stop("`pt` has ", length(pt), " starting-state component(s), but state '",
         j, "' resolves to index ", ji, ". Check `pt` was built over the same ",
         "state space/order as `x$states`.", call. = FALSE)

  pt_j <- pt[[ji]]
  ## mstate::probtrans() names its predicted-probability and standard-error
  ## columns pstate<k>/se<k> by the *destination* state's index k -- here
  ## that's li, the position of l in the (shared) state space.
  pcol <- paste0("pstate", li); scol <- paste0("se", li)
  miss <- setdiff(c("time", pcol, scol), names(pt_j))
  if (length(miss))
    stop("`pt[[j]]` is missing expected mstate::probtrans column(s): ",
         paste(miss, collapse = ", "),
         ". Was `pt` computed with predt = 1, i.e. mstate::probtrans(msf, predt = 1)?",
         call. = FALSE)

  if (!nrow(pt_j) || anyNA(pt_j$time))
    stop("`pt[[j]]` has no valid `time` values.", call. = FALSE)
  ## `predt`, the landmark time pt was computed from, anchors the mapping
  ## between ckequations()'s step count n (n = 1, 2, ...) and pt's own
  ## absolute time scale: step n corresponds to calendar time predt + n.
  predt <- pt_j$time[1L]
  target_t <- predt + seq_len(nsteps)
  ## Bail out clearly if the caller asked for more steps than `pt` actually
  ## covers, rather than silently reusing (flat-lining) whatever the last
  ## available row happens to be for every step beyond pt's real range.
  if (max(pt_j$time) < max(target_t))
    stop("`pt` does not extend far enough in time to cover nsteps = ", nsteps,
         " step(s) from predt = ", predt, " (`pt` only reaches time = ",
         max(pt_j$time), ").", call. = FALSE)
  ## mstate::probtrans() reports a right-continuous step function of
  ## cumulative transition probabilities, so within the range just checked,
  ## carrying the last observed row forward to each target time is correct:
  ## for each target time, find the most recent pt row at or before it.
  idx <- vapply(target_t, function(tt) max(which(pt_j$time <= tt)), integer(1))

  ## The first-order baseline curve, in the same (h, n, estimate, lower,
  ## upper) shape compare2() produces, so it can be row-bound onto the
  ## second-order curves below and handled identically by every downstream
  ## msm2pred method (print/summary/plot/overlap_step).
  first <- data.frame(h = base_lab, n = seq_len(nsteps),
                      estimate = pt_j[[pcol]][idx], se = pt_j[[scol]][idx])
  first$lower <- pmax(0, first$estimate - z * first$se)
  first$upper <- pmin(1, first$estimate + z * first$se)
  first$se <- NULL

  ## One second-order evolution curve per requested preceding state h, each
  ## computed by ckequations() exactly as compare2() would.
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
  ## Same "msm2pred" class as compare2()'s output (so plot()/print() etc.
  ## work unchanged), plus a `baseline` attribute recording which group is
  ## the first-order curve -- that's what lets divergence() automatically
  ## pair every second-order group against the baseline rather than the
  ## caller having to say so again.
  ## Paired bootstrap: when the first-order fit was built from the same
  ## P2boot replicates (P1est(<P2boot>)), keep the replicate curves of every
  ## group, so that divergence() can test the difference itself rather than
  ## only the overlap of two separate intervals.
  curves <- NULL
  if (inherits(x, "P2boot") && !is.null(p1fit$boot) &&
      identical(dim(p1fit$boot), dim(x$boot1)) &&
      isTRUE(all.equal(p1fit$boot, x$boot1, check.attributes = FALSE))) {
    curves <- list()
    curves[[base_lab]] <- .boot_curves_first(x$boot1, ji, li, nsteps)
    for (hh in h) {
      hlab <- if (is.numeric(hh)) states[hh] else as.character(hh)
      curves[[hlab]] <- .boot_curves(x$boot, .resolve(hh, states), ji, li, nsteps)
    }
  }
  structure(df, class = c("msm2pred", "data.frame"),
            j = jlab, l = llab, bounds = TRUE, estimator = x$estimator,
            conf.level = conf.level, baseline = base_lab,
            bands = if (inherits(x, "P2boot")) "bootstrap" else "evolution",
            boot_curves = curves)
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
#'   \code{h, n, s, separated_steps, diff_steps, max_abs_diff}:
#'   \code{separated_steps} counts the leading steps whose two intervals do not
#'   overlap; \code{diff_steps} counts the leading steps whose percentile
#'   bootstrap interval of the paired difference excludes 0 (available when
#'   \code{x} is a \code{\link{P2boot}} fit and the baseline is
#'   \code{P1est(x)}, i.e. both curves come from the same resamples;
#'   \code{NA} otherwise); \code{max_abs_diff} is the largest absolute
#'   difference in point estimates over the whole horizon. The paired
#'   difference is the direct test of whether the two predictions differ; the
#'   overlap criterion is more conservative (Schenker and Gentleman, 2001).
#' @examples
#' st   <- c("A", "B", "C")                                 # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' fit <- P2est(prep2(panel))
#'
#' # see ?compare_order for what this stand-in `pt` represents
#' est   <- fit$estimate
#' denom <- sum(unique(est[est$j == "B", c("h", "at.risk")])$at.risk)
#' p_BB  <- sum(est$n.trans[est$j == "B" & est$l == "B"]) / denom
#' p_BC  <- sum(est$n.trans[est$j == "B" & est$l == "C"]) / denom
#' se_BB <- sqrt(p_BB * (1 - p_BB) / denom)
#' se_BC <- sqrt(p_BC * (1 - p_BC) / denom)
#'
#' P1 <- matrix(0, 3, 3, dimnames = list(st, st))
#' P1["A", "A"] <- 1; P1["B", "B"] <- p_BB; P1["B", "C"] <- p_BC; P1["C", "C"] <- 1
#' P1_tensor <- array(rep(as.vector(P1), 3), c(3, 3, 3), dimnames = list(st, st, st))
#'
#' nsteps <- 6
#' curve <- ckequations(P1_tensor, h = "A", j = "B", l = NULL, nsteps = nsteps)
#' pt_B  <- data.frame(time = c(1, 1 + seq_len(nsteps)))
#' mat   <- rbind(c(A = 0, B = 1, C = 0), curve[, st])
#' pt_B$pstate1 <- mat[, "A"]; pt_B$se1 <- 0
#' pt_B$pstate2 <- mat[, "B"]; pt_B$se2 <- se_BB
#' pt_B$pstate3 <- mat[, "C"]; pt_B$se3 <- se_BC
#' pt <- vector("list", 3); pt[[2]] <- pt_B; class(pt) <- "probtrans"
#'
#' cmp2 <- compare_order(fit, pt, h = c("A", "B"), j = "B", l = "B", nsteps = nsteps)
#' divergence(cmp2)
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

  ## Run overlap_step() once per (baseline, h) pair -- it only ever compares
  ## exactly two groups, so with more than one second-order h present in `x`
  ## they have to be peeled off and paired against the baseline one at a time.
  rows <- lapply(groups, function(g) {
    sub <- x[x$h %in% c(base_lab, g), ]
    ## Subsetting a data.frame with `[` drops the custom attributes compare2()/
    ## compare_order() attach (j, l, bounds, estimator, conf.level), so they
    ## have to be reattached by hand before overlap_step() -- which requires
    ## a proper "msm2pred" object -- can be called on this two-group subset.
    sub <- structure(sub, class = c("msm2pred", "data.frame"),
                     j = attr(x, "j"), l = attr(x, "l"), bounds = TRUE,
                     estimator = attr(x, "estimator"), conf.level = attr(x, "conf.level"),
                     boot_curves = attr(x, "boot_curves"))
    os <- overlap_step(sub)
    data.frame(h = g, n = os$n, s = os$s, separated_steps = os$separated_steps,
              diff_steps = os$diff_steps,
              max_abs_diff = max(abs(sub$estimate[sub$h == g] -
                                     sub$estimate[sub$h == base_lab])))
  })
  out <- do.call(rbind, rows)
  ## Longest-diverging preceding state first: the one(s) knowing h changes
  ## the forecast for the longest are the most informative to report first.
  out <- out[order(-out$separated_steps), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Systematic first- versus second-order comparison across all transitions
#'
#' Applies \code{\link{compare_order}} and \code{\link{divergence}} to every
#' current state \eqn{j} reached from at least \code{min.h} different preceding
#' states and to every target state \eqn{\ell}, and stacks the results: for
#' each \eqn{(j, \ell, h)}, how many leading steps the second-order prediction
#' is separated from the first-order one (the memory horizon) and the largest
#' absolute difference between them. It answers, for the whole model at once,
#' for which transitions and up to which horizon conditioning on the preceding
#' state matters.
#'
#' @param x A "P2est" or "P2boot" object (with a "P2boot", the second-order
#'   bands are percentile bootstrap intervals).
#' @param pt First-order baseline: a \code{\link{P1est}} fit with bootstrap
#'   replicates (recommended: same data and time grid), or an
#'   \code{mstate::probtrans} object computed with \code{predt = 1}.
#' @param nsteps Horizon. Default 9.
#' @param min.h Minimum number of preceding states for a current state to be
#'   included. Default 2.
#' @param targets Optional target states (default: every state).
#' @return A data frame with one row per (j, l, h) and columns \code{j, l, h,
#'   n, s, separated_steps, diff_steps, max_abs_diff} (see \code{\link{divergence}}),
#'   sorted by decreasing \code{separated_steps} and \code{max_abs_diff}.
#'   Targets whose probability is zero for every history and the baseline
#'   over the whole horizon are omitted.
#' @examples
#' st   <- c("A", "B", "C")
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' d <- prep2(panel)
#' divergence_table(P2est(d), P1est(d, B = 100, seed = 1), nsteps = 5)
#' @export
divergence_table <- function(x, pt, nsteps = 9L, min.h = 2, targets = NULL) {
  stopifnot(inherits(x, "P2est"))
  est <- x$estimate
  states <- x$states
  if (is.null(targets)) targets <- states
  cur <- setdiff(unique(as.character(est$j)), x$absorbing)
  rows <- list()
  for (jj in cur) {
    hs <- sort(unique(as.character(est$h[as.character(est$j) == jj])))
    if (length(hs) < min.h) next
    for (ll in targets) {
      co <- compare_order(x, pt, h = hs, j = jj, l = ll, nsteps = nsteps)
      if (all(co$estimate == 0)) next                       # never reachable
      dv <- divergence(co)
      rows[[length(rows) + 1L]] <- data.frame(j = jj, l = ll, dv)
    }
  }
  if (!length(rows))
    stop("No current state is reached from at least `min.h` preceding states.",
         call. = FALSE)
  out <- do.call(rbind, rows)
  out <- out[order(-out$separated_steps, -out$max_abs_diff), , drop = FALSE]
  rownames(out) <- NULL
  out
}
