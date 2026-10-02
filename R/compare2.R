#' Compare n-step transitions across preceding states (evolution intervals)
#'
#' For a fixed current state \eqn{j} and target \eqn{\ell}, computes the n-step
#' transition probabilities for several preceding states \eqn{h}, with their
#' evolution intervals. The construction behind Figures 4-5: it assesses whether,
#' and for how long, the preceding state affects the future trajectory.
#'
#' @param object A "P2est" object.
#' @param h Vector of preceding states to compare (labels or indices).
#' @param j Current state shared across the compared paths.
#' @param l Target state.
#' @param nsteps Number of steps (default 9).
#' @param bounds If TRUE (default) also compute evolution-interval bounds; FALSE
#'   returns curves only (e.g. Figure-4 overlays) and skips two thirds of the work.
#'   If \code{object} is a \code{\link{P2boot}} fit, the bounds are percentile
#'   bootstrap intervals instead of evolution intervals.
#' @return An object of class "msm2pred" (a data frame): columns h, n, estimate
#'   and, when bounds = TRUE, lower and upper.
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
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' summary(cmp)
#' @export
compare2 <- function(object, h, j, l, nsteps = 9L, bounds = TRUE) {
  stopifnot(inherits(object, "P2est"))
  states <- object$states; M <- length(states)
  lab <- function(s) if (is.numeric(s)) states[s] else as.character(s)
  ji <- .resolve(j, states); li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)

  ## build pair-transition matrices ONCE, reuse across every h
  ## (.pair_matrix() is an O(M^2) build; doing it once here instead of once
  ## per h avoids repeating that work for every compared preceding state).
  Q  <- .pair_matrix(object$P, M)
  boot <- inherits(object, "P2boot")
  if (bounds && !boot) {
    Ql <- .pair_matrix(object$P.lower, M)
    Qu <- .pair_matrix(object$P.upper, M)
  }

  ## One n-step curve (and, if requested, its evolution-interval bounds) per
  ## preceding state h, stacked into a single long data frame with an `h`
  ## column identifying which group each row belongs to -- the layout
  ## overlap_step(), summary() and plot() all expect.
  curves <- list()
  parts <- lapply(h, function(hh) {
    hi <- .resolve(hh, states)
    if (is.na(hi)) stop("preceding state '", hh, "' not found.", call. = FALSE)
    out <- data.frame(h = lab(hh), n = seq_len(nsteps),
                      estimate = .propagate(Q, hi, ji, nsteps, M)[, li])
    if (bounds && boot) {
      ## percentile bootstrap intervals (P2boot): propagate every replicate.
      ## The replicate curves are kept, so that overlap_step() can also test
      ## the paired difference between groups (same resamples).
      cv <- .boot_curves(object$boot, hi, ji, li, nsteps)
      curves[[lab(hh)]] <<- cv
      bb <- .boot_bands(object$boot, hi, ji, li, nsteps, object$conf.level, curves = cv)
      out$lower <- bb$lower; out$upper <- bb$upper
    } else if (bounds) {
      ## clipped to [0, 1], as in ckequations(): the limit tensors' rows do
      ## not sum to 1, so unclipped bounds can leave the probability scale
      out$lower <- pmax(0, .propagate(Ql, hi, ji, nsteps, M)[, li])
      out$upper <- pmin(1, .propagate(Qu, hi, ji, nsteps, M)[, li])
    }
    out
  })

  df <- do.call(rbind, parts); rownames(df) <- NULL
  ## j, l, bounds, estimator and conf.level travel as attributes (not
  ## columns) because they describe the comparison as a whole, not any one
  ## row; print/summary/plot methods and overlap_step() all read them back
  ## via attr(). Subsetting a data.frame with `[` drops custom attributes,
  ## which is why divergence() has to reattach them explicitly after
  ## subsetting (see compare_order.R).
  structure(df, class = c("msm2pred", "data.frame"),
            j = lab(j), l = lab(l), bounds = bounds,
            bands = if (!bounds) "none" else if (boot) "bootstrap" else "evolution",
            estimator = object$estimator, conf.level = object$conf.level,
            boot_curves = if (length(curves)) curves else NULL)
}

#' First step at which two evolution intervals overlap
#'
#' @param x A two-group "msm2pred" with evolution-interval bounds.
#' Non-overlap of two 95\% intervals is the criterion of the methods paper
#' (Section 6.3). It is conservative: two estimates whose intervals overlap
#' can still differ significantly (Schenker and Gentleman, 2001). When the two
#' curves come from the same bootstrap replicates (\code{compare2()} on a
#' \code{\link{P2boot}} fit, or \code{\link{compare_order}} with a
#' \code{\link{P1est}} built from that fit), the paired difference is also
#' tested directly with a percentile bootstrap interval (\code{diff_steps}).
#'
#' @return A list: first overlapping step \code{n} and time \code{s = n + 1}; the
#'   leading separated-step count \code{separated_steps}; per-step \code{overlap}
#'   and signed \code{separation} (= max lower - min upper, > 0 when separated);
#'   the two \code{groups}; and, when bootstrap replicate curves are available,
#'   \code{diff_steps} (leading steps whose percentile interval of the paired
#'   difference excludes 0; \code{NA} otherwise) and \code{diff} (per-step
#'   interval of the difference).
#' @references Schenker, N. and Gentleman, J. F. (2001). On judging the
#'   significance of differences by examining the overlap between confidence
#'   intervals. \emph{The American Statistician}, 55(3), 182-186.
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
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' overlap_step(cmp)
#' @export
overlap_step <- function(x) {
  stopifnot(inherits(x, "msm2pred"))
  if (!all(c("lower", "upper") %in% names(x)))
    stop("No evolution-interval bounds present. Re-run compare2() with bounds = TRUE.",
         call. = FALSE)
  groups <- unique(x$h)
  if (length(groups) != 2L)
    stop("overlap_step() compares exactly two preceding states; got ",
         length(groups), ".", call. = FALSE)

  a <- x[x$h == groups[1], ]; a <- a[order(a$n), ]
  b <- x[x$h == groups[2], ]; b <- b[order(b$n), ]
  ## Two intervals [lo1, up1] and [lo2, up2] fail to overlap exactly when one
  ## sits entirely above the other, i.e. max(lo1, lo2) > min(up1, up2); that
  ## same quantity, `separation`, doubles as a signed "gap" measure: positive
  ## means separated (and by how much), <= 0 means they overlap.
  separation <- pmax(a$lower, b$lower) - pmin(a$upper, b$upper)  # >0  <=>  separated
  overlap    <- separation <= 0
  ## The *first* step, scanning forward from n = 1, at which the two
  ## intervals stop being separated -- i.e. the point after which the
  ## preceding state h no longer provably changes the forecast.
  first      <- if (any(overlap)) which(overlap)[1L] else NA_integer_

  ## With bootstrap replicate curves for both groups (compare2() on a P2boot,
  ## or compare_order() with a P1est built from the same P2boot), the
  ## difference between the two curves is tested directly: percentile
  ## interval of the paired difference at each step. This is the proper test
  ## of "the two predictions differ"; non-overlap of two separate intervals is
  ## a more conservative criterion (Schenker and Gentleman, 2001).
  bc <- attr(x, "boot_curves")
  pd <- if (!is.null(bc) && all(groups %in% names(bc)))
    .boot_diff(bc[[groups[1]]], bc[[groups[2]]], attr(x, "conf.level") %||% 0.95)

  list(n               = if (is.na(first)) NA_integer_ else a$n[first],
       ## `s = n + 1` converts a step count back to a calendar time: step n
       ## corresponds to X_{n+1}, so it first becomes indistinguishable at
       ## time s = n + 1 (matching prep2()'s own s = time-of-destination-state
       ## convention).
       s               = if (is.na(first)) NA_integer_ else a$n[first] + 1L,
       ## How many *leading* steps stayed significantly separated before the
       ## (possible) first overlap -- 0 if they overlapped from step 1, or
       ## every step in the horizon if they never overlapped at all.
       separated_steps = if (is.na(first)) length(a$n) else first - 1L,
       overlap         = stats::setNames(overlap,    a$n),
       separation      = stats::setNames(separation, a$n),
       groups          = groups,
       diff_steps      = if (is.null(pd)) NA_integer_ else pd$steps,
       diff            = if (is.null(pd)) NULL else pd$diff)
}

## Console representation of a comparison: header (target/current state,
## which preceding states were compared) followed by the full table of
## per-step estimates (and bounds, if present).
#' @export
print.msm2pred <- function(x, ...) {
  cat(sprintf("<msm2pred>  %s n-step transitions to '%s' via current state '%s'\n",
              attr(x, "estimator"), attr(x, "l"), attr(x, "j")))
  cat(sprintf("  preceding states compared: %s\n", paste(unique(x$h), collapse = ", ")))
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}

## Adds a plain-language verdict on top of print.msm2pred()'s table: for a
## two-group comparison with bounds, runs overlap_step() and reports where
## (if anywhere) the preceding state stops mattering. Invisibly returns the
## overlap_step() result itself (so `os <- summary(cmp)` is a convenient way
## to grab it without a second call), falling back to `object` when there's
## no overlap analysis to report (fewer/more than 2 groups, or no bounds).
#' @export
summary.msm2pred <- function(object, ...) {
  has_ci <- all(c("lower", "upper") %in% names(object))
  groups <- unique(object$h)
  cat(sprintf("Trajectory comparison (%s%s)\n",
              attr(object, "estimator"),
              if (has_ci) sprintf(", %.0f%% %s intervals",
                                  100 * attr(object, "conf.level"),
                                  attr(object, "bands") %||% "evolution")
              else ", curves only"))
  cat(sprintf("  target '%s' via current state '%s'; preceding: %s\n",
              attr(object, "l"), attr(object, "j"), paste(groups, collapse = ", ")))

  os <- NULL
  if (has_ci && length(groups) == 2L) {
    os <- overlap_step(object)
    if (is.na(os$n))
      cat(sprintf("  intervals never overlap in %d steps: '%s' vs '%s' differ throughout.\n",
                  max(object$n), groups[1], groups[2]))
    else
      cat(sprintf("  intervals first overlap at step %d (time s = %d); significant for the first %d step(s).\n",
                  os$n, os$s, os$separated_steps))
    if (!is.na(os$diff_steps))
      cat(sprintf("  paired bootstrap test of the difference: significant for the first %d step(s).\n",
                  os$diff_steps))
  } else if (!has_ci) {
    cat("  (no bounds: run compare2(bounds = TRUE) for the overlap analysis)\n")
  }
  invisible(os %||% object)
}

#' Plot n-step evolution curves or evolution intervals
#'
#' @param x A "msm2pred" object.
#' @param type "interval" (shaded bands) or "curve" (lines only). Defaults to
#'   "interval" when bounds exist, else "curve".
#' @param col,lty,lwd,alpha Appearance of lines and interval shading.
#' @param add Overlay on the current plot (e.g. bootstrap over evolution bands). Default FALSE.
#' @param legend,dualaxis Draw the legend / the dual step-and-time axis.
#' @param mark.overlap For a two-group interval plot, draw a vertical line at the
#'   first overlap step. Default TRUE.
#' @param ylim,main,xlab,ylab Standard graphical parameters.
#' @param ... Passed to the initial plot().
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
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' plot(cmp, type = "interval")
#' plot(cmp, type = "curve")
#' @export
plot.msm2pred <- function(x, type = NULL,
                          col = NULL, lty = 1, lwd = 2, alpha = 0.2,
                          add = FALSE, legend = TRUE, dualaxis = TRUE,
                          mark.overlap = TRUE, ylim = NULL, main = NULL,
                          xlab = "", ylab = "probability", ...) {
  has_ci <- all(c("lower", "upper") %in% names(x))
  if (is.null(type)) type <- if (has_ci) "interval" else "curve"
  if (type == "interval" && !has_ci) {
    warning("No bounds available; drawing curves instead.", call. = FALSE)
    type <- "curve"
  }
  groups <- unique(x$h)
  if (is.null(col)) col <- c("red", "blue", "darkgreen", "purple",
                             "orange")[seq_along(groups)]
  if (is.null(ylim)) ylim <- c(0, max(if (type == "interval") x$upper else x$estimate))
  if (is.null(main)) main <- sprintf("To %s via %s", attr(x, "l"), attr(x, "j"))
  ns <- sort(unique(x$n))
  ## A translucent version of a line colour, used for the interval-band
  ## fill so overlapping bands from different groups stay visually
  ## distinguishable instead of one opaque polygon hiding another.
  fade <- function(cl) { v <- grDevices::col2rgb(cl) / 255
                         grDevices::rgb(v[1], v[2], v[3], alpha) }

  if (!add) {
    ## An empty (type = "n") plot just to establish the axis ranges; the
    ## actual data goes on in the per-group loop below. x-axis labelling is
    ## built by hand (xaxt = "n" here, then explicit axis() calls) so both a
    ## step count and its corresponding calendar time (dualaxis) can be shown.
    plot(range(ns), ylim, type = "n", xaxt = "n", las = 2,
         xlab = xlab, ylab = ylab, main = main, ...)
    axis(1, ns, ns)
    if (dualaxis) {
      ## Second copy of the x-axis, offset below the first, labelled with
      ## s = n + 1 (calendar time of X_{n+1}) instead of the step count n.
      mtext("steps (n)", 1, line = 1, at = min(ns) - 0.7)
      axis(1, ns, ns + 1L, line = 2.5)
      mtext("time (s)", 1, line = 2.5, at = min(ns) - 0.7)
    }
  }
  for (k in seq_along(groups)) {
    g <- x[x$h == groups[k], ]; g <- g[order(g$n), ]
    if (type == "interval") {
      ## Shaded band between the lower and upper bounds: polygon() needs its
      ## points in a single closed loop, hence tracing forward along `lower`
      ## and back along the reverse of `upper`. The band's own outline is
      ## then re-drawn as two thin solid lines so the boundary stays crisp
      ## even though the fill itself is translucent.
      polygon(c(g$n, rev(g$n)), c(g$lower, rev(g$upper)), col = fade(col[k]), border = NA)
      lines(g$n, g$lower, col = col[k], lwd = 0.5)
      lines(g$n, g$upper, col = col[k], lwd = 0.5)
    }
    lines(g$n, g$estimate, col = col[k], lwd = lwd, lty = lty)
    points(g$n, g$estimate, col = col[k], pch = if (lty == 1) 1 else 4)
  }
  ## Visually flag where the two compared trajectories stop being
  ## significantly different, reusing the very same calculation
  ## summary.msm2pred() reports in words.
  if (mark.overlap && type == "interval" && length(groups) == 2L && !add) {
    os <- overlap_step(x)
    if (!is.na(os$n)) abline(v = os$n, lty = 3, col = "grey40")
  }
  if (legend && !add)
    legend("topright", legend = groups, col = col, lty = 1, lwd = lwd, bty = "n")
  invisible(x)
}

## --- internal: null-coalescing helper ---------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
