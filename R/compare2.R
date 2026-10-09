#' Compare n-step transitions across preceding states (evolution intervals)
#'
#' For a fixed current state \eqn{j} and target \eqn{\ell}, computes the n-step transition probabilities for several preceding states \eqn{h}, with their evolution intervals. This is the construction of the RPE curves: it assesses whether, and for how long, the state at the previous time affects the future trajectory.
#'
#' @param object A "P2est" object.
#' @param h Vector of preceding states to compare (labels or indices).
#' @param j Current state shared across the compared paths.
#' @param l Target state.
#' @param nsteps Number of steps (default 9).
#' @param bounds If TRUE (default) also compute evolution-interval bounds; FALSE returns curves only and skips two thirds of the work. If \code{object} is a \code{\link{P2boot}} fit, the bounds are percentile bootstrap intervals instead of evolution intervals.
#' @return An object of class "msm2pred" (a data frame): columns h, n, estimate and, when bounds = TRUE, lower and upper; one block of \code{nsteps} rows per previous state, in the order given in \code{h}.
#' @examples
#' st <- c("A", "B", "C") # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6
#' tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3
#' tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ] <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st))
#' first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' fit <- P2est(prep2(panel))
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' summary(cmp)
#' @export
compare2 <- function(object, h, j, l, nsteps = 9L, bounds = TRUE) {

  # Check the arguments and resolve the current and target states.

  stopifnot(inherits(object, "P2est"))
  states <- object$states
  M <- length(states)
  lab <- function(s) if (is.numeric(s)) states[s] else as.character(s)
  ji <- .resolve(j, states)
  li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)

  # Build the pair-transition matrices once and reuse them for every h: the estimates and, for evolution intervals, the lower and upper one-step limits.
  Q <- .pair_matrix(object$P, M)
  boot <- inherits(object, "P2boot")
  if (bounds && !boot) {
    Ql <- .pair_matrix(object$P.lower, M)
    Qu <- .pair_matrix(object$P.upper, M)
  }

  # One curve per previous state h: the n-step probability of l from the pair (h, j) and, if requested, its interval
  one_curve <- function(hh) {
    hi <- .resolve(hh, states)
    if (is.na(hi)) stop("preceding state '", hh, "' not found.", call. = FALSE)
    out <- tibble::tibble(h = lab(hh), n = seq_len(nsteps),
                          estimate = .propagate(Q, hi, ji, nsteps, M)[, li])
    if (bounds && boot) {
      bb <- .boot_bands(object$boot, hi, ji, li, nsteps, object$conf.level)
      out <- dplyr::mutate(out, lower = bb$lower, upper = bb$upper)
    } else if (bounds) {
      out <- dplyr::mutate(out,
                           lower = pmax(0, .propagate(Ql, hi, ji, nsteps, M)[, li]),
                           upper = pmin(1, .propagate(Qu, hi, ji, nsteps, M)[, li]))
    }
    out
  }
  curves <- purrr::map(h, one_curve) |>
    purrr::list_rbind()

  # Return a data frame of class "msm2pred"
  structure(as.data.frame(curves), class = c("msm2pred", "data.frame"),
            j = lab(j), l = lab(l), bounds = bounds,
            bands = if (!bounds) "none" else if (boot) "bootstrap" else "evolution",
            estimator = object$estimator, conf.level = object$conf.level)
}

#' First step at which two evolution intervals overlap
#'
#' Scans the two curves of a \code{\link{compare2}} comparison forward from \eqn{n = 1} and returns the first step at which their intervals overlap: before it, the state at the previous time changes the prediction.
#'
#' @param x A two-group "msm2pred" with interval bounds.
#' @return A list: first overlapping step \code{n} and time \code{s = n + 1}; the leading separated-step count \code{separated_steps}; per-step \code{overlap} and signed \code{separation} (= max lower - min upper, > 0 when separated); and the two \code{groups}.
#' @examples
#' st <- c("A", "B", "C") # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6
#' tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3
#' tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ] <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st))
#' first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' fit <- P2est(prep2(panel))
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' overlap_step(cmp)
#' @export
overlap_step <- function(x) {

  # Check that there are interval bounds and exactly two curves.

  stopifnot(inherits(x, "msm2pred"))
  if (!all(c("lower", "upper") %in% names(x)))
    stop("No evolution-interval bounds present. Re-run compare2() with bounds = TRUE.",
         call. = FALSE)
  groups <- unique(x$h)
  if (length(groups) != 2L)
    stop("overlap_step() compares exactly two preceding states; got ",
         length(groups), ".", call. = FALSE)

  # The two curves, each ordered by step
  curve <- function(g) tibble::as_tibble(x) |> dplyr::filter(h == g) |> dplyr::arrange(n)
  a <- curve(groups[1])
  b <- curve(groups[2])

  # Two intervals are separated when one lies entirely above the other, i.e. when max(lower) > min(upper). That difference is also a signed gap: > 0 separated (by that much), <= 0 overlapping.
  separation <- pmax(a$lower, b$lower) - pmin(a$upper, b$upper)
  overlap <- separation <= 0

  # First step, scanning forward from n = 1, at which they overlap: from then on the previous state no longer changes the prediction significantly. NA if they never overlap within the horizon.
  first <- if (any(overlap)) which(overlap)[1L] else NA_integer_

  # Return the step n, its time s = n + 1 (step n predicts X_{n+1}), the number of leading separated steps, and the per-step detail.
  list(n = if (is.na(first)) NA_integer_ else a$n[first],
       s = if (is.na(first)) NA_integer_ else a$n[first] + 1L,
       separated_steps = if (is.na(first)) length(a$n) else first - 1L,
       overlap = stats::setNames(overlap,    a$n),
       separation = stats::setNames(separation, a$n),
       groups = groups)
}

# Printed form of a comparison
#' @export
print.msm2pred <- function(x, ...) {
  cat(sprintf("<msm2pred>  %s n-step transitions to '%s' via current state '%s'\n", attr(x, "estimator"), attr(x, "l"), attr(x, "j")))
  cat(sprintf("  preceding states compared: %s\n", paste(unique(x$h), collapse = ", ")))
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}

# Verdict
#' @export
summary.msm2pred <- function(object, ...) {
  
  # Header: estimator, type and level of the intervals, compared paths.
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

  # Overlap analysis, when there are bounds and exactly two curves.
  os <- NULL
  if (has_ci && length(groups) == 2L) {
    os <- overlap_step(object)
    if (is.na(os$n))
      cat(sprintf("  intervals never overlap in %d steps: '%s' vs '%s' differ throughout.\n",
                  max(object$n), groups[1], groups[2]))
    else
      cat(sprintf("  intervals first overlap at step %d (time s = %d); significant for the first %d step(s).\n",
                  os$n, os$s, os$separated_steps))
  } else if (!has_ci) {
    cat("  (no bounds: run compare2(bounds = TRUE) for the overlap analysis)\n")
  }
  invisible(os %||% object)
}

#' Plot n-step evolution curves or evolution intervals
#'
#' @param x A "msm2pred" object.
#' @param type "interval" (shaded bands) or "curve" (lines only). Defaults to "interval" when bounds exist, else "curve".
#' @param col,lty,lwd,alpha Appearance of lines and interval shading. The colours are used in the order of \code{h} in \code{\link{compare2}}
#' @param add Overlay on the current plot. Default FALSE.
#' @param legend,dualaxis Draw the legend / the dual step-and-time axis.
#' @param mark.overlap For a two-group interval plot, draw a vertical line at the first overlap step. Default TRUE.
#' @param ylim,main,xlab,ylab Standard graphical parameters.
#' @param ... Passed to the initial plot().
#' @examples
#' st <- c("A", "B", "C") # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6
#' tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3
#' tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ] <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st))
#' first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' fit <- P2est(prep2(panel))
#' cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
#' plot(cmp, type = "interval")
#' plot(cmp, type = "curve")
#' @export
plot.msm2pred <- function(x, type = NULL, col = NULL, lty = 1, lwd = 2, alpha = 0.2, add = FALSE, legend = TRUE, dualaxis = TRUE, mark.overlap = TRUE, ylim = NULL, main = NULL, xlab = "", ylab = "probability", ...) {

  # Plot type: shaded intervals when there are bounds, lines otherwise.
  has_ci <- all(c("lower", "upper") %in% names(x))
  if (is.null(type)) type <- if (has_ci) "interval" else "curve"
  if (type == "interval" && !has_ci) {
    warning("No bounds available; drawing curves instead.", call. = FALSE)
    type <- "curve"
  }

  # Defaults: one colour per previous state, a y range that contains every curve (or interval), a title naming the transition.
  groups <- unique(x$h)
  if (is.null(col)) col <- c("red", "blue", "darkgreen", "purple",
                             "orange")[seq_along(groups)]
  if (is.null(ylim)) ylim <- c(0, max(if (type == "interval") x$upper else x$estimate))
  if (is.null(main)) main <- sprintf("To %s via %s", attr(x, "l"), attr(x, "j"))
  ns <- sort(unique(x$n))

  # Translucent version of a colour, for the interval fill, so that  overlapping intervals stay visible.
  fade <- function(cl) { 
    v <- grDevices::col2rgb(cl) / 255 
    grDevices::rgb(v[1], v[2], v[3], alpha) 
  }

  # Empty frame with the axes. With dualaxis = TRUE a second axis gives the time s = n + 1 under the step n
  if (!add) {
    plot(range(ns), ylim, type = "n", xaxt = "n", las = 2,
         xlab = xlab, ylab = ylab, main = main, ...)
    axis(1, ns, ns)
    if (dualaxis) {
      mtext("steps (n)", 1, line = 1, at = min(ns) - 0.7)
      axis(1, ns, ns + 1L, line = 2.5)
      mtext("time (s)", 1, line = 2.5, at = min(ns) - 0.7)
    }
  }

  # One curve per previous state
  for (k in seq_along(groups)) {
    g <- tibble::as_tibble(x) |>
      dplyr::filter(h == groups[k]) |>
      dplyr::arrange(n)
    if (type == "interval") {
      polygon(c(g$n, rev(g$n)), c(g$lower, rev(g$upper)), col = fade(col[k]), border = NA)
      lines(g$n, g$lower, col = col[k], lwd = 0.5)
      lines(g$n, g$upper, col = col[k], lwd = 0.5)
    }
    lines(g$n, g$estimate, col = col[k], lwd = lwd, lty = lty)
    points(g$n, g$estimate, col = col[k], pch = if (lty == 1) 1 else 4)
  }

  # Dotted vertical line at the first step at which two intervals overlap and the legend
  if (mark.overlap && type == "interval" && length(groups) == 2L && !add) {
    os <- overlap_step(x)
    if (!is.na(os$n)) abline(v = os$n, lty = 3, col = "grey40")
  }
  if (legend && !add)
    legend("topright", legend = groups, col = col, lty = 1, lwd = lwd, bty = "n")
  invisible(x)
}

# internal: `a` unless it is NULL, then `b`
`%||%` <- function(a, b) if (is.null(a)) b else a
