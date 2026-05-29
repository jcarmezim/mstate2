#' Compare n-step transitions across preceding states (evolution intervals)
#'
#' For a fixed current state \eqn{j} and target state \eqn{\ell}, computes the
#' n-step second-order transition probabilities \eqn{P_{hj\ell}(1,n)} for several
#' preceding states \eqn{h}, together with their evolution intervals. This is the
#' construction behind Figures 4 and 5 of the paper: it lets one assess whether,
#' and for how long, the state occupied at the \emph{preceding} time affects the
#' future trajectory (i.e. whether the first-order Markov assumption holds).
#'
#' Evolution intervals are obtained by propagating the lower and upper confidence
#' bounds of the 1-step estimates through the extended Chapman-Kolmogorov relation
#' (\code{\link{ckequations}}). They are not themselves confidence intervals, but
#' they encompass them, so \emph{non-overlapping} evolution intervals imply a
#' statistically significant difference between the trajectories.
#'
#' @param object A \code{"P2est"} object.
#' @param h Vector of preceding states to compare (labels or integer indices).
#' @param j The current state (shared across the compared paths).
#' @param l The target state.
#' @param nsteps Number of steps to compute (default 9).
#'
#' @return An object of class \code{"msm2pred"} (a data frame) with columns
#'   \code{h}, \code{n}, \code{estimate}, \code{lower}, \code{upper}, and
#'   attributes recording \code{j}, \code{l}, \code{estimator} and
#'   \code{conf.level}.
#' @export
compare2 <- function(object, h, j, l, nsteps = 9) {
  stopifnot(inherits(object, "P2est"))
  states <- object$states
  lab <- function(s) if (is.numeric(s)) states[s] else as.character(s)

  parts <- lapply(h, function(hh) {
    ck <- ckequations(object, hh, j, l, nsteps = nsteps, bounds = TRUE)
    data.frame(h = lab(hh), ck, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, parts)
  df <- df[, c("h", "n", "estimate", "lower", "upper")]
  rownames(df) <- NULL

  structure(df, class = c("msm2pred", "data.frame"),
            j = lab(j), l = lab(l),
            estimator = object$estimator, conf.level = object$conf.level)
}

#' @export
print.msm2pred <- function(x, ...) {
  cat(sprintf("<msm2pred>  %s n-step transitions to '%s' via current state '%s'\n",
              attr(x, "estimator"), attr(x, "l"), attr(x, "j")))
  cat(sprintf("  preceding states compared: %s\n",
              paste(unique(x$h), collapse = ", ")))
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}

#' First step at which two evolution intervals overlap
#'
#' For a two-group \code{"msm2pred"}, returns the smallest step \eqn{n} at which
#' the evolution intervals of the two preceding states overlap. Steps before this
#' have non-overlapping intervals, i.e. a statistically significant effect of the
#' preceding state.
#'
#' @param x A \code{"msm2pred"} object with exactly two preceding states.
#' @return A list with the first overlapping step \code{n}, the corresponding
#'   time index \code{s = n + 1}, and the number of leading separated steps.
#' @export
overlap_step <- function(x) {
  stopifnot(inherits(x, "msm2pred"))
  groups <- unique(x$h)
  if (length(groups) != 2L)
    stop("overlap_step() needs exactly two preceding states.", call. = FALSE)
  a <- x[x$h == groups[1], ]; b <- x[x$h == groups[2], ]
  ord <- order(a$n); a <- a[ord, ]; b <- b[order(b$n), ]
  ## intervals overlap iff max(lowers) <= min(uppers)
  overlap <- pmax(a$lower, b$lower) <= pmin(a$upper, b$upper)
  first <- which(overlap)
  first <- if (length(first)) min(first) else NA_integer_
  list(n = if (is.na(first)) NA_integer_ else a$n[first],
       s = if (is.na(first)) NA_integer_ else a$n[first] + 1L,
       separated_steps = if (is.na(first)) length(a$n) else first - 1L,
       groups = groups)
}

#' @export
summary.msm2pred <- function(object, ...) {
  cat(sprintf("Evolution-interval comparison (%s, %.0f%% bounds)\n",
              attr(object, "estimator"), 100 * attr(object, "conf.level")))
  cat(sprintf("  target '%s' via current state '%s'\n",
              attr(object, "l"), attr(object, "j")))
  groups <- unique(object$h)
  if (length(groups) == 2L) {
    os <- overlap_step(object)
    if (is.na(os$n)) {
      cat(sprintf("  evolution intervals never overlap within %d steps:\n",
                  max(object$n)))
      cat(sprintf("  preceding state ('%s' vs '%s') matters throughout.\n",
                  groups[1], groups[2]))
    } else {
      cat(sprintf("  intervals first overlap at step n = %d (time s = %d)\n",
                  os$n, os$s))
      cat(sprintf("  -> preceding state ('%s' vs '%s') has a significant effect\n",
                  groups[1], groups[2]))
      cat(sprintf("     for at least the first %d step(s).\n", os$separated_steps))
    }
  }
  invisible(object)
}

#' Plot n-step evolution curves or evolution intervals
#'
#' @param x A \code{"msm2pred"} object.
#' @param type \code{"interval"} (default) draws the estimate with a shaded
#'   evolution interval per preceding state (Figure 5 style); \code{"curve"}
#'   draws estimate lines only (Figure 4 style).
#' @param col Colours, one per preceding state. Defaults to red/blue/...
#' @param lty,lwd Line type and width for the estimate lines.
#' @param alpha Shading transparency for evolution intervals.
#' @param add If \code{TRUE}, overlay on the current plot (e.g. to compare two
#'   estimators). Default \code{FALSE}.
#' @param legend Draw a legend. Default \code{TRUE}.
#' @param dualaxis Draw both the step axis (n) and the time axis (s = n + 1),
#'   as in the paper. Default \code{TRUE}.
#' @param ylim,main,xlab,ylab Standard graphical parameters; sensible defaults.
#' @param ... Passed to the initial \code{plot}.
#' @export
plot.msm2pred <- function(x, type = c("interval", "curve"),
                          col = NULL, lty = 1, lwd = 2, alpha = 0.2,
                          add = FALSE, legend = TRUE, dualaxis = TRUE,
                          ylim = NULL, main = NULL, xlab = "",
                          ylab = "probability", ...) {
  type <- match.arg(type)
  groups <- unique(x$h)
  if (is.null(col)) col <- c("red", "blue", "darkgreen", "purple",
                             "orange")[seq_along(groups)]
  if (is.null(ylim))
    ylim <- c(0, max(if (type == "interval") x$upper else x$estimate))
  if (is.null(main))
    main <- sprintf("To %s via %s", attr(x, "l"), attr(x, "j"))
  ns <- sort(unique(x$n))

  rgb_a <- function(cl) {
    v <- grDevices::col2rgb(cl) / 255
    grDevices::rgb(v[1], v[2], v[3], alpha = alpha)
  }

  if (!add) {
    plot(range(ns), ylim, type = "n", xaxt = "n", las = 2,
         xlab = xlab, ylab = ylab, main = main, ...)
    axis(1, ns, ns)
    if (dualaxis) {
      mtext("steps (n)", side = 1, line = 1, at = min(ns) - 0.7)
      axis(1, ns, ns + 1L, line = 2.5)
      mtext("time (s)", side = 1, line = 2.5, at = min(ns) - 0.7)
    }
  }

  for (k in seq_along(groups)) {
    g <- x[x$h == groups[k], ]; g <- g[order(g$n), ]
    if (type == "interval") {
      polygon(c(g$n, rev(g$n)), c(g$lower, rev(g$upper)),
              col = rgb_a(col[k]), border = NA)
      lines(g$n, g$lower, col = col[k], lwd = 0.5)
      lines(g$n, g$upper, col = col[k], lwd = 0.5)
    }
    lines(g$n, g$estimate, col = col[k], lwd = lwd, lty = lty)
    points(g$n, g$estimate, col = col[k], pch = if (lty == 1) 1 else 4)
  }
  if (legend && !add)
    legend("topright", legend = groups, col = col, lty = 1, lwd = lwd, bty = "n")
  invisible(x)
}
