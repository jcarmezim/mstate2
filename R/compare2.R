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
#' @return An object of class "msm2pred" (a data frame): columns h, n, estimate
#'   and, when bounds = TRUE, lower and upper.
#' @export
compare2 <- function(object, h, j, l, nsteps = 9L, bounds = TRUE) {
  stopifnot(inherits(object, "P2est"))
  states <- object$states; M <- length(states)
  lab <- function(s) if (is.numeric(s)) states[s] else as.character(s)
  ji <- .resolve(j, states); li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)

  ## build pair-transition matrices ONCE, reuse across every h
  Q  <- .pair_matrix(object$P, M)
  if (bounds) {
    Ql <- .pair_matrix(object$P.lower, M)
    Qu <- .pair_matrix(object$P.upper, M)
  }

  parts <- lapply(h, function(hh) {
    hi <- .resolve(hh, states)
    if (is.na(hi)) stop("preceding state '", hh, "' not found.", call. = FALSE)
    out <- data.frame(h = lab(hh), n = seq_len(nsteps),
                      estimate = .propagate(Q, hi, ji, nsteps, M)[, li])
    if (bounds) {
      out$lower <- .propagate(Ql, hi, ji, nsteps, M)[, li]
      out$upper <- .propagate(Qu, hi, ji, nsteps, M)[, li]
    }
    out
  })

  df <- do.call(rbind, parts); rownames(df) <- NULL
  structure(df, class = c("msm2pred", "data.frame"),
            j = lab(j), l = lab(l), bounds = bounds,
            estimator = object$estimator, conf.level = object$conf.level)
}

#' First step at which two evolution intervals overlap
#'
#' @param x A two-group "msm2pred" with evolution-interval bounds.
#' @return A list: first overlapping step \code{n} and time \code{s = n + 1}; the
#'   leading separated-step count \code{separated_steps}; per-step \code{overlap}
#'   and signed \code{separation} (= max lower - min upper, > 0 when separated);
#'   and the two \code{groups}.
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
  separation <- pmax(a$lower, b$lower) - pmin(a$upper, b$upper)  # >0  <=>  separated
  overlap    <- separation <= 0
  first      <- if (any(overlap)) which(overlap)[1L] else NA_integer_

  list(n               = if (is.na(first)) NA_integer_ else a$n[first],
       s               = if (is.na(first)) NA_integer_ else a$n[first] + 1L,
       separated_steps = if (is.na(first)) length(a$n) else first - 1L,
       overlap         = stats::setNames(overlap,    a$n),
       separation      = stats::setNames(separation, a$n),
       groups          = groups)
}

#' @export
print.msm2pred <- function(x, ...) {
  cat(sprintf("<msm2pred>  %s n-step transitions to '%s' via current state '%s'\n",
              attr(x, "estimator"), attr(x, "l"), attr(x, "j")))
  cat(sprintf("  preceding states compared: %s\n", paste(unique(x$h), collapse = ", ")))
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}

#' @export
summary.msm2pred <- function(object, ...) {
  has_ci <- all(c("lower", "upper") %in% names(object))
  groups <- unique(object$h)
  cat(sprintf("Trajectory comparison (%s%s)\n",
              attr(object, "estimator"),
              if (has_ci) sprintf(", %.0f%% evolution intervals",
                                  100 * attr(object, "conf.level")) else ", curves only"))
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
#' @param add Overlay on the current plot (e.g. RPE over CPE). Default FALSE.
#' @param legend,dualaxis Draw the legend / the dual step-and-time axis.
#' @param mark.overlap For a two-group interval plot, draw a vertical line at the
#'   first overlap step. Default TRUE.
#' @param ylim,main,xlab,ylab Standard graphical parameters.
#' @param ... Passed to the initial plot().
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
  fade <- function(cl) { v <- grDevices::col2rgb(cl) / 255
                         grDevices::rgb(v[1], v[2], v[3], alpha) }

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
  for (k in seq_along(groups)) {
    g <- x[x$h == groups[k], ]; g <- g[order(g$n), ]
    if (type == "interval") {
      polygon(c(g$n, rev(g$n)), c(g$lower, rev(g$upper)), col = fade(col[k]), border = NA)
      lines(g$n, g$lower, col = col[k], lwd = 0.5)
      lines(g$n, g$upper, col = col[k], lwd = 0.5)
    }
    lines(g$n, g$estimate, col = col[k], lwd = lwd, lty = lty)
    points(g$n, g$estimate, col = col[k], pch = if (lty == 1) 1 else 4)
  }
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
