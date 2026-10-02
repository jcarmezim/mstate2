#' n-step state-occupation probabilities in mstate's probtrans format
#'
#' The second-order analogue of \code{mstate::probtrans()}: for a history
#' \eqn{(X_0, X_1) = (h, j)}, returns the probability of occupying every state
#' at times \eqn{s = 1, \dots, n + 1} (with \code{predt = 1}), computed with the
#' extended Chapman-Kolmogorov relation (\code{\link{ckequations}}). Given a
#' first-order \code{\link{P1est}} fit, it returns the first-order predictions
#' from every starting state instead.
#'
#' The result has exactly the structure of an \code{mstate} "probtrans"
#' object: a list with one data frame per starting state (columns \code{time},
#' \code{pstate1}, ..., \code{pstateM}, \code{se1}, ..., \code{seM}), followed
#' by \code{trans}, \code{method}, \code{predt = 1} and
#' \code{direction = "forward"}. It can therefore be passed to
#' \code{\link{compare_order}} as \code{pt}, plotted with
#' \code{plot()} (types \code{"filled"}, \code{"stacked"}, \code{"single"},
#' \code{"separate"}, as in \code{mstate}), or used with tools written for
#' \code{probtrans} objects. For a second-order fit only the component of the
#' current state \code{j} is filled (the others are \code{NULL}), since the
#' prediction is specific to the history \eqn{(h, j)}.
#'
#' Standard errors (\code{se} columns) are the bootstrap standard deviations of
#' the curves when \code{x} is a \code{\link{P2boot}} object or a
#' \code{\link{P1est}} with bootstrap replicates, and \code{NA} otherwise.
#'
#' @param x A "P2est" (or "P2boot") object, or a "P1est" object.
#' @param h Previous state (time 0); required for a second-order fit, ignored
#'   for a first-order one.
#' @param j Current state (time 1). For a "P1est", \code{NULL} (default) fills
#'   the components of every non-absorbing state.
#' @param nsteps Number of steps. Default 9.
#' @return An object of class \code{c("probtrans2", "probtrans")}, with the
#'   extra element \code{h} (\code{NA} for first-order predictions).
#' @references
#' de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package
#' for the analysis of competing risks and multi-state models. \emph{Journal
#' of Statistical Software}, 38(7), 1-30.
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
#' pt2 <- probtrans2(fit, h = "A", j = "B", nsteps = 5)
#' pt2[[2]]
#' plot(pt2, type = "filled")
#' @export
probtrans2 <- function(x, h = NULL, j = NULL, nsteps = 9L) {
  if (nsteps < 1L) stop("`nsteps` must be >= 1.", call. = FALSE)
  states <- x$states; M <- length(states)
  out <- vector("list", M)

  if (inherits(x, "P1est")) {
    from <- if (is.null(j)) setdiff(seq_len(M), match(x$absorbing, states))
            else .resolve(j, states)
    if (anyNA(from)) stop("`j` not found in the state space.", call. = FALSE)
    T1 <- array(x$P, c(M, M, M), dimnames = list(states, states, states))
    Q  <- .pair_matrix(T1, M)
    for (jj in from) {
      prob <- .propagate(Q, 1L, jj, nsteps, M)
      se <- if (is.null(x$boot)) matrix(NA_real_, nsteps, M)
            else .boot_sd_first(x$boot, jj, nsteps)
      out[[jj]] <- .pt_frame(prob, se, jj, M)
    }
    method <- "first-order discrete RPE"
    hlab <- NA_character_
  } else if (inherits(x, "P2est")) {
    if (is.null(h) || is.null(j))
      stop("`h` and `j` are required for a second-order fit.", call. = FALSE)
    hi <- .resolve(h, states); ji <- .resolve(j, states)
    if (anyNA(c(hi, ji))) stop("`h` or `j` not found in the state space.", call. = FALSE)
    prob <- .propagate(.pair_matrix(x$P, M), hi, ji, nsteps, M)
    se <- if (inherits(x, "P2boot")) .boot_sd_second(x$boot, hi, ji, nsteps)
          else matrix(NA_real_, nsteps, M)
    out[[ji]] <- .pt_frame(prob, se, ji, M)
    method <- sprintf("second-order %s", x$estimator)
    hlab <- states[hi]
  } else {
    stop("`x` must be a 'P2est', 'P2boot' or 'P1est' object.", call. = FALSE)
  }

  out$trans     <- .trans_matrix(x)
  out$method    <- method
  out$predt     <- 1
  out$direction <- "forward"
  out$h         <- hlab
  class(out) <- c("probtrans2", "probtrans")
  out
}

## --- internal: one probtrans data frame (row 1 = time 1, the start) ----------
.pt_frame <- function(prob, se, ji, M) {
  p0 <- numeric(M); p0[ji] <- 1
  P  <- rbind(p0, prob); S <- rbind(numeric(M), se)
  df <- data.frame(time = seq_len(nrow(P)))
  for (k in seq_len(M)) df[[paste0("pstate", k)]] <- P[, k]
  for (k in seq_len(M)) df[[paste0("se", k)]]     <- S[, k]
  rownames(df) <- NULL
  df
}

## --- internal: bootstrap SD of the n-step curves -----------------------------
.boot_sd_second <- function(boot, hi, ji, nsteps) {
  M <- dim(boot)[1L]; B <- dim(boot)[4L]
  arr <- vapply(seq_len(B), function(b)
    .propagate(.pair_matrix(boot[, , , b], M), hi, ji, nsteps, M), matrix(0, nsteps, M))
  apply(arr, c(1L, 2L), stats::sd)
}
.boot_sd_first <- function(boot, ji, nsteps) {
  M <- dim(boot)[1L]; B <- dim(boot)[3L]
  arr <- vapply(seq_len(B), function(b)
    .propagate(.pair_matrix(array(boot[, , b], c(M, M, M)), M), 1L, ji, nsteps, M),
    matrix(0, nsteps, M))
  apply(arr, c(1L, 2L), stats::sd)
}

## --- internal: mstate-style transition matrix of the observed moves ----------
## Same layout as mstate::transMat(): integers numbering the transitions row
## by row, NA elsewhere, dimnames list(from = states, to = states).
.trans_matrix <- function(x) {
  states <- x$states; M <- length(states)
  e <- x$estimate
  mv <- unique(data.frame(j = match(as.character(e$j), states),
                          l = match(as.character(e$l), states)))
  mv <- mv[mv$j != mv$l, , drop = FALSE]
  mv <- mv[order(mv$j, mv$l), , drop = FALSE]
  tm <- matrix(NA_integer_, M, M, dimnames = list(from = states, to = states))
  if (nrow(mv)) tm[cbind(mv$j, mv$l)] <- seq_len(nrow(mv))
  tm
}

#' Plot n-step state-occupation probabilities
#'
#' Plot method for \code{\link{probtrans2}} objects, following the conventions
#' of \code{mstate}'s \code{plot.probtrans()}: \code{"filled"} (default) and
#' \code{"stacked"} stack the state-occupation probabilities, \code{"single"}
#' draws one curve per state, \code{"separate"} one panel per state.
#'
#' @param x A "probtrans2" object.
#' @param from Starting state whose predictions are plotted (index or label).
#'   Default: the first filled component.
#' @param type \code{"filled"}, \code{"stacked"}, \code{"single"} or \code{"separate"}.
#' @param ord Plotting order of the states for \code{"filled"}/\code{"stacked"}.
#' @param cols Colours, one per state.
#' @param xlab,ylab,xlim,ylim,lwd,lty,cex Usual graphical parameters.
#' @param legend State labels (default: the state names).
#' @param legend.pos Legend position for \code{type = "single"}.
#' @param bty Legend box type.
#' @param ... Passed to the initial \code{plot()}.
#' @return \code{x}, invisibly.
#' @export
plot.probtrans2 <- function(x, from = NULL,
                            type = c("filled", "stacked", "single", "separate"),
                            ord = NULL, cols = NULL, xlab = "Time (s)",
                            ylab = "Probability", xlim = NULL, ylim = NULL,
                            lwd = 1, lty = 1, cex = 0.9, legend = NULL,
                            legend.pos = "topright", bty = "n", ...) {
  type <- match.arg(type)
  trans <- x$trans; S <- nrow(trans)
  filled <- which(!vapply(x[seq_len(S)], is.null, logical(1)))
  from <- if (is.null(from)) filled[1L]
          else if (is.numeric(from)) as.integer(from)
          else match(from, rownames(trans))
  if (is.na(from) || is.null(x[[from]]))
    stop("No predictions stored for that starting state.", call. = FALSE)
  pt  <- x[[from]]
  tt  <- pt$time
  pp  <- as.matrix(pt[, paste0("pstate", seq_len(S))])
  if (is.null(legend)) legend <- rownames(trans)
  if (is.null(ord)) ord <- seq_len(S)
  if (is.null(cols)) cols <- if (type %in% c("filled", "single"))
    rev(grDevices::hcl.colors(S, "RdYlGn")) else rep("black", S)
  cols <- rep(cols, length.out = S); lwd <- rep(lwd, length.out = S)
  lty  <- rep(lty, length.out = S)
  if (is.null(xlim)) xlim <- range(tt)
  main_txt <- if (!is.na(x$h)) sprintf("From (%s, %s)", x$h, legend[from])
              else sprintf("From %s (first order)", legend[from])

  if (type == "single") {
    if (is.null(ylim)) ylim <- c(0, max(pp))
    plot(tt, pp[, 1], type = "s", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
         col = cols[1], lwd = lwd[1], lty = lty[1], main = main_txt, ...)
    for (s in seq_len(S)[-1]) lines(tt, pp[, s], type = "s", col = cols[s], lwd = lwd[s], lty = lty[s])
    graphics::legend(legend.pos, legend = legend, col = cols, lwd = lwd, lty = lty, bty = bty)
  } else if (type == "separate") {
    for (s in seq_len(S)) {
      plot(tt, pp[, s], type = "s", xlim = xlim,
           ylim = if (is.null(ylim)) range(pp[, s]) else ylim,
           xlab = xlab, ylab = ylab, col = cols[s], lwd = lwd[s], main = legend[s], ...)
    }
  } else {
    if (is.null(ylim)) ylim <- c(0, 1)
    plot(tt, rep(0, length(tt)), type = "n", xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main_txt, xaxs = "i", yaxs = "i", ...)
    low <- rep(0, length(tt)); y0 <- 0
    nt <- length(tt); eps <- diff(xlim) / 50
    for (s in ord) {
      up <- low + pp[, s]
      if (type == "filled") {
        ## step-function polygon between the lower and upper stacked curves
        xs <- c(rep(tt, each = 2)[-1], rev(rep(tt, each = 2)[-1]))
        ys <- c(rep(up, each = 2)[-(2 * nt)], rev(rep(low, each = 2)[-(2 * nt)]))
        polygon(xs, ys, col = cols[s], border = NA)
      }
      lines(tt, up, type = "s", col = if (type == "filled") "grey30" else cols[s], lwd = lwd[s])
      dy <- pp[nt, s]
      if (dy > 0.02) graphics::text(xlim[2] - eps, y0 + dy / 2, legend[s], adj = 1, cex = cex)
      y0 <- y0 + dy; low <- up
    }
    graphics::box()
  }
  invisible(x)
}
