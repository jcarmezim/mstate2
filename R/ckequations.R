#' Extended Chapman-Kolmogorov n-step second-order transition probabilities
#'
#' Computes \eqn{P_{hj\ell}(1,n) = P(X_{n+1} = \ell \mid X_1 = j, X_0 = h)} for
#' \eqn{n = 1, \dots, nsteps} via the extension of the Chapman-Kolmogorov
#' relation (Eq. 6). The second-order chain is lifted to a first-order chain on
#' ordered pairs of states, \eqn{Q_{(a,b)\to(b,c)} = P_{abc}}, and the initial
#' pair distribution is propagated; this is exact, works for any number of
#' states and any horizon, and is far cheaper than expanding the path sum.
#'
#' @param x A "P2est" object or an \eqn{M \times M \times M} tensor in
#'   P[j, l, h] layout.
#' @param h,j Starting states: h at time 0, j at time 1 (label or index).
#' @param l Target state(s). If NULL (default), the full distribution over all
#'   states is returned for each step.
#' @param nsteps Number of steps n (default 9), giving
#'   \eqn{X_2, \ldots, X_\{n+1\}}.
#' @param bounds If x is a "P2est" object and a single l is given, also propagate
#'   the CI tensors to return evolution-interval bounds.
#' @return If l is a single state and bounds = FALSE: a numeric vector of length
#'   nsteps. If l is a vector: an \eqn{nsteps \times length(l)} matrix. If l is
#'   NULL: an \eqn{nsteps \times M} matrix (one column per state). If
#'   bounds = TRUE: a data.frame with columns n, estimate, lower, upper.
#' @export
ckequations <- function(x, h, j, l = NULL, nsteps = 9L, bounds = FALSE) {
  if (inherits(x, "P2est")) {
    states <- x$states; P <- x$P
  } else if (is.array(x) && length(dim(x)) == 3L) {
    P <- x
    states <- dimnames(P)[[1]]
    if (is.null(states)) states <- as.character(seq_len(dim(P)[1L]))
    bounds <- FALSE
  } else {
    stop("`x` must be a 'P2est' object or an M x M x M tensor.", call. = FALSE)
  }
  if (nsteps < 1L) stop("`nsteps` must be >= 1.", call. = FALSE)

  dist <- .ck_distribution(P, h, j, nsteps, states)   # nsteps x M, one propagation

  if (is.null(l)) {
    if (bounds)
      warning("`bounds` needs an explicit `l`; returning estimates only.", call. = FALSE)
    return(dist)
  }

  li  <- .resolve(l, states)
  if (anyNA(li)) stop("Target state(s) `l` not found.", call. = FALSE)
  est <- dist[, li]                                   # vector if one l, else matrix
  if (!bounds) return(est)

  lo <- .ck_distribution(x$P.lower, h, j, nsteps, states)[, li]
  up <- .ck_distribution(x$P.upper, h, j, nsteps, states)[, li]
  data.frame(n = seq_len(nsteps), estimate = est, lower = lo, upper = up)
}

## --- internal: full per-step destination distribution -----------------------
.ck_distribution <- function(P, h, j, nsteps, states) {
  M  <- dim(P)[1L]
  hi <- .resolve(h, states); ji <- .resolve(j, states)
  if (anyNA(c(hi, ji))) stop("Starting state(s) not found.", call. = FALSE)
  res <- .propagate(.pair_matrix(P, M), hi, ji, nsteps, M)
  dimnames(res) <- list(NULL, states)
  res
}

## --- internal: propagate a prebuilt pair-matrix Q (reused by compare2) -------
.propagate <- function(Q, hi, ji, nsteps, M) {
  v <- numeric(M * M); v[(hi - 1L) * M + ji] <- 1     # start in pair (h, j)
  out <- matrix(0, nsteps, M)
  for (m in seq_len(nsteps)) {
    v <- as.numeric(v %*% Q)                          # advance to pair (X_m, X_{m+1})
    out[m, ] <- colSums(matrix(v, M, M, byrow = TRUE))  # marginal P(X_{m+1} = .)
  }
  out
}

## --- internal: pair-transition matrix Q[(a,b),(b,c)] = P_{abc} ---------------
.pair_matrix <- function(P, M) {
  Q <- matrix(0, M * M, M * M)
  for (a in seq_len(M)) for (b in seq_len(M))
    Q[(a - 1L) * M + b, (b - 1L) * M + seq_len(M)] <- P[b, , a]
  Q
}

## --- internal: resolve a label or index (vector) to state indices -----------
.resolve <- function(s, states) {
  if (is.numeric(s)) as.integer(s) else match(as.character(s), states)
}
