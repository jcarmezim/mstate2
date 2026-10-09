#' Extended Chapman-Kolmogorov n-step second-order transition probabilities
#'
#' Computes \eqn{P_{hj\ell}(1,n) = P(X_{n+1} = \ell \mid X_1 = j, X_0 = h)} for \eqn{n = 1, \dots, nsteps} via the extension of the Chapman-Kolmogorov relation. The second-order chain is lifted to a first-order chain on ordered pairs of states, \eqn{Q_{(a,b)\to(b,c)} = P_{abc}}, and the initial pair distribution is propagated; this is exact, works for any number of states and any horizon, and is far cheaper than expanding the path sum.
#'
#' Step \eqn{n} is the probability of being in \eqn{\ell} at time \eqn{s = n + 1}; \eqn{n = 1} is the one-step probability \eqn{P_{hj\ell}}. The evolution intervals are obtained by propagating the tensors of lower and upper confidence limits in the same way.
#'
#' @param x A "P2est" object or an \eqn{M \times M \times M} tensor in P[j, l, h] layout.
#' @param h,j Starting states: h at time 0, j at time 1 (label or index).
#' @param l Target state(s). If NULL (default), the full distribution over all states is returned for each step.
#' @param nsteps Number of steps n (default 9), giving \eqn{X_2, \ldots, X_{n+1}}.
#' @param bounds If x is a "P2est" object and a single l is given, also the CI tensors to return evolution-interval bounds. If x is a \code{\link{P2boot}} object, the bounds are instead bootstrap intervals of the n-step probability (every replicate tensor propagated), which have the nominal coverage the evolution intervals lack.
#' @return If l is a single state and bounds = FALSE: a numeric vector of length nsteps. If l is a vector: an \eqn{nsteps \times length(l)} matrix. If l is NULL: an \eqn{nsteps \times M} matrix (one column per state). If bounds = TRUE: a tibble with columns n, estimate, lower, upper.
#' @section Why the chain on pairs:
#' A second-order chain on \eqn{M} states is a first-order chain on the \eqn{M^2} pairs of consecutive states, a standard representation of higher-order Markov chains (e.g. Benson, Gleich and Lim, 2017). Multiplying the pair distribution by its transition matrix once per step gives the extended Chapman-Kolmogorov relation of Najera-Zuloaga, Besalu and Gomez Melis (2025) exactly, at a cost linear in the horizon, without enumerating the \eqn{M^{n-1}} paths.
#'
#' @references
#' Benson, A. R., Gleich, D. F. and Lim, L.-H. (2017). The spacey random walk:
#' a stochastic process for higher-order data. \emph{SIAM Review}, 59(2),
#' 321-345.
#' @examples
#' st <- c("A", "B", "C")                                 # C is absorbing
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
#' ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6)
#' ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6, bounds = TRUE)
#' @export

ckequations <- function(x, h, j, l = NULL, nsteps = 9L, bounds = FALSE) {

  # Accept a "P2est" fit (estimates and confidence-limit tensors) or a bare M x M x M tensor (estimates only, so no bounds can be given).

  if (inherits(x, "P2est")) {
    states <- x$states
    P <- x$P
  } else if (is.array(x) && length(dim(x)) == 3L) {
    P <- x
    states <- dimnames(P)[[1]]
    if (is.null(states)) states <- as.character(seq_len(dim(P)[1L]))
    bounds <- FALSE
  } else {
    stop("`x` must be a 'P2est' object or an M x M x M tensor.", call. = FALSE)
  }
  if (nsteps < 1L) stop("`nsteps` must be >= 1.", call. = FALSE)

  # Propagate once from the starting pair (h, j): an nsteps x M matrix with the distribution of X_{n+1} over all states at every step n.
  dist <- .ck_distribution(P, h, j, nsteps, states)

  # Without a target state, return that whole matrix.
  if (is.null(l)) {
    if (bounds)
      warning("`bounds` needs an explicit `l`; returning estimates only.", call. = FALSE)
    return(dist)
  }

  # Keep the column(s) of the target state(s) l.
  li <- .resolve(l, states)
  if (anyNA(li)) stop("Target state(s) `l` not found.", call. = FALSE)
  est <- dist[, li] # vector if one l, else matrix
  if (!bounds) return(est)

  # Interval bounds need a single target state.
  if (length(li) != 1L)
    stop("`bounds = TRUE` needs a single target state `l`.", call. = FALSE)

  # Bootstrap fit (P2boot): percentile intervals of the n-step probability, from the propagation of every replicate tensor.
  if (inherits(x, "P2boot")) {
    bb <- .boot_bands(x$boot, .resolve(h, states), .resolve(j, states), li,
                      nsteps, x$conf.level)
    return(tibble::tibble(n = seq_len(nsteps), estimate = est,
                          lower = bb$lower, upper = bb$upper))
  }

  # Evolution intervals: propagate the tensors of lower and upper one-step limits exactly as the estimates. Their rows do not sum to 1, so the result can leave [0, 1] (e.g. an upper limit that accumulates in an absorbing state); clip it.
  lo <- pmax(0, .ck_distribution(x$P.lower, h, j, nsteps, states)[, li])
  up <- pmin(1, .ck_distribution(x$P.upper, h, j, nsteps, states)[, li])
  tibble::tibble(n = seq_len(nsteps), estimate = est, lower = lo, upper = up)
}

# internal: distribution of X_{n+1} at every step -> Resolves the starting states, builds the pair-transition matrix and  propagates it. Shared by ckequations(), compare2() and the bootstrap bands.
.ck_distribution <- function(P, h, j, nsteps, states) {
  M <- dim(P)[1L]
  hi <- .resolve(h, states)
  ji <- .resolve(j, states)
  if (anyNA(c(hi, ji))) stop("Starting state(s) not found.", call. = FALSE)
  res <- .propagate(.pair_matrix(P, M), hi, ji, nsteps, M)
  dimnames(res) <- list(NULL, states)
  res
}

# internal: propagate the chain on pairs -> The extended Chapman-Kolmogorov relation written as a first-order  chain on pairs of states: 1. start with probability 1 on the pair (X_0, X_1) = (h, j); 2. at each step multiply by Q, which moves the distribution from the pairs (X_{m-1}, X_m) to the pairs (X_m, X_{m+1}); 3. add over the first state of the pair to get the distribution of X_{m+1}. Each step costs one (M^2 x M^2) matrix-vector product.
.propagate <- function(Q, hi, ji, nsteps, M) {
  v <- numeric(M * M)
  v[(hi - 1L) * M + ji] <- 1 # 1. point mass on pair (h, j)
  out <- matrix(0, nsteps, M)
  for (m in seq_len(nsteps)) {
    v <- as.numeric(v %*% Q) # 2. advance one step
    # 3. v is indexed by the pair code (a - 1) M + b; as an M x M matrix filled by row, rows are X_m and columns X_{m+1}, so the column sums are the distribution of X_{m+1}.
    out[m, ] <- colSums(matrix(v, M, M, byrow = TRUE))
  }
  out
}

#  internal: pair-transition matrix Q[(a,b),(b,c)] = P_{abc} -> From the pair (a, b) ("was in a, now in b") the chain can only move to a pair (b, c), with probability P_{abc} = P[b, c, a]. Each pair is coded as the integer (a - 1) M + b, so Q is an ordinary M^2 x M^2 matrix.
.pair_matrix <- function(P, M) {
  Q <- matrix(0, M * M, M * M)
  for (a in seq_len(M)) for (b in seq_len(M))
    Q[(a - 1L) * M + b, (b - 1L) * M + seq_len(M)] <- P[b, , a]
  Q
}

# internal: state label or index -> Every state argument (h, j, l) can be a label or an integer position. An unknown label gives NA, which the callers turn into a "not found" error.
.resolve <- function(s, states) {
  if (!is.numeric(s)) return(match(as.character(s), states))
  i <- as.integer(s)
  i[i < 1L | i > length(states)] <- NA_integer_ # an index outside the state space is "not found"
  i
}
