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
#'   the CI tensors to return evolution-interval bounds (clipped to [0, 1]).
#'   If x is a \code{\link{P2boot}} object, the bounds are instead percentile
#'   bootstrap intervals of the n-step probability (every replicate tensor is
#'   propagated), which have the nominal coverage the evolution intervals lack.
#' @return If l is a single state and bounds = FALSE: a numeric vector of length
#'   nsteps. If l is a vector: an \eqn{nsteps \times length(l)} matrix. If l is
#'   NULL: an \eqn{nsteps \times M} matrix (one column per state). If
#'   bounds = TRUE: a data.frame with columns n, estimate, lower, upper.
#' @section Why the chain on pairs:
#' A second-order chain on \eqn{M} states is a first-order chain on the
#' \eqn{M^2} pairs of consecutive states, a standard representation of
#' higher-order Markov chains (e.g. Benson, Gleich and Lim, 2017). Multiplying
#' the pair distribution by its transition matrix once per step gives the
#' extended Chapman-Kolmogorov relation of Najera-Zuloaga, Besalu and Gomez
#' Melis (2025) exactly, at a cost linear in the horizon, without enumerating
#' the \eqn{M^{n-1}} paths.
#'
#' @references
#' Benson, A. R., Gleich, D. F. and Lim, L.-H. (2017). The spacey random walk:
#' a stochastic process for higher-order data. \emph{SIAM Review}, 59(2),
#' 321-345.
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
#' ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6)
#' ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6, bounds = TRUE)
#' @export
ckequations <- function(x, h, j, l = NULL, nsteps = 9L, bounds = FALSE) {
  ## Accept either a full "P2est" object (point estimate + CI tensors, so
  ## `bounds = TRUE` can also propagate the confidence limits) or a bare
  ## M x M x M tensor (point estimates only -- e.g. a hand-built or
  ## simulation "ground truth" tensor with no notion of a CI, hence
  ## `bounds` is forced off in that branch regardless of what was asked for).
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

  ## One propagation gives every step's *full* destination distribution at
  ## once (an nsteps x M matrix); slicing out just the column(s) asked for
  ## below is cheap, so there's no need to re-propagate per requested `l`.
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

  ## Evolution-interval bounds: propagate the *lower* and *upper* CI tensors
  ## through the exact same Chapman-Kolmogorov machinery used for the point
  ## estimate. These are conservative bands (they compound each step's own
  ## uncertainty forward), not exact confidence intervals for the n-step
  ## probability -- but they're what markov_test()'s informal cousin,
  ## overlap_step(), needs to compare two trajectories.
  ## The rows of P.lower / P.upper do not sum to 1, so their propagation is
  ## not confined to [0, 1] (e.g. an upper limit can accumulate past 1 when it
  ## feeds an absorbing state); clip to the probability scale.
  if (length(li) != 1L)
    stop("`bounds = TRUE` needs a single target state `l`.", call. = FALSE)
  if (inherits(x, "P2boot")) {
    ## percentile bootstrap intervals: propagate every replicate tensor
    bb <- .boot_bands(x$boot, .resolve(h, states), .resolve(j, states), li,
                      nsteps, x$conf.level)
    return(data.frame(n = seq_len(nsteps), estimate = est,
                      lower = bb$lower, upper = bb$upper))
  }
  lo <- pmax(0, .ck_distribution(x$P.lower, h, j, nsteps, states)[, li])
  up <- pmin(1, .ck_distribution(x$P.upper, h, j, nsteps, states)[, li])
  data.frame(n = seq_len(nsteps), estimate = est, lower = lo, upper = up)
}

## --- internal: full per-step destination distribution -----------------------
## Resolves the starting states, lifts the tensor to a pair-transition
## matrix, and propagates it -- the shared machinery behind ckequations(),
## compare2() and compare_order()'s second-order curves.
.ck_distribution <- function(P, h, j, nsteps, states) {
  M  <- dim(P)[1L]
  hi <- .resolve(h, states); ji <- .resolve(j, states)
  if (anyNA(c(hi, ji))) stop("Starting state(s) not found.", call. = FALSE)
  res <- .propagate(.pair_matrix(P, M), hi, ji, nsteps, M)
  dimnames(res) <- list(NULL, states)
  res
}

## --- internal: propagate a prebuilt pair-matrix Q (reused by compare2) -------
## This is the extended Chapman-Kolmogorov relation (Eq. 6), implemented as
## ordinary first-order Markov-chain propagation on the *pair* state space:
## starting from a point mass on pair (h, j) = (X_0, X_1), each multiplication
## by Q advances one step, so after m iterations `v` holds the distribution
## over pairs (X_m, X_{m+1}). This is exact (no approximation) and, unlike
## summing over every length-m path h -> j -> ... -> l by hand, its cost
## grows only linearly in nsteps (one M^2 x M^2 matrix-vector product per
## step) rather than exponentially in the number of states.
.propagate <- function(Q, hi, ji, nsteps, M) {
  v <- numeric(M * M); v[(hi - 1L) * M + ji] <- 1     # start in pair (h, j)
  out <- matrix(0, nsteps, M)
  for (m in seq_len(nsteps)) {
    v <- as.numeric(v %*% Q)                          # advance to pair (X_m, X_{m+1})
    ## `v` is still indexed by the same linear pair-code (a-1)*M+b used to
    ## build Q; reshaping it byrow into an M x M matrix puts the *earlier*
    ## state of the pair (X_m) on the rows and the *later* one (X_{m+1}) on
    ## the columns, so summing over rows (colSums) marginalizes out X_m and
    ## leaves exactly the distribution of X_{m+1} -- the quantity this
    ## function returns for step m.
    out[m, ] <- colSums(matrix(v, M, M, byrow = TRUE))  # marginal P(X_{m+1} = .)
  }
  out
}

## --- internal: pair-transition matrix Q[(a,b),(b,c)] = P_{abc} ---------------
## Lifts the second-order tensor P[j, l, h] = P_{hjl} to an ordinary
## first-order transition matrix Q on ordered pairs of states: from pair
## (a, b) (meaning "was in a, now in b"), the only pairs reachable in one
## more step are (b, c) for each c, with probability P_{abc} (i.e. what P
## says about going b -> c given the preceding state was a). Pairs are
## encoded as a single integer (a-1)*M + b so Q can be an ordinary
## M^2 x M^2 matrix rather than a 4-dimensional array.
.pair_matrix <- function(P, M) {
  Q <- matrix(0, M * M, M * M)
  for (a in seq_len(M)) for (b in seq_len(M))
    Q[(a - 1L) * M + b, (b - 1L) * M + seq_len(M)] <- P[b, , a]
  Q
}

## --- internal: resolve a label or index (vector) to state indices -----------
## Every state-space argument throughout the package (h, j, l, ...) can be
## given either as a state label (matched against `states`) or as a bare
## integer position -- this is the single place that distinction is resolved,
## used consistently by ckequations(), compare2(), compare_order() and
## P2reg(). An unmatched label (or an index that isn't itself validated by
## the caller) comes back as NA, which callers check with anyNA() to raise a
## clear "not found" error rather than an obscure indexing failure later.
.resolve <- function(s, states) {
  if (is.numeric(s)) as.integer(s) else match(as.character(s), states)
}
