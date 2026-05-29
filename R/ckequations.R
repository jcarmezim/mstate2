#' Extended Chapman-Kolmogorov n-step second-order transition probabilities
#'
#' Computes the n-step second-order homogeneous transition probabilities
#' \eqn{P_{hj\ell}(1,n) = P(X_{n+1} = \ell \mid X_1 = j, X_0 = h)} for
#' \eqn{n = 1, \dots, nsteps}, by applying the extension of the
#' Chapman-Kolmogorov relation (Eq. 6 in the paper).
#'
#' Rather than expanding the path sum of Eq. (6) directly (which is
#' \eqn{O(M^{n})}), the second-order chain is lifted to an equivalent
#' first-order Markov chain on ordered pairs of states: the pair
#' \eqn{(X_{t-1}, X_t)} evolves with transition
#' \eqn{Q_{(a,b)\to(b,c)} = P_{abc}}. The n-step probabilities then follow from
#' propagating the initial pair distribution, which is exact, valid for any
#' number of states and any horizon, and far cheaper.
#'
#' @param x Either a \code{"P2est"} object or a numeric \eqn{M\times M\times M}
#'   tensor in the \code{P[j, l, h]} \eqn{= P_{hj\ell}} layout.
#' @param h,j,l Starting states \eqn{h} (at time 0), \eqn{j} (at time 1), and the
#'   target state \eqn{\ell}. May be given by label or by integer index.
#' @param nsteps Number of steps \eqn{n} to compute (default 9, giving
#'   \eqn{X_2, \dots, X_{n+1}}).
#' @param bounds If \code{x} is a \code{"P2est"} object, also propagate the lower
#'   and upper CI tensors to produce evolution-interval bounds. Ignored for a
#'   raw tensor.
#'
#' @return If \code{bounds = FALSE} (or \code{x} is a tensor), a numeric vector of
#'   length \code{nsteps} giving \eqn{P_{hj\ell}(1, n)} for \eqn{n = 1, \dots, nsteps}.
#'   If \code{bounds = TRUE}, a data frame with columns \code{n}, \code{estimate},
#'   \code{lower}, \code{upper}.
#' @export
ckequations <- function(x, h, j, l, nsteps = 9, bounds = FALSE) {
  if (inherits(x, "P2est")) {
    states <- x$states
    P <- x$P
  } else if (is.array(x) && length(dim(x)) == 3L) {
    states <- dimnames(x)[[1]]
    if (is.null(states)) states <- as.character(seq_len(dim(x)[1]))
    P <- x
    bounds <- FALSE
  } else {
    stop("`x` must be a 'P2est' object or an M x M x M tensor.", call. = FALSE)
  }

  est <- .nstep(P, h, j, l, nsteps, states)
  if (!bounds) return(est)

  lo <- .nstep(x$P.lower, h, j, l, nsteps, states)
  up <- .nstep(x$P.upper, h, j, l, nsteps, states)
  data.frame(n = seq_len(nsteps), estimate = est, lower = lo, upper = up)
}

## --- internal: pair-lift propagation on one tensor ---------------------------
.nstep <- function(P, h, j, l, nsteps, states) {
  M <- dim(P)[1]
  resolve <- function(s) if (is.numeric(s)) as.integer(s) else match(as.character(s), states)
  hi <- resolve(h); ji <- resolve(j); li <- resolve(l)
  if (anyNA(c(hi, ji, li)))
    stop("State(s) not found in the tensor.", call. = FALSE)

  Q <- .pairQ(P, M)                       # (M^2) x (M^2) pair-transition matrix
  v <- numeric(M * M)
  v[(hi - 1L) * M + ji] <- 1              # start in pair (h, j)
  target <- (seq_len(M) - 1L) * M + li    # all pairs (a, l)

  out <- numeric(nsteps)
  for (m in seq_len(nsteps)) {
    v <- as.numeric(v %*% Q)              # advance one step: pair (X_{m}, X_{m+1})
    out[m] <- sum(v[target])              # P(X_{m+1} = l)
  }
  out
}

## --- internal: build the pair-transition matrix Q ----------------------------
## Q[(a,b), (b,c)] = P_{abc} = P[b, c, a]; all other entries 0.
.pairQ <- function(P, M) {
  Q <- matrix(0, M * M, M * M)
  for (a in seq_len(M)) {
    for (b in seq_len(M)) {
      from <- (a - 1L) * M + b
      to   <- (b - 1L) * M + seq_len(M)   # pairs (b, c), c = 1..M
      Q[from, to] <- P[b, , a]            # (P_{a b c})_c
    }
  }
  Q
}
