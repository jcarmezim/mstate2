#' Subject-level bootstrap of the second-order estimates
#'
#' Resamples subjects with replacement and recomputes the RPE of every 1-step second-order transition probability in each replicate. The resulting object behaves like a \code{\link{P2est}} fit, but carries the bootstrap tensors, so that \code{\link{ckequations}} and \code{\link{compare2}} report \strong{percentile bootstrap intervals} for the \eqn{n}-step probabilities instead of the evolution intervals obtained by propagating the one-step confidence limits.
#'
#' Resampling whole subjects keeps the dependence between the several instants a subject contributes to the same or to different \eqn{(h, j)} risk sets. The implementation never re-reads the data: the per-subject counts of every observed \eqn{(h, j, \ell)} triple are tabulated once, and each replicate is a weighted sum of those counts, so thousands of replicates are cheap even for large cohorts.
#'
#' In a replicate where an observed \eqn{(h, j)} pair happens to have no subject at risk, its row is taken from the point estimate, so that the replicate tensor stays a proper set of transition probabilities.
#'
#' @param object An "msm2data" object from \code{\link{prep2}}.
#' @param B Number of bootstrap replicates. Default 200.
#' @param conf.level Confidence level of the intervals. Default 0.95.
#' @param seed Optional integer seed, for reproducible replicates.
#' @return An object of class \code{c("P2boot", "P2est")}: all the components of \code{P2est(object, conf.level)} plus \code{boot}, an \eqn{M \times M \times M \times B} array of replicate tensors (layout \code{[j, l, h, b]}), \code{boot1}, the \eqn{M \times M \times B} array of first-order (history-pooled) matrices of the same resamples (used by \code{\link{P1est}} and \code{\link{compare_order}}), and \code{B}. The \code{estimate} table gains a column \code{se.boot} (bootstrap standard deviation).
#' @seealso \code{\link{ckequations}}, \code{\link{compare2}}
#' @section Why a bootstrap of subjects:
#' Patients are the independent units; the days of one patient are not. Resampling whole patients keeps that within-patient dependence, which is the standard non-parametric bootstrap for clustered data (Davison and Hinkley, 1997; Field and Welsh, 2007). The n-step predictions are non-linear (polynomial) functions of all the estimated probabilities, so a delta-method variance would be cumbersome; propagating each replicate and taking percentiles (Efron and Tibshirani, 1993) gives their intervals directly.
#'
#' @references
#' Davison, A. C. and Hinkley, D. V. (1997). \emph{Bootstrap Methods and their
#' Application}. Cambridge University Press.
#'
#' Efron, B. and Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman and Hall.
#'
#' Field, C. A. and Welsh, A. H. (2007). Bootstrapping clustered data.
#' \emph{Journal of the Royal Statistical Society: Series B}, 69(3), 369-390.
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
#' bt <- P2boot(prep2(panel), B = 100, seed = 1)
#' bt
#' ckequations(bt, h = "A", j = "B", l = "C", nsteps = 4, bounds = TRUE)
#' @export
P2boot <- function(object, B = 200, conf.level = 0.95, seed = NULL) {

  # Check the arguments and set the seed, if given.
  stopifnot(inherits(object, "msm2data"))
  .check_conf_level(conf.level)
  if (!is.numeric(B) || length(B) != 1L || B < 2)
    stop("`B` must be a single number >= 2.", call. = FALSE)
  B <- as.integer(B)
  if (!is.null(seed)) set.seed(seed)

  # Point estimates, and the per-subject counts of every observed (h, j, l) triple (one row per subject, one column per triple).
  fit <- P2est(object, conf.level = conf.level)
  ct <- .id_counts(object)
  n <- nrow(ct$mat)

  # First-order (history-pooled) point estimate, the fallback of a replicate in which a state has nobody at risk.
  P1_point <- .first_from_counts(colSums(ct$mat), ct, object, NULL)

  # Replicates. Each one draws n subjects with replacement; w is the number of times each subject is drawn, so the replicate's count of every triple is the weighted sum w %*% mat, without rebuilding the data. The counts are then turned into a tensor of probabilities and, from the same resample, into the first-order matrix used by P1est() and compare_order().
  one_replicate <- function(b) {
    w <- tabulate(sample.int(n, n, replace = TRUE), n)
    tot <- as.vector(w %*% ct$mat)
    list(second = .tensor_from_counts(tot, ct, object, fit$P),
         first = .first_from_counts(tot, ct, object, P1_point))
  }
  reps <- purrr::map(seq_len(B), one_replicate)
  boot <- simplify2array(purrr::map(reps, "second"))  # M x M x M x B array
  boot1 <- simplify2array(purrr::map(reps, "first"))  # M x M x B array

  # Bootstrap standard deviation of every estimated probability, added to the estimates table next to the analytic standard error.
  key <- cbind(match(as.character(fit$estimate$j), object$states),
               match(as.character(fit$estimate$l), object$states),
               match(as.character(fit$estimate$h), object$states))
  fit$estimate <- fit$estimate |>
    dplyr::mutate(se.boot = apply(boot, 1:3, stats::sd)[key])

  # Return the fit with the replicate tensors (and first-order matrices) with class "P2boot"
  fit$boot <- boot
  fit$boot1 <- boot1
  fit$B <- B
  class(fit) <- c("P2boot", "P2est")
  fit
}

# Short description of a bootstrap fit
#' @export
print.P2boot <- function(x, ...) {
  cat(sprintf("<P2boot>  RPE estimates with %d subject-level bootstrap replicates\n", x$B))
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %.0f%% percentile intervals for n-step predictions; %d subjects\n", 100 * x$conf.level, x$n))
  cat(sprintf("  %d estimated transition probabilities (h -> j -> l)\n", nrow(x$estimate)))
  invisible(x)
}

# internal: per-subject counts of every observed (h, j, l) triple -> Returns the matrix `mat` (n subjects x K observed triples) with the number of times each subject made each triple, the triples themselves as state indices (h, j, l), and the (h, j) pair of each triple (`pair`), so that replicate totals can be turned into probabilities by grouped sums.
.id_counts <- function(object) {
  st <- object$states
  M <- length(st)
  # Count of every (subject, triple) combination.
  cnt <- object$triples |>
    dplyr::count(id, h, j, l, name = "N")
  # Code each triple as one integer, (h-1) M^2 + (j-1) M + l, and give each subject and each distinct triple a row and a column of `mat`.
  ids <- unique(object$triples$id)
  code <- (match(as.character(cnt$h), st) - 1L) * M^2 + (match(as.character(cnt$j), st) - 1L) * M + match(as.character(cnt$l), st)
  ucode <- sort(unique(code))
  mat <- matrix(0, length(ids), length(ucode))
  mat[cbind(match(cnt$id, ids), match(code, ucode))] <- cnt$N
  # Decode the triples back to state indices.
  th <- (ucode - 1L) %/% M^2 + 1L
  tj <- ((ucode - 1L) %/% M) %% M + 1L
  tl <- (ucode - 1L) %% M + 1L
  list(mat = mat, h = th, j = tj, l = tl,
       pair = (th - 1L) * M + tj,  # (h, j) pair of each triple
       first = (tj - 1L) * M + tl) # (j, l) first-order transition of each triple
}

# internal: replicate counts -> tensor of probabilities -> The RPE of one replicate: each triple's count divided by the count of its (h, j) pair.
.tensor_from_counts <- function(tot, ct, object, P_point) {
  M <- length(object$states)
  # At-risk count of the (h, j) pair of every triple (sum over l).
  Y <- stats::ave(tot, ct$pair, FUN = sum)
  # Probabilities of the pairs with somebody at risk.
  P <- array(0, c(M, M, M), dimnames = dimnames(P_point))
  ok <- Y > 0
  P[cbind(ct$j, ct$l, ct$h)[ok, , drop = FALSE]] <- tot[ok] / Y[ok]
  # A pair observed in the data but with nobody at risk in this replicate keeps the point estimate, so that every row of the replicate tensor is still a probability distribution.
  for (p in unique(ct$pair[!ok])) {
    hh <- (p - 1L) %/% M + 1L
    jj <- (p - 1L) %% M + 1L
    P[jj, , hh] <- P_point[jj, , hh]
  }
  # Absorbing states stay absorbing.
  for (a in match(object$absorbing, object$states)) P[a, a, ] <- 1
  P
}

# internal: replicate counts -> first-order matrix -> The first-order RPE of one replicate: the counts of every (j, l) transition pooled over the previous state h, divided by the subject-instants at risk in j. A state with nobody at risk in the replicate keeps the point estimate (`P1_point`); absorbing states keep probability 1.
.first_from_counts <- function(tot, ct, object, P1_point) {
  M <- length(object$states)
  N1 <- tapply(tot, ct$first, sum)                    # pooled over h
  code <- as.integer(names(N1))
  jj <- (code - 1L) %/% M + 1L
  ll <- (code - 1L) %% M + 1L
  Y1 <- tapply(as.numeric(N1), jj, sum)               # at risk in j
  den <- as.numeric(Y1[as.character(jj)])
  P1 <- matrix(0, M, M, dimnames = list(object$states, object$states))
  P1[cbind(jj, ll)] <- ifelse(den > 0, as.numeric(N1) / den, 0)
  if (!is.null(P1_point))
    for (j0 in unique(jj[den == 0])) P1[j0, ] <- P1_point[j0, ]
  for (a in match(object$absorbing, object$states)) {
    P1[a, ] <- 0
    P1[a, a] <- 1
  }
  P1
}

#  internal: n-step curves of every bootstrap replicate -> nsteps x B matrix: column b is P(X_{n+1} = l | X_1 = j, X_0 = h) under the b-th replicate tensor.
.boot_curves <- function(boot, hi, ji, li, nsteps) {
  M <- dim(boot)[1L]
  purrr::map(seq_len(dim(boot)[4L]), \(b)
    .propagate(.pair_matrix(boot[, , , b], M), hi, ji, nsteps, M)[, li]) |>
    unlist() |>
    matrix(nrow = nsteps)
}

# internal: percentile intervals of the n-step curves -> For every step, the alpha/2 and 1 - alpha/2 quantiles (and the standard deviation) of the replicate curves.
.boot_bands <- function(boot, hi, ji, li, nsteps, conf.level) {
  curves <- .boot_curves(boot, hi, ji, li, nsteps)
  a <- (1 - conf.level) / 2
  list(lower = apply(curves, 1L, stats::quantile, probs = a, names = FALSE),
       upper = apply(curves, 1L, stats::quantile, probs = 1 - a, names = FALSE),
       se    = apply(curves, 1L, stats::sd))
}
