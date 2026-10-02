#' Subject-level bootstrap of the second-order (and first-order) estimates
#'
#' Resamples subjects with replacement and recomputes the RPE of every
#' 1-step second-order transition probability (and of its first-order,
#' history-pooled counterpart) in each replicate. The resulting object behaves
#' like a \code{\link{P2est}} fit, but carries the bootstrap tensors, so that
#' \code{\link{ckequations}}, \code{\link{compare2}}, \code{\link{probtrans2}}
#' and \code{\link{compare_order}} can report \strong{percentile bootstrap
#' intervals} for the \eqn{n}-step probabilities instead of the heuristic
#' evolution intervals obtained by propagating the one-step confidence limits.
#'
#' Resampling whole subjects keeps the dependence between the several
#' instants a subject contributes to the same or to different \eqn{(h, j)}
#' risk sets. The implementation never re-reads the data: the per-subject
#' counts of every observed \eqn{(h, j, \ell)} triple are tabulated once, and
#' each replicate is a weighted sum of those counts, so thousands of
#' replicates are cheap even for large cohorts.
#'
#' In a replicate where an observed \eqn{(h, j)} pair happens to have no
#' subject at risk, its row is taken from the point estimate, so that the
#' replicate tensor stays a proper set of transition probabilities.
#'
#' @param object An "msm2data" object from \code{\link{prep2}}.
#' @param B Number of bootstrap replicates. Default 200.
#' @param conf.level Confidence level of the intervals. Default 0.95.
#' @param seed Optional integer seed, for reproducible replicates.
#' @return An object of class \code{c("P2boot", "P2est")}: all the components
#'   of \code{P2est(object, conf.level)} plus \code{boot}, an
#'   \eqn{M \times M \times M \times B} array of replicate tensors (layout
#'   \code{[j, l, h, b]}), \code{boot1}, an \eqn{M \times M \times B} array of
#'   replicate first-order matrices (see \code{\link{P1est}}), and \code{B}.
#'   The \code{estimate} table gains a column \code{se.boot} (bootstrap
#'   standard deviation).
#' @seealso \code{\link{ckequations}}, \code{\link{compare2}}, \code{\link{P1est}}
#' @section Why a bootstrap of subjects:
#' Patients are the independent units; the days of one patient are not.
#' Resampling whole patients keeps that within-patient dependence, which is
#' the standard non-parametric bootstrap for clustered data (Davison and
#' Hinkley, 1997; Field and Welsh, 2007). The n-step predictions are
#' non-linear (polynomial) functions of all the estimated probabilities, so a
#' delta-method variance would be cumbersome; propagating each replicate and
#' taking percentiles (Efron and Tibshirani, 1993) gives their intervals
#' directly. Because the same replicates serve every curve, the difference
#' between two curves can be tested on the paired replicates
#' (\code{\link{overlap_step}}, \code{\link{divergence}}). In the package's
#' simulation study the 95\% percentile intervals of n-step predictions had
#' coverage 0.91-0.94, while the evolution intervals had coverage 1 and were
#' up to five times wider.
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
#' st   <- c("A", "B", "C")                                 # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' bt <- P2boot(prep2(panel), B = 100, seed = 1)
#' bt
#' ckequations(bt, h = "A", j = "B", l = "C", nsteps = 4, bounds = TRUE)
#' @export
P2boot <- function(object, B = 200, conf.level = 0.95, seed = NULL) {
  stopifnot(inherits(object, "msm2data"))
  .check_conf_level(conf.level)
  if (!is.numeric(B) || length(B) != 1L || B < 2)
    stop("`B` must be a single number >= 2.", call. = FALSE)
  B <- as.integer(B)
  if (!is.null(seed)) set.seed(seed)

  fit <- P2est(object, conf.level = conf.level)
  ct  <- .id_counts(object)
  n   <- nrow(ct$mat)
  M   <- length(object$states)

  boot  <- array(0, c(M, M, M, B), dimnames = c(dimnames(fit$P), list(NULL)))
  boot1 <- array(0, c(M, M, B),    dimnames = c(dimnames(fit$P)[1:2], list(NULL)))
  for (b in seq_len(B)) {
    ## multiplicity of each subject in this replicate, then the replicate's
    ## total count of every (h, j, l) triple as a weighted column sum
    w   <- tabulate(sample.int(n, n, replace = TRUE), n)
    tot <- as.vector(w %*% ct$mat)
    boot[, , , b] <- .tensor_from_counts(tot, ct, object, fit$P)
    boot1[, , b]  <- .first_from_counts(tot, ct, object)
  }

  ## bootstrap SD of each estimated probability, next to the analytic se
  key <- cbind(match(as.character(fit$estimate$j), object$states),
               match(as.character(fit$estimate$l), object$states),
               match(as.character(fit$estimate$h), object$states))
  fit$estimate$se.boot <- apply(boot, 1:3, stats::sd)[key]

  fit$boot <- boot; fit$boot1 <- boot1; fit$B <- B
  class(fit) <- c("P2boot", "P2est")
  fit
}

#' @export
print.P2boot <- function(x, ...) {
  cat(sprintf("<P2boot>  RPE estimates with %d subject-level bootstrap replicates\n", x$B))
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %.0f%% percentile intervals for n-step predictions; %d subjects\n",
              100 * x$conf.level, x$n))
  cat(sprintf("  %d estimated transition probabilities (h -> j -> l)\n", nrow(x$estimate)))
  invisible(x)
}

## --- internal: per-subject counts of every observed (h, j, l) triple ---------
## Returns the n_subjects x K matrix `mat` (K = number of distinct observed
## triples), the K triples themselves (`types`, as state indices) and, for
## each triple, the index of its (h, j) pair and of its (j, l) first-order
## transition, so that replicate totals can be turned into probabilities by
## simple grouped sums.
.id_counts <- function(object) {
  tr <- object$triples
  st <- object$states
  cnt <- tr[, .(N = .N), by = .(id, h, j, l)]
  ids <- unique(tr$id)
  hi <- match(as.character(cnt$h), st); ji <- match(as.character(cnt$j), st)
  li <- match(as.character(cnt$l), st)
  code <- (hi - 1L) * length(st)^2 + (ji - 1L) * length(st) + li
  ucode <- sort(unique(code))
  mat <- matrix(0, length(ids), length(ucode))
  mat[cbind(match(cnt$id, ids), match(code, ucode))] <- cnt$N
  M <- length(st)
  th <- (ucode - 1L) %/% M^2 + 1L
  tj <- ((ucode - 1L) %/% M) %% M + 1L
  tl <- (ucode - 1L) %% M + 1L
  list(mat = mat, h = th, j = tj, l = tl,
       pair = (th - 1L) * M + tj,              # (h, j) group of each triple
       first = (tj - 1L) * M + tl)             # (j, l) group of each triple
}

## --- internal: replicate totals -> second-order tensor -----------------------
.tensor_from_counts <- function(tot, ct, object, P_point) {
  M  <- length(object$states)
  Y  <- stats::ave(tot, ct$pair, FUN = sum)               # at-risk count of the pair
  P  <- array(0, c(M, M, M), dimnames = dimnames(P_point))
  ok <- Y > 0
  P[cbind(ct$j, ct$l, ct$h)[ok, , drop = FALSE]] <- tot[ok] / Y[ok]
  ## pairs with nobody at risk in this replicate: fall back to the point estimate
  for (p in unique(ct$pair[!ok])) {
    hh <- (p - 1L) %/% M + 1L; jj <- (p - 1L) %% M + 1L
    P[jj, , hh] <- P_point[jj, , hh]
  }
  for (a in match(object$absorbing, object$states)) P[a, a, ] <- 1
  P
}

## --- internal: replicate totals -> first-order matrix (pooled over h) --------
.first_from_counts <- function(tot, ct, object) {
  M  <- length(object$states)
  N1 <- tapply(tot, ct$first, sum)                        # pooled over h
  jj <- (as.integer(names(N1)) - 1L) %/% M + 1L
  ll <- (as.integer(names(N1)) - 1L) %% M + 1L
  Y1 <- tapply(as.numeric(N1), jj, sum)
  P1 <- matrix(0, M, M, dimnames = list(object$states, object$states))
  den <- as.numeric(Y1[as.character(jj)])
  P1[cbind(jj, ll)] <- ifelse(den > 0, as.numeric(N1) / den, 0)
  for (a in match(object$absorbing, object$states)) { P1[a, ] <- 0; P1[a, a] <- 1 }
  P1
}

## --- internal: n-step curves of every bootstrap replicate --------------------
## nsteps x B matrix: column b = P(X_{n+1} = l | X_1 = j, X_0 = h) under the
## b-th replicate tensor (second order) or matrix (first order).
.boot_curves <- function(boot, hi, ji, li, nsteps) {
  M <- dim(boot)[1L]
  matrix(vapply(seq_len(dim(boot)[4L]), function(b)
    .propagate(.pair_matrix(boot[, , , b], M), hi, ji, nsteps, M)[, li],
    numeric(nsteps)), nrow = nsteps)
}
.boot_curves_first <- function(boot1, ji, li, nsteps) {
  M <- dim(boot1)[1L]
  matrix(vapply(seq_len(dim(boot1)[3L]), function(b)
    .propagate(.pair_matrix(array(boot1[, , b], c(M, M, M)), M), 1L, ji, nsteps, M)[, li],
    numeric(nsteps)), nrow = nsteps)
}

## --- internal: percentile interval of a paired difference ---------------------
## Given the replicate curves of two groups computed on the SAME bootstrap
## resamples, the percentile interval of their difference at each step, and
## the number of leading steps whose interval excludes 0.
.boot_diff <- function(c1, c2, conf.level) {
  dd <- c1 - c2
  a  <- (1 - conf.level) / 2
  lo <- apply(dd, 1L, stats::quantile, probs = a, names = FALSE)
  up <- apply(dd, 1L, stats::quantile, probs = 1 - a, names = FALSE)
  sep <- lo > 0 | up < 0
  list(diff = data.frame(n = seq_len(nrow(dd)), lower = lo, upper = up),
       steps = if (all(sep)) length(sep) else which(!sep)[1L] - 1L)
}

## --- internal: percentile bands of n-step curves over bootstrap replicates ---
## Propagates every replicate tensor from (h, j) and returns, per step, the
## (alpha/2, 1 - alpha/2) quantiles of the probability of being in state li.
.boot_bands <- function(boot, hi, ji, li, nsteps, conf.level,
                        curves = .boot_curves(boot, hi, ji, li, nsteps)) {
  a <- (1 - conf.level) / 2
  list(lower = apply(curves, 1L, stats::quantile, probs = a, names = FALSE),
       upper = apply(curves, 1L, stats::quantile, probs = 1 - a, names = FALSE),
       se    = apply(curves, 1L, stats::sd))
}
