#' Estimate 1-step second-order transition probabilities
#'
#' Estimates the 1-step second-order homogeneous transition probabilities
#' \eqn{P_{hj\ell}} from an \code{\link{prep2}} object, using either the
#' relative probability estimator (RPE) or the conditional probability
#' estimator (CPE), and returns standard errors and Wald confidence intervals.
#'
#' Let \eqn{\tilde N_{hj\ell}(s)} and \eqn{\tilde Y_{hj}(s-1)} be the triple
#' counts and at-risk counts (see \code{\link{prep2}}), and let
#' \eqn{t_{hj}} be the number of time points at which \eqn{\tilde Y_{hj} > 0}.
#'
#' \strong{RPE} (Definition 6, Eq. 9), with variance from Theorem 5:
#' \deqn{\tilde P_{hj\ell} = \frac{\sum_s \tilde N_{hj\ell}(s)}{\sum_s \tilde Y_{hj}(s-1)},
#'   \quad \widehat{\mathrm{se}} = \sqrt{\frac{\tilde P_{hj\ell}(1-\tilde P_{hj\ell})}{\sum_s \tilde Y_{hj}(s-1)}}.}
#'
#' \strong{CPE} (Definition 5, Eq. 8), with variance from Theorem 3:
#' \deqn{\hat P_{hj\ell} = \frac{1}{t_{hj}} \sum_s \frac{\tilde N_{hj\ell}(s)}{\tilde Y_{hj}(s-1)},
#'   \quad \widehat{\mathrm{se}} = \sqrt{\hat P_{hj\ell}(1-\hat P_{hj\ell}) \frac{\sum_s 1/\tilde Y_{hj}(s-1)}{t_{hj}^2}}.}
#'
#' Confidence intervals are Wald intervals (Corollaries 1 and 2) clipped to
#' \eqn{[0,1]}.
#'
#' @param object An object of class \code{"msm2data"} from \code{\link{prep2}}.
#' @param estimator Either \code{"RPE"} (default) or \code{"CPE"}. The RPE is
#'   recommended: it is at least as efficient as the CPE (Corollary 3).
#' @param conf.level Confidence level for the intervals. Default \code{0.95}.
#'
#' @return An object of class \code{"P2est"}: a list with
#'   \item{estimate}{data.frame of \eqn{(h, j, l, p, se, lower, upper)} for every
#'     observed triple.}
#'   \item{P}{the \eqn{M \times M \times M} point-estimate tensor with layout
#'     \code{P[j, l, h]} \eqn{= P_{hj\ell}} (the layout consumed by
#'     \code{ckequations}); absorbing states carry \eqn{P_{aaa}=1}.}
#'   \item{P.lower, P.upper}{tensors of CI bounds in the same layout, clipped to [0,1].}
#'   \item{P.se}{tensor of standard errors.}
#'   \item{states, absorbing, estimator, conf.level, n}{metadata.}
#' @export
P2est <- function(object, estimator = c("RPE", "CPE"), conf.level = 0.95) {
  stopifnot(inherits(object, "msm2data"))
  estimator <- match.arg(estimator)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  N <- object$N; Y <- object$Y
  states <- object$states
  M <- length(states)

  ## aggregate at-risk per (h, j): total exposure, sum of inverses, number of times
  Yagg <- Y[, list(sumY  = sum(Y),
                   sumInv = sum(1 / Y),
                   thj    = .N),               # times with Y > 0
            by = list(h, j)]

  ## total triple count per (h, j, l) summed over time s
  Nagg <- N[, list(sumN = sum(N)), by = list(h, j, l)]

  est <- merge(Nagg, Yagg, by = c("h", "j"), all.x = TRUE)

  if (estimator == "RPE") {
    est[, p  := sumN / sumY]
    est[, se := sqrt(p * (1 - p) / sumY)]
  } else {
    ## CPE needs the per-time ratio, so aggregate N_{hjl}(s)/Y_{hj}(s-1)
    NY <- merge(N, Y, by = c("h", "j", "s"))
    cpe <- NY[, list(p = sum(N / Y)), by = list(h, j, l)]
    cpe <- merge(cpe, Yagg, by = c("h", "j"))
    cpe[, p  := p / thj]
    cpe[, se := sqrt(p * (1 - p) * sumInv / thj^2)]
    est <- cpe
  }

  est[, lower := pmax(0, p - z * se)]
  est[, upper := pmin(1, p + z * se)]
  est <- est[order(h, j, l)]

  ## ---- build tensors in the P[j, l, h] = P_{hjl} layout used by ckequations ----
  P  <- array(0, dim = c(M, M, M), dimnames = list(states, states, states))
  Pl <- Pu <- Pse <- P
  idx <- function(v) as.integer(factor(as.character(v), levels = states))
  ih <- idx(est$h); ij <- idx(est$j); il <- idx(est$l)
  for (r in seq_len(nrow(est))) {
    P [ij[r], il[r], ih[r]] <- est$p[r]
    Pl[ij[r], il[r], ih[r]] <- est$lower[r]
    Pu[ij[r], il[r], ih[r]] <- est$upper[r]
    Pse[ij[r], il[r], ih[r]] <- est$se[r]
  }
  ## absorbing states: once in state a you stay, regardless of the previous
  ## state, so P_{x a a} = 1 for every x. In the P[j,l,h] layout that is
  ## P[a, a, h] = 1 for all h (and all other P[a, ., h] = 0).
  for (a in object$absorbing) {
    ia <- match(a, states)
    P [ia, ia, ] <- 1
    Pl[ia, ia, ] <- 1
    Pu[ia, ia, ] <- 1
  }

  out <- list(
    estimate  = as.data.frame(est[, list(h, j, l, p, se, lower, upper)]),
    P = P, P.lower = Pl, P.upper = Pu, P.se = Pse,
    states = states, absorbing = object$absorbing,
    estimator = estimator, conf.level = conf.level, n = object$n
  )
  class(out) <- "P2est"
  out
}

#' @export
print.P2est <- function(x, ...) {
  cat(sprintf("<P2est>  %s estimates of 1-step second-order transition probabilities\n",
              x$estimator))
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %.0f%% confidence intervals; %d subjects\n",
              100 * x$conf.level, x$n))
  cat(sprintf("  %d estimated transition probabilities (h -> j -> l)\n",
              nrow(x$estimate)))
  invisible(x)
}
