#' Estimate 1-step second-order transition probabilities
#'
#' Estimates \eqn{P_{hj\ell} = P(X_s = \ell \mid X_{s-1} = j, X_{s-2} = h)} from
#' an \code{\link{prep2}} object with the relative probability estimator (RPE,
#' Eq. 9 / Theorem 5), with standard errors and confidence intervals
#' (Corollaries 1-2). The RPE pools all transitions and all subject-instants at
#' risk of each pair \eqn{(h, j)}:
#' \deqn{\hat P_{hj\ell} = \sum_s \tilde N_{hj\ell}(s) / \sum_s \tilde Y_{hj}(s-1).}
#'
#' @param object An "msm2data" object from \code{\link{prep2}}.
#' @param conf.level Confidence level. Default 0.95.
#' @param ci Confidence-interval type: "wald" (default) or "logit". The logit
#'   interval is built on the log-odds scale (delta method) and stays inside
#'   (0, 1), behaving better for probabilities near 0 or 1. At the boundary
#'   (p = 0 or p = 1, e.g. an absorbing-state row) the logit transform is
#'   undefined and the interval degenerates to the point estimate [p, p]
#'   rather than [0, 1].
#' @param clip For ci = "wald", clip the interval to [0, 1] (default TRUE).
#' @return An object of class "P2est": the \code{estimate} data frame, the
#'   point/CI/se tensors \code{P}, \code{P.lower}, \code{P.upper}, \code{P.se}
#'   (layout \code{P[j, l, h]}), and metadata.
#' @section Interval choice:
#' Wald intervals for a proportion are known to behave poorly when the
#' probability is close to 0 or 1 or the number at risk is small (Brown, Cai
#' and DasGupta, 2001); \code{ci = "logit"} is then preferable. The estimator
#' itself is the relative probability estimator of the methods paper, which is
#' more efficient than the conditional probability estimator studied there and
#' is therefore the only one implemented.
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for a
#' binomial proportion. \emph{Statistical Science}, 16(2), 101-133.
#'
#' Najera-Zuloaga, J., Besalu, M. and Gomez Melis, G. (2025). Second-order
#' Markov multistate models: nonparametric estimation and inference.
#' Manuscript submitted for publication.
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
#' fit$estimate
#' @export
P2est <- function(object, conf.level = 0.95, ci = c("wald", "logit"), clip = TRUE) {
  stopifnot(inherits(object, "msm2data"))
  ci <- match.arg(ci)
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)   # e.g. 1.96 for a 95% Wald/logit CI

  N <- object$N; Y <- object$Y
  states <- object$states; M <- length(states)

  ## --- RPE point estimate + standard error (Eq. 9 / Theorem 5) ---
  ## at.risk = sum_s Y_hj(s-1), the total subject-instants at risk of (h, j);
  ## n.trans = sum_s N_hjl(s), the total observed h -> j -> l transitions.
  ## The estimate is their ratio, with the binomial standard error
  ## sqrt(p(1-p)/at.risk).
  Yagg <- Y[, .(at.risk = sum(Y)), keyby = .(h, j)]
  Nagg <- N[, .(n.trans = sum(N)), keyby = .(h, j, l)]
  est <- merge(Nagg, Yagg, by = c("h", "j"), all.x = TRUE)
  est[, p  := n.trans / at.risk]
  est[, se := sqrt(p * (1 - p) / at.risk)]

  ## --- confidence interval (Corollaries 1-2) ---
  if (ci == "wald") {
    ## Standard symmetric normal-approximation interval on the probability
    ## scale; `clip` truncates it to stay within [0, 1] since the raw
    ## interval can otherwise poke outside that range when p is close to a
    ## boundary and se isn't tiny.
    est[, lower := p - z * se]
    est[, upper := p + z * se]
    if (clip) { est[, lower := pmax(0, lower)]; est[, upper := pmin(1, upper)] }
  } else {                                                   # logit / delta method
    ## Build the interval on the log-odds (logit) scale, where a normal
    ## approximation behaves much better near 0/1, then map both endpoints
    ## back to the probability scale with the inverse-logit (plogis). The
    ## delta method turns se(p) into se(logit(p)) = se(p) / (p(1-p)) (the
    ## derivative of logit at p is 1/(p(1-p))).
    se_logit <- est$se / (est$p * (1 - est$p))
    ## At p = 0 or 1 (e.g. any absorbing-state row) logit(p) is +-Inf and
    ## se_logit is Inf/0 (undefined either way), so those rows fall back to
    ## the degenerate point interval [p, p] instead of computing garbage.
    inside   <- est$p > 0 & est$p < 1
    est[, lower := data.table::fifelse(inside,
            stats::plogis(stats::qlogis(p) - z * se_logit), p)]
    est[, upper := data.table::fifelse(inside,
            stats::plogis(stats::qlogis(p) + z * se_logit), p)]
  }

  data.table::setorder(est, h, j, l)

  ## --- build the M x M x M tensors by vectorized array indexing ---
  ## `key` is a 3-column integer matrix, one row per estimated (h, j, l)
  ## triple; P[key] <- est$p (etc.) is R's vectorized array-assignment form,
  ## equivalent to looping "P[j_i, l_i, h_i] <- p_i" over every row of `key`
  ## but without an explicit loop. This is what turns the long "one row per
  ## (h, j, l)" `est` table into the P[j, l, h] tensor layout ckequations()
  ## and compare2() build their pair-transition matrix from.
  key <- cbind(match(as.character(est$j), states),   # row    = current  state j
               match(as.character(est$l), states),   # column = next     state l
               match(as.character(est$h), states))   # slice  = previous state h
  blank <- array(0, c(M, M, M), dimnames = list(states, states, states))
  P <- Pl <- Pu <- Pse <- blank
  P[key] <- est$p; Pl[key] <- est$lower; Pu[key] <- est$upper; Pse[key] <- est$se

  ## An absorbing state a has, by construction, no possible destination other
  ## than itself, so P_{h,a,a} = 1 for every h -- including any (h, a) pair
  ## with zero observed data, which the loop above would otherwise leave at
  ## the tensor's blank default of 0 and which would then incorrectly break
  ## the chain in ckequations()'s propagation (a "0% chance of staying put"
  ## absorbing state isn't absorbing at all). For an (h, a) pair that *was*
  ## observed, this is a no-op: an absorbing state's self-transition
  ## proportion is already exactly 1 in the raw counts (nobody at risk can go
  ## anywhere else), so p = 1 and se = 0 there regardless.
  for (a in object$absorbing) {                       # P_{x a a} = 1 for all x
    ia <- match(a, states)
    P[ia, ia, ] <- Pl[ia, ia, ] <- Pu[ia, ia, ] <- 1
  }

  structure(
    list(estimate  = as.data.frame(est[, .(h, j, l, p, se, lower, upper,
                                           n.trans, at.risk)]),
         P = P, P.lower = Pl, P.upper = Pu, P.se = Pse,
         states = states, absorbing = object$absorbing,
         estimator = "RPE", conf.level = conf.level, ci = ci, n = object$n),
    class = "P2est"
  )
}

#' @export
print.P2est <- function(x, ...) {
  cat(sprintf("<P2est>  %s estimates of 1-step second-order transition probabilities\n",
              x$estimator))
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %.0f%% %s confidence intervals; %d subjects\n",
              100 * x$conf.level, x$ci, x$n))
  cat(sprintf("  %d estimated transition probabilities (h -> j -> l)\n",
              nrow(x$estimate)))
  invisible(x)
}

## --- internal: shared conf.level guard, used by P2est(), compare_order() and
## summary.P2reg() -- without it, an out-of-range conf.level (e.g. > 1) feeds
## qnorm() a probability outside [0, 1] and silently returns NaN bounds, with
## only R's generic "NaNs produced" warning (no mention of which argument, or
## which function, is at fault) to go on.
.check_conf_level <- function(conf.level) {
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
}
