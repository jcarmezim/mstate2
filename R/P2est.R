#' Estimate 1-step second-order transition probabilities
#'
#' Estimates \eqn{P_{hj\ell} = P(X_s = \ell \mid X_{s-1} = j, X_{s-2} = h)} from an \code{\link{prep2}} object with the relative probability estimator (RPE) of the methods paper (Najera-Zuloaga, Besalu and Gomez Melis, 2025), with its standard error and confidence interval. The RPE (Eq. 9) pools all transitions and all subject-instants at risk of each pair \eqn{(h, j)}: \deqn{\tilde P_{hj\ell} = \sum_s \tilde N_{hj\ell}(s) / \sum_s \tilde Y_{hj}(s-1).} Its estimated asymptotic variance (Theorems 4-5, Eqs. 13-14) gives the standard error \deqn{se = \sqrt{\tilde P_{hj\ell}(1 - \tilde P_{hj\ell}) / \sum_s \tilde Y_{hj}(s-1)},} and the Wald interval \eqn{\tilde P_{hj\ell} \pm z_{1-\alpha/2}\, se} is that of Corollary 2.
#'
#'
#' @param object An "msm2data" object from \code{\link{prep2}}.
#' @param conf.level Confidence level. Default 0.95.
#' @param ci Confidence-interval type: "wald" (default) or "logit". The logit interval is built on the log-odds scale (delta method) and stays inside (0, 1), behaving better for probabilities near 0 or 1. At the boundary (p = 0 or p = 1, e.g. an absorbing-state row) the logit transform is undefined and the interval degenerates to the point estimate [p, p] rather than [0, 1].
#' @param clip For ci = "wald", clip the interval to [0, 1] (default TRUE).
#' @return An object of class "P2est": the \code{estimate} tibble (one row per observed transition \eqn{h \to j \to \ell}, with columns \code{h, j, l, p, se, lower, upper, n.trans, at.risk}), the point/CI/se tensors \code{P}, \code{P.lower}, \code{P.upper}, \code{P.se} (layout \code{P[j, l, h]}), and metadata.
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
#' fit <- P2est(prep2(panel))
#' fit$estimate
#' @export
P2est <- function(object, conf.level = 0.95, ci = c("wald", "logit"), clip = TRUE) {

  # Check the arguments and get the normal quantile z_{1 - alpha/2}.

  stopifnot(inherits(object, "msm2data"))
  ci <- match.arg(ci)
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2) # 1.959964 for a 95% interval
  states <- object$states
  M <- length(states)

  # Denominator of the RPE: at.risk = sum_s Y_hj(s-1), the number of subject-instants at risk of each pair (h, j).
  at_risk <- object$Y |>
    dplyr::summarise(at.risk = sum(Y), .by = c("h", "j"))

  # Numerator and estimate: n.trans = sum_s N_hjl(s), the number of observed h -> j -> l transitions; p = n.trans / at.risk, with the standard error of Theorem 5, sqrt(p (1 - p) / at.risk).
  est <- object$N |>
    dplyr::summarise(n.trans = sum(N), .by = c("h", "j", "l")) |>
    dplyr::left_join(at_risk, by = c("h", "j")) |>
    dplyr::mutate(p  = n.trans / at.risk,
                  se = sqrt(p * (1 - p) / at.risk))

  # Confidence interval.
  if (ci == "wald") {
    # Wald interval clipped to [0, 1]
    est <- est |>
      dplyr::mutate(lower = p - z * se,
                    upper = p + z * se)
    if (clip)
      est <- est |>
        dplyr::mutate(lower = pmax(0, lower),
                      upper = pmin(1, upper))
  } else {
    # Logit interval: the Wald interval on the log-odds scale, mapped back with the inverse logit, so it always stays inside (0, 1). By the delta method, se(logit p) = se(p) / (p (1 - p)). At p = 0 or 1 the log-odds are infinite, and the interval is the point [p, p].
    est <- est |>
      dplyr::mutate(se_logit = se / (p * (1 - p)),
                    inside   = p > 0 & p < 1,
                    lower = dplyr::if_else(inside,
                              stats::plogis(stats::qlogis(p) - z * se_logit), p),
                    upper = dplyr::if_else(inside,
                              stats::plogis(stats::qlogis(p) + z * se_logit), p)) |>
      dplyr::select(-"se_logit", -"inside")
  }

  est <- est |>
    dplyr::arrange(h, j, l) |>
    dplyr::select("h", "j", "l", "p", "se", "lower", "upper", "n.trans", "at.risk")

  # Tensors P[j, l, h] = P_hjl: one M x M matrix per previous state h. Each row of `key` gives the (j, l, h) position of one estimated transition, so a single indexed assignment fills the tensor. Pairs (h, j) that are never observed stay at 0.
  key <- cbind(match(as.character(est$j), states), # row = current state j
               match(as.character(est$l), states), # column = next state l
               match(as.character(est$h), states)) # slice = previous state h
  
  blank <- array(0, c(M, M, M), dimnames = list(states, states, states))
  P <- Pl <- Pu <- Pse <- blank
  P[key] <- est$p
  Pl[key] <- est$lower
  Pu[key] <- est$upper
  Pse[key] <- est$se

  # Absorbing states: nobody leaves an absorbing state a, so P_{h a a} = 1 for every previous state h. Setting it also for pairs (h, a) without data keeps the chain of ckequations() from losing probability; for observed pairs the estimate is already exactly 1.
  for (a in object$absorbing) {
    ia <- match(a, states)
    P[ia, ia, ] <- Pl[ia, ia, ] <- Pu[ia, ia, ] <- 1
  }

  # Return the estimates table, the tensors and the settings used.
  structure(
    list(estimate  = est, 
        P = P, 
        P.lower = Pl, 
        P.upper = Pu, 
        P.se = Pse, 
        states = states, 
        absorbing = object$absorbing, 
        estimator = "RPE", 
        conf.level = conf.level, 
        ci = ci, 
        n = object$n),
    class = "P2est"
  )
}

# Short description of a fit, shown when it is printed.
#' @export
print.P2est <- function(x, ...) {
  cat(sprintf("<P2est>  %s estimates of 1-step second-order transition probabilities\n", x$estimator))
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %.0f%% %s confidence intervals; %d subjects\n", 100 * x$conf.level, x$ci, x$n))
  cat(sprintf("  %d estimated transition probabilities (h -> j -> l)\n", nrow(x$estimate)))
  invisible(x)
}

# internal: check conf.level (used by P2est() and P2boot()) -> An out-of-range conf.level would give qnorm() a probability outside [0, 1] and NaN limits with only a generic warning; stop with a clear message.
.check_conf_level <- function(conf.level) {
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
}
