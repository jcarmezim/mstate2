#' Estimate 1-step second-order transition probabilities
#'
#' Estimates \eqn{P_{hj\ell} = P(X_s = \ell \mid X_{s-1} = j, X_{s-2} = h)} from
#' an \code{\link{prep2}} object, using the relative probability estimator (RPE,
#' Eq. 9 / Theorem 5) or the conditional probability estimator (CPE, Eq. 8 /
#' Theorem 3), with standard errors and confidence intervals (Corollaries 1-2).
#'
#' @param object An "msm2data" object from \code{\link{prep2}}.
#' @param estimator "RPE" (default, recommended) or "CPE".
#' @param conf.level Confidence level. Default 0.95.
#' @param ci Confidence-interval type: "wald" (default) or "logit". The logit
#'   interval is built on the log-odds scale (delta method) and stays inside
#'   (0, 1), behaving better for probabilities near 0 or 1.
#' @param clip For ci = "wald", clip the interval to [0, 1] (default TRUE).
#' @return An object of class "P2est": the \code{estimate} data frame, the
#'   point/CI/se tensors \code{P}, \code{P.lower}, \code{P.upper}, \code{P.se}
#'   (layout \code{P[j, l, h]}), and metadata.
#' @export
P2est <- function(object, estimator = c("RPE", "CPE"), conf.level = 0.95,
                  ci = c("wald", "logit"), clip = TRUE) {
  stopifnot(inherits(object, "msm2data"))
  estimator <- match.arg(estimator)
  ci        <- match.arg(ci)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  N <- object$N; Y <- object$Y
  states <- object$states; M <- length(states)

  ## --- per-pair exposure, computed once ---
  Yagg <- Y[, .(at.risk = sum(Y), inv = sum(1 / Y), t.hj = .N), keyby = .(h, j)]
  Nagg <- N[, .(n.trans = sum(N)),                              keyby = .(h, j, l)]

  ## --- point estimate + standard error ---
  if (estimator == "RPE") {                                  # Eq. 9 / Theorem 5
    est <- merge(Nagg, Yagg, by = c("h", "j"), all.x = TRUE)
    est[, p  := n.trans / at.risk]
    est[, se := sqrt(p * (1 - p) / at.risk)]
  } else {                                                   # Eq. 8 / Theorem 3
    ratio <- merge(N, Y, by = c("h", "j", "s"))[, .(p = sum(N / Y)), by = .(h, j, l)]
    est <- merge(merge(ratio, Nagg, by = c("h", "j", "l")), Yagg, by = c("h", "j"))
    est[, p  := p / t.hj]
    est[, se := sqrt(p * (1 - p) * inv / t.hj^2)]
  }

  ## --- confidence interval (Corollaries 1-2) ---
  if (ci == "wald") {
    est[, lower := p - z * se]
    est[, upper := p + z * se]
    if (clip) { est[, lower := pmax(0, lower)]; est[, upper := pmin(1, upper)] }
  } else {                                                   # logit / delta method
    se_logit <- est$se / (est$p * (1 - est$p))
    inside   <- est$p > 0 & est$p < 1
    est[, lower := data.table::fifelse(inside,
            stats::plogis(stats::qlogis(p) - z * se_logit), p)]
    est[, upper := data.table::fifelse(inside,
            stats::plogis(stats::qlogis(p) + z * se_logit), p)]
  }

  data.table::setorder(est, h, j, l)

  ## --- build the M x M x M tensors by vectorized array indexing ---
  key <- cbind(match(as.character(est$j), states),   # row    = current  state j
               match(as.character(est$l), states),   # column = next     state l
               match(as.character(est$h), states))   # slice  = previous state h
  blank <- array(0, c(M, M, M), dimnames = list(states, states, states))
  P <- Pl <- Pu <- Pse <- blank
  P[key] <- est$p; Pl[key] <- est$lower; Pu[key] <- est$upper; Pse[key] <- est$se

  for (a in object$absorbing) {                       # P_{x a a} = 1 for all x
    ia <- match(a, states)
    P[ia, ia, ] <- Pl[ia, ia, ] <- Pu[ia, ia, ] <- 1
  }

  structure(
    list(estimate  = as.data.frame(est[, .(h, j, l, p, se, lower, upper,
                                           n.trans, at.risk, t.hj)]),
         P = P, P.lower = Pl, P.upper = Pu, P.se = Pse,
         states = states, absorbing = object$absorbing,
         estimator = estimator, conf.level = conf.level, ci = ci, n = object$n),
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
