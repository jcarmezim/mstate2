#' Test the first-order Markov assumption, transition by transition
#'
#' For each observed pair (current state \code{j}, target state \code{l}), tests
#' whether the 1-step second-order transition probability \eqn{P_{hj\ell}}
#' depends on the preceding state \eqn{h}. Under the null hypothesis of
#' first-order Markov behaviour for that transition, \eqn{P_{hj\ell}} is the
#' same for every observed preceding state \eqn{h}.
#'
#' The test statistic is an inverse-variance-weighted heterogeneity (Cochran's
#' Q-type) statistic built directly on the RPE/CPE estimates and their
#' consistent standard errors already computed by \code{\link{P2est}}
#' (Theorems 2-5 of Najera-Zuloaga, Besalu and Gomez Melis, 2025): writing
#' \eqn{w_h = 1/\mathrm{se}_{hj\ell}^2} and \eqn{\bar P_{j\ell} = \sum_h w_h
#' P_{hj\ell} / \sum_h w_h} for the pooled (history-blind) estimate,
#' \deqn{Q_{j\ell} = \sum_h w_h (P_{hj\ell} - \bar P_{j\ell})^2}
#' is asymptotically \eqn{\chi^2} with \eqn{K_{j\ell} - 1} degrees of freedom
#' under the null, where \eqn{K_{j\ell}} is the number of preceding states with
#' data for that (j, l) pair. This generalises, into a single p-value per
#' transition, the confidence-interval-overlap check used informally in the
#' paper's illustration (its Table 2) and formalised pairwise (for exactly two
#' preceding states) by \code{\link{overlap_step}}.
#'
#' This is \emph{not} the score/landmarking test of Titman and Putter (2020);
#' it is a Wald-type test built on the second-order model's own asymptotics,
#' in the same spirit and requiring no separate model fit. It treats the
#' subject-groups defined by different \eqn{h} as independent samples, which
#' holds for a non-reversible process in which each subject visits state
#' \eqn{j} with a single preceding state per visit (the typical multistate
#' disease-progression setting).
#'
#' @param object A "P2est" object.
#' @param min.h Minimum number of preceding states \eqn{h} with a well-defined
#'   (non-degenerate, se > 0) estimate required to test a given (j, l) pair.
#'   Default 2.
#' @return A data frame of class "markov_test", one row per tested (j, l) pair,
#'   sorted by ascending p-value, with columns \code{j, l, K, df, statistic,
#'   p.value, pooled} (the inverse-variance-weighted pooled estimate under the
#'   null).
#' @export
markov_test <- function(object, min.h = 2) {
  stopifnot(inherits(object, "P2est"))
  if (min.h < 2) stop("`min.h` must be >= 2: a heterogeneity test needs at ",
                      "least two preceding-state groups.", call. = FALSE)
  est <- data.table::as.data.table(object$estimate)
  est[, `:=`(h = as.character(h), j = as.character(j), l = as.character(l))]

  run_one <- function(sub) {
    sub <- sub[se > 0]                          # drop degenerate (se = 0) groups
    K <- nrow(sub)
    if (K < min.h) return(NULL)
    w    <- 1 / sub$se^2
    pbar <- sum(w * sub$p) / sum(w)
    Q    <- sum(w * (sub$p - pbar)^2)
    df   <- K - 1L
    data.table::data.table(K = K, df = df, statistic = Q,
                           p.value = stats::pchisq(Q, df, lower.tail = FALSE),
                           pooled = pbar)
  }

  out <- est[, run_one(.SD), by = .(j, l)]
  if (is.null(out) || nrow(out) == 0L) {
    warning("No (j, l) transition had >= min.h preceding states with a ",
            "non-degenerate estimate; nothing to test.", call. = FALSE)
    out <- data.table::data.table(j = character(), l = character(), K = integer(),
                                  df = integer(), statistic = numeric(),
                                  p.value = numeric(), pooled = numeric())
  }
  data.table::setorder(out, p.value)
  structure(as.data.frame(out), class = c("markov_test", "data.frame"),
            estimator = object$estimator, conf.level = object$conf.level)
}

#' @export
print.markov_test <- function(x, digits = 4, ...) {
  cat(sprintf("<markov_test>  %s-based test of the first-order Markov assumption\n",
              attr(x, "estimator")))
  cat("  H0 (per transition): P(X_s = l | X_{s-1} = j) does not depend on X_{s-2} = h\n")
  cat(sprintf("  %d transition(s) tested; sorted by ascending p-value\n\n", nrow(x)))
  if (nrow(x)) {
    y <- x
    for (nm in c("statistic", "p.value", "pooled")) y[[nm]] <- round(y[[nm]], digits)
    print(as.data.frame(y), row.names = FALSE)
  }
  invisible(x)
}
