#' Test the first-order Markov assumption, transition by transition
#'
#' For each observed pair (current state \code{j}, target state \code{l}), tests
#' whether the 1-step second-order transition probability \eqn{P_{hj\ell}}
#' depends on the preceding state \eqn{h}. Under the null hypothesis of
#' first-order Markov behaviour for that transition, \eqn{P_{hj\ell}} is the
#' same for every observed preceding state \eqn{h}.
#'
#' Two statistics are available.
#'
#' \code{method = "wald"} (default) is an inverse-variance-weighted
#' heterogeneity (Cochran's Q-type) statistic built directly on the RPE
#' estimates and their consistent standard errors already computed by
#' \code{\link{P2est}} (Najera-Zuloaga, Besalu and Gomez Melis, 2025): writing
#' \eqn{w_h = 1/\mathrm{se}_{hj\ell}^2} and \eqn{\bar P_{j\ell} = \sum_h w_h
#' P_{hj\ell} / \sum_h w_h} for the pooled (history-blind) estimate,
#' \deqn{Q_{j\ell} = \sum_h w_h (P_{hj\ell} - \bar P_{j\ell})^2}
#' is asymptotically \eqn{\chi^2_{K_{j\ell} - 1}} under the null, where
#' \eqn{K_{j\ell}} is the number of preceding states with a non-degenerate
#' estimate (se > 0).
#'
#' \code{method = "lrt"} is the likelihood-ratio statistic of the classical
#' test of Markov order (Anderson and Goodman, 1957), computed from the RPE
#' counts: with \eqn{N_h = \sum_s \tilde N_{hj\ell}(s)},
#' \eqn{Y_h = \sum_s \tilde Y_{hj}(s-1)} and the pooled first-order estimate
#' \eqn{\bar P_{j\ell} = \sum_h N_h / \sum_h Y_h},
#' \deqn{G^2_{j\ell} = 2 \sum_h \left[ N_h \log\frac{N_h}{Y_h \bar P_{j\ell}} +
#'   (Y_h - N_h) \log\frac{Y_h - N_h}{Y_h (1 - \bar P_{j\ell})} \right],}
#' asymptotically \eqn{\chi^2_{K-1}}, where \eqn{K} is now the number of
#' preceding states with subjects at risk in \eqn{j} (histories with no event
#' towards \eqn{\ell} count too: unlike the Wald version, nothing is dropped
#' for a zero standard error). With \code{by = "state"}, the
#' destinations are tested jointly for each current state \eqn{j}:
#' \deqn{G^2_j = 2 \sum_h \sum_\ell N_{h\ell} \log\frac{N_{h\ell}}{Y_h \bar P_{j\ell}},}
#' with \eqn{(K - 1)(L - 1)} degrees of freedom (\eqn{L} = number of observed
#' destinations from \eqn{j}), which is exactly the Anderson-Goodman test of a
#' first- against a second-order chain, state by state.
#'
#' Both versions rely on the likelihood of the second-order chain factorising
#' over \eqn{(h, j)} pairs, which makes the estimates for different preceding
#' states asymptotically independent even when the same subjects contribute to
#' several of them. This is \emph{not} the landmark-based test of Titman and
#' Putter (2022).
#'
#' @param object A "P2est" object.
#' @param min.h Minimum number of preceding states needed to test a given
#'   (j, l) pair (or state j). Default 2.
#' @param method \code{"wald"} (default) or \code{"lrt"}; see Details.
#' @param by \code{"transition"} (default): one test per (j, l) pair;
#'   \code{"state"}: one joint test per current state j (requires
#'   \code{method = "lrt"}).
#' @param p.adjust Method for the multiplicity-adjusted p-values (column
#'   \code{p.adj}), passed to \code{\link[stats]{p.adjust}}. Default
#'   \code{"holm"} (Holm, 1979), which controls the family-wise error rate
#'   without assuming independence between the tests; \code{"none"} keeps the
#'   unadjusted p-values.
#' @return A data frame of class "markov_test", sorted by ascending p-value.
#'   With \code{by = "transition"}, one row per tested (j, l) pair with columns
#'   \code{j, l, K, df, statistic, p.value, pooled, p.adj} (\code{pooled} is the
#'   history-blind estimate under the null: inverse-variance weighted for
#'   \code{"wald"}, the first-order RPE for \code{"lrt"}). With
#'   \code{by = "state"}, one row per current state with columns
#'   \code{j, K, L, df, statistic, p.value, p.adj}. Attributes \code{method}, \code{by},
#'   \code{estimator}, \code{conf.level}.
#' @section Why these tests:
#' The likelihood ratio test is the classical test of the
#' order of a Markov chain observed on many subjects (Anderson and Goodman,
#' 1957; Billingsley, 1961). The Wald version is the heterogeneity statistic
#' used to compare independent estimates of a common quantity (Cochran, 1954),
#' applied to the estimates of \eqn{P_{hj\ell}} for the different \eqn{h}. Both
#' test the order at the level of each transition or state, and are therefore
#' more specific than global tests of the Markov property in continuous-time
#' multistate models (Titman and Putter, 2022).
#'
#' @references Anderson, T. W. and Goodman, L. A. (1957). Statistical inference
#'   about Markov chains. \emph{Annals of Mathematical Statistics}, 28(1), 89-110.
#'
#' Billingsley, P. (1961). Statistical methods in Markov chains. \emph{Annals
#' of Mathematical Statistics}, 32(1), 12-40.
#'
#' Cochran, W. G. (1954). The combination of estimates from different
#' experiments. \emph{Biometrics}, 10(1), 101-129.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65-70.
#'
#' Titman, A. C. and Putter, H. (2022). General tests of the Markov property
#' in multi-state models. \emph{Biostatistics}, 23(2), 380-396.
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
#' markov_test(fit)                                # Wald, one test per transition
#' markov_test(fit, method = "lrt")                # likelihood ratio
#' markov_test(fit, method = "lrt", by = "state")  # joint test per current state
#' @export
markov_test <- function(object, min.h = 2, method = c("wald", "lrt"),
                        by = c("transition", "state"), p.adjust = "holm") {
  stopifnot(inherits(object, "P2est"))
  method <- match.arg(method)
  by     <- match.arg(by)
  if (min.h < 2) stop("`min.h` must be >= 2: a heterogeneity test needs at ",
                      "least two preceding-state groups.", call. = FALSE)
  if (by == "state" && method != "lrt")
    stop("`by = \"state\"` is only available with `method = \"lrt\"`.", call. = FALSE)

  est <- data.table::as.data.table(object$estimate)
  est[, `:=`(h = as.character(h), j = as.character(j), l = as.character(l))]

  out <- if (method == "wald") .mt_wald(est, min.h)
         else if (by == "transition") .mt_lrt_transition(est, object$absorbing, min.h)
         else .mt_lrt_state(est, object$absorbing, min.h)

  if (is.null(out) || nrow(out) == 0L) {
    ## Nothing testable (e.g. every destination is reached from a single
    ## history): return a well-formed, zero-row result rather than NULL, so
    ## callers can always rely on the documented columns.
    warning("No (j, l) transition had >= min.h preceding states with a ",
            "usable estimate; nothing to test.", call. = FALSE)
    out <- if (by == "state")
      data.table::data.table(j = character(), K = integer(), L = integer(),
                             df = integer(), statistic = numeric(), p.value = numeric(),
                             p.adj = numeric())
    else
      data.table::data.table(j = character(), l = character(), K = integer(),
                             df = integer(), statistic = numeric(),
                             p.value = numeric(), pooled = numeric(), p.adj = numeric())
  }
  ## Several transitions (or states) are tested at once: adjusted p-values
  ## control the family-wise error rate (Holm, 1979, by default).
  out[, p.adj := stats::p.adjust(p.value, method = p.adjust)]
  data.table::setorder(out, p.value)   # most significant (smallest p) first
  structure(as.data.frame(out), class = c("markov_test", "data.frame"),
            method = method, by = by, p.adjust = p.adjust,
            estimator = object$estimator, conf.level = object$conf.level)
}

## --- internal: Wald / Cochran's Q, one test per (j, l) -----------------------
.mt_wald <- function(est, min.h) {
  ## Runs the heterogeneity test on one (j, l) transition's slice of `est` --
  ## one row per preceding state h observed for this (j, l).
  run_one <- function(sub) {
    sub <- sub[se > 0]                          # drop degenerate (se = 0) groups
    K <- nrow(sub)
    if (K < min.h) return(NULL)                 # omitted from the grouped result
    ## Weight each history's estimate by its precision, pool them into a
    ## single history-blind estimate pbar, and measure how far each history
    ## sits from it in weighted-squared-distance units.
    w    <- 1 / sub$se^2
    pbar <- sum(w * sub$p) / sum(w)
    Q    <- sum(w * (sub$p - pbar)^2)
    df   <- K - 1L
    data.table::data.table(K = K, df = df, statistic = Q,
                           p.value = stats::pchisq(Q, df, lower.tail = FALSE),
                           pooled = pbar)
  }
  est[, run_one(.SD), by = .(j, l)]
}

## --- internal: the at-risk table of every (h, j) pair ------------------------
## One row per (h, j) with Y_h = sum_s Y_hj(s-1); taken from `estimate`, where
## at.risk is repeated on every (h, j, l) row of the pair. Absorbing current
## states are left out (nothing can be tested about leaving them).
.mt_risk <- function(est, absorbing) {
  unique(est[!j %in% absorbing, .(h, j, Y = at.risk)])
}

## x log(x / y) with the convention 0 log 0 = 0.
.xlogx <- function(x, y) ifelse(x > 0, x * log(x / y), 0)

## --- internal: likelihood ratio, one test per (j, l) -------------------------
.mt_lrt_transition <- function(est, absorbing, min.h) {
  risk <- .mt_risk(est, absorbing)
  pairs <- unique(est[!j %in% absorbing, .(j, l)])
  rows <- lapply(seq_len(nrow(pairs)), function(r) {
    jj <- pairs$j[r]; ll <- pairs$l[r]
    rk <- risk[j == jj]                                   # every history of j
    K  <- nrow(rk)
    if (K < min.h) return(NULL)
    ## events towards l for each history (0 when that history never moved to l)
    nt <- est[j == jj & l == ll, .(h, N = n.trans)]
    N  <- nt$N[match(rk$h, nt$h)]; N[is.na(N)] <- 0
    Y  <- rk$Y
    pbar <- sum(N) / sum(Y)
    if (pbar <= 0 || pbar >= 1) return(NULL)             # no information
    G2 <- 2 * sum(.xlogx(N, Y * pbar) + .xlogx(Y - N, Y * (1 - pbar)))
    df <- K - 1L
    data.table::data.table(j = jj, l = ll, K = K, df = df, statistic = G2,
                           p.value = stats::pchisq(G2, df, lower.tail = FALSE),
                           pooled = pbar)
  })
  data.table::rbindlist(rows)
}

## --- internal: likelihood ratio, one joint test per current state j ----------
.mt_lrt_state <- function(est, absorbing, min.h) {
  risk <- .mt_risk(est, absorbing)
  rows <- lapply(unique(risk$j), function(jj) {
    rk <- risk[j == jj]; K <- nrow(rk)
    if (K < min.h) return(NULL)
    dests <- sort(unique(est[j == jj, l]))
    L <- length(dests)
    if (L < 2L) return(NULL)                              # nothing to compare
    ## K x L matrix of event counts; rows sum to the at-risk counts Y_h
    Nm <- matrix(0, K, L, dimnames = list(rk$h, dests))
    sub <- est[j == jj]
    Nm[cbind(match(sub$h, rk$h), match(sub$l, dests))] <- sub$n.trans
    Y    <- rk$Y
    pbar <- colSums(Nm) / sum(Y)                          # first-order estimate
    G2 <- 2 * sum(.xlogx(Nm, outer(Y, pbar)))
    df <- (K - 1L) * (L - 1L)
    data.table::data.table(j = jj, K = K, L = L, df = df, statistic = G2,
                           p.value = stats::pchisq(G2, df, lower.tail = FALSE))
  })
  data.table::rbindlist(rows)
}

## Console representation: header (which statistic, which hypothesis) followed
## by the results table, with the numeric columns rounded for readability
## (the underlying object itself keeps full precision).
#' @export
print.markov_test <- function(x, digits = 4, ...) {
  method <- attr(x, "method") %||% "wald"
  by     <- attr(x, "by") %||% "transition"
  cat(sprintf("<markov_test>  %s test of the first-order Markov assumption (%s-based)\n",
              if (method == "lrt") "Likelihood-ratio" else "Wald (Cochran Q)",
              attr(x, "estimator")))
  if (by == "state")
    cat("  H0 (per current state j): P(X_s = . | X_{s-1} = j) does not depend on X_{s-2} = h\n")
  else
    cat("  H0 (per transition): P(X_s = l | X_{s-1} = j) does not depend on X_{s-2} = h\n")
  cat(sprintf("  %d test(s); sorted by ascending p-value\n\n", nrow(x)))
  if (nrow(x)) {
    y <- x
    ## intersect(): only round the columns actually present, so printing a
    ## column-subset of a markov_test object doesn't error on a missing one.
    for (nm in intersect(c("statistic", "p.value", "pooled", "p.adj"), names(y)))
      y[[nm]] <- round(y[[nm]], digits)
    print(as.data.frame(y), row.names = FALSE)
  }
  invisible(x)
}
