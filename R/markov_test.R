#' Test the first-order Markov assumption, transition by transition
#'
#' For each current state \eqn{j} and target \eqn{\ell}, tests whether the 1-step second-order transition probability \eqn{P_{hj\ell}} depends on the state at the previous time \eqn{h}. Under the null hypothesis of first-order Markov behaviour, \eqn{P_{hj\ell}} is the same for every observed \eqn{h}.
#'
#' \code{method = "wald"} (default) is an inverse-variance weighted heterogeneity statistic (Cochran's Q) built on the RPE estimates and standard errors of \code{\link{P2est}}: with \eqn{w_h = 1/\mathrm{se}_{hj\ell}^2} and the pooled estimate \eqn{\bar P_{j\ell} = \sum_h w_h P_{hj\ell} / \sum_h w_h},
#' \deqn{Q_{j\ell} = \sum_h w_h (P_{hj\ell} - \bar P_{j\ell})^2,}
#' asymptotically \eqn{\chi^2_{K-1}} under the null, where \eqn{K} is the number of previous states with a non-degenerate estimate (se > 0).
#'
#' \code{method = "lrt"} is the likelihood-ratio statistic of the classical test of Markov order (Anderson and Goodman, 1957), computed from the RPE counts: with \eqn{N_h = \sum_s \tilde N_{hj\ell}(s)}, \eqn{Y_h = \sum_s \tilde Y_{hj}(s-1)} and the first-order estimate \eqn{\bar P_{j\ell} = \sum_h N_h / \sum_h Y_h},
#' \deqn{G^2_{j\ell} = 2 \sum_h \left[ N_h \log\frac{N_h}{Y_h \bar P_{j\ell}} + (Y_h - N_h) \log\frac{Y_h - N_h}{Y_h (1 - \bar P_{j\ell})} \right],}
#' asymptotically \eqn{\chi^2_{K-1}}, where \eqn{K} is now the number of previous states with subjects at risk in \eqn{j} (histories with no move towards \eqn{\ell} count too). With \code{by = "state"} the destinations are tested jointly for each current state,
#' \deqn{G^2_j = 2 \sum_h \sum_\ell N_{h\ell} \log\frac{N_{h\ell}}{Y_h \bar P_{j\ell}},}
#' with \eqn{(K - 1)(L - 1)} degrees of freedom (\eqn{L} = number of observed destinations from \eqn{j}): the Anderson-Goodman test of a first- against a second-order chain, state by state.
#'
#' Both versions rely on the likelihood of the second-order chain factorising over the \eqn{(h, j)} pairs, which makes the estimates for different previous states asymptotically independent even when the same subjects contribute to several of them.
#'
#' @param object A "P2est" object.
#' @param min.h Minimum number of previous states needed to test a transition (or a state). Default 2.
#' @param method \code{"wald"} (default) or \code{"lrt"}; see Details.
#' @param by \code{"transition"} (default): one test per \eqn{(j, \ell)}; \code{"state"}: one joint test per current state \eqn{j} (requires \code{method = "lrt"}).
#' @param p.adjust Method for the multiplicity-adjusted p-values (column \code{p.adj}), passed to \code{\link[stats]{p.adjust}}. Default \code{"holm"}, which controls the family-wise error rate without assuming independent tests; \code{"none"} keeps the unadjusted p-values.
#' @return A data frame of class "markov_test", sorted by ascending p-value. With \code{by = "transition"}, one row per tested \eqn{(j, \ell)} with columns \code{j, l, K, df, statistic, p.value, pooled, p.adj} (\code{pooled} is the history-blind estimate under the null: inverse-variance weighted for \code{"wald"}, the first-order RPE for \code{"lrt"}). With \code{by = "state"}, one row per current state with columns \code{j, K, L, df, statistic, p.value, p.adj}. Attributes \code{method}, \code{by}, \code{p.adjust}, \code{estimator} and \code{conf.level}.
#' @section Why these tests:
#' The likelihood-ratio test is the classical test of the order of a Markov chain observed on many subjects (Anderson and Goodman, 1957; Billingsley, 1961). The Wald version is the heterogeneity statistic used to compare independent estimates of a common quantity (Cochran, 1954). Both test the order at the level of each transition or state, and are therefore more specific than global tests of the Markov property in continuous-time multistate models (Titman and Putter, 2022).
#' @seealso \code{\link{P1est}}, \code{\link{compare_order}}
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about Markov chains. \emph{Annals of Mathematical Statistics}, 28(1), 89-110.
#'
#' Billingsley, P. (1961). Statistical methods in Markov chains. \emph{Annals of Mathematical Statistics}, 32(1), 12-40.
#'
#' Cochran, W. G. (1954). The combination of estimates from different experiments. \emph{Biometrics}, 10(1), 101-129.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics}, 6(2), 65-70.
#'
#' Titman, A. C. and Putter, H. (2022). General tests of the Markov property in multi-state models. \emph{Biostatistics}, 23(2), 380-396.
#' @examples
#' st <- c("A", "B", "C")   # C is absorbing
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
#' markov_test(fit)                                 # Wald, one test per transition
#' markov_test(fit, method = "lrt")                 # likelihood ratio
#' markov_test(fit, method = "lrt", by = "state")   # joint test per current state
#' @export
markov_test <- function(object, min.h = 2, method = c("wald", "lrt"), by = c("transition", "state"), p.adjust = "holm") {

  # Check the arguments.
  stopifnot(inherits(object, "P2est"))
  method <- match.arg(method)
  by <- match.arg(by)
  if (min.h < 2)
    stop("`min.h` must be >= 2: a heterogeneity test needs at least two previous states.", call. = FALSE)
  if (by == "state" && method != "lrt")
    stop("`by = \"state\"` is only available with `method = \"lrt\"`.", call. = FALSE)

  # The estimates with the states as text, one row per (h, j, l).
  est <- object$estimate |>
    dplyr::mutate(dplyr::across(c("h", "j", "l"), as.character))

  # One test per transition (Wald or LRT) or per current state (LRT).
  out <- if (method == "wald") .mt_wald(est, min.h)
         else if (by == "transition") .mt_lrt_transition(est, object$absorbing, min.h)
         else .mt_lrt_state(est, object$absorbing, min.h)

  # Nothing testable (e.g. every destination reached from a single history): a well-formed, zero-row result with the documented columns.
  if (!nrow(out)) {
    warning("No transition had >= min.h previous states with a usable estimate; nothing to test.", call. = FALSE)
  }

  # Several transitions (or states) are tested at once: adjusted p-values control the family-wise error rate; most significant first.
  out <- out |>
    dplyr::mutate(p.adj = stats::p.adjust(p.value, method = p.adjust)) |>
    dplyr::arrange(p.value)

  structure(as.data.frame(out), class = c("markov_test", "data.frame"),
            method = method, by = by, p.adjust = p.adjust,
            estimator = object$estimator, conf.level = object$conf.level)
}

# internal: Wald (Cochran's Q), one test per (j, l) -> weight each history's estimate by its precision, pool them, and measure the weighted squared distance of each history to the pooled estimate. Histories with se = 0 carry no information and are dropped.
.mt_wald <- function(est, min.h) {
  est |>
    dplyr::filter(se > 0) |>
    dplyr::mutate(w = 1 / se^2) |>
    dplyr::summarise(K = dplyr::n(),
                     pooled = sum(w * p) / sum(w),
                     statistic = sum(w * (p - pooled)^2),
                     .by = c("j", "l")) |>
    dplyr::filter(K >= min.h) |>
    dplyr::mutate(df = K - 1L,
                  p.value = stats::pchisq(statistic, df, lower.tail = FALSE)) |>
    dplyr::select("j", "l", "K", "df", "statistic", "p.value", "pooled")
}

# internal: at-risk count of every (h, j) pair, Y_h = sum_s Y_hj(s-1) (repeated on every (h, j, l) row of the estimates). Absorbing current states are left out: nothing can be tested about leaving them.
.mt_risk <- function(est, absorbing) {
  est |>
    dplyr::filter(!j %in% absorbing) |>
    dplyr::distinct(h, j, Yh = at.risk)
}

# internal: x log(x / y), with the convention 0 log 0 = 0.
.xlogx <- function(x, y) ifelse(x > 0, x * log(x / y), 0)

# internal: likelihood ratio, one test per (j, l) -> for every history of j, the moves towards l (0 for a history that never moved there) and the subjects at risk; the first-order estimate pools them.
.mt_lrt_transition <- function(est, absorbing, min.h) {
  risk <- .mt_risk(est, absorbing)
  targets <- est |>
    dplyr::filter(!j %in% absorbing) |>
    dplyr::distinct(j, l)
  risk |>
    dplyr::inner_join(targets, by = "j", relationship = "many-to-many") |>
    dplyr::left_join(dplyr::select(est, "h", "j", "l", Nh = "n.trans"), by = c("h", "j", "l")) |>
    dplyr::mutate(Nh = dplyr::coalesce(as.numeric(Nh), 0)) |>
    dplyr::summarise(K = dplyr::n(),
                     pooled = sum(Nh) / sum(Yh),
                     statistic = 2 * sum(.xlogx(Nh, Yh * pooled) + .xlogx(Yh - Nh, Yh * (1 - pooled))),
                     .by = c("j", "l")) |>
    dplyr::filter(K >= min.h, pooled > 0, pooled < 1) |>       # a pooled 0 or 1 carries no information
    dplyr::mutate(df = K - 1L,
                  p.value = stats::pchisq(statistic, df, lower.tail = FALSE)) |>
    dplyr::select("j", "l", "K", "df", "statistic", "p.value", "pooled")
}

# internal: likelihood ratio, one joint test per current state j -> the K x L table of moves from every history to every destination, against the first-order (pooled) distribution of the next state.
.mt_lrt_state <- function(est, absorbing, min.h) {
  risk <- .mt_risk(est, absorbing)
  at_risk_j <- dplyr::summarise(risk, Yj = sum(Yh), K = dplyr::n(), .by = "j")
  cells <- est |>
    dplyr::filter(!j %in% absorbing) |>
    dplyr::select("h", "j", "l", Nh = "n.trans") |>
    dplyr::left_join(risk, by = c("h", "j")) |>
    dplyr::left_join(at_risk_j, by = "j") |>
    dplyr::mutate(pbar = sum(Nh) / Yj, .by = c("j", "l"))   # first-order estimate of j -> l
  cells |>
    dplyr::summarise(K = dplyr::first(K), L = dplyr::n_distinct(l),
                     statistic = 2 * sum(.xlogx(Nh, Yh * pbar)),
                     .by = "j") |>
    dplyr::filter(K >= min.h, L >= 2L) |>
    dplyr::mutate(df = (K - 1L) * (L - 1L),
                  p.value = stats::pchisq(statistic, df, lower.tail = FALSE)) |>
    dplyr::select("j", "K", "L", "df", "statistic", "p.value")
}

# Printed form of a test: which statistic and which hypothesis, then the table with the numbers rounded (the object keeps full precision).
#' @export
print.markov_test <- function(x, digits = 4, ...) {
  method <- attr(x, "method") %||% "wald"
  by <- attr(x, "by") %||% "transition"
  cat(sprintf("<markov_test>  %s test of the first-order Markov assumption (%s-based)\n",
              if (method == "lrt") "Likelihood-ratio" else "Wald (Cochran Q)", attr(x, "estimator")))
  if (by == "state")
    cat("  H0 (per current state j): P(X_s = . | X_{s-1} = j) does not depend on X_{s-2} = h\n")
  else
    cat("  H0 (per transition): P(X_s = l | X_{s-1} = j) does not depend on X_{s-2} = h\n")
  cat(sprintf("  %d test(s); sorted by ascending p-value\n\n", nrow(x)))
  if (nrow(x)) {
    y <- as.data.frame(x)
    for (nm in intersect(c("statistic", "p.value", "pooled", "p.adj"), names(y)))
      y[[nm]] <- round(y[[nm]], digits)
    print(y, row.names = FALSE)
  }
  invisible(x)
}
