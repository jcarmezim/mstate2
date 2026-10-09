#' First-order transition probabilities from the same data
#'
#' Estimates the first-order discrete-time transition probabilities \eqn{P_{j\ell} = P(X_s = \ell \mid X_{s-1} = j)} from the \emph{same} second-order counting processes used by \code{\link{P2est}}, by pooling over the state at the previous time:
#' \deqn{\hat P^{(1)}_{j\ell} = \frac{\sum_h \sum_s \tilde N_{hj\ell}(s)}{\sum_h \sum_s \tilde Y_{hj}(s-1)}.}
#' This is the RPE of a first-order homogeneous chain, the estimate a second-order analysis reduces to when the first-order Markov property holds. Because it uses exactly the same subjects, instants and time grid as the second-order fit, a comparison between the two isolates the effect of conditioning on the previous state: it is the natural first-order baseline of \code{\link{compare_order}} and \code{\link{divergence_table}}.
#'
#' Standard errors are binomial, \eqn{\sqrt{\hat p(1-\hat p)/Y_j}}, with Wald intervals clipped to \eqn{[0, 1]}. With \code{B > 0} the matrix is also bootstrapped by resampling subjects; the replicates give the standard errors of the \eqn{n}-step predictions in \code{\link{probtrans2}} and \code{\link{compare_order}}. Given a \code{\link{P2boot}} fit, its first-order replicates (from the same resamples as the second-order ones) are reused, which allows paired comparisons.
#'
#' @param object An "msm2data" object from \code{\link{prep2}}, or a "P2boot" object (whose first-order replicates are then reused).
#' @param B Number of bootstrap replicates (ignored for a "P2boot" input). Default 200; 0 skips the bootstrap.
#' @param conf.level Confidence level. Default 0.95.
#' @param seed Optional integer seed for the bootstrap.
#' @return An object of class "P1est": a list with \code{estimate} (tibble \code{j, l, p, se, lower, upper, n.trans, at.risk}), the \eqn{M \times M} matrices \code{P}, \code{P.lower}, \code{P.upper} and \code{P.se} (rows = current state, columns = next state), the bootstrap array \code{boot} (\eqn{M \times M \times B}, or \code{NULL}), \code{states}, \code{absorbing}, \code{conf.level} and \code{n}.
#' @section Why this first-order model:
#' Under the first-order hypothesis the transition probabilities do not depend on the state at the previous time, and their maximum likelihood estimator pools the second-order counts over that state (Anderson and Goodman, 1957). It is therefore the null model of \code{\link{markov_test}} and, being estimated from the same patients, days and time scale, the like-for-like first-order comparator: any difference with the second-order model is due to the order alone, not to a different time scale or estimator.
#' @seealso \code{\link{P2est}}, \code{\link{probtrans2}}, \code{\link{compare_order}}, \code{\link{markov_test}}
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about Markov chains. \emph{Annals of Mathematical Statistics}, 28(1), 89-110.
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
#' p1 <- P1est(prep2(panel), B = 50, seed = 1)
#' p1
#' p1$P   # B -> B pools the 0.6 (from A) and 0.3 (from B) histories
#' @export
P1est <- function(object, B = 200, conf.level = 0.95, seed = NULL) {

  # Check the arguments and get the normal quantile z_{1 - alpha/2}.
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # Pooled counts. From a P2boot fit: its estimate table (one row per (h, j, l), with the at-risk count of its (h, j) pair) and its first-order replicates. From an msm2data object: its counting processes, and a new bootstrap if B > 0.
  if (inherits(object, "P2boot")) {
    boot <- object$boot1
    N <- object$estimate |>
      dplyr::summarise(n.trans = sum(n.trans), .by = c("j", "l"))
    Y <- object$estimate |>
      dplyr::distinct(h, j, at.risk) |>
      dplyr::summarise(at.risk = sum(at.risk), .by = "j")
  } else {
    stopifnot(inherits(object, "msm2data"))
    if (!is.numeric(B) || length(B) != 1L || B < 0)
      stop("`B` must be a single number >= 0.", call. = FALSE)
    N <- object$N |>
      dplyr::summarise(n.trans = sum(N), .by = c("j", "l"))
    Y <- object$Y |>
      dplyr::summarise(at.risk = sum(Y), .by = "j")
    boot <- NULL
    if (B > 0) {
      # Bootstrap of subjects, as in P2boot(): per-subject counts of every triple, reweighted by how many times each subject is drawn, then pooled over h.
      if (!is.null(seed)) set.seed(seed)
      ct <- .id_counts(object)
      n <- nrow(ct$mat)
      P1_point <- .first_from_counts(colSums(ct$mat), ct, object, NULL)
      boot <- purrr::map(seq_len(B), \(b) {
        w <- tabulate(sample.int(n, n, replace = TRUE), n)
        .first_from_counts(as.vector(w %*% ct$mat), ct, object, P1_point)
      }) |>
        simplify2array()   # M x M x B array
    }
  }

  # Estimates: p = n.trans / at.risk, binomial standard error and Wald interval clipped to [0, 1].
  est <- dplyr::left_join(N, Y, by = "j") |>
    dplyr::mutate(p = n.trans / at.risk,
                  se = sqrt(p * (1 - p) / at.risk),
                  lower = pmax(0, p - z * se),
                  upper = pmin(1, p + z * se)) |>
    dplyr::arrange(j, l) |>
    dplyr::select("j", "l", "p", "se", "lower", "upper", "n.trans", "at.risk")

  # Matrices P[j, l]: one row per current state. Absorbing states keep probability 1.
  states <- object$states
  M <- length(states)
  key <- cbind(match(as.character(est$j), states), match(as.character(est$l), states))
  blank <- matrix(0, M, M, dimnames = list(states, states))
  P <- Pl <- Pu <- Pse <- blank
  P[key] <- est$p
  Pl[key] <- est$lower
  Pu[key] <- est$upper
  Pse[key] <- est$se
  for (a in match(object$absorbing, states)) {
    P[a, ] <- Pl[a, ] <- Pu[a, ] <- 0
    P[a, a] <- Pl[a, a] <- Pu[a, a] <- 1
  }

  structure(list(estimate = est, P = P, P.lower = Pl, P.upper = Pu, P.se = Pse, boot = boot,
                 states = states, absorbing = object$absorbing,
                 conf.level = conf.level, n = object$n),
            class = "P1est")
}

# Short description of a first-order fit, shown when it is printed.
#' @export
print.P1est <- function(x, ...) {
  cat("<P1est>  first-order (history-pooled) RPE transition probabilities\n")
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %d estimated transition probabilities (j -> l); %d subjects; %s\n",
              nrow(x$estimate), x$n,
              if (is.null(x$boot)) "no bootstrap" else sprintf("%d bootstrap replicates", dim(x$boot)[3L])))
  invisible(x)
}
