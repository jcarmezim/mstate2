#' First-order (history-pooled) transition probabilities from the same data
#'
#' Estimates the first-order discrete-time transition probabilities
#' \eqn{P_{j\ell} = P(X_s = \ell \mid X_{s-1} = j)} from the \emph{same}
#' second-order counting processes used by \code{\link{P2est}}, by pooling
#' over the preceding state:
#' \deqn{\hat P^{(1)}_{j\ell} = \frac{\sum_h \sum_s \tilde N_{hj\ell}(s)}{\sum_h \sum_s \tilde Y_{hj}(s-1)}.}
#' This is the RPE of a first-order homogeneous chain, and the estimate a
#' second-order analysis collapses to when the first-order Markov property
#' holds (it is the \code{pooled} column of \code{markov_test(method = "lrt")}).
#' Because it uses exactly the same subjects, instants and time grid as the
#' second-order fit, a comparison between the two isolates the effect of
#' conditioning on the preceding state: it is the natural first-order baseline
#' for \code{\link{compare_order}} and \code{\link{divergence_table}} when no
#' (continuous-time) \code{mstate} fit is wanted.
#'
#' Standard errors are binomial, \eqn{\sqrt{\hat p(1-\hat p)/Y_j}}. With
#' \code{B > 0}, the matrix is also bootstrapped by resampling subjects, and
#' the replicates are used by \code{\link{probtrans2}} (and therefore by
#' \code{compare_order()}) for the standard errors of the \eqn{n}-step
#' predictions.
#'
#' @param object An "msm2data" object, or a "P2boot" object (whose first-order
#'   replicates are then reused).
#' @param B Number of bootstrap replicates (ignored for a "P2boot" input).
#'   Default 200; 0 skips the bootstrap.
#' @param conf.level Confidence level. Default 0.95.
#' @param seed Optional integer seed for the bootstrap.
#' @return An object of class "P1est": a list with \code{estimate} (data frame
#'   \code{j, l, p, se, lower, upper, n.trans, at.risk}), the \eqn{M \times M}
#'   matrices \code{P}, \code{P.lower}, \code{P.upper}, \code{P.se}, the
#'   bootstrap array \code{boot} (\eqn{M \times M \times B}, or \code{NULL}),
#'   \code{states}, \code{absorbing}, \code{conf.level} and \code{n}.
#' @seealso \code{\link{probtrans2}}, \code{\link{compare_order}},
#'   \code{\link{divergence_table}}
#' @section Why this first-order model:
#' Under the first-order hypothesis the transition probabilities do not depend
#' on the state at the previous time, and their maximum likelihood estimator
#' is obtained by pooling the second-order counts over that state (Anderson
#' and Goodman, 1957). It is therefore the null model of
#' \code{\link{markov_test}} and, being estimated from exactly the same
#' patients, days and time scale, the like-for-like first-order comparator of
#' \code{\link{compare_order}}: any difference with the second-order model is
#' due to the order alone, not to a different time scale (discrete versus
#' continuous) or estimator.
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. \emph{Annals of Mathematical Statistics}, 28(1), 89-110.
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
#' p1 <- P1est(prep2(panel), B = 50, seed = 1)
#' p1$P            # B -> B pools the 0.6 (from A) and 0.3 (from B) histories
#' @export
P1est <- function(object, B = 200, conf.level = 0.95, seed = NULL) {
  .check_conf_level(conf.level)
  if (inherits(object, "P2boot")) {
    boot <- object$boot1
    states <- object$states; absorbing <- object$absorbing; n <- object$n
    est <- data.table::as.data.table(object$estimate)
  } else {
    stopifnot(inherits(object, "msm2data"))
    states <- object$states; absorbing <- object$absorbing; n <- object$n
    est <- NULL
    boot <- NULL
    if (B > 0) {
      if (!is.null(seed)) set.seed(seed)
      ct <- .id_counts(object)
      nn <- nrow(ct$mat)
      boot <- array(0, c(length(states), length(states), B),
                    dimnames = list(states, states, NULL))
      for (b in seq_len(B)) {
        w <- tabulate(sample.int(nn, nn, replace = TRUE), nn)
        boot[, , b] <- .first_from_counts(as.vector(w %*% ct$mat), ct, object)
      }
    }
  }
  ## pooled counts: N_jl = sum_h n.trans, Y_j = sum_h at.risk over the pairs (h, j)
  if (is.null(est)) {
    N <- object$N[, .(n.trans = sum(N)), keyby = .(j, l)]
    Y <- object$Y[, .(at.risk = sum(Y)), keyby = .(j)]
  } else {
    N <- est[, .(n.trans = sum(n.trans)), keyby = .(j, l)]
    Y <- unique(est[, .(h, j, at.risk)])[, .(at.risk = sum(at.risk)), keyby = .(j)]
  }
  tab <- merge(N, Y, by = "j")
  tab[, p := n.trans / at.risk]
  tab[, se := sqrt(p * (1 - p) / at.risk)]
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  tab[, lower := pmax(0, p - z * se)]
  tab[, upper := pmin(1, p + z * se)]

  M <- length(states)
  blank <- matrix(0, M, M, dimnames = list(states, states))
  P <- Pl <- Pu <- Pse <- blank
  key <- cbind(match(as.character(tab$j), states), match(as.character(tab$l), states))
  P[key] <- tab$p; Pl[key] <- tab$lower; Pu[key] <- tab$upper; Pse[key] <- tab$se
  for (a in match(absorbing, states)) {
    P[a, ] <- Pl[a, ] <- Pu[a, ] <- 0
    P[a, a] <- Pl[a, a] <- Pu[a, a] <- 1
  }
  structure(list(estimate = as.data.frame(tab[, .(j, l, p, se, lower, upper, n.trans, at.risk)]),
                 P = P, P.lower = Pl, P.upper = Pu, P.se = Pse, boot = boot,
                 states = states, absorbing = absorbing,
                 conf.level = conf.level, n = n),
            class = "P1est")
}

#' @export
print.P1est <- function(x, ...) {
  cat("<P1est>  first-order (history-pooled) RPE transition probabilities\n")
  cat(sprintf("  states: %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  %d estimated transition probabilities (j -> l); %d subjects; %s\n",
              nrow(x$estimate), x$n,
              if (is.null(x$boot)) "no bootstrap"
              else sprintf("%d bootstrap replicates", dim(x$boot)[3L])))
  invisible(x)
}
