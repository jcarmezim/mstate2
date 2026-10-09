#' Discrete-time cause-specific hazard regression for 1-step second-order
#' transitions
#'
#' Fits the covariate-adjusted discrete-time cause-specific hazard model
#' \deqn{\mathrm{cloglog}(\lambda_{hj\ell}(Z)) = \gamma_{hj\ell} + \beta_{hj\ell}^\top Z}
#' for a fixed history \eqn{(h, j)} (state at the previous and at the current
#' time), one binary GLM per destination \eqn{\ell \neq j}. Here
#' \eqn{\lambda_{hj\ell}(Z) = P(X_s = \ell \mid X_{s-1} = j, X_{s-2} = h, Z)}
#' is the discrete-time cause-specific hazard of moving to \eqn{\ell} (Tutz and
#' Schmid, 2016, ch. 8); the probability of staying in \eqn{j} is
#' \eqn{1 - \sum_{\ell \neq j} \lambda_{hj\ell}}. \code{\link{P2est}} is
#' non-parametric and covariate-free; \code{P2reg} adds covariates and time
#' scales.
#'
#' \strong{Interpretation of \eqn{\exp(\beta)}.} With the (default)
#' complementary log-log link, \eqn{\exp(\beta)} is the ratio of
#' \eqn{-\log(1 - \lambda_{hj\ell})} between two covariate values. For a single
#' type of event this is exactly the hazard ratio of a continuous-time
#' proportional-hazards model observed in grouped time (Prentice and
#' Gloeckler, 1978). With several destinations it approximates the
#' continuous-time cause-specific hazard ratio when the daily probabilities
#' of moving are small, because then \eqn{-\log(1 - \lambda) \approx \lambda}.
#'
#' \strong{Why one binary model per destination.} Each indicator "moves to
#' \eqn{\ell}" is, given the history, a Bernoulli variable with probability
#' \eqn{\lambda_{hj\ell}(Z)}: the marginal of the multinomial transition. The
#' per-destination binomial likelihoods are therefore valid (marginal)
#' likelihoods and give consistent estimators, at the price of some
#' efficiency and of not forcing the fitted probabilities to sum to at most
#' one (a composite-likelihood argument; Varin, Reid and Firth, 2011). A
#' multinomial model would impose that constraint; it is a possible
#' extension. \code{\link{P2reg_all}} obtains staying as one minus the moves.
#'
#' \strong{Standard errors.} Given the second-order model, the days of a
#' patient are conditionally independent and the likelihood factorises, so
#' the model-based GLM standard errors are valid if the model is correct.
#' \code{cluster = TRUE} gives cluster-robust (sandwich) standard errors by
#' patient (Liang and Zeger, 1986; Zeileis, Koll and Graham, 2020), which
#' remain valid if days of the same patient are correlated beyond what the
#' model captures (e.g. unobserved frailty). It is the recommended choice for
#' reporting.
#'
#' The risk set for (h, j) is \code{object$triples} filtered to that pair: one
#' row per subject per at-risk instant, exactly \code{\link{prep2}}'s
#' \eqn{\tilde Y_{hj}(s-1)} at the individual level. All destinations observed
#' from that risk set share it, so a single data preparation step fits every
#' destination-specific model. With \code{formula = ~1} (no covariates), the
#' fitted probability for a given l is, by construction, the sample proportion
#' of the risk set that moved to l -- which is exactly the RPE computed
#' by \code{P2est()}, regardless of the link function
#' used. This identity is a useful internal check that \code{P2reg} and
#' \code{P2est} agree, and is verified in the package's test suite.
#'
#' @param object An "msm2data" object, built with \code{prep2(..., covariates = )}
#'   when the formula uses baseline covariates.
#' @param h,j The preceding and current state defining the risk set (label or index).
#' @param l Optional target state(s) (label or index). If NULL (default), a
#'   separate cause-specific model is fit for every move \eqn{\ell \neq j}
#'   observed from (h, j). \code{l = j} (staying) can be requested
#'   explicitly; its \eqn{\exp(\beta)} is not a hazard ratio.
#' @param formula A one-sided formula giving the covariate predictors, e.g.
#'   \code{~ age + sex}. Default \code{~1} (baseline hazard only, no covariates).
#'   Variables must be \code{covariates} supplied to \code{\link{prep2}}, or
#'   the two time scales stored in \code{object$triples}: \code{s} (time of
#'   the destination) and \code{d} (time already spent in the current state).
#'   They give a flexible baseline: \code{~ age + splines::ns(s, 3)} relaxes
#'   time homogeneity (\eqn{\gamma_{hj\ell}(s)}), \code{~ age + factor(s)}
#'   gives one parameter per time unit, and \code{~ age + log(d)} or
#'   \code{~ age + splines::ns(d, 3)} gives the exploratory semi-Markov
#'   extension \eqn{\gamma_{hj\ell}(d)} (duration dependence). A baseline
#'   that varies with time is the usual discrete-time hazard model (Allison,
#'   1982; Tutz and Schmid, 2016); dependence on the time already spent in the
#'   current state is the semi-Markov ("clock-reset") assumption of multistate
#'   models (Putter, Fiocco and Geskus, 2007).
#' @param family The GLM family/link, as a family object. Default
#'   \code{stats::binomial("cloglog")} (see Details); \code{exp(coef)} is
#'   labelled a hazard ratio only when \code{family$link == "cloglog"}.
#' @param cluster If TRUE, replace the model-based standard errors with a
#'   cluster-robust ("sandwich") estimate clustered by subject id, via
#'   \code{sandwich::vcovCL()} (requires the \pkg{sandwich} package). Relevant
#'   when a subject can contribute more than one at-risk instant to the same
#'   (h, j) risk set. Default FALSE (model-based standard errors).
#' @param ... Passed to \code{\link[stats]{glm}}.
#' @references
#' Allison, P. D. (1982). Discrete-time methods for the analysis of event
#' histories. \emph{Sociological Methodology}, 13, 61-98.
#'
#' Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, 73(1), 13-22.
#'
#' Prentice, R. L. and Gloeckler, L. A. (1978). Regression analysis of grouped
#' survival data with application to breast cancer data. \emph{Biometrics},
#' 34(1), 57-67.
#'
#' Putter, H., Fiocco, M. and Geskus, R. B. (2007). Tutorial in biostatistics:
#' competing risks and multi-state models. \emph{Statistics in Medicine},
#' 26(11), 2389-2430.
#'
#' Tutz, G. and Schmid, M. (2016). \emph{Modeling Discrete Time-to-Event Data}.
#' Springer.
#'
#' Varin, C., Reid, N. and Firth, D. (2011). An overview of composite
#' likelihood methods. \emph{Statistica Sinica}, 21(1), 5-42.
#'
#' Zeileis, A., Koll, S. and Graham, N. (2020). Various versatile variances:
#' an object-oriented implementation of clustered covariances in R.
#' \emph{Journal of Statistical Software}, 95(1), 1-36.
#' @return An object of class "P2reg": a list with \code{models} (one
#'   \code{stats::glm} fit per destination state, named by state), \code{h},
#'   \code{j}, \code{states}, \code{formula}, \code{family}, \code{cluster},
#'   \code{n} (subjects at risk) and \code{n.risk} (at-risk instants).
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
#'
#' # attach a synthetic, subject-level covariate (plays no role in the
#' # simulation, so its estimated effect should be negligible)
#' ids <- unique(panel$id)
#' set.seed(2)
#' panel$age <- stats::setNames(sample(40:80, length(ids), TRUE),
#'                              ids)[as.character(panel$id)]
#'
#' d2  <- prep2(panel, covariates = "age", states = st)
#' reg <- P2reg(d2, h = "A", j = "B", formula = ~age)
#' summary(reg)
#' @export
P2reg <- function(object, h, j, l = NULL, formula = ~1,
                  family = stats::binomial("cloglog"), cluster = FALSE, ...) {

  # Check the arguments.
  stopifnot(inherits(object, "msm2data"))
  if (cluster && !requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for cluster = TRUE.", call. = FALSE)
  states <- object$states
  hi <- .resolve(h, states)
  ji <- .resolve(j, states)
  if (anyNA(c(hi, ji))) stop("`h` or `j` not found in the state space.", call. = FALSE)
  hlab <- states[hi]
  jlab <- states[ji]

  # The risk set of (h, j): one row per subject and instant in j after h, i.e. prep2()'s Y_hj(s-1) at the individual level, with the covariates of each row.
  risk <- object$triples |>
    dplyr::filter(as.character(h) == hlab, as.character(j) == jlab)
  if (!nrow(risk))
    stop("No observed (h, j) = (", hlab, ", ", jlab, ") triples in `object`.", call. = FALSE)

  # Destinations: by default the moves l != j observed from the risk set (a hazard is the probability of leaving towards l; staying is its complement). Staying can be requested with l = j.
  targets <- if (is.null(l)) setdiff(sort(unique(as.character(risk$l))), jlab) else {
    li <- .resolve(l, states)
    if (anyNA(li)) stop("`l` not found in the state space.", call. = FALSE)
    states[li]
  }
  if (!length(targets))
    stop("No move out of '", jlab, "' is observed from (h, j) = (", hlab, ", ", jlab,
         "); request destinations explicitly with `l`.", call. = FALSE)

  # The formula: one-sided (the response is built internally for each l), with baseline covariates and the time scales s (time of the destination) and d (time already spent in j).
  if (length(formula) != 2L)
    stop("`formula` must be one-sided (e.g. ~ age + sex): the response is the indicator of moving to each destination, built internally.", call. = FALSE)
  rhs <- formula[[length(formula)]]
  miss <- setdiff(all.vars(rhs), c(object$covariates, "s", "d"))
  if (length(miss))
    stop(if (!length(object$covariates)) "`object` carries no covariates; build it with prep2(..., covariates = ). " else "",
         "`formula` uses variable(s) that are not covariates of `object` (or s, d): ",
         paste(miss, collapse = ", "), call. = FALSE)

  # One binary GLM per destination: the response is 1 for the rows that moved to l and 0 for all the others in the same risk set (staying or moving elsewhere), the usual cause-specific hazard construction. The formula is evaluated in its own environment, so that functions such as splines::ns() are found. With cluster = TRUE, the cluster-robust covariance by subject is stored with the fit.
  fm <- stats::as.formula(paste0(".y ~ ", deparse1(rhs)), env = environment(formula))
  fit_one <- function(tgt) {
    dat <- dplyr::mutate(risk, .y = as.integer(as.character(l) == tgt))
    fit <- stats::glm(fm, data = dat, family = family, ...)
    if (cluster) fit$P2reg.vcov <- sandwich::vcovCL(fit, cluster = dat$id)
    fit
  }
  fits <- stats::setNames(purrr::map(targets, fit_one), targets)

  structure(list(models = fits, h = hlab, j = jlab, states = states,
                 formula = formula, family = family, cluster = cluster,
                 n = dplyr::n_distinct(risk$id), n.risk = nrow(risk)),
            class = "P2reg")
}

#' Predict covariate-adjusted 1-step second-order transition probabilities
#'
#' @param object A "P2reg" object.
#' @param newdata A data frame with one row per covariate profile, columns
#'   matching the variables in \code{object$formula}.
#' @param type Passed to \code{stats::predict.glm}; \code{"response"} (default)
#'   gives probabilities.
#' @param ... Passed to \code{stats::predict.glm}.
#' @return A numeric matrix, one row per row of \code{newdata} and one column
#'   per destination state (named), giving \eqn{P_{hj\ell}(Z)} for each l.
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
#' ids <- unique(panel$id)
#' set.seed(2)
#' panel$age <- stats::setNames(sample(40:80, length(ids), TRUE),
#'                              ids)[as.character(panel$id)]
#'
#' d2  <- prep2(panel, covariates = "age", states = st)
#' reg <- P2reg(d2, h = "A", j = "B", formula = ~age)
#' predict(reg, newdata = data.frame(age = c(40, 60, 80)))
#' @export
predict.P2reg <- function(object, newdata, type = "response", ...) {

  # One column per destination model, evaluated on the same newdata (vapply and matrix() keep the destination names even with one row).
  out <- vapply(object$models, \(fit) unname(stats::predict(fit, newdata = newdata, type = type, ...)),
                numeric(nrow(newdata)))
  matrix(out, nrow = nrow(newdata),
         dimnames = list(rownames(newdata), names(object$models)))
}

# Short description of a fit: risk set, destinations, subjects and instants, covariates and type of standard errors.
#' @export
print.P2reg <- function(x, ...) {
  cat(sprintf("<P2reg>  %d destination-specific %s hazard model(s) from (h, j) = (%s, %s)\n",
              length(x$models), x$family$link, x$h, x$j))
  cat(sprintf("  destinations: %s\n", paste(names(x$models), collapse = ", ")))
  cat(sprintf("  %d subject(s), %d at-risk instant(s)\n", x$n, x$n.risk))
  cat(sprintf("  covariates   : %s\n", deparse1(x$formula[[length(x$formula)]])))
  cat(sprintf("  std. errors  : %s\n", if (x$cluster) "cluster-robust (subject id)" else "model-based"))
  invisible(x)
}

# Coefficient table of every destination model: estimate, standard error (cluster-robust if requested), exponentiated coefficient with its Wald interval and p-value. exp(coef) is labelled a hazard ratio only under the cloglog link.
#' @export
summary.P2reg <- function(object, conf.level = 0.95, ...) {
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  hr_ok <- identical(object$family$link, "cloglog")

  one_model <- function(fit, tgt) {
    V <- if (!is.null(fit$P2reg.vcov)) fit$P2reg.vcov else stats::vcov(fit)
    b <- stats::coef(fit)
    se <- sqrt(diag(V))[names(b)]   # terms not estimable (NA coefficient) have no variance: NA
    tibble::tibble(l = tgt, term = names(b), estimate = unname(b), se = unname(se)) |>
      dplyr::mutate(exp.coef = exp(estimate),
                    exp.lower = exp(estimate - z * se),
                    exp.upper = exp(estimate + z * se),
                    p.value = 2 * stats::pnorm(-abs(estimate / se)))
  }
  out <- purrr::imap(object$models, one_model) |>
    purrr::list_rbind() |>
    as.data.frame()
  if (hr_ok) names(out)[names(out) %in% c("exp.coef", "exp.lower", "exp.upper")] <- c("HR", "HR.lower", "HR.upper")

  structure(out, class = c("summary.P2reg", "data.frame"),
            h = object$h, j = object$j, cluster = object$cluster,
            link = object$family$link, conf.level = conf.level)
}

# Printed summary: risk set, link and type of standard errors, then the rounded table.
#' @export
print.summary.P2reg <- function(x, digits = 4, ...) {
  cat(sprintf("<P2reg>  discrete-time cause-specific hazard regression, (h, j) = (%s, %s)\n",
              attr(x, "h"), attr(x, "j")))
  cat(sprintf("  link: %s%s\n", attr(x, "link"),
              if (attr(x, "link") == "cloglog") " (exp(coef) is a hazard ratio)" else ""))
  cat(sprintf("  %.0f%% CIs; %s standard errors\n", 100 * attr(x, "conf.level"),
              if (isTRUE(attr(x, "cluster"))) "cluster-robust (subject id)" else "model-based"))
  y <- x
  for (nm in intersect(names(y), c("estimate", "se", "HR", "HR.lower", "HR.upper",
                                   "exp.coef", "exp.lower", "exp.upper", "p.value")))
    y[[nm]] <- round(y[[nm]], digits)
  print(as.data.frame(y), row.names = FALSE)
  invisible(x)
}

# Subsetting keeps the header attributes, so that e.g. summary(fit)[summary(fit)$term != "(Intercept)", ] still prints correctly.
#' @export
`[.summary.P2reg` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    for (a in c("h", "j", "cluster", "link", "conf.level")) attr(out, a) <- attr(x, a)
    class(out) <- class(x)
  }
  out
}

#' Forest plot of the effects estimated by P2reg()
#'
#' Draws, for every destination and covariate term, the exponentiated
#' coefficient (a hazard ratio under the cloglog link) with its confidence
#' interval, on a log scale with a reference line at 1. Intercepts are
#' omitted.
#'
#' @param x A "P2reg" object.
#' @param conf.level Confidence level of the intervals. Default 0.95.
#' @param terms Optional character vector of terms to show (default: all
#'   non-intercept terms).
#' @param col Colour of points and intervals.
#' @param xlab,main Axis label and title.
#' @param ... Passed to \code{plot()}.
#' @return The \code{summary()} table of the plotted terms, invisibly.
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
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' ids <- unique(panel$id)
#' panel$age <- stats::setNames(sample(40:80, length(ids), TRUE), ids)[as.character(panel$id)]
#' reg <- P2reg(prep2(panel, covariates = "age", states = st), h = "B", j = "B",
#'              formula = ~ age + log(d))
#' plot(reg)
#' @export
plot.P2reg <- function(x, conf.level = 0.95, terms = NULL, col = "#1E2761",
                       xlab = NULL, main = NULL, ...) {

  # The non-intercept terms to draw, with their exponentiated estimates and intervals.
  sm <- summary(x, conf.level = conf.level)
  sm <- sm[sm$term != "(Intercept)", , drop = FALSE]
  if (!is.null(terms)) sm <- sm[sm$term %in% terms, , drop = FALSE]
  if (!nrow(sm)) stop("No covariate terms to plot.", call. = FALSE)
  hr <- identical(attr(sm, "link"), "cloglog")
  est <- if (hr) sm$HR else sm$exp.coef
  lo  <- if (hr) sm$HR.lower else sm$exp.lower
  up  <- if (hr) sm$HR.upper else sm$exp.upper
  k <- nrow(sm)
  y <- rev(seq_len(k))

  # Points and intervals on a log scale, with a reference line at 1 and one labelled row per destination and term.
  if (is.null(xlab)) xlab <- if (hr) "Hazard ratio (log scale)" else "exp(coef) (log scale)"
  if (is.null(main)) main <- sprintf("(h, j) = (%s, %s)", x$h, x$j)
  op <- graphics::par(mar = c(4.5, max(8, max(nchar(paste(sm$l, sm$term))) * 0.55), 2.5, 1))
  on.exit(graphics::par(op))
  plot(est, y, log = "x", xlim = range(c(lo, up, 1), finite = TRUE), ylim = c(0.5, k + 0.5),
       yaxt = "n", ylab = "", xlab = xlab, main = main, pch = 15, col = col, ...)
  graphics::segments(lo, y, up, y, col = col, lwd = 2)
  graphics::abline(v = 1, lty = 2, col = "grey50")
  graphics::axis(2, at = y, labels = paste0(sm$l, ": ", sm$term), las = 1, cex.axis = 0.8)
  invisible(sm)
}

#' Covariate-adjusted second-order model for every history
#'
#' Fits \code{\link{P2reg}} for \strong{every} observed history \eqn{(h, j)}
#' (with \eqn{j} transient), modelling the moves \eqn{j \to \ell}, \eqn{\ell \ne j};
#' staying in \eqn{j} is obtained by difference. The \code{predict()} method
#' then assembles, for any covariate profile, the full covariate-adjusted
#' tensor \eqn{\hat P_{hj\ell}(Z)}, which \code{\link{ckequations}},
#' \code{\link{probtrans2}}-style propagation or \code{\link{simulate2}} accept
#' directly: this gives covariate-adjusted \eqn{n}-step predictions through
#' the extended Chapman-Kolmogorov relation.
#'
#' The tensor is time-homogeneous when the formula uses only baseline
#' covariates. If it uses \code{s} or \code{d}, \code{predict()} needs their
#' values in \code{newdata} and returns the one-step tensor at those values;
#' propagating it with \code{ckequations()} then assumes those values stay
#' fixed, which is only an approximation.
#'
#' @inheritParams P2reg
#' @param min.events Minimum number of observed moves from a history for its
#'   destination models to be fitted; histories below it keep the
#'   nonparametric RPE row (as in \code{\link{P2est}}). Default 5.
#' @return An object of class "P2reg_all": a list with \code{fits} (one
#'   \code{P2reg} per fitted history, named \code{"h|j"}), \code{rpe} (the
#'   \code{P2est} fit used for the rows that are not modelled), \code{pairs}
#'   (data frame of histories, whether each was modelled, and which terms were
#'   not estimable in it), \code{formula},
#'   \code{family}, \code{states}, \code{absorbing}.
#' @section Why this construction:
#' Each move is modelled as in \code{\link{P2reg}} and staying in \eqn{j} is
#' obtained as one minus the moves, so every row of the predicted tensor is a
#' probability distribution and the tensor can be propagated by
#' \code{\link{ckequations}}, which gives covariate-adjusted predictions with
#' the same machinery as \code{\link{P2est}}. Histories with fewer than
#' \code{min.events} moves keep the non-parametric estimate: a regression
#' there would rest on too few events (the usual guidance asks for roughly 5 to
#' 10 events per parameter; Vittinghoff and McCulloch, 2007). If the predicted
#' moves add up to more than one (possible because the per-destination models
#' are fitted separately), they are rescaled with a warning; a multinomial
#' model would avoid this and is a natural extension.
#'
#' @references
#' Vittinghoff, E. and McCulloch, C. E. (2007). Relaxing the rule of ten events
#' per variable in logistic and Cox regression. \emph{American Journal of
#' Epidemiology}, 165(6), 710-718.
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
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' ids <- unique(panel$id)
#' panel$age <- stats::setNames(sample(40:80, length(ids), TRUE), ids)[as.character(panel$id)]
#' all_fit <- P2reg_all(prep2(panel, covariates = "age", states = st), formula = ~ age)
#' P70 <- predict(all_fit, newdata = data.frame(age = 70))   # tensor for a 70-year-old
#' ckequations(P70, h = "A", j = "B", l = "C", nsteps = 4)
#' @export
P2reg_all <- function(object, formula = ~1, family = stats::binomial("cloglog"),
                      cluster = FALSE, min.events = 5, ...) {

  # The non-parametric fit (kept for the histories that are not modelled) and the histories (h, j) with j transient.
  stopifnot(inherits(object, "msm2data"))
  rpe <- P2est(object)
  histories <- object$Y |>
    dplyr::distinct(h = as.character(h), j = as.character(j)) |>
    dplyr::filter(!j %in% object$absorbing)

  # Moves out of j observed in each history.
  moves <- object$triples |>
    dplyr::filter(as.character(l) != as.character(j)) |>
    dplyr::summarise(n.moves = dplyr::n(), dest = list(sort(unique(as.character(l)))),
                     .by = c("h", "j")) |>
    dplyr::mutate(h = as.character(h), j = as.character(j))
  histories <- histories |>
    dplyr::left_join(moves, by = c("h", "j")) |>
    dplyr::mutate(n.moves = dplyr::coalesce(n.moves, 0L),
                  modelled = n.moves >= min.events)

  # P2reg() for every history with enough moves; terms that cannot be estimated in a history (glm() returns NA, typically a covariate constant in its risk set) are recorded.
  fit_one <- function(hh, jj, dest) {
    P2reg(object, h = hh, j = jj, l = dest, formula = formula, family = family, cluster = cluster, ...)
  }
  modelled <- dplyr::filter(histories, modelled)
  fits <- stats::setNames(purrr::pmap(list(modelled$h, modelled$j, modelled$dest), fit_one),
                          paste(modelled$h, modelled$j, sep = "|"))
  aliased <- purrr::map_chr(fits, \(fit) {
    na_terms <- unique(unlist(purrr::map(fit$models, \(m) names(which(is.na(stats::coef(m)))))))
    paste(na_terms, collapse = ", ")
  })
  pairs <- histories |>
    dplyr::mutate(aliased = dplyr::coalesce(aliased[paste(h, j, sep = "|")], "")) |>
    dplyr::select("h", "j", "modelled", "aliased") |>
    as.data.frame()

  structure(list(fits = fits, rpe = rpe, pairs = pairs, formula = formula,
                 family = family, states = object$states, absorbing = object$absorbing),
            class = "P2reg_all")
}

# Short description: link, covariates, histories modelled and terms that could not be estimated.
#' @export
print.P2reg_all <- function(x, ...) {
  cat(sprintf("<P2reg_all>  covariate-adjusted second-order model (%s link)\n", x$family$link))
  cat(sprintf("  covariates  : %s\n", deparse1(x$formula[[length(x$formula)]])))
  cat(sprintf("  histories   : %d modelled, %d kept at the RPE (too few moves)\n",
              sum(x$pairs$modelled), sum(!x$pairs$modelled)))
  al <- x$pairs[nzchar(x$pairs$aliased), , drop = FALSE]
  if (nrow(al))
    cat(sprintf("  not estimable (constant in the risk set; set to 0 in predictions): %s\n",
                paste0("(", al$h, ", ", al$j, ") ", al$aliased, collapse = "; ")))
  invisible(x)
}

# Coefficient tables of every modelled history, stacked with their (h, j).
#' @export
summary.P2reg_all <- function(object, conf.level = 0.95, ...) {
  purrr::map(object$fits, \(fit)
    tibble::tibble(h = fit$h, j = fit$j, as.data.frame(summary(fit, conf.level = conf.level)))) |>
    purrr::list_rbind() |>
    as.data.frame()
}

#' Covariate-adjusted second-order tensor for given covariate profiles
#'
#' @param object A "P2reg_all" object.
#' @param newdata Data frame of covariate profiles (one row per profile), with
#'   the variables used in the formula.
#' @param ... Unused.
#' @return For a single-row \code{newdata}, an \eqn{M \times M \times M} tensor
#'   in \code{P[j, l, h]} layout; otherwise a list of such tensors, one per row.
#'   For every modelled history the predicted move probabilities are taken
#'   from the GLMs and the probability of staying is one minus their sum (if
#'   the predicted moves sum to more than 1, they are rescaled to 1 and a
#'   warning is given). Non-modelled histories keep their RPE row, and
#'   absorbing states keep probability 1.
#' @export
predict.P2reg_all <- function(object, newdata, ...) {

  # Check newdata.
  if (!is.data.frame(newdata) || !nrow(newdata))
    stop("`newdata` must be a data frame with at least one row.", call. = FALSE)
  rescaled <- FALSE

  # One tensor per profile: start from the RPE tensor and replace the row of every modelled history by its predicted moves, staying = 1 - moves. Terms not estimable in a history act as 0 (predict.lm drops them); its generic rank-deficiency warning is replaced by one explicit warning below.
  one <- function(r) {
    P <- object$rpe$P
    for (k in names(object$fits)) {
      fit <- object$fits[[k]]
      pm <- withCallingHandlers(
        stats::predict(fit, newdata = newdata[r, , drop = FALSE]),
        warning = function(w) if (grepl("rank-deficient", conditionMessage(w)))
          invokeRestart("muffleWarning"))
      pr <- stats::setNames(as.numeric(pm[1L, ]), colnames(pm))   # keeps names with one move
      tot <- sum(pr)
      if (tot > 1) {
        pr <- pr / tot
        tot <- 1
        rescaled <<- TRUE
      }
      P[fit$j, , fit$h] <- 0
      P[fit$j, names(pr), fit$h] <- pr
      P[fit$j, fit$j, fit$h] <- 1 - tot
    }
    for (a in object$absorbing) P[a, a, ] <- 1
    P
  }
  out <- purrr::map(seq_len(nrow(newdata)), one)
  if (rescaled)
    warning("Predicted move probabilities summed to more than 1 for some history; ",
            "they were rescaled.", call. = FALSE)
  if (any(nzchar(object$pairs$aliased)))
    warning("Some terms are not estimable in some histories (the covariate is constant ",
            "in their risk set) and were set to 0; see print(object).", call. = FALSE)
  if (length(out) == 1L) out[[1L]] else out
}
