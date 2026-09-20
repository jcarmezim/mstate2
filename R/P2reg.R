#' Discrete-time cause-specific hazard regression for 1-step second-order
#' transitions
#'
#' Fits the covariate-adjusted discrete-time cause-specific hazard model
#' \deqn{\mathrm{cloglog}(\lambda_{hj\ell}(Z)) = \gamma_{hj\ell} + \beta_{hj\ell}^\top Z}
#' for a fixed (preceding, current) risk set \eqn{(h, j)}, one binary GLM per
#' destination state \eqn{\ell}. This is the covariate extension of the
#' second-order framework: \code{\link{P2est}} is non-parametric and
#' covariate-free, so it cannot answer "does this change for a 70-year-old
#' with kidney disease?" -- \code{P2reg} can.
#'
#' Under the (default) complementary log-log link, this is the discrete-time
#' analogue of the Cox proportional-hazards model (Prentice & Gloeckler, 1978):
#' \eqn{\exp(\beta_{hj\ell})} retains a hazard-ratio interpretation, directly
#' comparable with hazard ratios from a first-order \code{coxph}-based
#' analysis of the same data.
#'
#' The risk set for (h, j) is \code{object$triples} filtered to that pair: one
#' row per subject per at-risk instant, exactly \code{\link{prep2}}'s
#' \eqn{\tilde Y_{hj}(s-1)} at the individual level. All destinations observed
#' from that risk set share it, so a single data preparation step fits every
#' destination-specific model. With \code{formula = ~1} (no covariates), the
#' fitted probability for a given l is, by construction, the sample proportion
#' of the risk set that moved to l -- which is exactly the RPE (Eq. 9) computed
#' by \code{P2est(..., estimator = "RPE")}, regardless of the link function
#' used. This identity is a useful internal check that \code{P2reg} and
#' \code{P2est} agree, and is verified in the package's test suite.
#'
#' @param object An "msm2data" object built with \code{prep2(..., covariates = )}.
#' @param h,j The preceding and current state defining the risk set (label or index).
#' @param l Optional target state(s) (label or index). If NULL (default), a
#'   separate cause-specific model is fit for every destination state observed
#'   from (h, j).
#' @param formula A one-sided formula giving the covariate predictors, e.g.
#'   \code{~ age + sex}. Default \code{~1} (baseline hazard only, no covariates).
#'   Variables must be columns of \code{object$triples} (ordinarily the
#'   \code{covariates} supplied to \code{\link{prep2}}).
#' @param family The GLM family/link, as a family object. Default
#'   \code{stats::binomial("cloglog")} (see Details); \code{exp(coef)} is
#'   labelled a hazard ratio only when \code{family$link == "cloglog"}.
#' @param cluster If TRUE, replace the model-based standard errors with a
#'   cluster-robust ("sandwich") estimate clustered by subject id, via
#'   \code{sandwich::vcovCL()} (requires the \pkg{sandwich} package). Relevant
#'   when a subject can contribute more than one at-risk instant to the same
#'   (h, j) risk set. Default FALSE (model-based standard errors).
#' @param ... Passed to \code{\link[stats]{glm}}.
#' @return An object of class "P2reg": a list with \code{models} (one
#'   \code{stats::glm} fit per destination state, named by state), \code{h},
#'   \code{j}, \code{states}, \code{formula}, \code{family}, \code{cluster},
#'   \code{n} (subjects at risk) and \code{n.risk} (at-risk instants).
#' @export
P2reg <- function(object, h, j, l = NULL, formula = ~1,
                  family = stats::binomial("cloglog"), cluster = FALSE, ...) {
  stopifnot(inherits(object, "msm2data"))
  if (!length(object$covariates))
    stop("`object` carries no covariates; build it with prep2(..., covariates = ).",
         call. = FALSE)
  if (cluster && !requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for cluster = TRUE.", call. = FALSE)

  states <- object$states
  hi <- .resolve(h, states); ji <- .resolve(j, states)
  if (anyNA(c(hi, ji))) stop("`h` or `j` not found in the state space.", call. = FALSE)
  hlab <- states[hi]; jlab <- states[ji]

  trip <- as.data.frame(object$triples)
  risk <- trip[as.character(trip$h) == hlab & as.character(trip$j) == jlab, , drop = FALSE]
  if (!nrow(risk))
    stop("No observed (h, j) = (", hlab, ", ", jlab, ") triples in `object`.", call. = FALSE)

  targets <- if (is.null(l)) sort(unique(as.character(risk$l))) else {
    li <- .resolve(l, states)
    if (anyNA(li)) stop("`l` not found in the state space.", call. = FALSE)
    states[li]
  }

  rhs  <- formula[[length(formula)]]
  vars <- all.vars(rhs)
  miss <- setdiff(vars, object$covariates)
  if (length(miss))
    stop("`formula` references variable(s) not in object$covariates: ",
         paste(miss, collapse = ", "), call. = FALSE)

  fits <- lapply(targets, function(tgt) {
    dat <- risk
    dat$.y <- as.integer(as.character(dat$l) == tgt)
    fm  <- stats::as.formula(paste0(".y ~ ", deparse(rhs)))
    fit <- stats::glm(fm, data = dat, family = family, ...)
    if (cluster) fit$P2reg.vcov <- sandwich::vcovCL(fit, cluster = dat$id)
    fit
  })
  names(fits) <- targets

  structure(list(models = fits, h = hlab, j = jlab, states = states,
                 formula = formula, family = family, cluster = cluster,
                 n = data.table::uniqueN(risk$id), n.risk = nrow(risk)),
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
#' @export
predict.P2reg <- function(object, newdata, type = "response", ...) {
  sapply(object$models, function(fit) stats::predict(fit, newdata = newdata, type = type, ...))
}

#' @export
print.P2reg <- function(x, ...) {
  cat(sprintf("<P2reg>  %d destination-specific %s hazard model(s) from (h, j) = (%s, %s)\n",
              length(x$models), x$family$link, x$h, x$j))
  cat(sprintf("  destinations: %s\n", paste(names(x$models), collapse = ", ")))
  cat(sprintf("  %d subject(s), %d at-risk instant(s)\n", x$n, x$n.risk))
  cat(sprintf("  covariates   : %s\n", deparse(x$formula[[length(x$formula)]])))
  cat(sprintf("  std. errors  : %s\n", if (x$cluster) "cluster-robust (subject id)" else "model-based"))
  invisible(x)
}

#' @export
summary.P2reg <- function(object, conf.level = 0.95, ...) {
  z     <- stats::qnorm(1 - (1 - conf.level) / 2)
  hr_ok <- identical(object$family$link, "cloglog")

  rows <- lapply(names(object$models), function(tgt) {
    fit <- object$models[[tgt]]
    V   <- if (!is.null(fit$P2reg.vcov)) fit$P2reg.vcov else stats::vcov(fit)
    b   <- stats::coef(fit)
    se  <- sqrt(diag(V))
    data.frame(l = tgt, term = names(b), estimate = unname(b), se = unname(se),
              exp.coef  = exp(unname(b)),
              exp.lower = exp(unname(b) - z * unname(se)),
              exp.upper = exp(unname(b) + z * unname(se)),
              p.value   = 2 * stats::pnorm(-abs(unname(b) / unname(se))),
              row.names = NULL)
  })
  out <- do.call(rbind, rows)
  if (hr_ok) names(out)[names(out) %in% c("exp.coef", "exp.lower", "exp.upper")] <-
    c("HR", "HR.lower", "HR.upper")

  structure(out, class = c("summary.P2reg", "data.frame"),
            h = object$h, j = object$j, cluster = object$cluster,
            link = object$family$link, conf.level = conf.level)
}

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
  print(y, row.names = FALSE)
  invisible(x)
}
