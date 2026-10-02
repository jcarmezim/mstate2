#' mstate2: Second-Order Markov Multistate Models
#'
#' Nonparametric estimation and inference for second-order Markov multistate
#' models in discrete time, in which the next state depends on the state
#' occupied at the current time \emph{and} at the previous time. Complements the first-order \pkg{mstate}
#' framework and follows its conventions.
#'
#' @section Workflow:
#' \enumerate{
#'   \item \strong{Data.} \code{\link{sojourn_to_panel}},
#'     \code{\link{from_msdata}} or \code{\link{msprep2}} build a daily panel;
#'     \code{\link{prep2}} computes the second-order counting processes.
#'   \item \strong{Estimation.} \code{\link{P2est}} (relative probability estimator with
#'     variances and confidence intervals) and \code{\link{P2boot}}
#'     (subject-level bootstrap).
#'   \item \strong{Prediction.} \code{\link{ckequations}} (extended
#'     Chapman-Kolmogorov relation) and \code{\link{probtrans2}}
#'     (\pkg{mstate} \code{probtrans} format, with \code{plot()}).
#'   \item \strong{Does the history matter?} \code{\link{compare2}},
#'     \code{\link{overlap_step}} and \code{\link{markov_test}} (Wald and
#'     likelihood-ratio tests of the first-order assumption).
#'   \item \strong{First versus second order.} \code{\link{P1est}},
#'     \code{\link{compare_order}}, \code{\link{divergence}},
#'     \code{\link{divergence_table}} and \code{\link{as_tmat}}.
#'   \item \strong{Covariates.} \code{\link{P2reg}} (discrete-time cloglog
#'     cause-specific hazard regression, with time-varying or
#'     duration-dependent baselines) and \code{\link{P2reg_all}}
#'     (covariate-adjusted multi-step prediction).
#'   \item \strong{Simulation.} \code{\link{simulate2}}.
#' }
#'
#' See \code{vignette("mstate2")} for a guided analysis and
#' \code{vignette("reference")} for the complete function reference.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats qnorm runif plogis qlogis setNames
#' @importFrom graphics axis legend lines mtext points polygon abline
#' @importFrom grDevices col2rgb rgb
#' @importFrom utils globalVariables
"_PACKAGE"

## Quiet R CMD check on data.table non-standard-evaluation symbols.
utils::globalVariables(c(
  "h", "j", "l", "s", "N", "Y", "state", "time", "id",
  "from", "to", "Tstart", "Tstop", "status",
  "n.trans", "at.risk", "p", "se", "lower", "upper",
  "ok", "total_at_risk", "s_min", "s_max", ".", "..cols", "..map",
  "..out_cols", "p.value", "d", "w", "p.adj"
))
