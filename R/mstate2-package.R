#' mstate2: Second-Order Markov Multistate Models
#'
#' Nonparametric estimation and inference for second-order Markov multistate
#' models in discrete time: the RPE and CPE estimators of the 1-step second-order
#' transition probabilities with variances and confidence intervals, an extended
#' Chapman-Kolmogorov relation for n-step prediction, evolution intervals, a
#' sojourn-time / mstate panel converter, and a simulator. Complements the
#' first-order \pkg{mstate} framework.
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
  "n.trans", "at.risk", "inv", "t.hj", "p", "se", "lower", "upper",
  "ok", "t_hj", "total_at_risk", "s_min", "s_max", ".", "..cols", "..map"
))
