#' mstate2: Second-Order Markov Multistate Models
#'
#' Nonparametric estimation and inference for second-order Markov multistate models in discrete time, in which the next state depends on the state occupied at the current time \emph{and} at the previous time. Complements the first-order \pkg{mstate} framework and follows its conventions.
#'
#' The package implements the methods of Najera-Zuloaga, Besalu and Gomez Melis (2025), so that their analysis of the DIVINE cohort can be reproduced, and adds subject-level bootstrap intervals for the multi-step predictions.
#'
#' @section Workflow:
#' \enumerate{
#'   \item \strong{Data.} \code{\link{msprep2}} turns sojourn data (as DIVINE), wide data with one \code{Surv(time, status)} per state (as for \code{mstate::msprep()}) or an \pkg{mstate} \code{msdata} object into a daily panel and reports every record it drops or changes (\code{\link{rnd}} discretises the times); \code{\link{prep2}} computes the second-order counting processes.
#'   \item \strong{Estimation.} \code{\link{P2est}} (relative probability estimator with variances and confidence intervals).
#'   \item \strong{Prediction.} \code{\link{ckequations}} (extended Chapman-Kolmogorov relation, with evolution intervals).
#'   \item \strong{Does the previous time matter, and for how long?} \code{\link{compare2}} and \code{\link{overlap_step}}.
#'   \item \strong{Bootstrap intervals.} \code{\link{P2boot}} resamples patients; \code{ckequations()} and \code{compare2()} then report percentile bootstrap intervals.
#'   \item \strong{Simulation.} \code{\link{simulate2}} (the simulation design of the paper).
#' }
#'
#' @section Implementation:
#'
#' See \code{vignette("mstate2")} for a guided analysis, \code{vignette("paper")} for the reproduction of the methods paper and \code{vignette("reference")} for the complete function reference.
#'
#' @keywords internal
#' @importFrom stats qnorm runif plogis qlogis setNames
#' @importFrom graphics axis legend lines mtext points polygon abline
#' @importFrom grDevices col2rgb rgb
#' @importFrom utils globalVariables
"_PACKAGE"

utils::globalVariables(c(
  "h", "j", "l", "s", "N", "Y", "state", "time", "id",
  "n.trans", "at.risk", "p", "se", "lower", "upper", "se_logit", "inside",
  "n", "ok", "total_at_risk", "s_min", "s_max",
  ".row", ".pos", "raw", "dur", "ind",
  "from", "to", "issue", "exit", "first", "last", "initial", "origin", "t0",
  "t_max", ".col", ".scol", ".ord", ".t", "tval", "step", "end_step",
  ".first", ".k", ".n", "prev", "kept_state", ".abs", ".next", ".len", ".i",
  "total", "Tstart", "Tstop", "status", ".seen", "end", ".nxt", ".has",
  ".run", ".pos", "d"
))
