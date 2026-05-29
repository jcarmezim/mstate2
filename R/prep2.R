#' Build second-order counting processes from panel data
#'
#' \code{prep2} takes discrete-time panel data (one row per subject-time with the
#' occupied state) and constructs the second-order counting processes needed to
#' estimate 1-step second-order homogeneous transition probabilities
#' \eqn{P_{hj\ell} = P(X_s = \ell \mid X_{s-1} = j, X_{s-2} = h)}.
#'
#' For each subject the function forms every consecutive triple of observations
#' \eqn{(X_{s-2}, X_{s-1}, X_s) = (h, j, \ell)} and tabulates:
#' \itemize{
#'   \item \eqn{\tilde N_{hj\ell}(s)}: number of subjects observing the path
#'         \eqn{h \to j \to \ell} ending at time \eqn{s};
#'   \item \eqn{\tilde Y_{hj}(s-1)}: number of subjects at risk, i.e. occupying
#'         \eqn{h} at \eqn{s-2} and \eqn{j} at \eqn{s-1}.
#' }
#' Consecutive \emph{observations} are treated as consecutive time steps. Supply
#' a regular integer time grid (e.g. days) for the counting-process indices to
#' match the paper's definitions.
#'
#' @param data A data frame in long/panel format.
#' @param id Name of the subject identifier column.
#' @param time Name of the (discrete) time column.
#' @param state Name of the state column (factor, integer, or character).
#' @param states Optional vector giving the full state space and its ordering.
#'   Defaults to the sorted unique observed states. Use this to fix the state
#'   ordering / include never-occupied states.
#' @param absorbing Optional vector of absorbing states. If \code{NULL} (default)
#'   absorbing states are inferred as those with no observed onward transition to
#'   a different state.
#'
#' @return An object of class \code{"msm2data"}: a list with
#'   \item{N}{data.table of \eqn{(h, j, l, s, N)} triple counts.}
#'   \item{Y}{data.table of \eqn{(h, j, s, Y)} at-risk counts.}
#'   \item{states}{the ordered state space.}
#'   \item{absorbing}{the absorbing states.}
#'   \item{n}{number of subjects.}
#'   \item{ntriples}{total number of observed triples.}
#' @export
prep2 <- function(data, id = "id", time = "time", state = "state",
                  states = NULL, absorbing = NULL) {

  stopifnot(is.data.frame(data))
  for (col in c(id, time, state)) {
    if (!col %in% names(data))
      stop(sprintf("Column '%s' not found in `data`.", col), call. = FALSE)
  }

  DT <- data.table::as.data.table(data)
  DT <- DT[, c(id, time, state), with = FALSE]
  data.table::setnames(DT, c("id", "time", "state"))

  if (anyNA(DT$state) || anyNA(DT$time) || anyNA(DT$id))
    stop("`id`, `time` and `state` must not contain NA.", call. = FALSE)

  ## fix the state space / ordering
  if (is.null(states)) {
    states <- if (is.factor(DT$state)) levels(droplevels(DT$state))
              else sort(unique(DT$state))
  }
  DT$state <- factor(DT$state, levels = states)
  if (anyNA(DT$state))
    stop("Some observed states are not in `states`.", call. = FALSE)

  n <- data.table::uniqueN(DT$id)

  ## order within subject and form lagged triples: h = X_{s-2}, j = X_{s-1}, l = X_s
  data.table::setorder(DT, id, time)
  DT[, c("h", "j") := list(data.table::shift(state, 2L),
                           data.table::shift(state, 1L)), by = "id"]
  trip <- DT[!is.na(h) & !is.na(j)]
  data.table::setnames(trip, "state", "l")

  if (nrow(trip) == 0L)
    stop("No second-order triples found: each subject needs >= 3 observations.",
         call. = FALSE)

  ## N_{hjl}(s): triple counts indexed by ending time s
  N <- trip[, list(N = .N), by = list(h, j, l, s = time)]
  ## Y_{hj}(s-1): at-risk = number of triples (h, j, .) ending at time s
  Y <- trip[, list(Y = .N), by = list(h, j, s = time)]

  ## infer absorbing states if not supplied: occupied as j but no j -> l, l != j
  if (is.null(absorbing)) {
    moved <- unique(trip[as.character(j) != as.character(l), as.character(j)])
    occupied <- unique(c(as.character(trip$j), as.character(trip$l)))
    absorbing <- setdiff(occupied, moved)
  }

  structure(
    list(N = N, Y = Y, states = states,
         absorbing = as.character(absorbing),
         n = n, ntriples = nrow(trip)),
    class = "msm2data"
  )
}

#' @export
print.msm2data <- function(x, ...) {
  cat("<msm2data>  second-order counting processes\n")
  cat(sprintf("  subjects        : %d\n", x$n))
  cat(sprintf("  observed triples: %d\n", x$ntriples))
  cat(sprintf("  states (%d)      : %s\n", length(x$states),
              paste(x$states, collapse = ", ")))
  cat(sprintf("  absorbing       : %s\n",
              if (length(x$absorbing)) paste(x$absorbing, collapse = ", ") else "none"))
  cat(sprintf("  distinct (h,j)  : %d\n",
              nrow(unique(x$Y[, list(h, j)]))))
  invisible(x)
}
