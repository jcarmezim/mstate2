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
#' @param data A data frame in long/panel format, or an mstate "msdata" object
#'   (converted automatically via from_msdata()).
#' @param id,time,state Column names for subject id, discrete time, and state.
#' @param states Optional state space / ordering. Defaults to sorted observed states.
#' @param absorbing Optional absorbing states; inferred if NULL.
#' @param drop.na If TRUE, drop rows with NA in id/time/state (with a warning)
#'   instead of erroring. Default FALSE.
#' @param check.consecutive If TRUE (default), warn when any subject's times are
#'   not consecutive integers.
#' @return An object of class "msm2data".
#' @export
prep2 <- function(data, id = "id", time = "time", state = "state",
                  states = NULL, absorbing = NULL,
                  drop.na = FALSE, check.consecutive = TRUE) {

  if (inherits(data, "msdata")) data <- from_msdata(data)
  stopifnot(is.data.frame(data))

  cols <- c(id = id, time = time, state = state)
  missing_cols <- cols[!cols %in% names(data)]
  if (length(missing_cols))
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  DT <- data.table::as.data.table(data)[, .SD, .SDcols = cols]
  data.table::setnames(DT, names(cols))

  bad <- is.na(DT$id) | is.na(DT$time) | is.na(DT$state)
  if (any(bad)) {
    if (!drop.na)
      stop("`id`, `time` and `state` must not contain NA ",
           "(set drop.na = TRUE to drop these rows).", call. = FALSE)
    warning(sprintf("Dropping %d row(s) with NA in id/time/state.", sum(bad)),
            call. = FALSE)
    DT <- DT[!bad]
  }

  if (is.null(states))
    states <- if (is.factor(DT$state)) levels(droplevels(DT$state))
              else sort(unique(DT$state))
  DT[, state := factor(state, levels = states)]
  if (anyNA(DT$state)) stop("Some observed states are not in `states`.", call. = FALSE)

  n <- data.table::uniqueN(DT$id)
  data.table::setorder(DT, id, time)

  if (check.consecutive) {
    gaps <- DT[, .(ok = all(diff(time) == 1L)), by = id]
    if (!all(gaps$ok))
      warning(sprintf(paste("%d subject(s) have non-consecutive time values;",
                            "triples are formed from consecutive observations."),
                      sum(!gaps$ok)), call. = FALSE)
  }

  DT[, c("h", "j") := .(data.table::shift(state, 2L),
                        data.table::shift(state, 1L)), by = id]
  trip <- DT[!is.na(h) & !is.na(j)]
  data.table::setnames(trip, "state", "l")
  if (nrow(trip) == 0L)
    stop("No second-order triples found: each subject needs >= 3 observations.",
         call. = FALSE)

  N <- trip[, .(N = .N),      keyby = .(h, j, l, s = time)]
  Y <- N[,    .(Y = sum(N)),  keyby = .(h, j, s)]            # Y = N summed over l

  if (is.null(absorbing)) {
    moved    <- unique(trip[as.integer(j) != as.integer(l), j])
    occupied <- union(levels(droplevels(trip$j)), levels(droplevels(trip$l)))
    absorbing <- setdiff(occupied, as.character(moved))
  }

  structure(
    list(N = N, Y = Y, states = states, absorbing = as.character(absorbing),
         n = n, ntriples = nrow(trip), time.range = range(DT$time)),
    class = "msm2data"
  )
}

#' @export
print.msm2data <- function(x, ...) {
  cat("<msm2data>  second-order counting processes\n")
  cat(sprintf("  subjects        : %d\n", x$n))
  cat(sprintf("  observed triples: %d\n", x$ntriples))
  cat(sprintf("  time range      : %d - %d\n", x$time.range[1], x$time.range[2]))
  cat(sprintf("  states (%d)      : %s\n", length(x$states),
              paste(x$states, collapse = ", ")))
  cat(sprintf("  absorbing       : %s\n",
              if (length(x$absorbing)) paste(x$absorbing, collapse = ", ") else "none"))
  cat(sprintf("  distinct (h,j)  : %d\n", nrow(unique(x$Y[, .(h, j)]))))
  invisible(x)
}

#' @export
summary.msm2data <- function(object, ...) {
  exposure <- object$Y[, .(t_hj = .N, total_at_risk = sum(Y),
                           s_min = min(s), s_max = max(s)), keyby = .(h, j)]
  cat("<msm2data summary>\n")
  cat(sprintf("  %d subjects, %d triples, time %d-%d\n",
              object$n, object$ntriples,
              object$time.range[1], object$time.range[2]))
  cat("  exposure per (h, j) pair:\n")
  print(exposure)
  invisible(exposure)
}

#' Convert an mstate 'msdata' object to discrete-time panel format
#'
#' @param x An mstate "msdata" object (or a data frame with the named columns).
#' @param id,from,to,Tstart,Tstop,status Column names.
#' @param round_fun Discretisation function for transition times. Default rnd().
#' @return A data.table in panel format (id, time, state).
#' @export
from_msdata <- function(x, id = "id", from = "from", to = "to",
                        Tstart = "Tstart", Tstop = "Tstop",
                        status = "status", round_fun = rnd) {
  DT  <- data.table::as.data.table(x)
  map <- c(id, from, to, Tstart, Tstop, status)
  miss <- setdiff(map, names(DT))
  if (length(miss))
    stop("msdata missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  DT <- DT[, .SD, .SDcols = map]
  data.table::setnames(DT, c("id", "from", "to", "Tstart", "Tstop", "status"))

  build <- function(s, gid) {
    data.table::setorder(s, Tstart, Tstop)
    ev     <- s[status == 1][order(Tstop)]            # realized jumps, in order
    states <- c(s$from[1L], ev$to)                    # start state, then destinations
    rt     <- round_fun(c(s$Tstart[1L], ev$Tstop))    # their absolute (rounded) times
    grid   <- seq.int(rt[1L], round_fun(max(s$Tstop)))
    data.table::data.table(id = gid, time = grid,
                           state = states[findInterval(grid, rt)])
  }
  DT[, build(.SD, .BY[[1]]), by = id][, .(id, time, state)]
}
