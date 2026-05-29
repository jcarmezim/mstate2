#' Round half up (matching the convention used in the DIVINE analysis)
#'
#' @param x numeric vector.
#' @return \code{trunc(x + sign(x) * 0.5)}.
#' @export
rnd <- function(x) trunc(x + sign(x) * 0.5)

#' Convert sojourn-time data to discrete-time panel format
#'
#' Many multistate datasets (including the DIVINE COVID-19 cohort) store each
#' subject's trajectory as the time spent in each state plus indicators for the
#' absorbing states reached. \code{sojourn_to_panel} reconstructs the implied
#' day-by-day state sequence so it can be passed to \code{\link{prep2}}.
#'
#' Each transient segment is included for its (rounded) duration, in the
#' canonical order given by \code{segments}; the absorbing state whose indicator
#' equals 1 is appended once at the end. Durations are discretised with
#' \code{round_fun} (default \code{\link{rnd}}, round-half-up), so a segment whose
#' rounded duration is 0 is skipped.
#'
#' @param data A data frame, one row per subject.
#' @param id Name of the subject id column.
#' @param segments Named character vector mapping transient state labels to their
#'   duration columns, in the order states are visited, e.g.
#'   \code{c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")}.
#' @param absorbing Named character vector mapping absorbing state labels to their
#'   0/1 indicator columns, e.g. \code{c(Disch = "disch.s", Death = "death.s")}.
#' @param round_fun Function used to discretise durations. Default \code{rnd}.
#'
#' @return A \code{data.table} in panel format with columns \code{id}, \code{time}
#'   (integer days from 0), and \code{state}.
#' @export
sojourn_to_panel <- function(data, id, segments, absorbing, round_fun = rnd) {
  stopifnot(is.data.frame(data))
  seg_lab <- names(segments); seg_col <- unname(segments)
  abs_lab <- names(absorbing); abs_col <- unname(absorbing)
  need <- c(id, seg_col, abs_col)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  DT <- data.table::as.data.table(data)
  out <- vector("list", nrow(DT))
  for (r in seq_len(nrow(DT))) {
    durs <- round_fun(as.numeric(DT[r, seg_col, with = FALSE]))
    durs[is.na(durs) | durs < 0] <- 0
    seq_states <- rep(seg_lab, durs)

    ind <- as.numeric(DT[r, abs_col, with = FALSE])
    hit <- which(ind == 1)
    if (length(hit)) seq_states <- c(seq_states, abs_lab[hit[1]])

    if (length(seq_states) == 0L) next
    out[[r]] <- data.table::data.table(
      id    = DT[[id]][r],
      time  = seq_along(seq_states) - 1L,
      state = seq_states
    )
  }
  data.table::rbindlist(out)
}

#' Coerce a panel / msm2data object to an mstate transition matrix
#'
#' Builds an \pkg{mstate}-style transition matrix (via
#' \code{mstate::transMat}) describing which one-step transitions are observed,
#' giving interoperability with first-order \pkg{mstate} workflows. The
#' second-order model of this package can then be read as a memory-augmented
#' layer over that first-order structure.
#'
#' @param object A \code{"msm2data"} object or a panel data frame.
#' @param ... Passed to \code{\link{prep2}} when \code{object} is a data frame.
#' @return A transition matrix as returned by \code{mstate::transMat}, with one
#'   row/column per state and integer transition numbers in the allowed cells.
#' @export
as_tmat <- function(object, ...) {
  if (!requireNamespace("mstate", quietly = TRUE))
    stop("Package 'mstate' is required for as_tmat().", call. = FALSE)
  if (!inherits(object, "msm2data")) object <- prep2(object, ...)

  states <- object$states
  ## observed one-step transitions j -> l (l != j) from the triple table
  tr <- unique(object$N[as.character(j) != as.character(l), list(j, l)])
  to_list <- lapply(states, function(s) {
    as.character(tr[as.character(j) == s, l])
  })
  names(to_list) <- states
  mstate::transMat(to_list, names = states)
}
