#' Round half up (the convention used in the DIVINE analysis)
#' @param x numeric vector.
#' @return \code{trunc(x + sign(x) * 0.5)}.
#' @export
rnd <- function(x) trunc(x + sign(x) * 0.5)

#' Convert sojourn-time data to discrete-time panel format
#'
#' @param data A data frame, one row per subject.
#' @param id Name of the subject id column.
#' @param segments Named character vector mapping transient state labels to their
#'   duration columns, in visiting order, e.g.
#'   \code{c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")}.
#' @param absorbing Named character vector mapping absorbing state labels to their
#'   0/1 indicator columns, e.g. \code{c(Disch = "disch.s", Death = "death.s")}.
#' @param round_fun Discretisation function for durations. Default \code{rnd}.
#' @return A data.table in panel format (id, time, state).
#' @export
sojourn_to_panel <- function(data, id, segments, absorbing, round_fun = rnd) {
  stopifnot(is.data.frame(data))
  seg_lab <- names(segments); seg_col <- unname(segments)
  abs_lab <- names(absorbing); abs_col <- unname(absorbing)
  miss <- setdiff(c(id, seg_col, abs_col), names(data))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

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
      id = DT[[id]][r], time = seq_along(seq_states) - 1L, state = seq_states)
  }
  data.table::rbindlist(out)
}

#' Coerce a panel / msm2data object to an mstate transition matrix
#'
#' @param object A "msm2data" object or a panel data frame.
#' @param ... Passed to \code{\link{prep2}} when \code{object} is a data frame.
#' @return A transition matrix as returned by \code{mstate::transMat}.
#' @export
as_tmat <- function(object, ...) {
  if (!requireNamespace("mstate", quietly = TRUE))
    stop("Package 'mstate' is required for as_tmat().", call. = FALSE)
  if (!inherits(object, "msm2data")) object <- prep2(object, ...)

  states <- object$states
  tr <- unique(object$N[as.integer(j) != as.integer(l), .(j, l)])
  to_list <- lapply(states, function(s)
    match(as.character(tr[as.character(j) == s, l]), states))   # <- integer indices
  names(to_list) <- states
  mstate::transMat(to_list, names = states)
}
