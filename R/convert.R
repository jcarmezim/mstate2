#' Round half away from zero
#'
#' Function to turn the sojourn times into whole days (the default discretisation of \code{\link{msprep2}}). Unlike \code{\link{round}}, which rounds halves to the even number (\code{round(0.5) = 0}, \code{round(2.5) = 2}), \code{rnd} always rounds halves away from zero (\code{rnd(0.5) = 1}, \code{rnd(2.5) = 3}).
#'
#' @param x numeric vector.
#' @return \code{trunc(x + sign(x) * 0.5)}.
#' @examples
#' rnd(c(-2.5, -0.4, 0.4, 0.5, 2.5))
#' @export
rnd <- function(x) trunc(x + sign(x) * 0.5)

#' Convert sojourn-time data to discrete-time panel format
#'
#' Turns one row per subject, with the time spent in each transient state and an indicator of the absorbing state that ended follow-up, into one row per subject and time unit with the state occupied.
#'
#' @param data A data frame, one row per subject.
#' @param id Name of the subject id column.
#' @param segments Named character vector mapping transient state labels to their duration columns, in visiting order.
#' @param absorbing Named character vector mapping absorbing state labels to their 0/1 indicator columns.
#' @param round_fun Discretisation function for durations. Default \code{rnd}. Each sojourn is rounded once and then expanded into that many time units. A visit with a positive duration that rounds to 0 (with \code{rnd}, a duration below 0.5) is dropped, together with the transitions into and out of it; a warning reports how many visits were lost, by state.
#' @return A tibble in panel format with columns \code{id}, \code{time} (0, 1, 2, ... within each subject) and \code{state} (character). A subject whose follow-up ended in an absorbing state gets one final row with that state; a subject with no indicator equal to 1 (still under observation) has no absorbing row.
#' @details Internal function, kept for reference and for the tests: users call \code{msprep2(data, durations = , outcome = )}, which gives the same panel and reports the visits it loses.
#' @noRd
sojourn_to_panel <- function(data, id, segments, absorbing, round_fun = rnd) {

  # Check the input: every column named in `id`, `segments` and `absorbing` must exist in `data`.
  stopifnot(is.data.frame(data))

  seg_lab <- names(segments)
  seg_col <- unname(segments)
  abs_lab <- names(absorbing)
  abs_col <- unname(absorbing)
  miss <- setdiff(c(id, seg_col, abs_col), names(data))

  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  # Number the subjects (.row) so that the output keeps the row order of `data`, and rename the id column to `id`.
  subjects <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::select(".row", id = dplyr::all_of(id), dplyr::all_of(c(seg_col, abs_col)))

  # Transient visits: one row per subject and transient state, with the recorded duration, its position in the visiting order and the rounded number of time units.
  visits <- subjects |>
    dplyr::select(".row", "id", dplyr::all_of(seg_col)) |>
    tidyr::pivot_longer(dplyr::all_of(seg_col), names_to = "state", values_to = "raw") |>
    dplyr::mutate(.pos  = match(state, seg_col), # visiting order
                  state = seg_lab[.pos], # column name -> state label
                  dur   = round_fun(raw))

  # Visits with a positive recorded duration that round to 0 time units disappear from the panel, and with them the transitions into and out of that state (A -> B -> C becomes A -> C). Count them by state and warn, so that the loss is never silent.
  lost <- visits |>
    dplyr::filter(!is.na(raw), raw > 0, is.na(dur) | dur <= 0) |>
    dplyr::count(state = factor(state, levels = seg_lab), .drop = TRUE)

  if (nrow(lost))
    warning(sprintf(paste("%d visit(s) with a positive duration rounded to 0 time units and were dropped (%s); the transitions into and out of them are lost. Use a finer time unit or another `round_fun`."),
                    sum(lost$n), paste0(lost$state, ": ", lost$n, collapse = ", ")),
            call. = FALSE)

  # Expand each visit into one row per time unit spent there. A missing or negative duration means "never visited this state" and, like a zero duration, contributes no rows.
  transient_rows <- visits |>
    dplyr::mutate(dur = as.integer(dplyr::if_else(is.na(dur) | dur < 0, 0, dur))) |>
    tidyr::uncount(dur) |>
    dplyr::select(".row", "id", "state", ".pos")

  # Absorbing state: the first indicator column equal to 1 says how follow-up ended and adds one final row. Subjects with no indicator equal to 1 (still under observation) get no absorbing row.
  absorbing_rows <- subjects |>
    dplyr::select(".row", "id", dplyr::all_of(abs_col)) |>
    tidyr::pivot_longer(dplyr::all_of(abs_col), names_to = "state", values_to = "ind") |>
    dplyr::mutate(.pos  = length(seg_lab) + match(state, abs_col),  # after every transient state
                  state = abs_lab[match(state, abs_col)]) |>
    dplyr::filter(as.numeric(ind) == 1) |>
    dplyr::filter(.pos == min(.pos), .by = ".row") |>
    dplyr::select(".row", "id", "state", ".pos")

  # Put each subject's rows in visiting order and number the time units 0, 1, 2, ... within the subject.
  dplyr::bind_rows(transient_rows, absorbing_rows) |>
    dplyr::arrange(.row, .pos) |>
    dplyr::mutate(time = dplyr::row_number() - 1L, .by = ".row") |>
    dplyr::select("id", "time", "state")
}

#' Transition matrix of the observed moves, in mstate format
#'
#' Builds the transition matrix of every move \eqn{j \to \ell} (\eqn{\ell \neq j}) observed in the data, including each subject's first move, in the format of \code{mstate::transMat()}: integers numbering the transitions row by row, \code{NA} elsewhere, dimnames \code{from}/\code{to}. A first-order \pkg{mstate} analysis built on it (\code{msprep()}, \code{msfit()}, \code{probtrans()}) uses the same state space, in the same order, as the second-order fit, which \code{\link{compare_order}} requires. The \pkg{mstate} package is not needed.
#'
#' @param object An "msm2data" object from \code{\link{prep2}}, an "msm2prep" object from \code{\link{msprep2}}, or a panel data frame.
#' @param ... Passed to \code{\link{prep2}} when \code{object} is not an "msm2data" object.
#' @return A transition matrix as returned by \code{mstate::transMat()}.
#' @seealso \code{\link{compare_order}}
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
#' panel <- simulate2(200, tens, first, init = c(A = 1, B = 0, C = 0))
#' as_tmat(prep2(panel))   # A -> B (first move only) and B -> C
#' @export
as_tmat <- function(object, ...) {

  # The counting processes, with every observed one-step move (pairs), including first moves that have no triple.
  if (!inherits(object, "msm2data")) object <- prep2(object, ...)
  states <- object$states
  moves <- object$pairs |>
    dplyr::filter(as.character(from) != as.character(to)) |>
    dplyr::transmute(from = match(as.character(from), states), to = match(as.character(to), states)) |>
    dplyr::arrange(from, to)

  # Number the transitions row by row, as mstate::transMat() does.
  tmat <- matrix(NA_integer_, length(states), length(states), dimnames = list(from = states, to = states))
  if (nrow(moves)) tmat[cbind(moves$from, moves$to)] <- seq_len(nrow(moves))
  tmat
}
