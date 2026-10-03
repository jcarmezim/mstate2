#' Round half away from zero (the rounding of the paper's DIVINE analysis)
#'
#' The function used in the code of the methods paper (Najera-Zuloaga, Besalu
#' and Gomez Melis, 2025) to turn the sojourn times of the DIVINE cohort into
#' whole days. Unlike \code{\link{round}}, which rounds halves to the even
#' number (\code{round(0.5) = 0}, \code{round(2.5) = 2}), \code{rnd} always
#' rounds halves away from zero (\code{rnd(0.5) = 1}, \code{rnd(2.5) = 3}).
#' Only this convention reproduces Table 2 of the paper.
#'
#' @param x numeric vector.
#' @return \code{trunc(x + sign(x) * 0.5)}.
#' @examples
#' rnd(c(-2.5, -0.4, 0.4, 0.5, 2.5))
#' round(c(-2.5, -0.4, 0.4, 0.5, 2.5))   # base R: halves to even
#' @export
rnd <- function(x) trunc(x + sign(x) * 0.5)

#' Convert sojourn-time data to discrete-time panel format
#'
#' Turns one row per subject, with the time spent in each transient state and
#' an indicator of the absorbing state that ended follow-up, into one row per
#' subject and time unit with the state occupied. This is the discretisation
#' used in the illustration of the methods paper: each sojourn is rounded with
#' \code{\link{rnd}} and the states are visited in the order given in
#' \code{segments}.
#'
#' @param data A data frame, one row per subject.
#' @param id Name of the subject id column.
#' @param segments Named character vector mapping transient state labels to their
#'   duration columns, in visiting order, e.g.
#'   \code{c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")}.
#' @param absorbing Named character vector mapping absorbing state labels to their
#'   0/1 indicator columns, e.g. \code{c(Disch = "disch.s", Death = "death.s")}.
#' @param round_fun Discretisation function for durations. Default \code{rnd}.
#'   Each sojourn is rounded once and then expanded into that many time units.
#'   A visit with a positive duration that rounds to 0 (with \code{rnd}, a
#'   duration below 0.5) is dropped, together with the transitions into and
#'   out of it; a warning reports how many visits were lost, by state.
#' @return A tibble in panel format with columns \code{id}, \code{time}
#'   (0, 1, 2, ... within each subject) and \code{state} (character). A
#'   subject whose follow-up ended in an absorbing state gets one final row
#'   with that state; a subject with no indicator equal to 1 (still under
#'   observation) has no absorbing row.
#' @examples
#' d <- data.frame(id = 1:2, t.a = c(2, 1), t.b = c(1, 3),
#'                 disch.s = c(1, 0), death.s = c(0, 1))
#' sojourn_to_panel(d, id = "id",
#'                  segments = c(A = "t.a", B = "t.b"),
#'                  absorbing = c(Disch = "disch.s", Death = "death.s"))
#' @export
sojourn_to_panel <- function(data, id, segments, absorbing, round_fun = rnd) {
  ## 1. Check the input: every column named in `id`, `segments` and
  ##    `absorbing` must exist in `data`.
  stopifnot(is.data.frame(data))
  seg_lab <- names(segments); seg_col <- unname(segments)
  abs_lab <- names(absorbing); abs_col <- unname(absorbing)
  miss <- setdiff(c(id, seg_col, abs_col), names(data))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  ## 2. Number the subjects (.row) so that the output keeps the row order of
  ##    `data`, and rename the id column to `id`.
  subjects <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::select(".row", id = dplyr::all_of(id), dplyr::all_of(c(seg_col, abs_col)))

  ## 3. Transient visits: one row per subject and transient state, with the
  ##    recorded duration (`raw`), its position in the visiting order (`.pos`)
  ##    and the rounded number of time units (`dur`).
  visits <- subjects |>
    dplyr::select(".row", "id", dplyr::all_of(seg_col)) |>
    tidyr::pivot_longer(dplyr::all_of(seg_col), names_to = "state", values_to = "raw") |>
    dplyr::mutate(.pos  = match(state, seg_col),          # visiting order
                  state = seg_lab[.pos],                  # column name -> state label
                  dur   = round_fun(raw))

  ## 4. Visits with a positive recorded duration that round to 0 time units
  ##    disappear from the panel, and with them the transitions into and out
  ##    of that state (A -> B -> C becomes A -> C). Count them by state and
  ##    warn, so that the loss is never silent.
  lost <- visits |>
    dplyr::filter(!is.na(raw), raw > 0, is.na(dur) | dur <= 0) |>
    dplyr::count(state = factor(state, levels = seg_lab), .drop = TRUE)
  if (nrow(lost))
    warning(sprintf(paste("%d visit(s) with a positive duration rounded to 0 time units",
                          "and were dropped (%s); the transitions into and out of them",
                          "are lost. Use a finer time unit or another `round_fun`."),
                    sum(lost$n), paste0(lost$state, ": ", lost$n, collapse = ", ")),
            call. = FALSE)

  ## 5. Expand each visit into one row per time unit spent there. A missing
  ##    or negative duration means "never visited this state" and, like a
  ##    zero duration, contributes no rows.
  transient_rows <- visits |>
    dplyr::mutate(dur = as.integer(dplyr::if_else(is.na(dur) | dur < 0, 0, dur))) |>
    tidyr::uncount(dur) |>
    dplyr::select(".row", "id", "state", ".pos")

  ## 6. Absorbing state: the first indicator column equal to 1 says how
  ##    follow-up ended and adds one final row. Subjects with no indicator
  ##    equal to 1 (still under observation) get no absorbing row.
  absorbing_rows <- subjects |>
    dplyr::select(".row", "id", dplyr::all_of(abs_col)) |>
    tidyr::pivot_longer(dplyr::all_of(abs_col), names_to = "state", values_to = "ind") |>
    dplyr::mutate(.pos  = length(seg_lab) + match(state, abs_col),  # after every transient state
                  state = abs_lab[match(state, abs_col)]) |>
    dplyr::filter(as.numeric(ind) == 1) |>
    dplyr::filter(.pos == min(.pos), .by = ".row") |>
    dplyr::select(".row", "id", "state", ".pos")

  ## 7. Put each subject's rows in visiting order and number the time units
  ##    0, 1, 2, ... within the subject.
  dplyr::bind_rows(transient_rows, absorbing_rows) |>
    dplyr::arrange(.row, .pos) |>
    dplyr::mutate(time = dplyr::row_number() - 1L, .by = ".row") |>
    dplyr::select("id", "time", "state")
}
