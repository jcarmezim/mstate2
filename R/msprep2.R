#' Prepare raw multistate data for the analysis
#'
#' One entry point that turns multistate data as they are usually collected into the discrete-time panel that \code{\link{prep2}} needs (one row per subject and time unit, with the state occupied), and reports every record it had to drop or change. It plays the role of \code{mstate::msprep()}, but works from four input layouts, understands dates, discretises time, cleans the data and never changes them silently.
#'
#' \strong{Input layouts} (\code{format}):
#' \describe{
#'   \item{\code{"events"}}{One row per subject and recorded state, with the time at which the state was entered or observed (\code{time}, \code{state}). Repeated records of the same state are merged, so the rows can be state changes or repeated observations of the same patient.}
#'   \item{\code{"wide"}}{One row per subject, with one time column per state (\code{times}), as in \code{mstate::msprep()}: the time the state was entered, or \code{NA} if it was never visited. With \code{status}, a state is visited when its status is 1, and the times with status 0 are censoring times, as in \code{msprep()}.}
#'   \item{\code{"sojourn"}}{One row per subject, with the time spent in each transient state (\code{durations}, in visiting order) and 0/1 indicators of the absorbing state that ended follow-up (\code{outcome}), as in the DIVINE data. Gives the same panel as \code{\link{sojourn_to_panel}}.}
#'   \item{\code{"msdata"}}{An object of class \code{"msdata"} made by \code{mstate::msprep()}; its transition matrix is used as \code{trans}.}
#' }
#' With \code{format = "auto"} (default) the layout is taken from the arguments given: an \code{msdata} object, \code{durations}, \code{times}, or else events.
#'
#' \strong{Time.} Times can be numbers or dates (\code{Date}, \code{POSIXct}, or text such as \code{"2020-03-15"} or \code{"15/03/2020"}). Each time is measured from the subject's origin (\code{start}: a column, a common value, or by default the subject's first record) in units of \code{unit} (a number in the units of the times, or \code{"hour"}, \code{"day"}, \code{"week"}, \code{"month"} or \code{"year"} for dates) and discretised with \code{round_fun}. A subject occupies each state from the time unit in which it was entered until the unit before the next state; follow-up ends at the entry into an absorbing state, at the end of follow-up (\code{end}), or, without \code{end}, at the last record (the subject is then censored in that unit; a repeated record of the same state counts as a record).
#'
#' \strong{Cleaning.} Every record that is dropped or changed is listed in the \code{issues} table of the result, and a single warning gives the counts:
#' \describe{
#'   \item{\code{missing}}{an id, time or state is missing (events);}
#'   \item{\code{before_start}}{the record is before the subject's origin;}
#'   \item{\code{after_end}}{the record is after the end of follow-up;}
#'   \item{\code{same_unit}}{two different states fall in the same time unit; \code{ties = "last"} keeps the state entered last (an absorbing state is always kept), so the other one occupies no unit and its transitions are lost;}
#'   \item{\code{rounded_to_zero}}{a visit with a positive duration rounds to 0 time units (sojourn layout);}
#'   \item{\code{after_absorbing}}{the record is after the entry into an absorbing state;}
#'   \item{\code{not_allowed}}{the transition is not allowed by \code{trans} (kept, unless \code{check = "error"}, which stops);}
#'   \item{\code{no_data}}{the subject has no usable record left.}
#' }
#'
#' @param data A data frame with the raw data (or an \code{msdata} object).
#' @param id Name of the subject id column.
#' @param format Input layout: \code{"auto"} (default), \code{"events"}, \code{"wide"}, \code{"sojourn"} or \code{"msdata"} (see Details).
#' @param time,state Events layout: names of the time and state columns.
#' @param times Wide layout: named character vector mapping each state to its time column, e.g. \code{c(SP = "date_sp", IMV = "date_imv", Death = "date_death")}.
#' @param status Wide layout: optional named character vector mapping each state (names as in \code{times}) to its 0/1 status column.
#' @param durations Sojourn layout: named character vector mapping each transient state to its duration column, in visiting order.
#' @param outcome Sojourn layout: named character vector mapping each absorbing state to its 0/1 indicator column.
#' @param initial Optional initial state, entered at the origin: a column of \code{data} or a single state label. Useful when the records only list the changes after the start (e.g. the state at admission is in a separate column).
#' @param start Origin of time for each subject: a column of \code{data} (e.g. the admission date) or a single value. Default: the subject's first record (events, dates in the wide layout) or 0 (numeric wide layout).
#' @param end End of follow-up for subjects not absorbed: a column of \code{data} or a single value (e.g. the date the data were extracted). Default: none (events), or the largest recorded time of the subject when \code{status} is given (wide layout, as in \code{msprep()}).
#' @param unit Length of one time unit: a number in the units of the times, or a name (\code{"hour"}, \code{"day"}, \code{"week"}, \code{"month"}, \code{"year"}) when the times are dates. Default 1 (one day for dates).
#' @param round_fun Discretisation of the elapsed times. Default \code{\link{rnd}}.
#' @param recode Optional named character vector to relabel the recorded states, \code{c(old = "new")}, e.g. \code{c("1" = "NSP", "2" = "SP")}. Applied to \code{state} and \code{initial}.
#' @param states Optional state space and order. Default: the row names of \code{trans}, the order of \code{times} or \code{durations} and \code{outcome}, or the sorted observed states.
#' @param absorbing Optional absorbing states. Default: the states with no allowed exit in \code{trans}, the states of \code{outcome}, or the states nobody is seen leaving.
#' @param trans Optional matrix of allowed transitions, with the states as row and column names: an \code{mstate::transMat()} matrix (\code{NA} = not allowed) or a logical or 0/1 matrix.
#' @param ties Which state to keep when several fall in the same time unit: \code{"last"} (default, the state at the end of the unit) or \code{"first"}.
#' @param check What to do with transitions not allowed by \code{trans}: \code{"warn"} (default; they are kept and listed) or \code{"error"}.
#' @param keep Optional names of baseline covariate columns to carry into the panel.
#' @return An object of class \code{"msm2prep"}, which \code{\link{prep2}} accepts directly: a list with
#' \describe{
#'   \item{\code{panel}}{tibble \code{(id, time, state, keep...)}, one row per subject and time unit, \code{time} counted from the origin, \code{state} a factor with levels \code{states};}
#'   \item{\code{states}, \code{absorbing}, \code{trans}}{state space, absorbing states and the allowed-transition matrix (logical, or \code{NULL});}
#'   \item{\code{transitions}}{tibble \code{(from, to, n, allowed)} with the number of each observed transition;}
#'   \item{\code{subjects}}{tibble with one row per subject: first and last time, entry and exit state, \code{status} (\code{"absorbed"} or \code{"censored"}) and number of rows;}
#'   \item{\code{issues}}{tibble \code{(id, issue, state, time, detail)} with every record dropped or changed;}
#'   \item{\code{settings}}{the layout, unit and tie rule used.}
#' }
#' @section Differences with mstate::msprep():
#' \code{msprep()} takes only the wide layout, needs the transition matrix, works with numeric times and returns one row per subject and possible transition (the counting-process format of the Cox model). \code{msprep2()} returns the discrete-time panel of second-order models, accepts four layouts and dates, infers the states and absorbing states when no matrix is given, orders the states by the recorded times, and checks and reports the data instead of failing or changing them silently.
#' @seealso \code{\link{prep2}}, \code{\link{sojourn_to_panel}}
#' @references
#' de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package for the analysis of competing risks and multi-state models. \emph{Journal of Statistical Software}, 38(7), 1-30.
#' @examples
#' # Events layout: one row per state change, with dates
#' raw <- data.frame(
#'   id    = c(1, 1, 1, 2, 2, 3, 3, 3),
#'   date  = c("2020-03-01", "2020-03-03", "2020-03-06",
#'             "2020-03-02", "2020-03-09",
#'             "2020-03-05", "2020-03-06", "2020-03-12"),
#'   state = c("NSP", "SP", "Disch", "SP", "Death", "NSP", "SP", "SP"))
#' x <- msprep2(raw, time = "date", state = "state", absorbing = c("Disch", "Death"))
#' x
#' x$panel
#' prep2(x)
#'
#' # Wide layout, as in mstate::msprep(): entry time and status per state
#' wide <- data.frame(id = 1:3,
#'                    ill.t = c(2, 5, 4), ill.s = c(1, 0, 1),
#'                    dth.t = c(6, 5, 4.5), dth.s = c(1, 0, 1))
#' msprep2(wide, times = c(ill = "ill.t", dead = "dth.t"),
#'         status = c(ill = "ill.s", dead = "dth.s"), initial = "healthy",
#'         states = c("healthy", "ill", "dead"), absorbing = "dead")$panel
#' @export
msprep2 <- function(data, id = "id", format = c("auto", "events", "wide", "sojourn", "msdata"),
                    time = "time", state = "state", times = NULL, status = NULL,
                    durations = NULL, outcome = NULL, initial = NULL, start = NULL, end = NULL,
                    unit = 1, round_fun = rnd, recode = NULL, states = NULL, absorbing = NULL,
                    trans = NULL, ties = c("last", "first"), check = c("warn", "error"), keep = NULL) {

  # Check the arguments and choose the input layout from the arguments given.
  stopifnot(is.data.frame(data))
  format <- match.arg(format)
  ties <- match.arg(ties)
  check <- match.arg(check)
  if (format == "auto")
    format <- if (inherits(data, "msdata")) "msdata"
              else if (!is.null(durations)) "sojourn"
              else if (!is.null(times)) "wide"
              else "events"

  # Allowed transitions. An msdata object carries its own matrix. The matrix also gives the default state space and absorbing states (no allowed exit).
  if (format == "msdata" && is.null(trans)) trans <- attr(data, "trans")
  allowed <- .allowed_matrix(trans)
  if (is.null(states) && !is.null(allowed)) states <- rownames(allowed)
  if (is.null(absorbing) && !is.null(allowed)) absorbing <- rownames(allowed)[rowSums(allowed) == 0]

  # Columns that must exist in `data` for the chosen layout.
  needed <- c(id, keep, switch(format,
                               events  = c(time, state),
                               wide    = c(unname(times), unname(status)),
                               sojourn = c(unname(durations), unname(outcome)),
                               msdata  = c("from", "to", "Tstart", "Tstop", "status")))
  needed <- c(needed, .column_arg(initial, data), .column_arg(start, data), .column_arg(end, data))
  miss <- setdiff(needed, names(data))
  if (length(miss)) stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)

  # Turn each layout into the same intermediate form: `ev`, one row per recorded state (id, step, state, with the original time and the order of the record), and `subj`, one row per subject with the last time unit of follow-up (end_step, NA if unknown). `issues` collects the records dropped on the way.
  built <- switch(format,
                  events  = .from_events(data, id, time, state, initial, start, end, unit, round_fun, recode),
                  wide    = .from_wide(data, id, times, status, initial, start, end, unit, round_fun, recode),
                  sojourn = .from_sojourn(data, id, durations, outcome, round_fun),
                  msdata  = .from_msdata(data, id, rownames(allowed), start, end, unit, round_fun))

  # State space: every recorded state must belong to it.
  observed <- unique(built$ev$state)
  if (is.null(states))
    states <- switch(format,
                     wide    = union(c(.value_arg(initial, data), names(times)), observed),
                     sojourn = c(names(durations), names(outcome)),
                     c(sort(setdiff(observed, absorbing)), absorbing))
  unknown <- setdiff(observed, states)
  if (length(unknown))
    stop("State(s) not in `states`: ", paste(unknown, collapse = ", "), " (use `recode` to relabel them).", call. = FALSE)
  if (is.null(absorbing) && format == "sojourn") absorbing <- names(outcome)

  # Build the panel: drop records out of follow-up, stop at the absorbing state, resolve several states in one time unit, merge repeated records and expand each state over the units it occupies.
  res <- .build_panel(built$ev, built$subj, absorbing, ties)
  issues <- dplyr::bind_rows(built$issues, res$issues)
  panel <- res$panel

  # Absorbing states, if still unknown: the states nobody is seen leaving.
  if (is.null(absorbing)) {
    left <- res$transitions |>
      dplyr::pull("from") |>
      unique()
    absorbing <- setdiff(unique(as.character(panel$state)), left)
  }

  # Transitions observed and, with a transition matrix, whether each one is allowed.
  transitions <- res$transitions |>
    dplyr::count(from, to, name = "n") |>
    dplyr::mutate(allowed = if (is.null(allowed)) NA else allowed[cbind(from, to)])
  if (!is.null(allowed)) {
    bad <- res$transitions |>
      dplyr::filter(!allowed[cbind(from, to)])
    if (nrow(bad)) {
      msg <- sprintf("%d transition(s) not allowed by `trans`: %s.", nrow(bad),
                     paste(unique(paste(bad$from, "->", bad$to)), collapse = ", "))
      if (check == "error") stop(msg, call. = FALSE)
      issues <- dplyr::bind_rows(issues, bad |>
        dplyr::transmute(id, issue = "not_allowed", state = to, time = NA_character_,
                         detail = paste(from, "->", to, "at unit", step)))
    }
  }

  # Subject summary: first and last unit, entry and exit state, and whether follow-up ended in an absorbing state.
  subjects <- panel |>
    dplyr::summarise(first = min(time), last = max(time),
                     entry = dplyr::first(state), exit = dplyr::last(state), rows = dplyr::n(),
                     .by = "id") |>
    dplyr::mutate(status = dplyr::if_else(exit %in% absorbing, "absorbed", "censored"))

  # Baseline covariates (`keep`), one value per subject, carried into the panel.
  if (length(keep)) {
    covs <- data |>
      tibble::as_tibble() |>
      dplyr::select(id = dplyr::all_of(id), dplyr::all_of(keep)) |>
      dplyr::distinct(id, .keep_all = TRUE)
    panel <- dplyr::left_join(panel, covs, by = "id")
  }

  # One warning with the counts of records dropped or changed; the details are in `issues`.
  if (nrow(issues)) {
    counts <- dplyr::count(issues, issue)
    warning(sprintf("%d record(s) were dropped or changed while building the panel (%s); see the `issues` table of the result.",
                    nrow(issues), paste0(counts$issue, ": ", counts$n, collapse = ", ")), call. = FALSE)
  }

  structure(
    list(panel = dplyr::mutate(panel, state = factor(state, levels = states)),
         states = states,
         absorbing = as.character(absorbing),
         trans = allowed,
         transitions = transitions,
         subjects = subjects,
         issues = issues,
         settings = list(format = format, unit = unit, ties = ties)),
    class = "msm2prep")
}

# Short report of a prepared data set, shown when it is printed.
#' @export
print.msm2prep <- function(x, ...) {
  st <- x$subjects
  cat("<msm2prep>  discrete-time panel ready for prep2()\n")
  cat(sprintf("  layout          : %s\n", x$settings$format))
  cat(sprintf("  subjects        : %d (%d absorbed, %d censored)\n", nrow(st), sum(st$status == "absorbed"), sum(st$status == "censored")))
  cat(sprintf("  panel rows      : %d (time %d - %d)\n", nrow(x$panel), min(x$panel$time), max(x$panel$time)))
  cat(sprintf("  states (%d)      : %s\n", length(x$states), paste(x$states, collapse = ", ")))
  cat(sprintf("  absorbing       : %s\n", if (length(x$absorbing)) paste(x$absorbing, collapse = ", ") else "none"))
  cat(sprintf("  transitions     : %d types, %d in total%s\n", nrow(x$transitions), sum(x$transitions$n),
              if (is.null(x$trans)) "" else sprintf(" (%d not allowed by `trans`)", sum(!x$transitions$allowed))))
  if (nrow(x$issues)) {
    counts <- dplyr::count(x$issues, issue)
    cat(sprintf("  issues          : %d record(s) dropped or changed (%s); see x$issues\n",
                nrow(x$issues), paste0(counts$issue, ": ", counts$n, collapse = ", ")))
  } else {
    cat("  issues          : none\n")
  }
  invisible(x)
}

# Detailed report: the table of observed transitions (from x to, as mstate::events()), the issues by type and the follow-up per subject. Returned invisibly as a list.
#' @export
summary.msm2prep <- function(object, ...) {
  from_to <- object$transitions |>
    dplyr::mutate(from = factor(from, levels = object$states), to = factor(to, levels = object$states)) |>
    stats::xtabs(formula = n ~ from + to)
  issues <- dplyr::count(object$issues, issue)
  follow_up <- object$subjects |>
    dplyr::summarise(subjects = dplyr::n(), min = min(last - first + 1), median = stats::median(last - first + 1),
                     max = max(last - first + 1), .by = "status")
  cat("<msm2prep summary>\n\nObserved transitions (from rows to columns):\n")
  print(from_to)
  cat("\nTime units of follow-up by status:\n")
  print(follow_up)
  cat("\nRecords dropped or changed:\n")
  if (nrow(issues)) print(issues) else cat("  none\n")
  invisible(list(transitions = from_to, follow_up = follow_up, issues = issues))
}

# --- internal: layout-specific readers ------------------------------------------
# Each one returns list(ev, subj, issues): `ev` has one row per recorded state (id, step, state, tval = original time as text, .t = elapsed time before rounding, .ord = order of the record), `subj` has one row per subject (id, end_step, .first = order of appearance).

# Events layout: one row per recorded state.
.from_events <- function(data, id, time, state, initial, start, end, unit, round_fun, recode) {

  # Keep the needed columns with canonical names, convert the times and relabel the states.
  raw <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.ord = dplyr::row_number())
  ev <- raw |>
    dplyr::select(id = dplyr::all_of(id), time = dplyr::all_of(time), state = dplyr::all_of(state), ".ord") |>
    dplyr::mutate(time = .as_time(time), state = .recode(state, recode))

  # Records with a missing id, time or state cannot be placed: drop and log them.
  miss <- is.na(ev$id) | is.na(ev$time) | is.na(ev$state)
  issues <- ev[miss, ] |>
    dplyr::transmute(id, issue = "missing", state, time = format(time), detail = "missing id, time or state")
  ev <- ev[!miss, ]

  # Subjects in order of appearance, with their origin and end of follow-up (column, common value, or default).
  subj <- tibble::tibble(id = raw[[id]],
                         origin = .subject_value(start, data, raw),
                         end = .subject_value(end, data, raw),
                         initial = .subject_value(initial, data, raw)) |>
    dplyr::filter(!is.na(id)) |>
    dplyr::summarise(origin = dplyr::first(stats::na.omit(origin)),
                     end = dplyr::first(stats::na.omit(end)),
                     initial = dplyr::first(stats::na.omit(initial)),
                     .by = "id") |>
    dplyr::mutate(.first = dplyr::row_number())
  first_time <- ev |>
    dplyr::summarise(t0 = min(time), .by = "id")
  subj <- subj |>
    dplyr::left_join(first_time, by = "id") |>
    dplyr::mutate(origin = .as_time(origin), end = .as_time(end))
  subj$origin <- .coalesce_time(subj$origin, subj$t0)

  # The initial state, if given, is entered at the origin, before any record of the same unit.
  ev <- .add_initial(ev, subj, recode)
  .finish_reader(ev, subj, unit, round_fun, issues)
}

# Wide layout: one row per subject, one time column (and optionally one status column) per state.
.from_wide <- function(data, id, times, status, initial, start, end, unit, round_fun, recode) {
  if (is.null(names(times)) || any(names(times) == ""))
    stop("`times` must be named by state, e.g. c(SP = \"date_sp\").", call. = FALSE)
  raw <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.row = dplyr::row_number())
  if (anyDuplicated(raw[[id]])) stop("The wide layout needs one row per subject; `id` has duplicates.", call. = FALSE)

  # One row per subject and state with its time and, optionally, its status.
  ev <- raw |>
    dplyr::select(".row", id = dplyr::all_of(id), dplyr::all_of(unname(times))) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(unname(times)), .as_time)) |>
    tidyr::pivot_longer(dplyr::all_of(unname(times)), names_to = ".col", values_to = "time") |>
    dplyr::mutate(.ord = match(.col, times), state = names(times)[.ord])
  if (!is.null(status)) {
    # Status columns in the order of `times` (kept in a variable: inside the pivot, `status` is the new column).
    status_cols <- unname(status[names(times)])
    stat <- raw |>
      dplyr::select(".row", dplyr::all_of(status_cols)) |>
      tidyr::pivot_longer(-".row", names_to = ".scol", values_to = "status") |>
      dplyr::mutate(.ord = match(.scol, status_cols))
    ev <- dplyr::left_join(ev, dplyr::select(stat, ".row", ".ord", "status"), by = c(".row", ".ord"))
  } else {
    ev <- dplyr::mutate(ev, status = as.numeric(!is.na(time)))
  }

  # Follow-up ends at `end` or, with status columns, at the largest recorded time (entry or censoring), as in msprep().
  last_time <- ev |>
    dplyr::filter(!is.na(time)) |>
    dplyr::summarise(t_max = max(time), .by = ".row")
  subj <- tibble::tibble(.row = raw$.row, id = raw[[id]],
                         origin = .as_time(.subject_value(start, data, raw)),
                         end = .as_time(.subject_value(end, data, raw)),
                         initial = .subject_value(initial, data, raw),
                         .first = raw$.row) |>
    dplyr::left_join(last_time, by = ".row")
  if (is.null(end) && !is.null(status)) subj$end <- subj$t_max

  # Only visited states are records; their times give the order of the visits.
  ev <- ev |>
    dplyr::filter(!is.na(time), status == 1) |>
    dplyr::select("id", "time", "state", ".ord")
  t0 <- ev |>
    dplyr::summarise(t0 = min(time), .by = "id")
  subj <- dplyr::left_join(subj, t0, by = "id")
  default_origin <- if (inherits(ev$time, c("Date", "POSIXt"))) subj$t0 else rep(0, nrow(subj))
  subj$origin <- .coalesce_time(subj$origin, default_origin)

  ev <- .add_initial(ev, subj, recode)
  .finish_reader(ev, subj, unit, round_fun, NULL)
}

# Sojourn layout (as DIVINE): durations per transient state in visiting order and 0/1 indicators of the absorbing state.
.from_sojourn <- function(data, id, durations, outcome, round_fun) {
  raw <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.first = dplyr::row_number())
  if (anyDuplicated(raw[[id]])) stop("The sojourn layout needs one row per subject; `id` has duplicates.", call. = FALSE)

  # Rounded number of units of every visit; missing or negative durations mean "not visited".
  visits <- raw |>
    dplyr::select(".first", id = dplyr::all_of(id), dplyr::all_of(unname(durations))) |>
    tidyr::pivot_longer(dplyr::all_of(unname(durations)), names_to = ".col", values_to = "raw") |>
    dplyr::mutate(.ord = match(.col, durations), state = names(durations)[.ord],
                  dur = round_fun(raw), dur = dplyr::if_else(is.na(dur) | dur < 0, 0, dur))

  # Visits with a positive duration that round to 0 units are lost: log them.
  issues <- visits |>
    dplyr::filter(!is.na(raw), raw > 0, dur == 0) |>
    dplyr::transmute(id, issue = "rounded_to_zero", state, time = NA_character_,
                     detail = paste("duration", raw, "rounds to 0 units"))

  # Entry unit of each visit = units spent in the previous visits.
  ev <- visits |>
    dplyr::filter(dur > 0) |>
    dplyr::arrange(.first, .ord) |>
    dplyr::mutate(step = cumsum(dur) - dur, .last = cumsum(dur) - 1, .by = ".first")

  # The absorbing state (first indicator equal to 1) is entered after the last visit.
  total <- visits |>
    dplyr::summarise(total = sum(dur), .by = c(".first", "id"))
  absorbed <- raw |>
    dplyr::select(".first", dplyr::all_of(unname(outcome))) |>
    tidyr::pivot_longer(-".first", names_to = ".col", values_to = "ind") |>
    dplyr::filter(as.numeric(ind) == 1) |>
    dplyr::mutate(.ord = length(durations) + match(.col, outcome), state = names(outcome)[match(.col, outcome)]) |>
    dplyr::filter(.ord == min(.ord), .by = ".first") |>
    dplyr::left_join(total, by = ".first") |>
    dplyr::mutate(step = total)

  # Without an absorbing state, follow-up ends with the last unit of the last visit.
  subj <- total |>
    dplyr::mutate(end_step = dplyr::if_else(.first %in% absorbed$.first, NA_real_, total - 1)) |>
    dplyr::select("id", "end_step", ".first")
  ev <- dplyr::bind_rows(dplyr::select(ev, "id", "step", "state", ".ord"),
                         dplyr::select(absorbed, "id", "step", "state", ".ord")) |>
    dplyr::mutate(.t = step, tval = NA_character_)
  list(ev = ev, subj = subj, issues = issues)
}

# msdata layout (mstate::msprep()): one row per subject and possible transition.
.from_msdata <- function(data, id, state_names, start, end, unit, round_fun) {
  if (is.null(state_names)) stop("The msdata layout needs its transition matrix (`trans`).", call. = FALSE)
  ms <- data |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::rename(id = dplyr::all_of(id))

  # Entries: the transitions that happened (status 1) at their stop time; the first state is the origin state of the subject's first row.
  entries <- ms |>
    dplyr::filter(status == 1) |>
    dplyr::transmute(id, time = Tstop, state = state_names[to], .ord = 1L)
  firsts <- ms |>
    dplyr::filter(Tstart == min(Tstart), .by = "id") |>
    dplyr::distinct(id, .keep_all = TRUE) |>
    dplyr::transmute(id, time = Tstart, state = state_names[from], .ord = 0L)
  subj <- ms |>
    dplyr::summarise(origin = min(Tstart), end = max(Tstop), .by = "id") |>
    dplyr::mutate(.first = dplyr::row_number(), initial = NA)
  ev <- dplyr::bind_rows(firsts, entries)
  .finish_reader(ev, subj, unit, round_fun, NULL)
}

# --- internal: common steps of the readers ----------------------------------------

# Elapsed time from the origin in units, its rounded value (step) and the end of follow-up in units.
.finish_reader <- function(ev, subj, unit, round_fun, issues) {
  ev <- ev |>
    dplyr::left_join(dplyr::select(subj, "id", "origin"), by = "id") |>
    dplyr::mutate(.t = .elapsed(time, origin, unit), step = round_fun(.t), tval = format(time)) |>
    dplyr::select("id", "step", "state", ".ord", ".t", "tval")
  subj <- subj |>
    dplyr::mutate(end_step = round_fun(.elapsed(end, origin, unit))) |>
    dplyr::select("id", "end_step", ".first")
  list(ev = ev, subj = subj, issues = issues)
}

# Adds the initial state of each subject (if given) as a record at the origin, ordered before the other records of that unit.
.add_initial <- function(ev, subj, recode) {
  init <- subj |>
    dplyr::filter(!is.na(initial)) |>
    dplyr::transmute(id, time = origin, state = .recode(initial, recode), .ord = 0L)
  if (!nrow(init)) return(ev)
  dplyr::bind_rows(init, ev)
}

# --- internal: the panel builder ---------------------------------------------------
# From the records (id, step, state, .ord, .t, tval) and the subjects (id, end_step, .first) to the panel, logging every record dropped or changed.
.build_panel <- function(ev, subj, absorbing, ties) {
  log_rows <- function(x, issue, detail) {
    dplyr::transmute(x, id, issue = issue, state, time = tval, detail = detail)
  }

  # Records before the origin or after the end of follow-up.
  ev <- dplyr::left_join(ev, subj, by = "id")
  before <- ev$step < 0
  after <- !is.na(ev$end_step) & ev$step > ev$end_step
  issues <- dplyr::bind_rows(log_rows(ev[before, ], "before_start", "before the origin"),
                             log_rows(ev[after & !before, ], "after_end", "after the end of follow-up"))
  ev <- ev[!before & !after, ]

  # Stop at the first absorbing state, in time order: later records are errors (or records after death) and are dropped. Done before the ties below, so that a record later in the same unit cannot replace the absorbing state.
  ev <- ev |>
    dplyr::arrange(.first, .t, .ord) |>
    dplyr::mutate(.abs = cumsum(cumsum(state %in% absorbing)), .by = "id")
  issues <- dplyr::bind_rows(issues, log_rows(ev[ev$.abs > 1, ], "after_absorbing", "after the entry into an absorbing state"))
  ev <- ev[ev$.abs <= 1, ]

  # Several records in the same unit: keep the last (or first) one in time; an absorbing state, if present, is always kept.
  ev <- ev |>
    dplyr::mutate(.k = dplyr::row_number(), .n = dplyr::n(), .by = c("id", "step"))
  keep <- if (ties == "last") ev$.k == ev$.n else ev$.k == 1L
  has_abs <- ev |>
    dplyr::mutate(.has = any(state %in% absorbing), .by = c("id", "step")) |>
    dplyr::pull(".has")
  keep <- dplyr::if_else(has_abs, ev$state %in% absorbing, keep)
  kept <- ev[keep, ]
  dropped <- ev[!keep, ]

  # Last unit with a record of each subject: a repeated observation of the same state extends follow-up even though it is merged below.
  seen <- kept |>
    dplyr::summarise(.seen = max(step), .by = "id")

  # Merge repeated records of the same state (repeated observations, not transitions).
  kept <- kept |>
    dplyr::filter(is.na(dplyr::lag(state)) | state != dplyr::lag(state), .by = "id")

  # A dropped record is a lost visit only if its state differs from the state kept in that unit and from the state before it.
  if (nrow(dropped)) {
    ctx <- kept |>
      dplyr::mutate(prev = dplyr::lag(state), .by = "id") |>
      dplyr::select("id", "step", kept_state = "state", "prev")
    lost <- dropped |>
      dplyr::left_join(ctx, by = c("id", "step")) |>
      dplyr::filter(state != kept_state, is.na(prev) | state != prev)
    issues <- dplyr::bind_rows(issues, log_rows(lost, "same_unit",
                                                paste("same time unit as", lost$kept_state, "(unit", lost$step, ")")))
  }

  # Subjects with no record left.
  gone <- dplyr::filter(subj, !id %in% kept$id)
  issues <- dplyr::bind_rows(issues, dplyr::transmute(gone, id, issue = "no_data", state = NA_character_,
                                                      time = NA_character_, detail = "no usable record"))

  # Expand every state over the units it occupies: until the unit before the next state; the last state until the end of follow-up (one unit if absorbing; the last record if the end is unknown).
  kept <- kept |>
    dplyr::left_join(seen, by = "id") |>
    dplyr::mutate(.next = dplyr::lead(step), .by = "id") |>
    dplyr::mutate(.len = dplyr::case_when(!is.na(.next) ~ .next - step,
                                          state %in% absorbing ~ 1,
                                          is.na(end_step) ~ .seen - step + 1,
                                          TRUE ~ end_step - step + 1))
  panel <- kept |>
    dplyr::arrange(.first, step) |>
    tidyr::uncount(as.integer(.len), .id = ".i") |>
    dplyr::mutate(time = as.integer(step + .i - 1)) |>
    dplyr::select("id", "time", "state")

  # Transitions between consecutive states of each subject.
  transitions <- kept |>
    dplyr::arrange(.first, step) |>
    dplyr::mutate(from = dplyr::lag(state), .by = "id") |>
    dplyr::filter(!is.na(from)) |>
    dplyr::transmute(id, from, to = state, step)

  list(panel = panel, issues = issues, transitions = transitions)
}

# --- internal: small helpers ---------------------------------------------------------

# Logical matrix of allowed transitions (diagonal excluded) from an mstate transMat, a logical or a 0/1 matrix.
.allowed_matrix <- function(trans) {
  if (is.null(trans)) return(NULL)
  trans <- as.matrix(trans)
  nm <- rownames(trans) %||% colnames(trans)
  if (nrow(trans) != ncol(trans) || is.null(nm))
    stop("`trans` must be a square matrix with the states as row (or column) names.", call. = FALSE)
  ok <- !is.na(trans) & trans != 0
  dimnames(ok) <- list(nm, nm)
  diag(ok) <- FALSE
  ok
}

# Times as numbers or dates; text is read as a date in the usual formats.
.as_time <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  d <- as.Date(x, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d"), optional = TRUE)
  if (any(!is.na(x) & is.na(d)))
    stop("Some times could not be read as dates (use the format 2020-03-15 or 15/03/2020, or convert them with as.Date()).", call. = FALSE)
  d
}

# Elapsed time from the origin, in units of `unit`.
.elapsed <- function(t, origin, unit) {
  if (all(is.na(t))) return(rep(NA_real_, length(t)))
  if (inherits(t, c("Date", "POSIXt"))) {
    days <- c(hour = 1 / 24, day = 1, week = 7, month = 365.25 / 12, year = 365.25)
    if (is.character(unit)) {
      key <- sub("s$", "", unit)
      if (!key %in% names(days)) stop("`unit` must be a number or one of: hour, day, week, month, year.", call. = FALSE)
      unit <- days[[key]]
    }
    as.numeric(difftime(t, origin, units = "days")) / unit
  } else {
    if (is.character(unit)) stop("A named `unit` ('day', 'week', ...) needs dates as times.", call. = FALSE)
    (as.numeric(t) - as.numeric(origin)) / unit
  }
}

# Is the argument the name of a column of `data`? Returns that name (or NULL).
.column_arg <- function(x, data) {
  if (is.character(x) && length(x) == 1L && x %in% names(data)) x else NULL
}

# The value of a column argument, if it is a single value rather than a column.
.value_arg <- function(x, data) {
  if (is.null(x) || !is.null(.column_arg(x, data))) NULL else x
}

# Per-row value of an argument that can be a column, a single value or NULL (NA).
.subject_value <- function(x, data, raw) {
  if (is.null(x)) return(rep(NA, nrow(raw)))
  if (!is.null(.column_arg(x, data))) return(raw[[x]])
  rep(x, nrow(raw))
}

# Use `x` where it is known and `y` elsewhere, keeping the date class.
.coalesce_time <- function(x, y) {
  if (all(is.na(x))) return(y)
  out <- x
  out[is.na(x)] <- y[is.na(x)]
  out
}

# Relabel recorded states with c(old = "new").
.recode <- function(x, recode) {
  x <- as.character(x)
  if (is.null(recode)) return(x)
  hit <- x %in% names(recode)
  x[hit] <- recode[x[hit]]
  x
}
