#' Prepare multistate data for the analysis
#'
#' Turns multistate data into the discrete-time panel that \code{\link{prep2}} needs (one row per subject and time unit, with the state occupied), and reports every record it had to drop or change.
#'
#' The times are divided by \code{unit} and rounded with \code{\link{rnd}}. A subject occupies each state from the unit in which it was entered until the unit before the next state; follow-up ends at the entry into an absorbing state or, otherwise, at the last recorded time (entry or censoring). The absorbing states are the states with no exit in \code{trans}, the states of \code{outcome}, or the states nobody is seen leaving. All the columns that are not used (\code{id}, the columns in \code{Surv()}, \code{durations}, \code{outcome}, or the columns of \code{msdata}) are kept as baseline covariates.
#'
#' Every record that is dropped or changed is listed in the \code{issues} table of the result, and a single warning gives the counts:
#' \describe{
#'   \item{\code{missing}}{a state with status 1 but no time, or a time without a valid status (wide);}
#'   \item{\code{rounded_to_zero}}{a stay with a positive duration rounds to 0 units (sojourn);}
#'   \item{\code{same_unit}}{two different states fall in the same time unit: the state entered last (or an absorbing state) is kept, so the other one occupies no unit and its transitions are lost;}
#'   \item{\code{after_absorbing}}{a state entered after the absorbing state;}
#'   \item{\code{before_start}}{a negative time;}
#'   \item{\code{not_allowed}}{a transition not allowed by \code{trans} (kept);}
#'   \item{\code{no_data}}{a subject with no usable record.}
#' }
#'
#' @param data A data frame, or an \code{msdata} object from \code{mstate::msprep()}.
#' @param states Wide data: a named list with the states in order, each one \code{Surv(time, status)} or \code{NULL} (the initial state), e.g. \code{list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))}. Sojourn data: optionally, the order of the states (a character vector; default the order of \code{durations} and \code{outcome}).
#' @param trans Optional allowed transitions, written as text, \code{c("healthy -> ill -> dead", "healthy -> dead")} (a chain \code{"A -> B -> C"} gives A -> B and B -> C), or a matrix from \code{mstate::transMat()}. It also gives the absorbing states (no exit). An \code{msdata} object carries its own.
#' @param durations Sojourn data: named character vector mapping each transient state to the column with the time spent in it, in visiting order, e.g. \code{c(healthy = "t_healthy", ill = "t_ill")}.
#' @param outcome Sojourn data: named character vector mapping each absorbing state to its 0/1 indicator column, e.g. \code{c(dead = "dead")}.
#' @param id Name of the subject id column (if missing, the rows are numbered).
#' @param keep Optional names of the baseline covariate columns to carry into the panel (default: all the columns that are not used).
#' @param unit Length of one time unit in the units of the data, e.g. \code{30.4375} to go from days to months. Default 1.
#' @return An object of class \code{"msm2prep"}, which \code{\link{prep2}} accepts directly: a list with
#' \describe{
#'   \item{\code{panel}}{tibble \code{(id, time, state, covariates...)}, one row per subject and time unit, \code{state} a factor with levels \code{states};}
#'   \item{\code{states}, \code{absorbing}, \code{trans}}{state space, absorbing states and the allowed-transition matrix (logical, or \code{NULL});}
#'   \item{\code{transitions}}{tibble \code{(from, to, n, allowed)} with the number of each observed transition;}
#'   \item{\code{subjects}}{tibble with one row per subject: first and last time, entry and exit state, number of rows and \code{status} (\code{"absorbed"} or \code{"censored"});}
#'   \item{\code{issues}}{tibble \code{(id, issue, state, time, detail)} with every record dropped or changed;}
#'   \item{\code{settings}}{the layout and the unit.}
#' }
#' @section Differences with mstate::msprep():
#' \code{msprep()} needs the transition matrix from state numbers, the time and status columns in two parallel vectors (\code{NA} for the initial state) and the covariates listed in \code{keep}, and returns one row per subject and possible transition (the counting-process format of the Cox model). \code{msprep2()} gives each state its own \code{Surv(time, status)}, reads the transitions as text, keeps the other columns as covariates, and returns the discrete-time panel of second-order models, with a report of the data.
#' @seealso \code{\link{prep2}}, \code{\link{rnd}}
#' @references
#' de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package for the analysis of competing risks and multi-state models. \emph{Journal of Statistical Software}, 38(7), 1-30.
#' @examples
#' # Illness-death model: healthy -> ill -> dead and healthy -> dead
#'
#' ## Sojourn: the structure of the DIVINE data
#' sojourn <- data.frame(id = 1:4, t_healthy = c(3, 5, 2, 7), t_ill = c(4, 0, 8, 0),
#'                       dead = c(1, 1, 0, 0))
#' p1 <- msprep2(sojourn, durations = c(healthy = "t_healthy", ill = "t_ill"),
#'               outcome = c(dead = "dead"))
#' p1
#'
#' ## Wide data
#' wide <- data.frame(id = 1:4, ill_time = c(3, 5, 2, 6), ill_status = c(1, 0, 1, 0),
#'                    dead_time = c(7, 5, 9, 6), dead_status = c(1, 1, 0, 0),
#'                    age = c(60, 72, 55, 49))
#' p2 <- msprep2(wide,
#'               states = list(healthy = NULL,
#'                             ill     = Surv(ill_time, ill_status),
#'                             dead    = Surv(dead_time, dead_status)),
#'               trans  = c("healthy -> ill -> dead", "healthy -> dead"))
#' p2$panel                       # `age` is kept as a covariate
#' summary(p2)
#'
#' @export
msprep2 <- function(data, states = NULL, trans = NULL, durations = NULL, outcome = NULL, id = "id", keep = NULL, unit = 1) {

  # Check the arguments
  stopifnot(is.data.frame(data))
  if (!nrow(data)) stop("`data` has no rows.", call. = FALSE)
  if (!is.numeric(unit) || length(unit) != 1L || !(unit > 0))
    stop("`unit` must be a positive number, e.g. 30.4375 to go from days to months.", call. = FALSE)

  # The layout, from the arguments given
  states_expr <- substitute(states)
  format <- if (inherits(data, "msdata")) "msdata"
            else if (!is.null(durations)) "sojourn"
            else if (is.call(states_expr) && identical(states_expr[[1]], as.name("list"))) "wide"
            else stop("Give wide data with `states = list(initial = NULL, state = Surv(time, status), ...)`, ",
                      "sojourn data with `durations` and `outcome`, or an msdata object from mstate::msprep().", call. = FALSE)
  if (format == "sojourn" && is.null(outcome))
    stop("Sojourn data need `outcome`: the 0/1 indicators of the absorbing states, e.g. c(dead = \"dead\").", call. = FALSE)

  # Wide data: evaluate each Surv(time, status) on the columns of `data`.
  if (format == "wide") {
    sv <- .surv_states(states_expr, data, parent.frame())
    states <- sv$states
  } else if (!is.null(states) && !is.character(states)) {
    stop("With ", format, " data, `states` can only give the order of the states (a character vector).", call. = FALSE)
  }

  # One row per subject without an id column: number the rows.
  if (format != "msdata" && !id %in% names(data)) {
    data <- tibble::as_tibble(data)
    data[[id]] <- seq_len(nrow(data))
  }

  # Columns used to build the states, the other ones are the covariates.
  used <- switch(format,
                 wide    = c(id, sv$used),
                 sojourn = c(id, unname(durations), unname(outcome)),
                 msdata  = c(id, "from", "to", "trans", "Tstart", "Tstop", "time", "status"))
  if (is.null(keep)) keep <- setdiff(names(data), used)
  miss <- setdiff(c(id, keep, if (format == "sojourn") c(unname(durations), unname(outcome))), names(data))
  if (length(miss)) stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  clash <- intersect(keep, c("time", "state"))
  if (length(clash))
    stop("Covariate(s) named ", paste(clash, collapse = " and "), " would clash with the columns of the panel: rename them, or choose the covariates with `keep`.", call. = FALSE)

  # Allowed transitions: text or a transMat matrix
  if (format == "msdata" && is.null(trans)) trans <- attr(data, "trans")
  if (format == "sojourn" && is.null(states)) states <- c(names(durations), names(outcome))
  allowed <- .allowed_matrix(trans, states)
  if (is.null(states)) states <- rownames(allowed)
  if (!is.null(allowed)) allowed <- allowed[states, states, drop = FALSE]
  absorbing <- if (!is.null(allowed)) states[rowSums(allowed) == 0]
               else if (format == "sojourn") names(outcome)

  # Turn the data into the same intermediate form: `ev`, one row per recorded state (id, step, state, with the time before rounding and the order of the record), and `subj`, one row per subject with the last time unit of follow-up. `issues` collects the records dropped on the way
  built <- switch(format,
                  wide    = .from_wide(data, id, sv, unit),
                  sojourn = .from_sojourn(data, id, durations, outcome, unit),
                  msdata  = .from_msdata(data, id, rownames(.allowed_matrix(attr(data, "trans"))), unit))

  # Every recorded state must belong to the state space
  observed <- unique(built$ev$state)
  unknown <- setdiff(observed, states)
  if (length(unknown))
    stop("State(s) not in `states`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)

  # Absorbing states, if not given by `trans` or `outcome`
  if (is.null(absorbing)) {
    left <- built$ev |>
      dplyr::arrange(id, .t, .ord) |>
      dplyr::mutate(.nxt = dplyr::lead(state), .by = "id") |>
      dplyr::filter(!is.na(.nxt), .nxt != state) |>
      dplyr::pull("state") |>
      unique()
    absorbing <- setdiff(observed, left)
  }

  # Build the panel: stop at the absorbing state, resolve several states in one time unit and expand each state over the units it occupies.
  res <- .build_panel(built$ev, built$subj, absorbing)
  issues <- dplyr::bind_rows(built$issues, res$issues)
  panel <- res$panel

  # Transitions observed and, with `trans`, whether each one is allowed; those not allowed are kept and reported.
  transitions <- res$transitions |>
    dplyr::count(from, to, name = "n") |>
    dplyr::mutate(allowed = if (is.null(allowed)) NA else allowed[cbind(from, to)])
  if (!is.null(allowed)) {
    bad <- dplyr::filter(res$transitions, !allowed[cbind(from, to)])
    issues <- dplyr::bind_rows(issues, bad |>
      dplyr::transmute(id, issue = "not_allowed", state = to, time = NA_character_,
                       detail = paste(from, "->", to, "at unit", step)))
  }

  # Subject summary: first and last unit, entry and exit state, and whether follow-up ended in an absorbing state.
  subjects <- panel |>
    dplyr::summarise(first = min(time), last = max(time),
                     entry = dplyr::first(state), exit = dplyr::last(state), rows = dplyr::n(),
                     .by = "id") |>
    dplyr::mutate(status = dplyr::if_else(exit %in% absorbing, "absorbed", "censored"))

  # Baseline covariates, one value per subject, carried into the panel
  if (length(keep)) {
    covs <- data |>
      as.data.frame() |>
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
         settings = list(format = format, unit = unit)),
    class = "msm2prep")
}

# Short report of a prepared data set, shown when it is printed
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

# Detailed report: the table of observed transitions (from x to, as mstate::events()), the follow-up per subject and the issues by type
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
  cat("\nTime units observed per subject (rows of the panel), by status:\n")
  print(follow_up)
  cat("\nRecords dropped or changed:\n")
  if (nrow(issues)) print(issues) else cat("  none\n")
  invisible(list(transitions = from_to, follow_up = follow_up, issues = issues))
}

# Wide data: the initial state at time 0, then each state with status 1 at its time
.from_wide <- function(data, id, sv, unit) {
  ids <- data[[id]]
  if (anyDuplicated(ids)) stop("Wide data need one row per subject; `id` has duplicates.", call. = FALSE)

  # One row per subject and state with its time and status, as given by Surv().
  ev <- purrr::imap(sv$values, \(v, st) tibble::tibble(.first = seq_along(ids), id = ids, state = st,
                                                       time = v$time, status = v$status)) |>
    dplyr::bind_rows() |>
    dplyr::mutate(.ord = match(state, sv$states))

  # A state marked as visited but without a time, or a time without a valid status, cannot be used
  issues <- dplyr::bind_rows(
    ev |>
      dplyr::filter(is.na(time), !is.na(status), status == 1) |>
      dplyr::transmute(id, issue = "missing", state, time = NA_character_, detail = "visited (status 1) but without a time"),
    ev |>
      dplyr::filter(!is.na(time), is.na(status)) |>
      dplyr::transmute(id, issue = "missing", state, time = as.character(time),
                       detail = "time recorded but status missing: not taken as a visit"))

  # Follow-up ends at the largest recorded time (entry or censoring)
  subj <- ev |>
    dplyr::summarise(end = suppressWarnings(max(time, na.rm = TRUE)), .by = c(".first", "id")) |>
    dplyr::mutate(end_step = dplyr::if_else(is.finite(end), rnd(end / unit), NA_real_)) |>
    dplyr::select("id", "end_step", ".first")

  # Records: the initial state at time 0 and the visited states (status 1) at their times
  ev <- dplyr::bind_rows(
    tibble::tibble(.first = seq_along(ids), id = ids, state = sv$initial, time = 0, .ord = 0L),
    dplyr::filter(ev, !is.na(time), !is.na(status), status == 1)) |>
    dplyr::mutate(.t = time / unit, step = rnd(.t), tval = as.character(time)) |>
    dplyr::select("id", "step", "state", ".ord", ".t", "tval")
  list(ev = ev, subj = subj, issues = issues)
}

# Sojourn data (as DIVINE): durations per transient state in visiting order and 0/1 indicators of the absorbing state
.from_sojourn <- function(data, id, durations, outcome, unit) {
  raw <- data |>
    tibble::as_tibble() |>
    dplyr::mutate(.first = dplyr::row_number())
  if (anyDuplicated(raw[[id]])) stop("Sojourn data need one row per subject; `id` has duplicates.", call. = FALSE)

  # Rounded number of units of every visit. Missing or negative durations mean "not visited".
  visits <- raw |>
    dplyr::select(".first", id = dplyr::all_of(id), dplyr::all_of(unname(durations))) |>
    tidyr::pivot_longer(dplyr::all_of(unname(durations)), names_to = ".col", values_to = "raw") |>
    dplyr::mutate(.ord = match(.col, durations), state = names(durations)[.ord],
                  dur = rnd(raw / unit), dur = dplyr::if_else(is.na(dur) | dur < 0, 0, dur))

  # Visits with a positive duration that round to 0 units are lost
  issues <- visits |>
    dplyr::filter(!is.na(raw), raw > 0, dur == 0) |>
    dplyr::transmute(id, issue = "rounded_to_zero", state, time = NA_character_,
                     detail = paste("duration", raw, "rounds to 0 units"))

  # Entry unit of each visit = units spent in the previous visits.
  ev <- visits |>
    dplyr::filter(dur > 0) |>
    dplyr::arrange(.first, .ord) |>
    dplyr::mutate(step = cumsum(dur) - dur, .by = ".first")

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

  # Without an absorbing state, follow-up ends with the last unit of the last visit
  subj <- total |>
    dplyr::mutate(end_step = dplyr::if_else(.first %in% absorbed$.first, NA_real_, total - 1)) |>
    dplyr::select("id", "end_step", ".first")
  ev <- dplyr::bind_rows(dplyr::select(ev, "id", "step", "state", ".ord"),
                         dplyr::select(absorbed, "id", "step", "state", ".ord")) |>
    dplyr::mutate(.t = step, tval = NA_character_)
  list(ev = ev, subj = subj, issues = issues)
}

# msdata (mstate::msprep())
.from_msdata <- function(data, id, state_names, unit) {
  if (is.null(state_names)) stop("The msdata object has no transition matrix (attribute \"trans\").", call. = FALSE)
  ms <- data |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::rename(id = dplyr::all_of(id))

  # Entries: the transitions that happened (status 1) at their stop time. The first state is the origin state of the subject's first row.
  entries <- ms |>
    dplyr::filter(status == 1) |>
    dplyr::transmute(id, time = Tstop, state = state_names[to], .ord = 1L)
  firsts <- ms |>
    dplyr::filter(Tstart == min(Tstart), .by = "id") |>
    dplyr::distinct(id, .keep_all = TRUE) |>
    dplyr::transmute(id, time = Tstart, state = state_names[from], .ord = 0L)

  # Time from each subject's first start. Follow-up ends at the last stop time.
  subj <- ms |>
    dplyr::summarise(origin = min(Tstart), end = max(Tstop), .by = "id") |>
    dplyr::mutate(.first = dplyr::row_number(), end_step = rnd((end - origin) / unit))
  ev <- dplyr::bind_rows(firsts, entries) |>
    dplyr::left_join(dplyr::select(subj, "id", "origin"), by = "id") |>
    dplyr::mutate(.t = (time - origin) / unit, step = rnd(.t), tval = as.character(time)) |>
    dplyr::select("id", "step", "state", ".ord", ".t", "tval")
  list(ev = ev, subj = dplyr::select(subj, "id", "end_step", ".first"), issues = NULL)
}

# From the records (id, step, state, .ord, .t, tval) and the subjects (id, end_step, .first) to the panel, logging every record dropped or changed
.build_panel <- function(ev, subj, absorbing) {
  log_rows <- function(x, issue, detail) {
    dplyr::transmute(x, id, issue = issue, state, time = tval, detail = detail)
  }

  # Records with a negative time
  ev <- dplyr::left_join(ev, subj, by = "id")
  before <- ev$step < 0
  issues <- log_rows(ev[before, ], "before_start", "negative time")
  ev <- ev[!before, ]

  # Stop at the first absorbing state, in time order: states entered later are dropped. Done before the ties below, so that a record later in the same unit cannot replace the absorbing state.
  ev <- ev |>
    dplyr::arrange(.first, .t, .ord) |>
    dplyr::mutate(.abs = cumsum(cumsum(state %in% absorbing)), .by = "id")
  issues <- dplyr::bind_rows(issues, log_rows(ev[ev$.abs > 1, ], "after_absorbing", "after the entry into an absorbing state"))
  ev <- ev[ev$.abs <= 1, ]

  # Several records in the same unit: keep the last one in time. An absorbing state, if present, is always kept.
  ev <- ev |>
    dplyr::mutate(.k = dplyr::row_number(), .n = dplyr::n(), .has = any(state %in% absorbing), .by = c("id", "step"))
  keep <- dplyr::if_else(ev$.has, ev$state %in% absorbing, ev$.k == ev$.n)
  kept <- ev[keep, ]
  dropped <- ev[!keep, ]

  # Last unit with a record of each subject: a repeated record of the same state extends follow-up even though it is merged below.
  seen <- kept |>
    dplyr::summarise(.seen = max(step), .by = "id")

  # Merge consecutive records of the same state (not transitions).
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
                                                paste0("same time unit as ", lost$kept_state, " (unit ", lost$step, ")")))
  }

  # Subjects with no record left
  gone <- dplyr::filter(subj, !id %in% kept$id)
  issues <- dplyr::bind_rows(issues, dplyr::transmute(gone, id, issue = "no_data", state = NA_character_,
                                                      time = NA_character_, detail = "no usable record"))

  # Expand every state over the units it occupies: until the unit before the next state. The last state until the end of follow-up (one unit if absorbing; the last record if the end is unknown).
  kept <- kept |>
    dplyr::left_join(seen, by = "id") |>
    dplyr::mutate(.next = dplyr::lead(step), .by = "id") |>
    dplyr::mutate(.len = dplyr::case_when(!is.na(.next) ~ .next - step,
                                          state %in% absorbing ~ 1,
                                          is.na(end_step) ~ .seen - step + 1,
                                          TRUE ~ pmax(end_step, .seen) - step + 1))
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


# Logical matrix of allowed transitions, from text such as c("A -> B -> C", "A -> C") or from an mstate::transMat() (or logical, or 0/1) matrix. Every state named in `trans` must be in `states`, if given.
.allowed_matrix <- function(trans, states = NULL) {
  if (is.null(trans)) return(NULL)
  if (is.character(trans) && is.null(dim(trans))) {
    # Text: each "A -> B -> C" gives A -> B and B -> C.
    steps <- purrr::map(strsplit(trans, "->", fixed = TRUE), trimws)
    if (any(lengths(steps) < 2L) || any(unlist(steps) == ""))
      stop("Write each transition as \"from -> to\", e.g. c(\"healthy -> ill -> dead\", \"healthy -> dead\").", call. = FALSE)
    pairs <- purrr::map(steps, \(x) cbind(utils::head(x, -1), x[-1])) |>
      purrr::reduce(rbind)
    nm <- unique(c(states, t(pairs)))
    ok <- matrix(FALSE, length(nm), length(nm), dimnames = list(nm, nm))
    ok[pairs] <- TRUE
  } else {
    # A matrix: mstate::transMat() (NA = not allowed), logical or 0/1.
    trans <- as.matrix(trans)
    nm <- if (!is.null(rownames(trans))) rownames(trans) else colnames(trans)
    if (nrow(trans) != ncol(trans) || is.null(nm))
      stop("`trans` must be text such as \"A -> B\" or a square matrix with the states as row names.", call. = FALSE)
    ok <- !is.na(trans) & trans != 0
    dimnames(ok) <- list(nm, nm)
    extra <- setdiff(states, nm)
    if (length(extra)) {
      ok <- rbind(cbind(ok, matrix(FALSE, nrow(ok), length(extra))), matrix(FALSE, length(extra), length(nm) + length(extra)))
      dimnames(ok) <- list(c(nm, extra), c(nm, extra))
    }
  }
  diag(ok) <- FALSE
  if (!is.null(states) && length(setdiff(rownames(ok), states)))
    stop("State(s) in `trans` but not in `states`: ", paste(setdiff(rownames(ok), states), collapse = ", "), ".", call. = FALSE)
  ok
}

# Evaluates the list of states of wide data
.surv_states <- function(expr, data, env) {
  els <- as.list(expr)[-1]
  if (!length(els) || is.null(names(els)) || any(names(els) == "") || anyDuplicated(names(els)))
    stop("Name every state in `states`: list(healthy = NULL, ill = Surv(ill_time, ill_status), ...).", call. = FALSE)

  # The state without columns is the initial state of everybody.
  initial <- names(els)[purrr::map_lgl(els, is.null)]
  if (length(initial) != 1L)
    stop("Give exactly one state as NULL in `states`: the initial state, e.g. list(healthy = NULL, ...).", call. = FALSE)

  # The other states
  surv <- purrr::compact(els)
  is_surv <- function(e) is.call(e) && (identical(e[[1]], as.name("Surv")) || identical(e[[1]], quote(survival::Surv)))
  bad <- names(surv)[!purrr::map_lgl(surv, is_surv)]
  if (length(bad))
    stop("Write state(s) ", paste(bad, collapse = ", "), " as Surv(time, status), or NULL for the initial state.", call. = FALSE)

  # The columns named inside Surv() must exist in `data`
  vars <- unique(unlist(purrr::map(surv, all.vars)))
  miss <- setdiff(vars, c(names(data), ls(env, all.names = TRUE)))
  if (length(miss)) stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)

  # Evaluate each Surv() on the columns of `data`.
  mask <- c(as.list(data), list(Surv = .surv_in_states))
  values <- lapply(surv, \(e) {
    e[[1]] <- as.name("Surv")
    eval(e, mask, env)
  })
  list(states = names(els), initial = initial, values = values, used = intersect(vars, names(data)))
}

# Surv() inside the list of states
.surv_in_states <- function(time, event, ...) {
  if (length(list(...))) stop("Use Surv(time, status) inside `states`: one time and one status.", call. = FALSE)
  if (missing(event)) stop("Surv() inside `states` needs a status: Surv(time, status).", call. = FALSE)
  if (!is.numeric(time)) stop("The time of Surv() must be numeric (e.g. days since the start).", call. = FALSE)
  if (is.numeric(event) && any(event == 2, na.rm = TRUE) && all(event %in% c(1, 2, NA)))
    warning("A status coded 1/2 is read as in survival::Surv(): 1 = censored, 2 = event.", call. = FALSE)
  # A status written as text or factor: Surv() refuses it, so say how to state which value is the event.
  if (is.character(event) || is.factor(event))
    stop("The status of Surv() must be 0/1, 1/2 or TRUE/FALSE, not text (values: ",
         paste(utils::head(unique(stats::na.omit(as.character(event))), 4), collapse = ", "),
         "). Say which value is the event, e.g. Surv(time, status == \"", stats::na.omit(as.character(event))[1], "\").",
         call. = FALSE)
  checked <- survival::Surv(time, event)
  list(time = as.numeric(time), status = unclass(checked)[, "status"])
}
