#' Prepare wide multistate data for second-order analysis (msprep-style)
#'
#' The second-order counterpart of \code{mstate::msprep()}: takes data in the
#' same \emph{wide} format (one row per subject, one time and one status column
#' per state) and returns the discrete-time panel \code{(id, time, state)} that
#' \code{\link{prep2}} needs, without going through the long \code{msdata}
#' format. The arguments \code{time}, \code{status}, \code{data}, \code{start},
#' \code{id} and \code{keep} have the same meaning as in \code{msprep()}.
#'
#' For each subject, the visited states are those with \code{status == 1};
#' their \code{time} is the time the state was entered. They are ordered by
#' entry time (ties broken by the order of the states), the subject starts in
#' \code{start$state} at \code{start$time}, and follow-up ends at the largest
#' time recorded for the subject (an entry or a censoring time), or at the
#' entry into an absorbing state, as in \code{msprep}. Times are discretised with \code{round_fun}, a daily grid
#' is built from the start to the end of follow-up, and each day gets the most
#' recently entered state (last observation carried forward). A state entered
#' and left within the same time unit therefore does not occupy a row, as with
#' \code{\link{from_msdata}}.
#'
#' Unlike \code{msprep()}, no transition matrix is required: the transition
#' structure is whatever the data show (\code{\link{as_tmat}} recovers it).
#' If \code{trans} is supplied, it only provides the state names and is used
#' to warn about observed transitions it does not allow.
#'
#' @param time Character vector, one entry per state, with the names of the
#'   columns of \code{data} holding the entry (or censoring) times; \code{NA}
#'   for the initial state(s), as in \code{msprep()}.
#' @param status Character vector of the same length with the 0/1 status
#'   columns (\code{NA} for the initial state(s)).
#' @param data Data frame in wide format.
#' @param states State labels, one per element of \code{time}. Defaults to the
#'   row names of \code{trans}, or to the names of \code{time}, or to
#'   \code{1, 2, ...}.
#' @param trans Optional \code{mstate} transition matrix (see Details).
#' @param start Optional list with elements \code{state} (index or label) and
#'   \code{time}, each of length 1 or \code{nrow(data)}. Default: every subject
#'   starts in the first state at time 0.
#' @param id Name of the id column in \code{data}; if missing, ids
#'   \code{1, ..., n} are created.
#' @param keep Character vector of baseline covariate columns to carry along.
#' @param absorbing Absorbing states (labels). Default: the states with no
#'   outgoing transition in \code{trans}. Follow-up stops when one of them is
#'   entered; without \code{trans} or \code{absorbing} a warning is given and
#'   no truncation is done.
#' @param round_fun Discretisation function for times. Default \code{\link{rnd}}.
#'   A state entered on the same rounded time as the next one occupies no
#'   grid point and is dropped, together with its transitions; a warning
#'   reports how many visits were lost.
#' @return A \code{data.table} with columns \code{id}, \code{time}, \code{state}
#'   (character) and the \code{keep} columns.
#' @references
#' de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package
#' for the analysis of competing risks and multi-state models. \emph{Journal
#' of Statistical Software}, 38(7), 1-30.
#' @examples
#' # illness-death data in the wide format of mstate::msprep()
#' wide <- data.frame(id = 1:3,
#'                    ill.t = c(2, 5, 4), ill.s = c(1, 0, 1),
#'                    dth.t = c(6, 5, 4.5), dth.s = c(1, 0, 1),
#'                    age = c(60, 72, 55))
#' msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
#'         data = wide, states = c("healthy", "ill", "dead"), absorbing = "dead",
#'         id = "id", keep = "age")
#' @export
msprep2 <- function(time, status, data, states = NULL, trans = NULL, start = NULL,
                    id = NULL, keep = NULL, absorbing = NULL, round_fun = rnd) {
  stopifnot(is.data.frame(data))
  S <- length(time)
  if (length(status) != S)
    stop("`time` and `status` must have the same length (one entry per state).", call. = FALSE)
  if (is.null(states))
    states <- if (!is.null(trans)) rownames(trans)
              else if (!is.null(names(time))) names(time)
              else as.character(seq_len(S))
  if (length(states) != S)
    stop("`states` must have one label per element of `time`.", call. = FALSE)
  ## absorbing states: given, or the states with no outgoing transition in trans
  if (is.null(absorbing) && !is.null(trans))
    absorbing <- states[rowSums(!is.na(trans)) == 0L]
  if (is.null(absorbing))
    warning("No `trans` or `absorbing` given: follow-up is not truncated at ",
            "absorbing states.", call. = FALSE)
  abs_idx <- match(absorbing, states)
  if (anyNA(abs_idx)) stop("`absorbing` state(s) not found in `states`.", call. = FALSE)
  cols <- c(stats::na.omit(time), stats::na.omit(status), id, keep)
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  n <- nrow(data)
  ids <- if (is.null(id)) seq_len(n) else data[[id]]

  st_state <- if (is.null(start)) rep(1L, n) else rep_len(start$state, n)
  if (!is.numeric(st_state)) st_state <- match(st_state, states)
  if (anyNA(st_state)) stop("`start$state` not found in `states`.", call. = FALSE)
  st_time <- if (is.null(start)) rep(0, n) else rep_len(start$time, n)

  tm <- vapply(time,   function(cl) if (is.na(cl)) rep(NA_real_, n) else as.numeric(data[[cl]]), numeric(n))
  sm <- vapply(status, function(cl) if (is.na(cl)) rep(NA_real_, n) else as.numeric(data[[cl]]), numeric(n))
  tm <- matrix(tm, n, S); sm <- matrix(sm, n, S)

  out <- vector("list", n)
  lost <- 0L   # states entered on the same discretised time as the next one
  for (i in seq_len(n)) {
    vis <- which(!is.na(sm[i, ]) & sm[i, ] == 1)
    vis <- vis[order(tm[i, vis], vis)]                  # by entry time, then state order
    ## the path stops at the first absorbing state entered (later entries or
    ## censoring times are ignored, as msprep() does)
    hit <- which(vis %in% abs_idx)
    if (length(hit)) vis <- vis[seq_len(hit[1L])]
    seq_st <- c(st_state[i], vis)
    rt <- round_fun(c(st_time[i], tm[i, vis]))
    lost <- lost + sum(diff(rt) == 0)
    end <- if (length(hit)) max(rt)
           else round_fun(max(c(st_time[i], tm[i, ]), na.rm = TRUE))
    grid <- seq.int(rt[1L], end)
    out[[i]] <- data.table::data.table(id = ids[i], time = grid,
                                       state = states[seq_st[findInterval(grid, rt)]])
    if (length(keep)) out[[i]][, (keep) := as.list(data[i, keep, drop = FALSE])]
  }
  panel <- data.table::rbindlist(out)
  if (lost)
    warning(sprintf(paste("%d state visit(s) were entered on the same discretised time",
                          "as the next state and were dropped; the transitions into and",
                          "out of them are lost. Use a finer time unit or another `round_fun`."),
                    lost), call. = FALSE)

  if (!is.null(trans)) {
    mv <- panel[, .(from = state[-.N], to = state[-1L]), by = id][from != to]
    bad <- unique(mv[is.na(trans[cbind(match(from, states), match(to, states))]), .(from, to)])
    if (nrow(bad))
      warning("Observed transition(s) not allowed by `trans`: ",
              paste(bad$from, bad$to, sep = " -> ", collapse = ", "), call. = FALSE)
  }
  panel
}
