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
#' @param covariates Optional character vector of baseline (time-fixed) covariate
#'   column names in \code{data} to carry through to the individual-level triples
#'   table (see Value). When \code{data} is an "msdata" object, these are also
#'   passed as the \code{keep} argument of \code{\link{from_msdata}}. Values are
#'   carried through as-is (a warning is issued if a covariate is not constant
#'   within subject); NA covariate values are not flagged here and are left for
#'   downstream modelling functions to handle.
#' @param states Optional state space / ordering. Defaults to sorted observed states.
#' @param absorbing Optional absorbing states; inferred if NULL.
#' @param drop.na If TRUE, drop rows with NA in id/time/state (with a warning)
#'   instead of erroring. Default FALSE.
#' @param check.consecutive If TRUE (default), warn when any subject's times are
#'   not consecutive integers.
#' @return An object of class "msm2data", with (among other components) a
#'   \code{triples} data.table: one row per subject per observed second-order
#'   triple, with columns \code{id, h, j, l, s, d} and, if requested, the
#'   \code{covariates} columns. \code{s} is the time of the destination state
#'   and \code{d} the time already spent in the current state \code{j} at
#'   \code{s - 1} (1 = first time point in \code{j}; counted from the start of
#'   follow-up), usable in \code{\link{P2reg}} for a time-varying or
#'   semi-Markov baseline. The names \code{id, time, state, h, j, l, s, d} are
#'   reserved and cannot be used as covariates. This is the individual-level, "long" analogue of
#'   \code{mstate::msprep()} for the second-order process, and is the table
#'   consumed by covariate-based extensions of the model. The component
#'   \code{pairs} (columns \code{from, to, N}) counts every observed one-step
#'   move, including each subject's first move (which has no preceding state
#'   and therefore no triple); it is the first-order transition structure used
#'   by \code{\link{as_tmat}}.
#' @examples
#' st   <- c("A", "B", "C")                                 # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' d <- prep2(panel, id = "id", time = "time", state = "state")
#' d
#' @export
prep2 <- function(data, id = "id", time = "time", state = "state",
                  covariates = NULL,
                  states = NULL, absorbing = NULL,
                  drop.na = FALSE, check.consecutive = TRUE) {

  ## Accept an mstate "msdata" object directly: convert it to the long
  ## (id, time, state) panel format this function actually works on, keeping
  ## any requested covariates along the way (see from_msdata()).
  if (inherits(data, "msdata")) data <- from_msdata(data, keep = covariates)
  stopifnot(is.data.frame(data))

  ## `cols` maps the package's fixed internal names (id/time/state) to
  ## whatever column names the caller's data actually uses, so the rest of
  ## the function can refer to them uniformly after the rename below.
  cols <- c(id = id, time = time, state = state)
  missing_cols <- cols[!cols %in% names(data)]
  if (length(missing_cols))
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  ## These names are created internally in the `triples` table (see Value),
  ## so a covariate with one of them would be silently overwritten.
  reserved <- intersect(covariates, c("id", "time", "state", "h", "j", "l", "s", "d"))
  if (length(reserved))
    stop("Covariate name(s) reserved for internal columns: ",
         paste(reserved, collapse = ", "), ". Rename them before calling prep2().",
         call. = FALSE)
  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov))
    stop("Covariate column(s) not found in `data`: ",
         paste(missing_cov, collapse = ", "), call. = FALSE)

  ## Copy only the columns actually needed (id/time/state + covariates) into
  ## a fresh data.table, then rename id/time/state to their canonical names.
  DT <- data.table::as.data.table(data)[, .SD, .SDcols = c(unname(cols), covariates)]
  data.table::setnames(DT, old = unname(cols), new = names(cols))

  ## A missing id/time/state makes a row impossible to place in the panel
  ## (which subject? which instant? which state?), so by default it's a hard
  ## error; drop.na = TRUE downgrades this to a warning-and-drop instead.
  bad <- is.na(DT$id) | is.na(DT$time) | is.na(DT$state)
  if (any(bad)) {
    if (!drop.na)
      stop("`id`, `time` and `state` must not contain NA ",
           "(set drop.na = TRUE to drop these rows).", call. = FALSE)
    warning(sprintf("Dropping %d row(s) with NA in id/time/state.", sum(bad)),
            call. = FALSE)
    DT <- DT[!bad]
  }

  ## Fix the state space and its ordering once, up front: everything
  ## downstream (the N/Y counting tables, the P2est() tensors, ckequations()'s
  ## pair-transition matrix) indexes states positionally against this vector,
  ## so it has to be decided before any of that. If the caller didn't supply
  ## one, use every state actually observed, sorted (alphabetically for
  ## character labels, numerically for numeric ones).
  if (is.null(states))
    states <- if (is.factor(DT$state)) levels(droplevels(DT$state))
              else sort(unique(DT$state))
  DT[, state := factor(state, levels = states)]
  ## factor() silently turns any value not in `levels` into NA, so this is
  ## the check that catches a caller-supplied `states` that omits an
  ## observed value (rather than letting it vanish unnoticed).
  if (anyNA(DT$state)) stop("Some observed states are not in `states`.", call. = FALSE)

  n <- data.table::uniqueN(DT$id)
  ## Sorting by (id, time) is required before shift() below, which looks at
  ## *row position* within each id-group to find "1 step back"/"2 steps
  ## back" -- it has no idea what `time` means, so the rows must already be
  ## in chronological order per subject.
  data.table::setorder(DT, id, time)

  if (check.consecutive) {
    ## shift() (below) treats consecutive *rows* as consecutive *time steps*
    ## regardless of the actual gap between their `time` values. That's fine
    ## for a genuinely regular grid, but silently wrong if a subject has
    ## gaps (e.g. observed on days 0, 1, 4): the triple formed from rows 2
    ## and 3 would treat day-1 and day-4 as adjacent. This is a warning, not
    ## an error, because occasional gaps (e.g. a missed visit) are common in
    ## real panel data and the alternative -- silently padding gaps with an
    ## assumed state -- would be a worse and less visible assumption.
    gaps <- DT[, .(ok = all(diff(time) == 1L)), by = id]
    if (!all(gaps$ok))
      warning(sprintf(paste("%d subject(s) have non-consecutive time values;",
                            "triples are formed from consecutive observations."),
                      sum(!gaps$ok)), call. = FALSE)
  }

  ## The heart of the second-order construction: for each subject, look back
  ## 2 rows (h = X_{s-2}) and 1 row (j = X_{s-1}) from the current row's own
  ## state (which becomes l = X_s once renamed below). shift() returns NA for
  ## the first 1-2 rows of each subject (nothing to look back to), which is
  ## exactly how rows with fewer than 3 observations get excluded next.
  ## Time already spent in the current state j at s - 1 (semi-Markov
  ## duration): rowid(rleid(state)) numbers the rows within each run of
  ## identical consecutive states (1, 2, 3, ...), so for the row of time s the
  ## *previous* row's run position is how many consecutive time points the
  ## subject had been in j up to and including s - 1. A subject whose first
  ## observation is already in j is counted from that first observation (the
  ## duration is left-truncated by the start of follow-up).
  DT[, d := data.table::shift(data.table::rowid(data.table::rleid(state)), 1L), by = id]
  DT[, c("h", "j") := .(data.table::shift(state, 2L),
                        data.table::shift(state, 1L)), by = id]
  ## Every observed one-step move (X_{s-1}, X_s), including each subject's
  ## first move, which has no preceding state and so never appears in the
  ## triples: this is the first-order transition structure of the data
  ## (used by as_tmat()).
  pairs <- DT[!is.na(j), .(N = .N), keyby = .(from = j, to = state)]
  trip <- DT[!is.na(h) & !is.na(j)]
  data.table::setnames(trip, "state", "l")
  if (nrow(trip) == 0L)
    stop("No second-order triples found: each subject needs >= 3 observations.",
         call. = FALSE)

  ## N_hjl(s): how many subjects were observed making the h -> j -> l jump
  ## ending at time s (one row of `trip` = one such subject-instant).
  ## Y_hj(s) = sum_l N_hjl(s): how many subjects shared the (h, j) history at
  ## that same instant, i.e. how many were "at risk" of any l. Grouping Y
  ## directly from N (rather than recomputing from `trip`) guarantees the two
  ## tables are perfectly aligned pointwise, which P2est() relies on when it
  ## merges them by (h, j, s).
  N <- trip[, .(N = .N),      keyby = .(h, j, l, s = time)]
  Y <- N[,    .(Y = sum(N)),  keyby = .(h, j, s)]            # Y = N summed over l

  ## Absorbing states aren't declared anywhere in a plain panel, so infer
  ## them: a state is absorbing if it was ever *occupied* (appeared as either
  ## the "current" state j or the destination l of some triple) but was
  ## *never* the origin of an actual j -> l move (j != l). Anything left over
  ## after removing every state seen moving away from itself is, by
  ## elimination, one that -- as far as the data shows -- nobody ever leaves.
  if (is.null(absorbing)) {
    moved    <- unique(trip[as.integer(j) != as.integer(l), j])
    occupied <- union(levels(droplevels(trip$j)), levels(droplevels(trip$l)))
    absorbing <- setdiff(occupied, as.character(moved))
  }

  data.table::setnames(trip, "time", "s")
  data.table::setcolorder(trip, c("id", "h", "j", "l", "s", "d"))
  if (length(covariates)) {
    ## Covariates are documented as baseline (time-fixed) values, so warn
    ## (rather than silently averaging or taking the first observation) if
    ## one actually varies within a subject across their triples -- that
    ## almost always means the wrong column was passed as a covariate.
    varying <- vapply(covariates, function(cv)
      any(trip[, data.table::uniqueN(get(cv), na.rm = TRUE) > 1L, by = id]$V1),
      logical(1))
    if (any(varying))
      warning("Covariate(s) are not constant within subject and are carried ",
              "through per-row rather than enforced as baseline: ",
              paste(covariates[varying], collapse = ", "), call. = FALSE)
  }

  structure(
    list(N = N, Y = Y, triples = trip, pairs = pairs, covariates = covariates,
         states = states, absorbing = as.character(absorbing),
         n = n, ntriples = nrow(trip), time.range = range(DT$time)),
    class = "msm2data"
  )
}

## Compact, human-readable summary of an "msm2data" object; also what gets
## shown automatically when such an object is printed at the console.
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
  cat(sprintf("  covariates      : %s\n",
              if (length(x$covariates)) paste(x$covariates, collapse = ", ") else "none"))
  invisible(x)
}

#' @export
summary.msm2data <- function(object, ...) {
  ## Per (h, j) pair: total at-risk exposure (sum of Y, the RPE denominator)
  ## and the earliest/latest time point it appeared at.
  exposure <- object$Y[, .(total_at_risk = sum(Y),
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
#' @param keep Optional character vector of baseline (time-fixed) covariate
#'   column names in \code{x} to carry through to the output panel, replicated
#'   on every row of a subject. Mirrors the \code{keep} argument of
#'   \code{mstate::msprep()}, so an msdata object built with covariates kept
#'   there converts directly.
#' @param round_fun Discretisation function for transition times. Default rnd().
#'   A state entered on the same rounded time as the next one occupies no
#'   grid point and is dropped, together with its transitions; a warning
#'   reports how many visits were lost.
#' @return A data.table in panel format (id, time, state, and, if requested,
#'   the \code{keep} covariate columns).
#' @examples
#' # A minimal mstate-style transition table (as mstate::msprep() would
#' # produce): subject 1 goes 1 -> 2 at Tstop = 2, then 2 -> 3 at Tstop = 5;
#' # subject 2 goes 1 -> 2 at Tstop = 1, then 2 -> 3 at Tstop = 4.
#' msdat <- data.frame(
#'   id     = c(1, 1, 2, 2),
#'   from   = c(1, 2, 1, 2),
#'   to     = c(2, 3, 2, 3),
#'   Tstart = c(0, 2, 0, 1),
#'   Tstop  = c(2, 5, 1, 4),
#'   status = c(1, 1, 1, 1)
#' )
#' from_msdata(msdat)
#' @export
from_msdata <- function(x, id = "id", from = "from", to = "to",
                        Tstart = "Tstart", Tstop = "Tstop",
                        status = "status", keep = NULL, round_fun = rnd) {
  DT  <- data.table::as.data.table(x)
  ## `map` holds the 6 required column names, in a fixed order; selecting
  ## .SDcols = c(map, keep) and then renaming the first 6 columns
  ## *positionally* to the canonical names is what lets this function accept
  ## any actual column-naming convention the caller's msdata object used.
  map <- c(id, from, to, Tstart, Tstop, status)
  miss <- setdiff(c(map, keep), names(DT))
  if (length(miss))
    stop("msdata missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  DT <- DT[, .SD, .SDcols = c(map, keep)]
  data.table::setnames(DT, old = 1:6, new = c("id", "from", "to", "Tstart", "Tstop", "status"))

  ## Reconstruct one subject's discrete-time trajectory from their mstate
  ## "long format" rows (one row per state they were ever at risk of moving
  ## to, from/to/Tstart/Tstop marking the interval, status = 1 marking which
  ## rows are the transitions that actually happened).
  build <- function(s, gid) {
    data.table::setorder(s, Tstart, Tstop)
    ev     <- s[status == 1][order(Tstop)]            # realized jumps, in order
    ## The trajectory is: start in `from` of the very first interval, then
    ## step through each realized jump's destination, in chronological order.
    states <- c(s$from[1L], ev$to)                    # start state, then destinations
    ## `rt` holds the (rounded, discretised) instant each of those states was
    ## *entered*: the start time, then each jump's own Tstop.
    rt     <- round_fun(c(s$Tstart[1L], ev$Tstop))    # their absolute (rounded) times
    ## A state entered on the same rounded time as the next one never
    ## occupies a grid point, so findInterval() below skips it.
    lost <<- lost + sum(diff(rt) == 0)
    ## A daily grid covering the subject's whole observed window, from their
    ## first Tstart to the very last Tstop on record (which may be a
    ## censoring time with no further transition, not just the last jump).
    grid   <- seq.int(rt[1L], round_fun(max(s$Tstop)))
    ## findInterval(grid, rt) is a right-continuous step function: for each
    ## day in `grid` it returns the index of the *last* entry time in `rt`
    ## that is <= that day, i.e. "which state was most recently entered by
    ## this day" -- exactly last-observation-carried-forward.
    out <- data.table::data.table(id = gid, time = grid,
                                  state = states[findInterval(grid, rt)])
    ## Baseline covariates don't vary over time, so just take the first row's
    ## value and replicate it onto every day of this subject's panel.
    if (length(keep)) out[, (keep) := as.list(s[1L, keep, with = FALSE])]  # baseline, replicated
    out
  }
  ## Apply build() to each subject's rows in turn (.SD = that subject's rows,
  ## .BY[[1]] = their id), then keep only the output columns actually wanted.
  out_cols <- c("id", "time", "state", keep)
  lost <- 0L
  out <- DT[, build(.SD, .BY[[1]]), by = id][, ..out_cols]
  if (lost)
    warning(sprintf(paste("%d state visit(s) were entered on the same discretised time",
                          "as the next state and were dropped; the transitions into and",
                          "out of them are lost. Use a finer time unit or another `round_fun`."),
                    lost), call. = FALSE)
  out
}
