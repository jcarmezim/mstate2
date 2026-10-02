#' Round half up (the convention used in the DIVINE analysis)
#' @param x numeric vector.
#' @return \code{trunc(x + sign(x) * 0.5)}.
#' @examples
#' rnd(c(-2.5, -0.4, 0.4, 2.5))
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
#'   Each sojourn is rounded once and then expanded into that many time units.
#'   A visit with a positive duration that rounds to 0 (with \code{rnd}, a
#'   duration below 0.5) is dropped, together with the transitions into and
#'   out of it; a warning reports how many visits were lost, by state.
#' @return A data.table in panel format (id, time, state).
#' @examples
#' d <- data.frame(id = 1:2, t.a = c(2, 1), t.b = c(1, 3),
#'                 disch.s = c(1, 0), death.s = c(0, 1))
#' sojourn_to_panel(d, id = "id",
#'                  segments = c(A = "t.a", B = "t.b"),
#'                  absorbing = c(Disch = "disch.s", Death = "death.s"))
#' @export
sojourn_to_panel <- function(data, id, segments, absorbing, round_fun = rnd) {
  stopifnot(is.data.frame(data))
  seg_lab <- names(segments); seg_col <- unname(segments)
  abs_lab <- names(absorbing); abs_col <- unname(absorbing)
  miss <- setdiff(c(id, seg_col, abs_col), names(data))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  DT <- data.table::as.data.table(data)
  out <- vector("list", nrow(DT))
  ## Visits with a positive recorded duration that round to 0 time units,
  ## per state: they vanish from the panel, and with them the transitions
  ## into and out of that state (A -> B -> C becomes A -> C).
  lost <- stats::setNames(integer(length(seg_lab)), seg_lab)
  ## One subject (row of `data`) at a time: turn their per-state total
  ## durations into a run of repeated state labels, one entry per discrete
  ## time unit spent there, in the visiting order `segments` was given in.
  for (r in seq_len(nrow(DT))) {
    raw  <- as.numeric(DT[r, seg_col, with = FALSE])
    durs <- round_fun(raw)
    lost <- lost + (!is.na(raw) & raw > 0 & (is.na(durs) | durs <= 0))
    ## A negative or missing duration means "never visited this state" --
    ## treated the same as a genuine zero-duration visit (contributes no
    ## rows), rather than propagating an NA or a nonsensical negative rep().
    durs[is.na(durs) | durs < 0] <- 0
    seq_states <- rep(seg_lab, durs)   # e.g. rep(c("A","B"), c(3,1)) = A,A,A,B

    ## Whichever absorbing-state indicator column is 1 marks how (or
    ## whether) this subject's follow-up ended; if none is (a subject still
    ## under observation, right-censored with no terminal event recorded),
    ## the trajectory simply has no final absorbing-state row.
    ind <- as.numeric(DT[r, abs_col, with = FALSE])
    hit <- which(ind == 1)
    if (length(hit)) seq_states <- c(seq_states, abs_lab[hit[1]])

    if (length(seq_states) == 0L) next   # no duration anywhere, no event: nothing to record
    out[[r]] <- data.table::data.table(
      id = DT[[id]][r], time = seq_along(seq_states) - 1L, state = seq_states)
  }
  if (sum(lost))
    warning(sprintf(paste("%d visit(s) with a positive duration rounded to 0 time units",
                          "and were dropped (%s); the transitions into and out of them",
                          "are lost. Use a finer time unit or another `round_fun`."),
                    sum(lost), paste0(names(lost)[lost > 0], ": ", lost[lost > 0],
                                      collapse = ", ")), call. = FALSE)
  data.table::rbindlist(out)
}

#' Coerce a panel / msm2data object to an mstate transition matrix
#'
#' Builds the transition matrix of every move \eqn{j \to \ell} (\eqn{\ell \neq j})
#' observed in the data, including each subject's first move, in exactly the
#' format of \code{mstate::transMat()} (integers numbering the transitions row
#' by row, \code{NA} elsewhere, dimnames \code{from}/\code{to}). A first-order
#' \code{mstate} analysis built on it (\code{msprep()}, \code{msfit()},
#' \code{probtrans()}) uses the same state space, in the same order, as the
#' second-order fit, which \code{\link{compare_order}} requires. The
#' \pkg{mstate} package is not needed.
#'
#' @param object A "msm2data" object or a panel data frame.
#' @param ... Passed to \code{\link{prep2}} when \code{object} is a data frame.
#' @return A transition matrix as returned by \code{mstate::transMat}.
#' @examples
#' st   <- c("A", "B", "C")                                 # C is absorbing
#' tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
#' tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
#' tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
#' tens["C", "C", ]    <- 1
#' first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
#'
#' set.seed(1)
#' panel <- simulate2(200, tens, first, init = c(A = 1, B = 0, C = 0))
#' as_tmat(prep2(panel))      # A -> B (first move only) and B -> C
#' @export
as_tmat <- function(object, ...) {
  if (!inherits(object, "msm2data")) object <- prep2(object, ...)
  states <- object$states
  ## Which j -> l moves (j != l) were ever observed. `pairs` includes the
  ## first move of each subject; objects built before it existed fall back
  ## on the triples' N table.
  tr <- if (!is.null(object$pairs))
          unique(object$pairs[as.character(from) != as.character(to), .(j = from, l = to)])
        else unique(object$N[as.integer(j) != as.integer(l), .(j, l)])
  ## mstate::transMat() wants, for each state, the positions of the states
  ## reachable from it; .trans_mat() builds the same matrix without mstate.
  to_list <- lapply(states, function(s)
    sort(match(as.character(tr[as.character(j) == s, l]), states)))
  names(to_list) <- states
  .trans_mat(to_list, states)
}

## --- internal: identical to mstate::transMat(x, names) -----------------------
.trans_mat <- function(x, names) {
  ns <- length(x)
  tmat <- matrix(NA, nrow = ns, ncol = ns)
  idx <- cbind(rep(seq_len(ns), lengths(x)), unlist(x))
  if (nrow(idx)) tmat[idx] <- seq_len(nrow(idx))
  dimnames(tmat) <- list(from = names, to = names)
  tmat
}
