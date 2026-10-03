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
#' These are the counting processes of Section 2.2 of the methods paper
#' (Najera-Zuloaga, Besalu and Gomez Melis, 2025), indexed by the time
#' \eqn{s \ge 2} of the destination state. The at-risk count is computed as
#' \eqn{\tilde Y_{hj}(s-1) = \sum_\ell \tilde N_{hj\ell}(s)}: with complete
#' follow-up, as in the paper and in DIVINE, this is exactly the number of
#' subjects in \eqn{h} at \eqn{s-2} and \eqn{j} at \eqn{s-1}; a subject
#' whose follow-up stops in \eqn{j} at \eqn{s-1} (right-censored) is not
#' counted at risk at \eqn{s}, because the next state is not observed.
#'
#' Consecutive \emph{observations} are treated as consecutive time steps. Supply
#' a regular integer time grid (e.g. days) for the counting-process indices to
#' match the paper's definitions.
#'
#' @param data A data frame in long/panel format (one row per subject and
#'   time point), e.g. the output of \code{\link{sojourn_to_panel}} or
#'   \code{\link{simulate2}}.
#' @param id,time,state Column names for subject id, discrete time, and state.
#' @param states Optional state space / ordering. Defaults to sorted observed states.
#' @param absorbing Optional absorbing states; inferred if NULL.
#' @param drop.na If TRUE, drop rows with NA in id/time/state (with a warning)
#'   instead of erroring. Default FALSE.
#' @param check.consecutive If TRUE (default), warn when any subject's times are
#'   not consecutive integers.
#' @return An object of class "msm2data": a list with the tibbles \code{N}
#'   (the counts \eqn{\tilde N_{hj\ell}(s)}, columns \code{h, j, l, s, N}),
#'   \code{Y} (the at-risk counts \eqn{\tilde Y_{hj}(s-1)}, columns
#'   \code{h, j, s, Y}) and \code{triples} (one row per subject per observed
#'   triple, columns \code{id, h, j, l, s}; used by \code{\link{P2boot}} to
#'   resample subjects), where \code{h, j, l} are factors with levels
#'   \code{states}; and \code{states}, \code{absorbing}, \code{n}
#'   (subjects), \code{ntriples} and \code{time.range}.
#' @section Why these counts:
#' The likelihood of a second-order Markov chain depends on the data only
#' through the transition counts \eqn{\tilde N_{hj\ell}} and the at-risk
#' counts \eqn{\tilde Y_{hj}} (Anderson and Goodman, 1957), so they are
#' computed once and every other function reuses them.
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. \emph{Annals of Mathematical Statistics}, 28(1), 89-110.
#'
#' Najera-Zuloaga, J., Besalu, M. and Gomez Melis, G. (2025). Second-order
#' Markov multistate models: nonparametric estimation and inference.
#' Manuscript submitted for publication.
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
                  states = NULL, absorbing = NULL,
                  drop.na = FALSE, check.consecutive = TRUE) {

  ## 1. Check the input and keep only the three columns that are needed,
  ##    renamed to the canonical names id, time and state.
  stopifnot(is.data.frame(data))
  cols <- c(id = id, time = time, state = state)
  missing_cols <- cols[!cols %in% names(data)]
  if (length(missing_cols))
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  panel <- data |>
    tibble::as_tibble() |>
    dplyr::select(dplyr::all_of(cols))

  ## 2. Missing values. A row with a missing id, time or state cannot be
  ##    placed in the panel, so it is an error unless drop.na = TRUE, in
  ##    which case those rows are dropped with a warning.
  bad <- is.na(panel$id) | is.na(panel$time) | is.na(panel$state)
  if (any(bad)) {
    if (!drop.na)
      stop("`id`, `time` and `state` must not contain NA ",
           "(set drop.na = TRUE to drop these rows).", call. = FALSE)
    warning(sprintf("Dropping %d row(s) with NA in id/time/state.", sum(bad)),
            call. = FALSE)
    panel <- dplyr::filter(panel, !bad)
  }

  ## 3. State space. Every later step (the counts, the P2est() tensors, the
  ##    pair chain of ckequations()) indexes the states by their position in
  ##    `states`, so it is fixed here: the value given by the user or, by
  ##    default, the observed states, sorted. Storing the state as a factor
  ##    with these levels makes an observed state that is not in `states`
  ##    turn into NA, which is then reported as an error.
  if (is.null(states))
    states <- if (is.factor(panel$state)) levels(droplevels(panel$state))
              else sort(unique(panel$state))
  panel <- dplyr::mutate(panel, state = factor(state, levels = states))
  if (anyNA(panel$state)) stop("Some observed states are not in `states`.", call. = FALSE)

  ## 4. Order each subject's observations in time: the triples below are
  ##    built from consecutive rows.
  n     <- dplyr::n_distinct(panel$id)
  panel <- dplyr::arrange(panel, id, time)

  ## 5. Regular time grid. Consecutive rows are taken as consecutive time
  ##    steps, so a subject observed on days 0, 1 and 4 would have days 1 and
  ##    4 treated as adjacent. This is a warning rather than an error, because
  ##    occasional gaps are common in real panels and filling them would mean
  ##    assuming the unobserved states.
  if (check.consecutive) {
    gaps <- panel |>
      dplyr::summarise(ok = all(diff(time) == 1L), .by = "id")
    if (!all(gaps$ok))
      warning(sprintf(paste("%d subject(s) have non-consecutive time values;",
                            "triples are formed from consecutive observations."),
                      sum(!gaps$ok)), call. = FALSE)
  }

  ## 6. Second-order triples. For every observation X_s = l of a subject,
  ##    the states two rows and one row earlier are h = X_{s-2} and
  ##    j = X_{s-1}. The first two observations of each subject have no such
  ##    history (lag() gives NA) and are dropped.
  trip <- panel |>
    dplyr::mutate(h = dplyr::lag(state, 2L),
                  j = dplyr::lag(state, 1L), .by = "id") |>
    dplyr::filter(!is.na(h), !is.na(j)) |>
    dplyr::select("id", "h", "j", l = "state", s = "time")
  if (nrow(trip) == 0L)
    stop("No second-order triples found: each subject needs >= 3 observations.",
         call. = FALSE)

  ## 7. Counting processes of the paper (Section 2.2).
  ##    N_hjl(s): number of subjects with h -> j -> l at times (s-2, s-1, s).
  ##    Y_hj(s-1): number of subjects at risk, i.e. with h -> j at
  ##    (s-2, s-1); computed as the sum of N over l, so that N and Y are
  ##    aligned at every (h, j, s), which P2est() relies on.
  N <- trip |>
    dplyr::count(h, j, l, s, name = "N")
  Y <- N |>
    dplyr::summarise(Y = sum(N), .by = c("h", "j", "s")) |>
    dplyr::arrange(h, j, s)

  ## 8. Absorbing states, if not given: a state that is occupied (appears as
  ##    j or l) but that nobody is ever seen leaving (it is never the j of a
  ##    triple with l != j).
  if (is.null(absorbing)) {
    moved <- trip |>
      dplyr::filter(as.integer(j) != as.integer(l)) |>
      dplyr::distinct(j) |>
      dplyr::pull(j)
    occupied  <- union(levels(droplevels(trip$j)), levels(droplevels(trip$l)))
    absorbing <- setdiff(occupied, as.character(moved))
  }

  ## 9. Return the counts together with the information later functions need.
  structure(
    list(N = N, Y = Y, triples = trip,
         states = states, absorbing = as.character(absorbing),
         n = n, ntriples = nrow(trip), time.range = range(panel$time)),
    class = "msm2data"
  )
}

## Compact summary of an "msm2data" object, shown when it is printed.
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
  cat(sprintf("  distinct (h,j)  : %d\n", nrow(dplyr::distinct(x$Y, h, j))))
  invisible(x)
}

## Exposure of every observed (h, j) pair: the total number of subject-instants
## at risk (the denominator of the RPE) and the first and last time at which
## the pair is observed. Printed and returned invisibly as a tibble.
#' @export
summary.msm2data <- function(object, ...) {
  exposure <- object$Y |>
    dplyr::summarise(total_at_risk = sum(Y), s_min = min(s), s_max = max(s),
                     .by = c("h", "j")) |>
    dplyr::arrange(h, j)
  cat("<msm2data summary>\n")
  cat(sprintf("  %d subjects, %d triples, time %d-%d\n",
              object$n, object$ntriples,
              object$time.range[1], object$time.range[2]))
  cat("  exposure per (h, j) pair:\n")
  print(exposure, n = Inf)
  invisible(exposure)
}
