#' Compare second-order and first-order n-step predictions
#'
#' Quantifies how much the second-order prediction \eqn{P_{hj\ell}(1, n)}, which depends on the state at the previous time \eqn{h}, diverges from a first-order prediction \eqn{P(X_{n+1} = \ell \mid X_1 = j)}, which ignores it. A second-order model needs \eqn{M} transition matrices instead of one; this comparison shows, transition by transition, how much and for how long knowing \eqn{h} changes the forecast.
#'
#' The first-order baseline \code{pt} can be:
#' \itemize{
#'   \item a \code{\link{P1est}} fit with bootstrap replicates: the discrete-time first-order chain estimated from the same data and time grid (recommended; labelled \code{"1st order (discrete)"}). If it was built from the same \code{\link{P2boot}} fit as \code{x} (\code{P1est(x)}), both curves come from the same resamples and \code{\link{divergence}} also tests their paired difference;
#'   \item an \code{mstate::probtrans} object of a first-order (continuous-time) fit, computed with \code{predt = 1} so that its curves start at the same time 1 as \code{n = 1}, with its components indexed by starting state in the order of \code{x$states} (see \code{\link{as_tmat}}); labelled \code{"1st order (mstate)"}. Only its structure is used, so \pkg{mstate} does not need to be installed.
#' }
#' The first-order curve gets a Wald interval from its standard errors; the second-order curves get the intervals of \code{\link{ckequations}} (evolution intervals for a \code{P2est} fit, percentile bootstrap intervals for a \code{P2boot} fit).
#'
#' @param x A "P2est" or "P2boot" object (second-order fit).
#' @param pt The first-order baseline: a \code{\link{P1est}} fit with bootstrap replicates, or an \code{mstate::probtrans} object (see Details).
#' @param h Previous states to compare against the first-order baseline (labels or indices).
#' @param j Current state shared by the compared paths.
#' @param l Target state.
#' @param nsteps Number of steps (default 9).
#' @param conf.level Confidence level of the first-order interval. Default \code{x$conf.level}.
#' @return An object of class "msm2pred" (as returned by \code{\link{compare2}}): one group per previous state \code{h} plus the first-order baseline group, so that \code{plot()} works directly. The attribute \code{baseline} names the baseline group, and \code{boot_curves} holds the replicate curves when the comparison is paired. See \code{\link{divergence}}.
#' @seealso \code{\link{divergence}}, \code{\link{divergence_table}}, \code{\link{P1est}}, \code{\link{probtrans2}}
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
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' d <- prep2(panel)
#' fit <- P2est(d)
#' p1 <- P1est(d, B = 100, seed = 1)
#' co <- compare_order(fit, p1, h = c("A", "B"), j = "B", l = "C", nsteps = 6)
#' co
#' divergence(co)
#' plot(co, dualaxis = FALSE)
#' @export
compare_order <- function(x, pt, h, j, l, nsteps = 9L, conf.level = NULL) {

  # Check the arguments. A P1est baseline is turned into the equivalent probtrans object (first-order predictions from j with bootstrap SEs) on the same time grid.
  stopifnot(inherits(x, "P2est"))
  p1fit <- if (inherits(pt, "P1est")) pt else NULL
  base_lab <- "1st order (mstate)"
  if (!is.null(p1fit)) {
    if (is.null(pt$boot))
      stop("`pt` (P1est) has no bootstrap replicates, so the first-order predictions have no standard errors; use P1est(..., B > 0).", call. = FALSE)
    pt <- probtrans2(pt, j = j, nsteps = nsteps)
    base_lab <- "1st order (discrete)"
  }
  if (!inherits(pt, "probtrans"))
    stop("`pt` must be an mstate::probtrans object or a P1est fit (see ?compare_order).", call. = FALSE)
  if (is.null(conf.level)) conf.level <- x$conf.level
  .check_conf_level(conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # The current and target states, as positions in the state space.
  states <- x$states
  lab <- function(s) if (is.numeric(s)) states[s] else as.character(s)
  ji <- .resolve(j, states)
  li <- .resolve(l, states)
  if (anyNA(c(ji, li))) stop("`j` or `l` not found in the state space.", call. = FALSE)

  # The first-order predictions from j: component ji of pt (mstate indexes probtrans objects by starting state), with the columns pstate<l>/se<l> of the target.
  if (ji > length(pt) || is.null(pt[[ji]]) || !is.data.frame(pt[[ji]]))
    stop("`pt` has no predictions from state '", lab(j), "'. Check that it was built over the same state space and order as `x$states`.", call. = FALSE)
  pt_j <- pt[[ji]]
  pcol <- paste0("pstate", li)
  scol <- paste0("se", li)
  miss <- setdiff(c("time", pcol, scol), names(pt_j))
  if (length(miss))
    stop("`pt[[j]]` is missing the probtrans column(s) ", paste(miss, collapse = ", "),
         ". Was `pt` computed with predt = 1, i.e. mstate::probtrans(msf, predt = 1)?", call. = FALSE)
  if (!nrow(pt_j) || anyNA(pt_j$time))
    stop("`pt[[j]]` has no valid `time` values.", call. = FALSE)

  # Step n corresponds to time predt + n. probtrans gives right-continuous step functions, so the last row at or before each target time is used; pt must reach the last one.
  predt <- pt_j$time[1L]
  target_t <- predt + seq_len(nsteps)
  if (max(pt_j$time) < max(target_t))
    stop("`pt` does not reach nsteps = ", nsteps, " step(s) from predt = ", predt,
         " (it ends at time ", max(pt_j$time), ").", call. = FALSE)
  idx <- purrr::map_int(target_t, \(tt) max(which(pt_j$time <= tt)))

  # The first-order baseline curve with its Wald interval, in the layout of compare2().
  baseline <- tibble::tibble(h = base_lab, n = seq_len(nsteps),
                             estimate = pt_j[[pcol]][idx], se = pt_j[[scol]][idx]) |>
    dplyr::mutate(lower = pmax(0, estimate - z * se),
                  upper = pmin(1, estimate + z * se)) |>
    dplyr::select(-"se")

  # One second-order curve per previous state h, with the intervals of ckequations().
  second <- purrr::map(h, \(hh)
    ckequations(x, hh, j, l, nsteps = nsteps, bounds = TRUE) |>
      dplyr::mutate(h = lab(hh), .before = 1)) |>
    purrr::list_rbind()

  # Paired bootstrap: when the first-order fit comes from the same P2boot replicates (P1est(x)), keep the replicate curves of every group, so that divergence() can test the difference itself rather than only the overlap of two intervals.
  curves <- NULL
  if (inherits(x, "P2boot") && !is.null(p1fit) && identical(dim(p1fit$boot), dim(x$boot1)) &&
      isTRUE(all.equal(p1fit$boot, x$boot1, check.attributes = FALSE))) {
    curves <- c(stats::setNames(list(.boot_curves_first(x$boot1, ji, li, nsteps)), base_lab),
                stats::setNames(purrr::map(h, \(hh) .boot_curves(x$boot, .resolve(hh, states), ji, li, nsteps)),
                                purrr::map_chr(h, lab)))
  }

  # Same class as compare2() (so print, summary and plot work unchanged), plus the name of the baseline group.
  structure(as.data.frame(dplyr::bind_rows(baseline, second)), class = c("msm2pred", "data.frame"),
            j = lab(j), l = lab(l), bounds = TRUE, estimator = x$estimator,
            conf.level = conf.level, baseline = base_lab,
            bands = if (inherits(x, "P2boot")) "bootstrap" else "evolution",
            boot_curves = curves)
}

#' Memory horizon: for how long second-order predictions diverge from first-order ones
#'
#' For each previous state in a \code{\link{compare_order}} comparison, reports for how many leading steps its second-order interval and the first-order interval do not overlap (the "memory horizon" over which knowing the previous state changes the forecast), by calling \code{\link{overlap_step}} on each (baseline, h) pair.
#'
#' @param x An object returned by \code{\link{compare_order}}.
#' @return A data frame with one row per previous state \code{h}, sorted by decreasing \code{separated_steps}, with columns \code{h, n, s, separated_steps, diff_steps, max_abs_diff}: \code{n} and \code{s} are the first step (and time) at which the intervals overlap; \code{separated_steps} counts the leading steps whose intervals do not overlap; \code{diff_steps} counts the leading steps whose percentile interval of the paired difference excludes 0 (when \code{x} is a \code{\link{P2boot}} fit and the baseline is \code{P1est(x)}; \code{NA} otherwise); \code{max_abs_diff} is the largest absolute difference between the point estimates over the horizon. The paired difference is the direct test of whether the two predictions differ; the overlap criterion is more conservative (Schenker and Gentleman, 2001).
#' @references Schenker, N. and Gentleman, J. F. (2001). On judging the significance of differences by examining the overlap between confidence intervals. \emph{The American Statistician}, 55(3), 182-186.
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
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' bt <- P2boot(prep2(panel), B = 100, seed = 1)
#' divergence(compare_order(bt, P1est(bt), h = c("A", "B"), j = "B", l = "C", nsteps = 6))
#' @export
divergence <- function(x) {

  # Check that x has a first-order baseline and at least one second-order group.
  stopifnot(inherits(x, "msm2pred"))
  base_lab <- attr(x, "baseline")
  if (is.null(base_lab))
    stop("`x` has no first-order baseline group; build it with compare_order().", call. = FALSE)
  groups <- setdiff(unique(x$h), base_lab)
  if (!length(groups))
    stop("No second-order group found besides the baseline.", call. = FALSE)

  # overlap_step() compares exactly two groups: pair each second-order group with the baseline in turn (subsetting drops the attributes, so they are put back).
  one_group <- function(g) {
    sub <- structure(x[x$h %in% c(base_lab, g), ], class = c("msm2pred", "data.frame"),
                     j = attr(x, "j"), l = attr(x, "l"), bounds = TRUE,
                     estimator = attr(x, "estimator"), conf.level = attr(x, "conf.level"),
                     boot_curves = attr(x, "boot_curves"))
    os <- overlap_step(sub)
    tibble::tibble(h = g, n = os$n, s = os$s, separated_steps = os$separated_steps,
                   diff_steps = os$diff_steps,
                   max_abs_diff = max(abs(sub$estimate[sub$h == g] - sub$estimate[sub$h == base_lab])))
  }

  # Longest-diverging previous state first.
  purrr::map(groups, one_group) |>
    purrr::list_rbind() |>
    dplyr::arrange(dplyr::desc(separated_steps)) |>
    as.data.frame()
}

#' First- versus second-order comparison across all transitions
#'
#' Applies \code{\link{compare_order}} and \code{\link{divergence}} to every current state \eqn{j} reached from at least \code{min.h} different previous states and to every target state \eqn{\ell}, and stacks the results: for each \eqn{(j, \ell, h)}, how many leading steps the second-order prediction is separated from the first-order one (the memory horizon) and the largest absolute difference between them. It answers, for the whole model at once, for which transitions and up to which horizon the previous state matters.
#'
#' @param x A "P2est" or "P2boot" object (with a "P2boot", the second-order intervals are percentile bootstrap intervals).
#' @param pt First-order baseline: a \code{\link{P1est}} fit with bootstrap replicates (recommended: same data and time grid), or an \code{mstate::probtrans} object computed with \code{predt = 1}.
#' @param nsteps Horizon. Default 9.
#' @param min.h Minimum number of previous states for a current state to be included. Default 2.
#' @param targets Optional target states (default: every state).
#' @return A data frame with one row per \eqn{(j, \ell, h)} and columns \code{j, l, h, n, s, separated_steps, diff_steps, max_abs_diff} (see \code{\link{divergence}}), sorted by decreasing \code{separated_steps} and \code{max_abs_diff}. Targets with probability zero for every history and the baseline over the whole horizon are omitted.
#' @examples
#' st <- c("A", "B", "C")
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
#' panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
#' d <- prep2(panel)
#' divergence_table(P2est(d), P1est(d, B = 100, seed = 1), nsteps = 5)
#' @export
divergence_table <- function(x, pt, nsteps = 9L, min.h = 2, targets = NULL) {

  # Current states (not absorbing) and the previous states observed with each of them.
  stopifnot(inherits(x, "P2est"))
  if (is.null(targets)) targets <- x$states
  histories <- x$estimate |>
    dplyr::distinct(h = as.character(h), j = as.character(j)) |>
    dplyr::filter(!j %in% x$absorbing) |>
    dplyr::summarise(hs = list(sort(h)), .by = "j") |>
    dplyr::filter(lengths(hs) >= min.h)
  if (!nrow(histories))
    stop("No current state is reached from at least `min.h` previous states.", call. = FALSE)

  # compare_order() and divergence() for every (j, l); targets never reached are skipped.
  grid <- tidyr::expand_grid(histories, l = targets)
  rows <- purrr::pmap(grid, \(j, hs, l) {
    co <- compare_order(x, pt, h = hs, j = j, l = l, nsteps = nsteps)
    if (all(co$estimate == 0)) return(NULL)
    tibble::tibble(j = j, l = l, divergence(co))
  })

  # Longest and largest divergences first.
  purrr::list_rbind(rows) |>
    dplyr::arrange(dplyr::desc(separated_steps), dplyr::desc(max_abs_diff)) |>
    as.data.frame()
}
