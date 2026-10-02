#' Simulate a second-order Markov multistate process
#'
#' Generates discrete-time panel data from a second-order Markov model defined by
#' a 1-step second-order transition tensor, a first-step (first-order)
#' transition matrix, and entry probabilities. This generalises the simulation in
#' Section 5 of the paper and produces output ready for \code{\link{prep2}}.
#'
#' Each individual receives an entry state \eqn{X_0} (its personal time origin),
#' then \eqn{X_1} from the first-step matrix, and subsequent states from the
#' second-order tensor \eqn{X_s \sim P_{X_{s-2}, X_{s-1}, \cdot}}. Sampling stops
#' when an absorbing state is reached.
#'
#' With \code{entry} supplied, individuals enter the process stochastically over
#' \emph{global} time (the staggered-entry mechanism of the paper's auxiliary
#' state 0): at each step every not-yet-entered individual enters state \eqn{h}
#' with probability \code{entry[h]}, otherwise waits, so that the number of
#' individuals at risk varies across global time. With
#' \code{entry = NULL} all individuals start at global time 0.
#'
#' @param n Number of individuals.
#' @param tensor An \eqn{M\times M\times M} tensor in the \code{P[j,l,h]}
#'   \eqn{= P_{hj\ell}} layout (e.g. the \code{P} component of a
#'   \code{\link{P2est}} object). Absorbing states must satisfy
#'   \code{tensor[a,a,h] = 1}.
#' @param first An \eqn{M\times M} first-order matrix with \code{first[h,l] =}
#'   \eqn{P(X_1 = \ell \mid X_0 = h)}.
#' @param init Probability vector over states for the entry state \eqn{X_0}
#'   (used when \code{entry = NULL}). Defaults to uniform over non-absorbing
#'   states.
#' @param entry Optional named probability vector over entry states; if its sum
#'   is below 1 the remainder is the per-step probability of \emph{not} yet
#'   entering, giving staggered entry. If \code{NULL}, all individuals enter at
#'   time 0 according to \code{init}.
#' @param states Optional state labels (defaults to the tensor's dimnames or
#'   \code{1:M}).
#' @param maxT Safety cap on the number of global time steps. Default 1000.
#'
#' @return A \code{data.table} in panel format (\code{id}, \code{time},
#'   \code{state}), suitable for \code{\link{prep2}}.
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
#' head(panel)
#' @export
simulate2 <- function(n, tensor, first, init = NULL, entry = NULL,
                      states = NULL, maxT = 1000) {
  M <- dim(tensor)[1]
  if (is.null(states)) {
    states <- dimnames(tensor)[[1]]
    if (is.null(states)) states <- as.character(seq_len(M))
  }
  ## A state is treated as absorbing purely by looking at `tensor`: one whose
  ## self-transition probability is 1 no matter the preceding state h. This
  ## mirrors P2est()'s own convention (see its absorbing-forcing loop) and is
  ## how the simulation knows when an individual is done, without the caller
  ## having to separately declare which states are absorbing.
  absorbing <- which(vapply(seq_len(M),
                            function(a) isTRUE(all(tensor[a, a, ] == 1)), logical(1)))
  transient <- setdiff(seq_len(M), absorbing)
  if (is.null(init)) { init <- numeric(M); init[transient] <- 1 / length(transient) }

  ## One categorical draw from an (unnormalized) probability vector -> a
  ## state index, or NA if there's nowhere to go at all (every probability
  ## zero) -- e.g. an under-specified row of `tensor`/`first` for an (h, j)
  ## or starting state that does occur in the simulation but was never given
  ## outgoing probabilities. sample.int() itself would error on an all-zero
  ## `prob`; returning NA instead lets the caller handle it as "this
  ## individual's path ends here" (see the `stuck` mechanism below) rather
  ## than crashing the whole simulation.
  draw <- function(probs) {                 # one categorical draw -> state index
    if (sum(probs) <= 0) return(NA_integer_)
    sample.int(M, 1L, prob = probs)
  }

  ## per-individual state, tracked in parallel vectors of length n (a
  ## data.table/data.frame of individuals would be more idiomatic but much
  ## slower to update in a tight loop like the one below):
  ##   entered  -- has this individual entered the process yet (always TRUE
  ##               up front unless staggered `entry` is used)
  ##   absorbed -- has this individual reached an absorbing state
  ##   stuck    -- did this individual draw into a zero-probability row (see
  ##               `draw()` above): if so, there is nowhere left for them to
  ##               go, and they must stop being re-drawn every step even
  ##               though they were never formally absorbed
  ##   prev/cur -- the individual's previous and current state (X_{s-1}, X_s)
  ##   nobs     -- how many observations this individual has recorded so far,
  ##               which decides whether their *next* draw uses `first`
  ##               (nobs == 1, no history yet) or `tensor` (nobs >= 2)
  entered  <- rep(is.null(entry), n)
  absorbed <- rep(FALSE, n)
  stuck    <- rep(FALSE, n)          # drew into a zero-probability row: no further draws
  prev <- cur <- rep(NA_integer_, n)
  nobs <- rep(0L, n)
  rec_id <- rec_t <- rec_s <- integer(0)   # the output panel, built up incrementally

  if (is.null(entry)) {                     # everyone in at t = 0
    cur <- sample.int(M, n, replace = TRUE, prob = init)
    nobs <- rep(1L, n)
    rec_id <- seq_len(n); rec_t <- rep(0L, n); rec_s <- cur
    absorbed <- cur %in% absorbing
  } else {
    ## `entry` must be named: an unnamed vector has names(entry) == NULL, so
    ## match(NULL, states) silently resolves to no positions at all, leaving
    ## entry_full entirely zero (p_enter = 0) -- nobody would ever enter, and
    ## the only visible symptom is an empty output panel after burning
    ## through the whole maxT loop with a warning that (wrongly) points at
    ## `tensor` instead of the real cause.
    if (is.null(names(entry)) || anyNA(match(names(entry), states)))
      stop("`entry` must be a named vector, named by entries of `states` ",
           "(or the tensor's dimnames), e.g. c(A = 0.1, B = 0.05).",
           call. = FALSE)
    entry_full <- numeric(M)
    entry_full[match(names(entry), states)] <- entry
    p_enter <- sum(entry_full)
  }

  ## `s` is global calendar time; the loop advances it one step at a time,
  ## handling every individual's step in a vectorized (grouped) fashion
  ## rather than looping over individuals one by one.
  s <- 0L
  repeat {
    s <- s + 1L
    ## --- staggered entries at this global step ---
    ## Each not-yet-entered individual enters this step with probability
    ## p_enter (the auxiliary "state 0" mechanism from the paper); the
    ## remaining 1 - p_enter probability mass is implicitly "keep waiting",
    ## which is why entry doesn't need its own explicit "not yet" category.
    ein <- integer(0)                     # who enters at this step (if anyone)
    if (!is.null(entry)) {
      cand <- which(!entered)
      if (length(cand)) {
        u <- stats::runif(length(cand))
        ein <- cand[u < p_enter]
        if (length(ein)) {
          st <- sample.int(M, length(ein), replace = TRUE, prob = entry_full)
          cur[ein] <- st; nobs[ein] <- 1L; entered[ein] <- TRUE
          rec_id <- c(rec_id, ein); rec_t <- c(rec_t, rep(s, length(ein))); rec_s <- c(rec_s, st)
          absorbed[ein[st %in% absorbing]] <- TRUE
        }
      }
    }
    ## --- advance active individuals ---
    ## Only individuals who have entered, aren't yet absorbed, and aren't
    ## stuck in a dead end get a new draw this step. Anyone who entered at
    ## *this very* step has just recorded X_0 at time s, so they wait until
    ## s + 1 for their first move -- advancing them now would record X_1 at
    ## the same time s as X_0 (two rows with one time value per subject).
    act <- setdiff(which(entered & !absorbed & !stuck), ein)
    if (length(act)) {
      newst <- integer(length(act))
      ## first step (nobs == 1): use first-order matrix, grouped by current state
      ## An individual's very first move (X_0 -> X_1) has no preceding state
      ## to condition on yet, so it uses the plain first-order matrix
      ## `first`, not the second-order `tensor`. Individuals sharing the
      ## same current state g get drawn together in one vectorized
      ## sample.int() call per distinct g, rather than one call per person.
      i1 <- act[nobs[act] == 1L]
      for (g in unique(cur[i1])) {
        idx <- i1[cur[i1] == g]
        newst_idx <- vapply(idx, function(z) draw(first[cur[z], ]), integer(1))
        newst[match(idx, act)] <- newst_idx
      }
      ## second order (nobs >= 2): use tensor, grouped by (prev, cur) pair
      ## From the second move onward, the model is genuinely second-order:
      ## the draw depends on both the previous state h = prev and the
      ## current one j = cur. `key` packs each (prev, cur) pair into a
      ## single integer so individuals who share the exact same (h, j) risk
      ## set -- and therefore the exact same draw distribution
      ## tensor[j, , h] -- can again be sampled together in one call.
      i2 <- act[nobs[act] >= 2L]
      if (length(i2)) {
        key <- prev[i2] * M + cur[i2]
        for (g in unique(key)) {
          idx <- i2[key == g]
          h <- prev[idx[1]]; j <- cur[idx[1]]
          probs <- tensor[j, , h]
          newst[match(idx, act)] <-
            if (sum(probs) > 0) sample.int(M, length(idx), replace = TRUE, prob = probs)
            else NA_integer_
        }
      }
      ## Record only the individuals who actually got a valid next state
      ## this step; an NA (from either branch above) means "nowhere to go",
      ## so that individual's trajectory simply stops here.
      keep <- !is.na(newst)
      ak <- act[keep]; ns <- newst[keep]
      rec_id <- c(rec_id, ak); rec_t <- c(rec_t, rep(s, length(ak))); rec_s <- c(rec_s, ns)
      prev[ak] <- cur[ak]; cur[ak] <- ns; nobs[ak] <- nobs[ak] + 1L
      absorbed[ak[ns %in% absorbing]] <- TRUE
      stuck[act[!keep]] <- TRUE       # zero-probability row for this (h, j): stop retrying
    }
    ## Everyone is either absorbed or permanently stuck, and (in staggered-
    ## entry mode) everyone who will ever enter has already entered: nothing
    ## left to simulate, so stop before reaching maxT.
    done <- all(absorbed | stuck) && (is.null(entry) || all(entered))
    if (done) break
    if (s >= maxT) {
      warning(sprintf(
        paste("simulate2() stopped at maxT = %d with %d/%d individual(s) not yet absorbed",
              "(check for a zero-probability row in `tensor` for an observed (h, j) pair)."),
        maxT, sum(!absorbed & !stuck), n), call. = FALSE)
      break
    }
  }

  ## Assemble the recorded (id, time, state) triples into the panel format
  ## prep2() expects, sorted so every subject's rows are chronological.
  data.table::data.table(id = rec_id, time = rec_t,
                         state = factor(states[rec_s], levels = states)
  )[order(id, time)]
}
