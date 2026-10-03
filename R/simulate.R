#' Simulate a second-order Markov multistate process
#'
#' Generates discrete-time panel data from a second-order Markov model defined by a 1-step second-order transition tensor, a first-step (first-order) transition matrix, and entry probabilities.
#'
#' Each individual receives an entry state \eqn{X_0} (its personal time origin), then \eqn{X_1} from the first-step matrix, and subsequent states from the second-order tensor \eqn{X_s \sim P_{X_{s-2}, X_{s-1}, \cdot}}. Sampling stops when an absorbing state is reached.
#'
#' With \code{entry} supplied, individuals enter the process stochastically over \emph{global} time: at each step every not-yet-entered individual enters state \eqn{h} with probability \code{entry[h]}, otherwise waits, so that the number of individuals at risk varies across global time. With \code{entry = NULL} all individuals start at global time 0.
#'
#'
#' @param n Number of individuals.
#' @param tensor An \eqn{M\times M\times M} tensor in the \code{P[j,l,h]} \eqn{= P_{hj\ell}} layout (e.g. the \code{P} component of a \code{\link{P2est}} object). Absorbing states must satisfy \code{tensor[a,a,h] = 1}.
#' @param first An \eqn{M\times M} first-order matrix with \code{first[h,l] =} \eqn{P(X_1 = \ell \mid X_0 = h)}.
#' @param init Probability vector over states for the entry state \eqn{X_0} (used when \code{entry = NULL}). Defaults to uniform over non-absorbing states.
#' @param entry Optional named probability vector over entry states; if its sum is below 1 the remainder is the per-step probability of \emph{not} yet entering, giving staggered entry. If \code{NULL}, all individuals enter at time 0 according to \code{init}.
#' @param states Optional state labels (defaults to the tensor's dimnames or \code{1:M}).
#' @param maxT Safety cap on the number of global time steps. Default 1000.
#'
#' @return A tibble in panel format (\code{id}, \code{time}, \code{state}, with \code{state} a factor with levels \code{states}), ordered by id and time and suitable for \code{\link{prep2}}.
#' @examples
#' st <- c("A", "B", "C") # C is absorbing
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
#' head(panel)
#' @export
simulate2 <- function(n, tensor, first, init = NULL, entry = NULL, states = NULL, maxT = 1000) {

  # State labels and absorbing states
  M <- dim(tensor)[1]
  if (is.null(states)) {
    states <- dimnames(tensor)[[1]]
    if (is.null(states)) states <- as.character(seq_len(M))
  }
  absorbing <- which(purrr::map_lgl(seq_len(M), \(a) isTRUE(all(tensor[a, a, ] == 1))))
  transient <- setdiff(seq_len(M), absorbing)
  if (is.null(init)) { 
    init <- numeric(M)
    init[transient] <- 1 / length(transient) 
  }

  # One categorical draw from a vector of probabilities; NA when they are all zero, which ends that individual's path instead of stopping the simulation.
  draw <- function(probs) {
    if (sum(probs) <= 0) return(NA_integer_)
    sample.int(M, 1L, prob = probs)
  }

  # State of every individual, kept in vectors of length n (faster to update in the loop below than a table)
  entered <- rep(is.null(entry), n)
  absorbed <- rep(FALSE, n)
  stuck <- rep(FALSE, n)
  prev <- cur <- rep(NA_integer_, n)
  nobs <- rep(0L, n)
  rec_id <- rec_t <- rec_s <- integer(0)

  # Entry. Without `entry`, everybody starts at time 0 in a state drawn from `init`. With `entry`, check that it is named by states and keep the probability of entering at each step, p_enter.
  if (is.null(entry)) {
    cur <- sample.int(M, n, replace = TRUE, prob = init)
    nobs <- rep(1L, n)
    rec_id <- seq_len(n)
    rec_t <- rep(0L, n)
    rec_s <- cur
    absorbed <- cur %in% absorbing
  } else {
    # An unnamed `entry` would match no state and nobody would ever enter.
    if (is.null(names(entry)) || anyNA(match(names(entry), states)))
      stop("`entry` must be a named vector, named by entries of `states` ",
           "(or the tensor's dimnames), e.g. c(A = 0.1, B = 0.05).",
           call. = FALSE)
    entry_full <- numeric(M)
    entry_full[match(names(entry), states)] <- entry
    p_enter <- sum(entry_full)
  }

  # Advance global time s one step at a time, drawing the next state of all active individuals together
  s <- 0L
  repeat {
    s <- s + 1L

    # Staggered entry: each individual still waiting enters with probability p_enter, in a state drawn from `entry`; X_0 is recorded at time s
    ein <- integer(0)
    if (!is.null(entry)) {
      cand <- which(!entered)
      if (length(cand)) {
        u <- stats::runif(length(cand))
        ein <- cand[u < p_enter]
        if (length(ein)) {
          st <- sample.int(M, length(ein), replace = TRUE, prob = entry_full)
          cur[ein] <- st
          nobs[ein] <- 1L
          entered[ein] <- TRUE
          rec_id <- c(rec_id, ein)
          rec_t <- c(rec_t, rep(s, length(ein)))
          rec_s <- c(rec_s, st)
          absorbed[ein[st %in% absorbing]] <- TRUE
        }
      }
    }

    # Active individuals: entered, not absorbed, not stuck, and not entered at this same step (their first move is at s + 1)
    act <- setdiff(which(entered & !absorbed & !stuck), ein)
    if (length(act)) {
      newst <- integer(length(act))

      # First move (X_0 -> X_1): there is no previous state yet, so the first-order matrix `first` is used, row cur.
      i1 <- act[nobs[act] == 1L]
      for (g in unique(cur[i1])) {
        idx <- i1[cur[i1] == g]
        newst[match(idx, act)] <- purrr::map_int(idx, \(z) draw(first[cur[z], ]))
      }

      # Later moves: second order. Individuals with the same (previous, current) pair (h, j) share the distribution tensor[j, , h] and are drawn together.
      i2 <- act[nobs[act] >= 2L]
      if (length(i2)) {
        key <- prev[i2] * M + cur[i2]
        for (g in unique(key)) {
          idx <- i2[key == g]
          h <- prev[idx[1]]
          j <- cur[idx[1]]
          probs <- tensor[j, , h]
          newst[match(idx, act)] <-
            if (sum(probs) > 0) sample.int(M, length(idx), replace = TRUE, prob = probs)
            else NA_integer_
        }
      }

      # Record the new states and update each individual; an NA (no outgoing probability) marks the individual as stuck.
      keep <- !is.na(newst)
      ak <- act[keep]
      ns <- newst[keep]
      rec_id <- c(rec_id, ak)
      rec_t <- c(rec_t, rep(s, length(ak)))
      rec_s <- c(rec_s, ns)
      prev[ak] <- cur[ak]; cur[ak] <- ns
      nobs[ak] <- nobs[ak] + 1L
      absorbed[ak[ns %in% absorbing]] <- TRUE
      stuck[act[!keep]] <- TRUE
    }

    # Stop when everybody has entered and is absorbed or stuck, or at maxT with a warning.
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

  # Return the recorded panel, ordered by individual and time.
  tibble::tibble(
    id = rec_id, 
    time = rec_t, 
    state = factor(states[rec_s], levels = states)
  ) |>
    dplyr::arrange(id, time)
}
