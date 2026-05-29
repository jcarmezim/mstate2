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
#' with probability \code{entry[h]}, otherwise waits. This induces the unbalanced
#' across-time exposure under which the RPE and CPE estimators differ. With
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
#' @export
simulate2 <- function(n, tensor, first, init = NULL, entry = NULL,
                      states = NULL, maxT = 1000) {
  M <- dim(tensor)[1]
  if (is.null(states)) {
    states <- dimnames(tensor)[[1]]
    if (is.null(states)) states <- as.character(seq_len(M))
  }
  absorbing <- which(vapply(seq_len(M),
                            function(a) isTRUE(all(tensor[a, a, ] == 1)), logical(1)))
  transient <- setdiff(seq_len(M), absorbing)
  if (is.null(init)) { init <- numeric(M); init[transient] <- 1 / length(transient) }

  draw <- function(probs) {                 # one categorical draw -> state index
    if (sum(probs) <= 0) return(NA_integer_)
    sample.int(M, 1L, prob = probs)
  }

  ## per-individual state
  entered  <- rep(is.null(entry), n)
  absorbed <- rep(FALSE, n)
  prev <- cur <- rep(NA_integer_, n)
  nobs <- rep(0L, n)
  rec_id <- rec_t <- rec_s <- integer(0)

  if (is.null(entry)) {                     # everyone in at t = 0
    cur <- sample.int(M, n, replace = TRUE, prob = init)
    nobs <- rep(1L, n)
    rec_id <- seq_len(n); rec_t <- rep(0L, n); rec_s <- cur
    absorbed <- cur %in% absorbing
  } else {
    entry_full <- numeric(M)
    entry_full[match(names(entry), states)] <- entry
    p_enter <- sum(entry_full)
  }

  s <- 0L
  repeat {
    s <- s + 1L
    ## --- staggered entries at this global step ---
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
    act <- which(entered & !absorbed)
    if (length(act)) {
      newst <- integer(length(act))
      ## first step (nobs == 1): use first-order matrix, grouped by current state
      i1 <- act[nobs[act] == 1L]
      for (g in unique(cur[i1])) {
        idx <- i1[cur[i1] == g]
        newst_idx <- vapply(idx, function(z) draw(first[cur[z], ]), integer(1))
        newst[match(idx, act)] <- newst_idx
      }
      ## second order (nobs >= 2): use tensor, grouped by (prev, cur) pair
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
      keep <- !is.na(newst)
      ak <- act[keep]; ns <- newst[keep]
      rec_id <- c(rec_id, ak); rec_t <- c(rec_t, rep(s, length(ak))); rec_s <- c(rec_s, ns)
      prev[ak] <- cur[ak]; cur[ak] <- ns; nobs[ak] <- nobs[ak] + 1L
      absorbed[ak[ns %in% absorbing]] <- TRUE
    }
    done <- all(absorbed) && (is.null(entry) || all(entered))
    if (done || s >= maxT) break
  }

  data.table::data.table(id = rec_id, time = rec_t,
                         state = factor(states[rec_s], levels = states)
  )[order(id, time)]
}
