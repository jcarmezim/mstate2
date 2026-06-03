## A known three-state second-order model: A -> B -> {B, C}, C absorbing.
##   from (A,B): stay B 0.6, exit C 0.4
##   from (B,B): stay B 0.3, exit C 0.7
toy_tensor <- function() {
  st <- c("A", "B", "C")
  P <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  P["B","B","A"] <- 0.6; P["B","C","A"] <- 0.4
  P["B","B","B"] <- 0.3; P["B","C","B"] <- 0.7
  P["C","C", ]   <- 1
  P
}

toy_first <- function() {
  st <- c("A", "B", "C")
  m <- matrix(0, 3, 3, dimnames = list(st, st)); m["A","B"] <- 1; m
}

## A simulated panel from that model (seeded, so tests are deterministic).
toy_panel <- function(n = 4000, seed = 1, ...) {
  set.seed(seed)
  simulate2(n, toy_tensor(), toy_first(), init = c(A = 1, B = 0, C = 0), ...)
}

## A fitted P2est on that panel.
toy_fit <- function(estimator = "RPE", ...) {
  d <- prep2(toy_panel(...), states = c("A", "B", "C"))
  P2est(d, estimator)
}
