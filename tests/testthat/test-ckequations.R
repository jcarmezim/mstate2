test_that("n = 1 equals the 1-step tensor entry", {
  P <- toy_tensor()
  expect_equal(ckequations(P, "A", "B", "B", nsteps = 3)[1], P["B","B","A"])
})

test_that("2-step probability matches the explicit path sum", {
  P <- toy_tensor(); st <- dimnames(P)[[1]]; M <- length(st)
  manual <- sum(vapply(seq_len(M), function(m) P["B", m, "A"] * P[m, "B", "B"], 0.0))
  expect_equal(ckequations(P, "A", "B", "B", nsteps = 2)[2], manual, tolerance = 1e-12)
})

test_that("each step's distribution is a valid (sub)distribution", {
  full <- ckequations(toy_tensor(), "A", "B", nsteps = 6)   # l = NULL -> matrix
  expect_equal(dim(full), c(6, 3))
  expect_true(all(abs(rowSums(full) - 1) < 1e-9))           # mass conserved here
})

test_that("vector l returns the matching columns of the full distribution", {
  full <- ckequations(toy_tensor(), "A", "B", nsteps = 5)
  two  <- ckequations(toy_tensor(), "A", "B", l = c("B","C"), nsteps = 5)
  expect_equal(unname(two), unname(full[, c("B","C")]))
})

test_that("label and integer state args agree", {
  expect_equal(ckequations(toy_tensor(), "A", "B", "C", 5),
               ckequations(toy_tensor(), 1,   2,   3,   5))
})

test_that("mass accumulates monotonically in an absorbing state", {
  v <- ckequations(toy_tensor(), "A", "B", "C", nsteps = 40)
  expect_true(all(diff(v) >= -1e-12))
  expect_equal(v[40], 1, tolerance = 1e-3)
})

test_that("bounds = TRUE returns ordered evolution-interval bounds", {
  out <- ckequations(toy_fit(), "A", "B", "B", nsteps = 6, bounds = TRUE)
  expect_s3_class(out, "data.frame")
  expect_true(all(out$lower <= out$estimate + 1e-9))
  expect_true(all(out$estimate <= out$upper + 1e-9))
})

test_that("ckequations errors on a bad `x` type", {
  expect_error(ckequations(list(1, 2, 3), "A", "B", "B"), "P2est.*tensor")
})

test_that("ckequations errors on nsteps < 1", {
  expect_error(ckequations(toy_tensor(), "A", "B", "B", nsteps = 0), "nsteps")
})

test_that("ckequations errors on an unknown target or starting state", {
  expect_error(ckequations(toy_tensor(), "A", "B", "Z"), "not found")
  expect_error(ckequations(toy_tensor(), "A", "Z", "B"), "not found")
})

test_that("ckequations warns when bounds = TRUE is requested without an explicit l", {
  expect_warning(ckequations(toy_fit(), "A", "B", l = NULL, nsteps = 3, bounds = TRUE),
                 "bounds")
})

test_that("a raw tensor without dimnames falls back to integer-string state labels", {
  P <- toy_tensor(); dimnames(P) <- NULL
  out <- ckequations(P, 1, 2, 3, nsteps = 3)
  expect_equal(out, ckequations(toy_tensor(), "A", "B", "C", nsteps = 3))
})

test_that("evolution-interval bounds stay within [0, 1]", {
  st   <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  set.seed(1)
  fit <- P2est(prep2(simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))))
  ## the absorbing target C accumulates the upper limits past 1 if unclipped
  b <- ckequations(fit, h = "A", j = "B", l = "C", nsteps = 6, bounds = TRUE)
  expect_true(all(b$lower >= 0 & b$upper <= 1))
  cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "C", nsteps = 6)
  expect_true(all(cmp$lower >= 0 & cmp$upper <= 1))
})

test_that("ckequations reproduces the trace formulas of the paper (Eq. 6) for n = 1-4", {
  ## A random 4-state second-order tensor P[j, l, h] (rows sum to 1).
  set.seed(3)
  M <- 4; P <- array(stats::runif(M^3), c(M, M, M))
  for (h in 1:M) for (j in 1:M) P[j, , h] <- P[j, , h] / sum(P[j, , h])
  tr <- function(m) sum(diag(m))
  Pc <- function(l) P[, l, ]                 # P(l): column l of every P(h)
  h <- 2; j <- 3; l <- 1
  ## The terms of the Chapman-Kolmogorov function of the paper's code.
  out1 <- P[j, l, h]
  out2 <- sum(P[j, , h] * P[, l, j])
  out3 <- tr(P[j, , h] * (P[, , j] %*% Pc(l)))
  out4 <- sum(vapply(1:M, function(i) P[j, i, h] * tr(P[i, , j] * (P[, , i] %*% Pc(l))), 0))
  expect_equal(ckequations(P, h, j, l, nsteps = 4), c(out1, out2, out3, out4),
               tolerance = 1e-12)
})
