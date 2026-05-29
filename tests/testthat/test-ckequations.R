make_tensor <- function() {
  set.seed(9); M <- 4
  st <- c("A", "B", "C", "D")
  P <- array(0, c(M, M, M), dimnames = list(st, st, st))
  for (h in 1:M) for (j in 1:M) { r <- runif(M); P[j, , h] <- r / sum(r) }
  P["D", , ] <- 0; P["D", "D", ] <- 1            # D absorbing
  P
}

test_that("n = 1 equals the 1-step tensor entry", {
  P <- make_tensor()
  v <- ckequations(P, "A", "B", "C", nsteps = 3)
  expect_equal(v[1], P["B", "C", "A"])
})

test_that("2-step probability matches the explicit path sum", {
  P <- make_tensor(); st <- dimnames(P)[[1]]; M <- length(st)
  h <- 1; j <- 2; l <- 3
  ## P(X3 = l | X0=h, X1=j) = sum_m P_{hjm} P_{jml}
  manual <- sum(sapply(seq_len(M), function(m) P[j, m, h] * P[m, l, j]))
  expect_equal(ckequations(P, h, j, l, nsteps = 2)[2], manual, tolerance = 1e-12)
})

test_that("mass accumulates monotonically in an absorbing state", {
  P <- make_tensor()
  v <- ckequations(P, "A", "B", "D", nsteps = 40)
  expect_true(all(diff(v) >= -1e-12))
  expect_equal(v[40], 1, tolerance = 1e-3)
})

test_that("works for label or integer state arguments", {
  P <- make_tensor()
  expect_equal(ckequations(P, "A", "B", "C", 5),
               ckequations(P, 1, 2, 3, 5))
})

test_that("ckequations on a P2est object returns evolution bounds", {
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  set.seed(4)
  d <- prep2(simulate2(2000, tens, first, init = c(A = 1, B = 0, C = 0)), states = st)
  fit <- P2est(d, "RPE")
  out <- ckequations(fit, "A", "B", "B", nsteps = 6, bounds = TRUE)
  expect_s3_class(out, "data.frame")
  expect_true(all(out$lower <= out$estimate + 1e-9))
  expect_true(all(out$estimate <= out$upper + 1e-9))
})
