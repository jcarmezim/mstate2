## A first-order chain written as a second-order tensor (identical slices).
first_order_tensor <- function() {
  st <- c("A", "B", "C")
  P1 <- matrix(c(.70, .25, .05,
                 .20, .60, .20,
                 0,   0,   1), 3, byrow = TRUE, dimnames = list(st, st))
  list(T = array(P1, c(3, 3, 3), dimnames = list(st, st, st)), P1 = P1)
}

test_that("P1est pools the second-order counts over the previous state", {
  d <- prep2(toy_panel(1500, seed = 4), states = c("A", "B", "C"))
  p1 <- P1est(d, B = 0)
  expect_s3_class(p1, "P1est")
  expect_null(p1$boot)
  ## by hand: transitions out of B pooled over h, over the subject-instants at risk in B
  tr <- d$triples
  expect_equal(p1$P["B", "C"], mean(tr$l[tr$j == "B"] == "C"))
  expect_equal(unname(rowSums(p1$P)), c(0, 1, 1))     # A only at time 0: no row
  expect_equal(p1$P["C", ], c(A = 0, B = 0, C = 1))   # absorbing
  e <- p1$estimate
  expect_equal(e$se, sqrt(e$p * (1 - e$p) / e$at.risk))
  expect_true(all(e$lower >= 0 & e$upper <= 1))
  expect_output(print(p1), "P1est")
})

test_that("with first-order data the second-order estimates reproduce the first-order ones", {
  ft <- first_order_tensor()
  set.seed(21)
  pan <- simulate2(3000, ft$T, ft$P1, init = c(A = 1, B = 0, C = 0))
  d <- prep2(pan, states = c("A", "B", "C"))
  fit <- P2est(d)
  p1 <- P1est(d, B = 0)
  for (jj in c("A", "B")) for (hh in c("A", "B")) {
    expect_equal(fit$P[jj, , hh], p1$P[jj, ], tolerance = 0.05)
    expect_equal(p1$P[jj, ], ft$P1[jj, ], tolerance = 0.05)
  }
})

test_that("P1est bootstraps subjects, and reuses the replicates of P2boot", {
  d <- prep2(toy_panel(600, seed = 12), states = c("A", "B", "C"))
  b1 <- P1est(d, B = 40, seed = 1)
  expect_equal(dim(b1$boot), c(3, 3, 40))
  expect_identical(b1$boot, P1est(d, B = 40, seed = 1)$boot)
  rs <- apply(b1$boot, c(1, 3), sum)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
  bt <- P2boot(d, B = 30, seed = 3)
  expect_equal(dim(bt$boot1), c(3, 3, 30))
  p1 <- P1est(bt)
  expect_identical(p1$boot, bt$boot1)
  expect_equal(p1$P, P1est(d, B = 0)$P)
  expect_error(P1est(P2est(d)), "msm2data")
  expect_error(P1est(d, B = -1), "B")
})
