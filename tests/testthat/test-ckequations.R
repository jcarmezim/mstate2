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
