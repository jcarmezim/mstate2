test_that("P2boot gives reproducible percentile bands inside [0, 1]", {
  d  <- prep2(toy_panel(800, seed = 5), states = c("A", "B", "C"))
  b1 <- P2boot(d, B = 60, seed = 9)
  b2 <- P2boot(d, B = 60, seed = 9)
  expect_s3_class(b1, "P2boot"); expect_s3_class(b1, "P2est")
  expect_identical(b1$boot, b2$boot)
  expect_equal(dim(b1$boot), c(3, 3, 3, 60))
  expect_equal(b1$P, P2est(d)$P)
  ## each replicate is a proper set of transition probabilities
  rs <- apply(b1$boot, c(1, 3, 4), sum)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
  ck <- ckequations(b1, "A", "B", "C", nsteps = 5, bounds = TRUE)
  expect_true(all(ck$lower >= 0 & ck$upper <= 1 & ck$lower <= ck$estimate + 1e-12))
  expect_true("se.boot" %in% names(b1$estimate))
  expect_output(print(b1), "bootstrap")
  expect_error(P2boot(d, B = 1), "B")
})

test_that("compare2 on a P2boot fit uses bootstrap bands, narrower than evolution intervals", {
  d   <- prep2(toy_panel(800, seed = 5), states = c("A", "B", "C"))
  bt  <- P2boot(d, B = 100, seed = 1)
  cb  <- compare2(bt, h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  ce  <- compare2(P2est(d), h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  expect_identical(attr(cb, "bands"), "bootstrap")
  expect_identical(attr(ce, "bands"), "evolution")
  expect_equal(cb$estimate, ce$estimate)
  expect_true(mean(cb$upper - cb$lower) < mean(ce$upper - ce$lower))
  os <- overlap_step(cb)
  expect_setequal(names(os), c("n", "s", "separated_steps", "overlap", "separation", "groups"))
  expect_output(summary(cb), "bootstrap intervals")
})

test_that("the bootstrap SD agrees with the analytic standard error", {
  d  <- prep2(toy_panel(1500, seed = 2), states = c("A", "B", "C"))
  bt <- P2boot(d, B = 400, seed = 4)
  e  <- bt$estimate[bt$estimate$se > 0, ]
  expect_true(all(abs(e$se.boot / e$se - 1) < 0.25))
})
