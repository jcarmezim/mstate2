test_that("P2reg(formula = ~1) recovers the RPE point estimate exactly", {
  panel <- toy_panel(n = 3000, seed = 31)
  panel$age <- 50                                   # unused by ~1, but required by prep2()
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  rpe <- P2est(d, "RPE")$estimate

  reg <- P2reg(d, h = "A", j = "B", formula = ~1)
  expect_s3_class(reg, "P2reg")
  expect_setequal(names(reg$models), c("B", "C"))

  p_hat    <- predict(reg, newdata = data.frame(age = 50))
  expect_true(is.matrix(p_hat))
  p_true_B <- rpe$p[rpe$h == "A" & rpe$j == "B" & rpe$l == "B"]
  p_true_C <- rpe$p[rpe$h == "A" & rpe$j == "B" & rpe$l == "C"]
  expect_equal(unname(p_hat[1, "B"]), p_true_B, tolerance = 1e-6)
  expect_equal(unname(p_hat[1, "C"]), p_true_C, tolerance = 1e-6)
})

test_that("P2reg errors when the object carries no covariates", {
  d <- prep2(toy_panel(), states = c("A", "B", "C"))
  expect_error(P2reg(d, h = "A", j = "B"), "covariates")
})

test_that("P2reg errors when formula references an unknown covariate", {
  panel <- toy_panel(); panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "A", j = "B", formula = ~sex), "object\\$covariates")
})

test_that("P2reg errors on an (h, j) pair with no data", {
  panel <- toy_panel(); panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "C", j = "C", formula = ~1), "No observed")
})

test_that("P2reg fits a single requested destination when l is given", {
  panel <- toy_panel(n = 2000, seed = 32); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", l = "C", formula = ~1)
  expect_identical(names(reg$models), "C")
})

test_that("summary.P2reg labels exp(coef) as HR under the default cloglog link", {
  panel <- toy_panel(n = 2000, seed = 33); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1)
  s   <- summary(reg)
  expect_s3_class(s, "summary.P2reg")
  expect_true(all(c("HR", "HR.lower", "HR.upper") %in% names(s)))
  expect_output(print(s), "hazard ratio")
})

test_that("summary.P2reg does not label exp(coef) as HR under a non-cloglog link", {
  panel <- toy_panel(n = 2000, seed = 34); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1, family = stats::binomial("logit"))
  s   <- summary(reg)
  expect_true(all(c("exp.coef", "exp.lower", "exp.upper") %in% names(s)))
})

test_that("print.P2reg reports the risk set and covariates", {
  panel <- toy_panel(n = 2000, seed = 35); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~age)
  expect_output(print(reg), "\\(h, j\\) = \\(A, B\\)")
  expect_output(print(reg), "age")
})

test_that("cluster = TRUE uses a cluster-robust vcov when sandwich is available", {
  skip_if_not_installed("sandwich")
  panel <- toy_panel(n = 2000, seed = 36); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1, cluster = TRUE)
  expect_true(reg$cluster)
  expect_false(is.null(reg$models[["B"]]$P2reg.vcov))
  expect_output(print(summary(reg)), "cluster-robust")
})
