test_that("P2reg(formula = ~1) recovers the RPE point estimate exactly", {
  panel <- toy_panel(n = 3000, seed = 31)
  panel$age <- 50                                   # unused by ~1, but required by prep2()
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  rpe <- P2est(d)$estimate

  expect_identical(names(P2reg(d, h = "A", j = "B")$models), "C")   # moves only
  reg <- P2reg(d, h = "A", j = "B", l = c("B", "C"), formula = ~1)      # staying on request
  expect_s3_class(reg, "P2reg")
  expect_setequal(names(reg$models), c("B", "C"))

  p_hat    <- predict(reg, newdata = data.frame(age = 50))
  expect_true(is.matrix(p_hat))
  p_true_B <- rpe$p[rpe$h == "A" & rpe$j == "B" & rpe$l == "B"]
  p_true_C <- rpe$p[rpe$h == "A" & rpe$j == "B" & rpe$l == "C"]
  expect_equal(unname(p_hat[1, "B"]), p_true_B, tolerance = 1e-6)
  expect_equal(unname(p_hat[1, "C"]), p_true_C, tolerance = 1e-6)
})

test_that("P2reg errors when the formula needs covariates the object lacks", {
  d <- prep2(toy_panel(), states = c("A", "B", "C"))
  expect_error(P2reg(d, h = "A", j = "B", formula = ~ age), "covariates")
  ## without baseline covariates, ~1 and the time scales s / d are still allowed
  expect_s3_class(P2reg(d, h = "A", j = "B"), "P2reg")
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

test_that("P2reg errors on a two-sided formula", {
  panel <- toy_panel(); panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "A", j = "B", formula = y ~ age), "one-sided")
})

test_that("P2reg handles right-hand sides long enough to line-wrap under deparse()", {
  panel <- toy_panel(n = 2000, seed = 37)
  ids <- unique(panel$id)
  set.seed(38)
  cov_names <- paste0("covariate_with_a_fairly_long_name_", 1:6)
  for (cv in cov_names) panel[[cv]] <- stats::setNames(sample(1:5, length(ids), TRUE), ids)[as.character(panel$id)]
  d <- prep2(panel, states = c("A", "B", "C"), covariates = cov_names)
  fm <- stats::as.formula(paste("~", paste(cov_names, collapse = " + ")))
  reg <- P2reg(d, h = "A", j = "B", formula = fm)
  ## the fitted model's formula must have exactly one "~" and one ".y" on the
  ## left; deparse() (vs. deparse1()) line-wrapping would duplicate both.
  fit_formula <- deparse1(formula(reg$models[[1]]))
  expect_equal(lengths(regmatches(fit_formula, gregexpr("~", fit_formula))), 1)
  expect_true(all(cov_names %in% all.vars(formula(reg$models[[1]]))))
})

test_that("summary.P2reg rejects an out-of-range conf.level", {
  panel <- toy_panel(n = 2000, seed = 39); panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1)
  expect_error(summary(reg, conf.level = 0), "conf.level")
})

test_that("P2reg errors on an unknown h, j or l", {
  panel <- toy_panel(); panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "Z", j = "B"), "not found")
  expect_error(P2reg(d, h = "A", j = "Z"), "not found")
  expect_error(P2reg(d, h = "A", j = "B", l = "Z"), "not found")
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
  expect_false(is.null(reg$models[["C"]]$P2reg.vcov))
  expect_output(print(summary(reg)), "cluster-robust")
})

test_that("subsetting summary.P2reg keeps its attributes and prints", {
  d <- prep2(toy_panel(), states = c("A", "B", "C"))
  sm <- summary(P2reg(d, "A", "B"))
  sub <- sm[sm$l == "C", c("l", "term", "HR")]
  expect_s3_class(sub, "summary.P2reg")
  expect_identical(attr(sub, "link"), "cloglog")
  expect_output(print(sub), "hazard ratio")
})
