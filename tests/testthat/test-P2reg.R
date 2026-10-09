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
  panel <- toy_panel()
  panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "A", j = "B", formula = ~sex), "not covariates")
})

test_that("P2reg errors on an (h, j) pair with no data", {
  panel <- toy_panel()
  panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "C", j = "C", formula = ~1), "No observed")
})

test_that("P2reg errors on a two-sided formula", {
  panel <- toy_panel()
  panel$age <- 50
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
  panel <- toy_panel(n = 2000, seed = 39)
  panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1)
  expect_error(summary(reg, conf.level = 0), "conf.level")
})

test_that("P2reg errors on an unknown h, j or l", {
  panel <- toy_panel()
  panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_error(P2reg(d, h = "Z", j = "B"), "not found")
  expect_error(P2reg(d, h = "A", j = "Z"), "not found")
  expect_error(P2reg(d, h = "A", j = "B", l = "Z"), "not found")
})

test_that("P2reg fits a single requested destination when l is given", {
  panel <- toy_panel(n = 2000, seed = 32)
  panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", l = "C", formula = ~1)
  expect_identical(names(reg$models), "C")
})

test_that("summary.P2reg labels exp(coef) as HR under the default cloglog link", {
  panel <- toy_panel(n = 2000, seed = 33)
  panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1)
  s   <- summary(reg)
  expect_s3_class(s, "summary.P2reg")
  expect_true(all(c("HR", "HR.lower", "HR.upper") %in% names(s)))
  expect_output(print(s), "hazard ratio")
})

test_that("summary.P2reg does not label exp(coef) as HR under a non-cloglog link", {
  panel <- toy_panel(n = 2000, seed = 34)
  panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~1, family = stats::binomial("logit"))
  s   <- summary(reg)
  expect_true(all(c("exp.coef", "exp.lower", "exp.upper") %in% names(s)))
})

test_that("print.P2reg reports the risk set and covariates", {
  panel <- toy_panel(n = 2000, seed = 35)
  panel$age <- 50
  d   <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  reg <- P2reg(d, h = "A", j = "B", formula = ~age)
  expect_output(print(reg), "\\(h, j\\) = \\(A, B\\)")
  expect_output(print(reg), "age")
})

test_that("cluster = TRUE uses a cluster-robust vcov when sandwich is available", {
  skip_if_not_installed("sandwich")
  panel <- toy_panel(n = 2000, seed = 36)
  panel$age <- 50
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

## A panel with a subject-level covariate that plays no role in the simulation.
panel_with_age <- function(n, seed) {
  pan <- toy_panel(n, seed = seed)
  ids <- unique(pan$id)
  set.seed(1)
  pan$age <- stats::setNames(sample(40:80, length(ids), TRUE), ids)[as.character(pan$id)]
  pan
}

test_that("P2reg accepts the time scales s and d, and P2reg_all reproduces P2est", {
  d <- prep2(panel_with_age(1500, 8), states = c("A", "B", "C"), covariates = "age")
  r <- P2reg(d, h = "B", j = "B", formula = ~ age + log(d) + s)
  expect_true(all(c("log(d)", "s") %in% summary(r)$term))
  expect_equal(dim(predict(r, newdata = data.frame(age = 60, d = 2, s = 4))), c(1L, 1L))
  ## splines in the formula are found through the formula's environment
  expect_s3_class(P2reg(d, h = "A", j = "B", formula = ~ splines::ns(age, 2)), "P2reg")

  all1 <- P2reg_all(d, formula = ~1)
  expect_equal(predict(all1, newdata = data.frame(age = 50)), P2est(d)$P)
  all2 <- P2reg_all(d, formula = ~ age)
  tl <- predict(all2, newdata = data.frame(age = c(40, 80)))
  expect_length(tl, 2L)
  for (P in tl) {
    rs <- apply(P, c(1, 3), sum)
    expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
  }
  expect_equal(nrow(summary(all2)), 4L)   # 2 histories x (intercept, age)
  expect_output(print(all2), "P2reg_all")
  ## covariate-adjusted n-step predictions through the extended Chapman-Kolmogorov relation
  expect_length(ckequations(tl[[1]], "A", "B", "C", nsteps = 4), 4L)
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(r))
  expect_error(plot(P2reg(d, h = "A", j = "B")), "No covariate")
})

test_that("P2reg_all keeps the RPE for histories with too few moves", {
  d <- prep2(panel_with_age(300, 9), states = c("A", "B", "C"), covariates = "age")
  a <- P2reg_all(d, formula = ~ age, min.events = 1e6)
  expect_false(any(a$pairs$modelled))
  expect_length(a$fits, 0L)
  expect_equal(predict(a, newdata = data.frame(age = 60)), P2est(d)$P)
  expect_error(predict(a, newdata = data.frame()), "at least one row")
})

test_that("P2reg_all reports terms that are constant within a history's risk set", {
  pan <- toy_panel(800, seed = 13)
  ## a covariate that is 1 only for some subjects absorbed straight from (A, B):
  ## it is constant (0) in the (B, B) risk set, so its effect is not estimable there
  set.seed(4)
  z <- pan |>
    dplyr::summarise(ex = dplyr::nth(state, 3) == "C", .by = "id") |>
    dplyr::mutate(z = as.integer(ex & stats::runif(dplyr::n()) < 0.5))
  pan <- dplyr::left_join(pan, dplyr::select(z, "id", "z"), by = "id")
  d <- prep2(pan, states = c("A", "B", "C"), covariates = "z")
  a <- P2reg_all(d, formula = ~ z)
  expect_true(any(nzchar(a$pairs$aliased)))
  expect_output(print(a), "not estimable")
  expect_warning(P <- predict(a, newdata = data.frame(z = 1)), "not estimable")
  rs <- apply(P, c(1, 3), sum)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
})

test_that("P2reg works on the covariates of an msprep2() result", {
  w <- data.frame(id = 1:200, ill_time = rep(c(2, 3, 1, 4), 50), ill_status = rep(c(1, 0, 1, 1), 50),
                  dead_time = rep(c(4, 3, 5, 6), 50), dead_status = rep(c(1, 0, 0, 1), 50),
                  age = rep(c(40, 50, 60, 70, 80), 40))
  x <- msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)))
  d <- prep2(x)
  expect_equal(d$covariates, "age")
  expect_s3_class(P2reg(d, h = "healthy", j = "healthy", formula = ~ age), "P2reg")
})
