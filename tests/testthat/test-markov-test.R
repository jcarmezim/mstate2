test_that("markov_test flags a genuinely history-dependent transition", {
  fit <- toy_fit(n = 6000, seed = 21)   # A -> B -> B 0.6 vs B -> B -> B 0.3
  mt <- markov_test(fit)
  expect_s3_class(mt, "markov_test")
  row <- mt[mt$j == "B" & mt$l == "B", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$K, 2)
  expect_equal(row$df, 1)
  expect_true(row$p.value < 0.001)
  ## Cochran's Q by hand
  e <- fit$estimate[fit$estimate$j == "B" & fit$estimate$l == "B", ]
  w <- 1 / e$se^2
  pbar <- sum(w * e$p) / sum(w)
  expect_equal(row$statistic, sum(w * (e$p - pbar)^2))
  expect_equal(row$pooled, pbar)
})

test_that("markov_test(method = 'lrt') matches the likelihood-ratio statistic by hand", {
  fit <- toy_fit(n = 1500, seed = 3)
  mt <- markov_test(fit, method = "lrt")
  e <- fit$estimate
  N <- e$n.trans[e$l == "C"]
  Y <- e$at.risk[e$l == "C"]
  pbar <- sum(N) / sum(Y)
  G2 <- 2 * sum(N * log(N / (Y * pbar)) + (Y - N) * log((Y - N) / (Y * (1 - pbar))))
  expect_equal(mt$statistic[mt$l == "C"], G2, tolerance = 1e-10)
  expect_equal(mt$pooled[mt$l == "C"], pbar)
  ## with two destinations the joint (per-state) test is the same statistic
  ms <- markov_test(fit, method = "lrt", by = "state")
  expect_equal(ms$statistic, G2, tolerance = 1e-10)
  expect_equal(ms$df, 1L)
  expect_s3_class(ms, "markov_test")
  expect_output(print(ms), "per current state")
  expect_output(print(mt), "Likelihood-ratio")
})

test_that("the LRT pooled estimate is the first-order P1est estimate", {
  d <- prep2(toy_panel(1500, seed = 4), states = c("A", "B", "C"))
  mt <- markov_test(P2est(d), method = "lrt")
  expect_equal(mt$pooled[mt$l == "C"], P1est(d, B = 0)$P["B", "C"])
})

test_that("markov_test does not reject a first-order chain", {
  st <- c("A", "B", "C")
  P1 <- matrix(c(.70, .25, .05, .20, .60, .20, 0, 0, 1), 3, byrow = TRUE, dimnames = list(st, st))
  set.seed(21)
  pan <- simulate2(3000, array(P1, c(3, 3, 3), dimnames = list(st, st, st)), P1, init = c(A = 1, B = 0, C = 0))
  fit <- P2est(prep2(pan, states = st))
  expect_true(all(markov_test(fit, method = "lrt", by = "state")$p.value > 0.001))
  expect_true(all(markov_test(fit)$p.value > 0.001))
})

test_that("markov_test adds Holm-adjusted p-values", {
  fit <- P2est(prep2(toy_panel(2000, seed = 4), states = c("A", "B", "C")))
  mt <- markov_test(fit, method = "lrt")
  expect_equal(mt$p.adj, stats::p.adjust(mt$p.value, "holm"))
  expect_equal(markov_test(fit, method = "lrt", p.adjust = "none")$p.adj, mt$p.value)
})

test_that("markov_test checks its arguments and returns an empty table when nothing is testable", {
  fit <- toy_fit()
  expect_error(markov_test(fit, min.h = 1), "min.h")
  expect_error(markov_test(fit, by = "state"), "lrt")
  expect_warning(mt <- markov_test(fit, min.h = 5), "min.h")
  expect_equal(nrow(mt), 0)
  expect_true(all(c("j", "l", "K", "df", "statistic", "p.value", "pooled", "p.adj") %in% names(mt)))
  expect_output(print(mt), "0 test")
})
