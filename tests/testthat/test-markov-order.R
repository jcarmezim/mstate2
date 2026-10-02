test_that("markov_test flags a genuinely history-dependent transition", {
  fit <- toy_fit(n = 6000, seed = 21)   # A->B->B ~ 0.6, B->B->B ~ 0.3: clearly different
  mt  <- markov_test(fit)
  expect_s3_class(mt, "markov_test")

  row <- mt[mt$j == "B" & mt$l == "B", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$K, 2)
  expect_equal(row$df, 1)
  expect_true(row$p.value < 0.001)
})

test_that("markov_test rejects min.h < 2", {
  expect_error(markov_test(toy_fit(), min.h = 1), "min.h")
})

test_that("markov_test warns and returns empty when min.h can't be met", {
  fit <- toy_fit()
  expect_warning(mt <- markov_test(fit, min.h = 5), "min.h")
  expect_equal(nrow(mt), 0)
})

test_that("print.markov_test runs without error", {
  fit <- toy_fit(n = 2000, seed = 22)
  expect_output(print(markov_test(fit)), "markov_test")
})

## --- compare_order() / divergence(): a hand-built stand-in for an
## mstate::probtrans object, so these tests don't require mstate to be
## installed (compare_order() itself never calls into the mstate namespace).
fake_probtrans <- function(states = c("A", "B", "C"), which = "B",
                           predt = 1, ntime = 12, p_start = 0.5, p_end = 0.3,
                           se = 0.02) {
  M   <- length(states)
  idx <- match(which, states)
  df  <- data.frame(time = seq(predt, predt + ntime - 1))
  p   <- seq(p_start, p_end, length.out = ntime)
  for (k in seq_len(M)) {
    df[[paste0("pstate", k)]] <- if (k == idx) p else (1 - p) / (M - 1)
    df[[paste0("se", k)]]     <- se
  }
  pt <- vector("list", M)
  pt[[idx]] <- df
  class(pt) <- "probtrans"
  pt
}

test_that("compare_order requires a probtrans object", {
  fit <- toy_fit()
  expect_error(compare_order(fit, list(1, 2, 3), h = "A", j = "B", l = "B"),
               "probtrans")
})

test_that("compare_order errors when pt is missing expected columns", {
  fit <- toy_fit()
  pt  <- fake_probtrans()
  pt[[2]]$se2 <- NULL
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B"), "probtrans column")
})

test_that("compare_order errors on an unknown j or l", {
  fit <- toy_fit()
  pt  <- fake_probtrans()
  expect_error(compare_order(fit, pt, h = "A", j = "Z", l = "B"), "not found")
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "Z"), "not found")
})

test_that("compare_order errors when pt does not extend far enough", {
  fit <- toy_fit()
  pt  <- fake_probtrans(ntime = 3)   # only 3 times from predt = 1
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B", nsteps = 20),
               "does not extend")
})

test_that("compare_order rejects an out-of-range conf.level", {
  fit <- toy_fit()
  pt  <- fake_probtrans()
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B", conf.level = 1.5),
               "conf.level")
})

test_that("compare_order errors when pt has fewer starting-state components than j needs", {
  fit <- toy_fit()
  pt  <- fake_probtrans()[1]        # only 1 component; j = "B" resolves to index 2
  class(pt) <- "probtrans"
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B"),
               "starting-state component")
})

test_that("compare_order errors when pt[[j]] has no valid time values", {
  fit <- toy_fit()
  pt  <- fake_probtrans()
  pt[[2]] <- pt[[2]][0, ]            # empty: no landmark row at all
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B"), "no valid")
})

test_that("compare_order builds a baseline + per-h msm2pred object", {
  fit <- toy_fit(n = 4000, seed = 23)
  pt  <- fake_probtrans()
  cmp <- compare_order(fit, pt, h = c("A", "B"), j = "B", l = "B", nsteps = 6)

  expect_s3_class(cmp, "msm2pred")
  expect_setequal(unique(cmp$h), c("1st order (mstate)", "A", "B"))
  expect_equal(attr(cmp, "baseline"), "1st order (mstate)")
  expect_equal(nrow(cmp), 3 * 6)
  expect_true(all(cmp$lower <= cmp$estimate + 1e-9))
  expect_true(all(cmp$estimate <= cmp$upper + 1e-9))
})

test_that("divergence() reports a memory horizon per preceding state", {
  fit <- toy_fit(n = 4000, seed = 23)
  pt  <- fake_probtrans()
  cmp <- compare_order(fit, pt, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
  div <- divergence(cmp)

  expect_equal(nrow(div), 2)
  expect_setequal(div$h, c("A", "B"))
  expect_true(all(div$separated_steps >= 0 & div$separated_steps <= 6))
  expect_true(all(div$max_abs_diff >= 0))
})

test_that("divergence() errors on an object without a baseline group", {
  cmp <- compare2(toy_fit(), c("A", "B"), "B", "B", nsteps = 6)
  expect_error(divergence(cmp), "baseline")
})

test_that("divergence() errors when no second-order group remains besides the baseline", {
  df <- data.frame(h = "1st order (mstate)", n = 1:3,
                   estimate = c(0.5, 0.4, 0.3),
                   lower = c(0.4, 0.3, 0.2), upper = c(0.6, 0.5, 0.4))
  x <- structure(df, class = c("msm2pred", "data.frame"), j = "B", l = "B",
                bounds = TRUE, estimator = "RPE", conf.level = 0.95,
                baseline = "1st order (mstate)")
  expect_error(divergence(x), "No second-order")
})
