## A hand-built stand-in for an mstate::probtrans object, so these tests do not
## need mstate (compare_order() only uses the structure of the object).
fake_probtrans <- function(states = c("A", "B", "C"), which = "B", predt = 1, ntime = 12,
                           p_start = 0.5, p_end = 0.3, se = 0.02) {
  M <- length(states)
  idx <- match(which, states)
  df <- data.frame(time = seq(predt, predt + ntime - 1))
  p <- seq(p_start, p_end, length.out = ntime)
  for (k in seq_len(M)) {
    df[[paste0("pstate", k)]] <- if (k == idx) p else (1 - p) / (M - 1)
    df[[paste0("se", k)]] <- se
  }
  pt <- vector("list", M)
  pt[[idx]] <- df
  class(pt) <- "probtrans"
  pt
}

test_that("probtrans2 has mstate's probtrans structure and plots", {
  d <- prep2(toy_panel(800, seed = 6), states = c("A", "B", "C"))
  fit <- P2est(d)
  pt <- probtrans2(fit, h = "A", j = "B", nsteps = 5)
  expect_s3_class(pt, "probtrans")
  expect_identical(names(pt[[2]]), c("time", paste0("pstate", 1:3), paste0("se", 1:3)))
  expect_equal(pt$predt, 1)
  expect_equal(pt$direction, "forward")
  expect_equal(unname(rowSums(pt[[2]][, 2:4])), rep(1, 6))
  expect_equal(pt[[2]]$pstate3[-1], ckequations(fit, "A", "B", "C", nsteps = 5))
  expect_null(pt[[1]])
  expect_true(all(is.na(pt[[2]]$se3[-1])))
  ## transMat layout: the estimated moves (here only B -> C) numbered row by row
  expect_equal(dimnames(pt$trans), list(from = c("A", "B", "C"), to = c("A", "B", "C")))
  expect_equal(pt$trans["B", "C"], 1L)
  expect_equal(sum(!is.na(pt$trans)), 1L)
  ## first-order predictions from every transient state, with bootstrap SEs
  p1 <- probtrans2(P1est(d, B = 40, seed = 1), nsteps = 4)
  expect_false(is.null(p1[[2]]))
  expect_true(all(p1[[2]]$se3[-1] > 0))
  expect_true(is.na(p1$h))
  ## second order from a P2boot: bootstrap SEs
  expect_true(all(probtrans2(P2boot(d, B = 20, seed = 1), "A", "B", nsteps = 3)[[2]]$se3[-1] > 0))
  pdf(NULL)
  on.exit(dev.off())
  for (ty in c("filled", "stacked", "single", "separate")) expect_silent(plot(pt, type = ty))
  expect_silent(plot(p1, from = "B", type = "stacked", ord = 3:1))
  expect_error(plot(pt, from = 1), "No predictions")
  expect_error(probtrans2(fit, j = "B"), "required")
  expect_error(probtrans2(P1est(d, B = 0), j = "Z"), "not found")
  expect_error(probtrans2(list(states = "A"), j = 1), "must be")
  expect_error(probtrans2(fit, "A", "B", nsteps = 0), "nsteps")
})

test_that("as_tmat gives the observed moves, first moves included, in transMat format", {
  d <- prep2(toy_panel(200, seed = 1), states = c("A", "B", "C"))
  tm <- as_tmat(d)
  expect_equal(dimnames(tm), list(from = c("A", "B", "C"), to = c("A", "B", "C")))
  expect_equal(tm["A", "B"], 1L)   # only ever a first move
  expect_equal(tm["B", "C"], 2L)
  expect_equal(sum(!is.na(tm)), 2L)
  expect_identical(as_tmat(toy_panel(200, seed = 1), states = c("A", "B", "C")), tm)
  skip_if_not_installed("mstate")
  expect_identical(unname(tm), unname(mstate::transMat(list(2, 3, integer(0)))))
})

test_that("compare_order with an mstate-like baseline builds a baseline + per-h comparison", {
  fit <- toy_fit(n = 4000, seed = 23)
  cmp <- compare_order(fit, fake_probtrans(), h = c("A", "B"), j = "B", l = "B", nsteps = 6)
  expect_s3_class(cmp, "msm2pred")
  expect_setequal(unique(cmp$h), c("1st order (mstate)", "A", "B"))
  expect_equal(attr(cmp, "baseline"), "1st order (mstate)")
  expect_equal(nrow(cmp), 3 * 6)
  expect_true(all(cmp$lower <= cmp$estimate + 1e-9 & cmp$estimate <= cmp$upper + 1e-9))
  div <- divergence(cmp)
  expect_equal(nrow(div), 2)
  expect_setequal(div$h, c("A", "B"))
  expect_true(all(div$separated_steps >= 0 & div$separated_steps <= 6))
  expect_true(all(div$max_abs_diff >= 0))
  expect_true(all(is.na(div$diff_steps)))
})

test_that("compare_order checks the first-order baseline", {
  fit <- toy_fit()
  pt <- fake_probtrans()
  expect_error(compare_order(fit, list(1, 2, 3), h = "A", j = "B", l = "B"), "probtrans")
  bad <- pt
  bad[[2]]$se2 <- NULL
  expect_error(compare_order(fit, bad, h = "A", j = "B", l = "B"), "probtrans column")
  expect_error(compare_order(fit, pt, h = "A", j = "Z", l = "B"), "not found")
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "Z"), "not found")
  expect_error(compare_order(fit, fake_probtrans(ntime = 3), h = "A", j = "B", l = "B", nsteps = 20), "does not reach")
  expect_error(compare_order(fit, pt, h = "A", j = "B", l = "B", conf.level = 1.5), "conf.level")
  short <- structure(pt[1], class = "probtrans")
  expect_error(compare_order(fit, short, h = "A", j = "B", l = "B"), "no predictions")
  empty <- pt
  empty[[2]] <- empty[[2]][0, ]
  expect_error(compare_order(fit, empty, h = "A", j = "B", l = "B"), "no valid")
  d <- prep2(toy_panel(300, seed = 2))
  expect_error(compare_order(P2est(d), P1est(d, B = 0), h = "A", j = "B", l = "C"), "bootstrap")
})

test_that("compare_order accepts a P1est baseline; divergence_table covers all transitions", {
  d <- prep2(toy_panel(800, seed = 7), states = c("A", "B", "C"))
  fit <- P2est(d)
  p1 <- P1est(d, B = 60, seed = 2)
  co <- compare_order(fit, p1, h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  expect_true("1st order (discrete)" %in% co$h)
  expect_equal(attr(co, "baseline"), "1st order (discrete)")
  expect_equal(co$estimate[co$h == "A"], ckequations(fit, "A", "B", "C", nsteps = 5))
  expect_equal(nrow(divergence(co)), 2L)
  dt <- divergence_table(fit, p1, nsteps = 5)
  expect_true(all(c("j", "l", "h", "separated_steps", "diff_steps", "max_abs_diff") %in% names(dt)))
  expect_setequal(unique(dt$l), c("B", "C"))   # A is never reached from B
  expect_true(all(diff(dt$separated_steps) <= 0))
  expect_error(divergence_table(fit, p1, min.h = 3), "min.h")
})

test_that("divergence() checks for a baseline and a second-order group", {
  expect_error(divergence(compare2(toy_fit(), c("A", "B"), "B", "B", nsteps = 6)), "baseline")
  x <- structure(data.frame(h = "1st order (mstate)", n = 1:3, estimate = c(.5, .4, .3),
                            lower = c(.4, .3, .2), upper = c(.6, .5, .4)),
                 class = c("msm2pred", "data.frame"), j = "B", l = "B", bounds = TRUE,
                 estimator = "RPE", conf.level = 0.95, baseline = "1st order (mstate)")
  expect_error(divergence(x), "No second-order")
})

test_that("paired bootstrap difference with shared replicates", {
  d <- prep2(toy_panel(2000, seed = 5), states = c("A", "B", "C"))
  bt <- P2boot(d, B = 100, seed = 1)
  cmp <- compare2(bt, h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  os <- overlap_step(cmp)
  expect_false(is.na(os$diff_steps))
  expect_equal(nrow(os$diff), 5L)
  ## the true curves differ a lot: the paired test detects it at least as long as the (conservative) overlap criterion
  expect_gte(os$diff_steps, os$separated_steps)
  expect_output(summary(cmp), "paired bootstrap")
  ## evolution intervals: no replicates, no paired test
  os0 <- overlap_step(compare2(P2est(d), c("A", "B"), "B", "C"))
  expect_true(is.na(os0$diff_steps))
  expect_null(os0$diff)
  ## P1est built from the same P2boot -> paired; from its own bootstrap -> not
  dv1 <- divergence(compare_order(bt, P1est(bt), h = c("A", "B"), j = "B", l = "C", nsteps = 5))
  expect_false(anyNA(dv1$diff_steps))
  dv2 <- divergence(compare_order(bt, P1est(d, B = 50, seed = 2), h = c("A", "B"), j = "B", l = "C", nsteps = 5))
  expect_true(all(is.na(dv2$diff_steps)))
})
