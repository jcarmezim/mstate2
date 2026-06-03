test_that("compare2 stacks groups and carries attributes", {
  cmp <- compare2(toy_fit(), c("A","B"), "B", "B", nsteps = 9)
  expect_s3_class(cmp, "msm2pred")
  expect_setequal(unique(cmp$h), c("A","B"))
  expect_equal(attr(cmp, "l"), "B")
  expect_true(attr(cmp, "bounds"))
})

test_that("bounds = FALSE omits interval columns and disables interval tools", {
  cmp <- compare2(toy_fit(), c("A","B"), "B", "B", bounds = FALSE)
  expect_false(any(c("lower","upper") %in% names(cmp)))
  expect_error(overlap_step(cmp), "bounds")
  expect_warning(plot(cmp, type = "interval"), "bounds")  # falls back to curves
})

test_that("compare2 errors on unknown states", {
  expect_error(compare2(toy_fit(), c("A","Z"), "B", "B"), "not found")
  expect_error(compare2(toy_fit(), "A", "B", "Z"), "not found")
})

test_that("overlap_step detects the early separation of 0.6 vs 0.3", {
  os <- overlap_step(compare2(toy_fit(), c("A","B"), "B", "B", nsteps = 9))
  expect_equal(length(os$overlap), 9)
  expect_gte(os$separated_steps, 1)
  expect_true(os$separation[1] > 0)               # separated at step 1
})

test_that("overlap_step requires exactly two groups", {
  expect_error(overlap_step(compare2(toy_fit(), "A", "B", "B")), "two")
})

test_that("summary.msm2pred returns the overlap result invisibly", {
  os <- summary(compare2(toy_fit(), c("A","B"), "B", "B"))
  expect_true(is.list(os) && !is.null(os$separated_steps))
})
