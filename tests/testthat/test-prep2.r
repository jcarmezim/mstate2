test_that("prep2 builds counting processes with the Y = sum(N) identity", {
  d <- prep2(toy_panel(), states = c("A", "B", "C"))
  expect_s3_class(d, "msm2data")
  expect_equal(sum(d$Y$Y), sum(d$N$N))            # Y is N marginalized over l
  expect_identical(d$absorbing, "C")
  expect_equal(d$time.range[1], 0L)        # panel starts at day 0
  expect_gte(min(d$N$s), 2L)               # first triple needs 3 observations
})

test_that("prep2 recovers known transition probabilities", {
  e <- P2est(prep2(toy_panel(), states = c("A","B","C")))$estimate
  g <- function(h, j, l) e$p[e$h == h & e$j == j & e$l == l]
  expect_equal(g("A","B","B"), 0.6, tolerance = 0.05)
  expect_equal(g("B","B","B"), 0.3, tolerance = 0.05)
  expect_equal(g("A","B","B") + g("A","B","C"), 1, tolerance = 1e-9)
})

test_that("prep2 errors on NA unless drop.na = TRUE", {
  p <- toy_panel(); p$state[5] <- NA
  expect_error(prep2(p, states = c("A","B","C")), "NA")
  expect_warning(prep2(p, states = c("A","B","C"), drop.na = TRUE), "Dropping")
})

test_that("prep2 warns on non-consecutive time grids", {
  p <- toy_panel()
  p <- p[p$time != 2]                              # punch a hole
  expect_warning(prep2(p, states = c("A","B","C")), "non-consecutive")
})

test_that("prep2 errors when a subject has < 3 observations", {
  p <- data.frame(id = 1, time = 0:1, state = c("A","B"))
  expect_error(prep2(p, states = c("A","B","C")), "3 observations")
})

test_that("summary.msm2data reports per-pair exposure", {
  d <- prep2(toy_panel(), states = c("A","B","C"))
  ex <- summary(d)
  expect_true(all(c("total_at_risk", "s_min", "s_max") %in% names(ex)))
  expect_true(any(ex$h == "A" & ex$j == "B"))
})

test_that("prep2 auto-detects `states` from the observed data when not supplied", {
  d <- prep2(toy_panel())    # no states = argument
  expect_setequal(d$states, c("A", "B", "C"))
})

test_that("prep2 errors when a required column is missing from `data`", {
  p <- toy_panel()
  expect_error(prep2(p, id = "subject"), "not found")
})

test_that("prep2 errors when an observed state is not in the supplied `states`", {
  p <- toy_panel()
  expect_error(prep2(p, states = c("A", "B")), "not in")   # 'C' is observed but omitted
})

test_that("from_msdata errors when a required column is missing", {
  md <- data.frame(id = 1, from = 1, to = 2, Tstart = 0, status = 1)  # no Tstop
  expect_error(from_msdata(md), "missing column")
})

test_that("from_msdata reconstructs the realized path and prep2 auto-dispatches", {
  md <- data.frame(
    id     = c(1,1,1, 2,2),
    from   = c(1,1,2, 1,1), to = c(2,3,3, 2,3),
    Tstart = c(0,0,3, 0,0), Tstop = c(3,3,5, 2,2),
    status = c(1,0,1, 0,1))
  class(md) <- c("msdata", "data.frame")
  pan <- from_msdata(md)
  expect_equal(as.character(pan[pan$id == 1]$state), c("1","1","1","2","2","3"))
  expect_equal(as.character(pan[pan$id == 2]$state), c("1","1","3"))
  expect_s3_class(prep2(md, states = c("1","2","3")), "msm2data")  # dispatch
})
