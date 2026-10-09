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
  p <- dplyr::filter(p, time != 2)                 # punch a hole
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

test_that("triples carry the time spent in the current state, d", {
  panel <- data.frame(id = 1, time = 0:5, state = c("A", "B", "B", "B", "C", "C"))
  d <- prep2(panel, check.consecutive = FALSE)
  ## l at s = 2..5; j at s - 1 = B, B, B, C; d = units in j up to s - 1
  expect_equal(d$triples$s, 2:5)
  expect_equal(d$triples$d, c(1L, 2L, 3L, 1L))
  ## a subject that starts in j is counted from its first observation
  d2 <- prep2(data.frame(id = 1, time = 0:3, state = c("B", "B", "B", "C")))
  expect_equal(d2$triples$d, c(2L, 3L))
})

test_that("pairs count every one-step move, including the first one", {
  panel <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2),
                      state = c("A", "B", "C", "A", "A", "B"))
  d <- prep2(panel)
  pr <- as.data.frame(d$pairs)
  expect_equal(nrow(pr), 3)   # A -> A, A -> B, B -> C
  expect_equal(pr$N[pr$from == "A" & pr$to == "B"], 2L)
  ## the first move A -> B of subject 1 has no triple, but it is a pair
  expect_false(any(d$triples$j == "A" & d$triples$l == "B" & d$triples$id == 1))
})

test_that("covariates are carried into the triples", {
  panel <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2),
                      state = c("A", "B", "C", "A", "B", "B"), age = rep(c(50, 70), each = 3))
  d <- prep2(panel, covariates = "age")
  expect_equal(d$covariates, "age")
  expect_equal(d$triples$age, c(50, 70))
  expect_output(print(d), "covariates      : age")
  expect_error(prep2(panel, covariates = "sex"), "sex")
  expect_error(prep2(transform(panel, s = 1), covariates = "s"), "reserved")
  panel$age[2] <- 51
  expect_warning(prep2(panel, covariates = "age"), "not constant")
  ## an msprep2() result passes its covariates by default
  w <- data.frame(id = 1:3, ill_time = c(2, 3, 1), ill_status = c(1, 0, 1),
                  dead_time = c(4, 3, 5), dead_status = c(1, 0, 0), age = c(40, 50, 60))
  x <- msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                dead = Surv(dead_time, dead_status)))
  expect_equal(prep2(x)$covariates, "age")
})
