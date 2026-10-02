test_that("RPE estimate and SE match the hand-rolled formulas", {
  panel <- toy_panel(n = 3000, seed = 2)
  d <- prep2(panel, states = c("A","B","C"))
  rpe <- P2est(d)

  P <- data.table::as.data.table(panel); data.table::setorder(P, id, time)
  P[, c("hh","jj") := .(data.table::shift(state,2), data.table::shift(state,1)), by = id]
  tr <- P[!is.na(hh) & !is.na(jj)]; M <- max(tr$time); Y <- Nn <- numeric(0)
  for (s in 2:M) {
    Y[s-1]  <- nrow(tr[time==s & hh=="B" & jj=="B"])
    Nn[s-1] <- nrow(tr[time==s & hh=="B" & jj=="B" & state=="B"])
  }
  k <- Y > 0; Y <- Y[k]; Nn <- Nn[k]
  rpe_p  <- sum(Nn)/sum(Y)
  rpe_se <- sqrt(rpe_p*(1-rpe_p)/sum(Y))

  er <- rpe$estimate
  expect_equal(er$p[er$h=="B"&er$j=="B"&er$l=="B"],  rpe_p,  tolerance = 1e-10)
  expect_equal(er$se[er$h=="B"&er$j=="B"&er$l=="B"], rpe_se, tolerance = 1e-10)
})

test_that("estimate table carries the transparency columns", {
  e <- toy_fit()$estimate
  expect_true(all(c("n.trans","at.risk") %in% names(e)))
  expect_true(all(e$n.trans <= e$at.risk))
})

test_that("wald CIs clip to [0,1]; logit CIs stay strictly inside", {
  w <- P2est(prep2(toy_panel(), states=c("A","B","C")), ci = "wald")
  lo <- P2est(prep2(toy_panel(), states=c("A","B","C")), ci = "logit")
  expect_true(all(w$estimate$lower >= 0 & w$estimate$upper <= 1))
  inside <- lo$estimate$p > 0 & lo$estimate$p < 1
  expect_true(all(lo$estimate$lower[inside] > 0 & lo$estimate$upper[inside] < 1))
  expect_identical(lo$ci, "logit")
})

test_that("P2est rejects an out-of-range conf.level instead of silently returning NaN CIs", {
  d <- prep2(toy_panel(), states = c("A", "B", "C"))
  expect_error(P2est(d, conf.level = 1.5), "conf.level")
  expect_error(P2est(d, conf.level = 0),   "conf.level")
  expect_error(P2est(d, conf.level = NA),  "conf.level")
})

test_that("absorbing states satisfy P[a,a,h] = 1 for every previous state", {
  fit <- toy_fit()
  expect_true(all(fit$P["C","C", ] == 1))
  expect_equal(sum(fit$P["C", -match("C", fit$states), ]), 0)
})

test_that("point estimates and tensor agree", {
  fit <- toy_fit(); e <- fit$estimate
  row <- e[e$h=="A" & e$j=="B" & e$l=="B", ]
  expect_equal(fit$P["B","B","A"], row$p)
})

test_that("print.P2est reports the estimator, states and sample size", {
  fit <- toy_fit()
  expect_output(print(fit), "P2est")
  expect_output(print(fit), "RPE")
  expect_output(print(fit), "A, B, C")
})
