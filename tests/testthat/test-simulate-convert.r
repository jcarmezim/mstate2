test_that("rnd rounds halves away from zero", {
  expect_equal(rnd(c(0.5, 1.5, 2.5, 0.4, -0.5)), c(1, 2, 3, 0, -1))
})

test_that("simulate2 returns a valid panel ending in absorption", {
  panel <- toy_panel(n = 500)
  expect_true(all(c("id","time","state") %in% names(panel)))
  expect_s3_class(panel, "tbl_df")
  last <- panel |> dplyr::summarise(s = dplyr::last(state), .by = id)
  expect_true(all(last$s == "C"))                 # everyone absorbed
})

test_that("staggered entry spreads exposure across time", {
  set.seed(7)
  p <- simulate2(800, toy_tensor(), toy_first(), entry = c(A = 0.1))
  d <- prep2(p, states = c("A","B","C"), check.consecutive = FALSE)   # <- silence
  expect_gt(nrow(dplyr::filter(d$Y, h == "B", j == "B")), 1)
})

test_that("sojourn_to_panel reconstructs the expected trajectory", {
  df <- data.frame(id = 1, t.nosp = 3, t.sp = 0.5, t.nimv = 2, t.mv = 0,
                   t.recov = 1, disch.s = 1, death.s = 0)
  segs <- c(NSP="t.nosp", SP="t.sp", NIMV="t.nimv", IMV="t.mv", Recov="t.recov")
  pan <- sojourn_to_panel(df, "id", segs, c(Disch="disch.s", Death="death.s"))
  expect_equal(as.character(pan$state),
               c("NSP","NSP","NSP","SP","NIMV","NIMV","Recov","Disch"))
})

test_that("a zero-duration segment is skipped", {
  df <- data.frame(id = 1, t.nosp = 0, t.sp = 2, t.nimv = 0, t.mv = 0,
                   t.recov = 0, disch.s = 0, death.s = 1)
  segs <- c(NSP="t.nosp", SP="t.sp", NIMV="t.nimv", IMV="t.mv", Recov="t.recov")
  pan <- sojourn_to_panel(df, "id", segs, c(Disch="disch.s", Death="death.s"))
  expect_equal(as.character(pan$state), c("SP","SP","Death"))
})

test_that("simulate -> estimate -> predict round-trips to the truth", {
  fit <- toy_fit(n = 6000, seed = 11)
  expect_equal(fit$P["B","B","A"], 0.6, tolerance = 0.05)
  v <- ckequations(fit, "A", "B", "C", nsteps = 30)   # eventual absorption
  expect_equal(v[30], 1, tolerance = 0.02)
})
