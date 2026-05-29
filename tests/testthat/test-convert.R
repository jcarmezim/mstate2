test_that("rnd rounds half up", {
  expect_equal(rnd(c(0.5, 1.5, 2.5, 0.4, -0.5)), c(1, 2, 3, 0, -1))
})

test_that("sojourn_to_panel reconstructs the expected trajectory", {
  df <- data.frame(id = 1,
                   t.nosp = 3, t.sp = 0.5, t.nimv = 2, t.mv = 0, t.recov = 1,
                   disch.s = 1, death.s = 0)
  segs <- c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")
  panel <- sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s"))
  ## NSP x3, SP x1 (0.5 -> 1), NIMV x2, Recov x1, then Disch
  expect_equal(as.character(panel$state),
               c("NSP","NSP","NSP","SP","NIMV","NIMV","Recov","Disch"))
  expect_equal(panel$time, 0:7)
})

test_that("a zero-duration segment is skipped", {
  df <- data.frame(id = 1, t.nosp = 0, t.sp = 2, t.nimv = 0, t.mv = 0, t.recov = 0,
                   disch.s = 0, death.s = 1)
  segs <- c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")
  panel <- sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s"))
  expect_equal(as.character(panel$state), c("SP", "SP", "Death"))
})

test_that("compare2 / overlap_step detect a separated then overlapping pair", {
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  set.seed(5)
  d   <- prep2(simulate2(4000, tens, first, init = c(A = 1, B = 0, C = 0)), states = st)
  cmp <- compare2(P2est(d, "RPE"), c("A", "B"), "B", "B", nsteps = 9)
  expect_s3_class(cmp, "msm2pred")
  os <- overlap_step(cmp)
  expect_true(os$separated_steps >= 1)          # the 0.6 vs 0.3 gap is detected early
})

test_that("as_tmat returns an mstate transition matrix", {
  skip_if_not_installed("mstate")
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  set.seed(6)
  d <- prep2(simulate2(1500, tens, first, init = c(A = 1, B = 0, C = 0)), states = st)
  tm <- as_tmat(d)
  expect_true(is.matrix(tm))
  expect_equal(dim(tm), c(3, 3))
  expect_false(is.na(tm["B", "C"]))             # B -> C is observed
})
