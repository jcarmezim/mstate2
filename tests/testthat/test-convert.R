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

test_that("sojourn_to_panel errors when a required column is missing", {
  df <- data.frame(id = 1, t.nosp = 1, disch.s = 1, death.s = 0)  # no t.sp
  segs <- c(NSP = "t.nosp", SP = "t.sp")
  expect_error(sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s")),
               "Missing column")
})

test_that("a subject with no duration and no absorbing event contributes no rows", {
  df <- data.frame(id = c(1, 2), t.nosp = c(0, 2), disch.s = c(0, 1), death.s = c(0, 0))
  segs <- c(NSP = "t.nosp")
  panel <- sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s"))
  expect_setequal(unique(panel$id), 2)   # subject 1 (zero duration, no event) is dropped
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
  cmp <- compare2(P2est(d), c("A", "B"), "B", "B", nsteps = 9)
  expect_s3_class(cmp, "msm2pred")
  os <- overlap_step(cmp)
  expect_true(os$separated_steps >= 1)          # the 0.6 vs 0.3 gap is detected early
})

test_that("simulate2 falls back to integer-string labels without tensor dimnames", {
  P <- toy_tensor(); dimnames(P) <- NULL
  first <- toy_first(); dimnames(first) <- NULL
  set.seed(7)
  panel <- simulate2(50, P, first, init = c(1, 0, 0))
  expect_setequal(as.character(unique(panel$state)), c("1", "2", "3"))
})

test_that("a first-order row with no outgoing mass ends that individual's path (no maxT wait)", {
  ## first["B", ] is all zero (toy_first only defines the "A" row): anyone
  ## starting at B has nowhere to go on their first step. They should get
  ## no further observations right away, without waiting out `maxT`.
  set.seed(9)
  expect_silent(
    panel <- simulate2(20, toy_tensor(), toy_first(), init = c(A = 0, B = 1, C = 0), maxT = 5)
  )
  expect_true(all(table(panel$id) == 1))
})

test_that("simulate2 warns and stops at maxT when absorption never occurs", {
  ## A <-> B oscillate forever (neither has a self-transition = 1, so neither
  ## is auto-detected as absorbing); C stays absorbing but is never reached.
  st   <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["A", "B", ] <- 1; tens["B", "A", ] <- 1; tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  expect_warning(
    simulate2(5, tens, first, init = c(A = 1, B = 0, C = 0), maxT = 3),
    "maxT"
  )
})

test_that("simulate2 errors on an unnamed `entry` vector instead of silently entering nobody", {
  st   <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  expect_error(simulate2(10, tens, first, entry = c(0.5, 0.3, 0), maxT = 5), "named")
})

test_that("simulate2 with staggered entry records one row per subject per time step", {
  st   <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  set.seed(1)
  panel <- simulate2(200, tens, first, entry = c(A = 0.2))
  ## the entry state X_0 and the first move X_1 must not share a time value
  gaps <- panel[, .(ok = all(diff(time) == 1L)), by = id]
  expect_true(all(gaps$ok))
  expect_true(all(panel[, time[1L] >= 1L, by = id]$V1))    # nobody enters at t = 0
  expect_silent(prep2(panel, states = st))                  # no non-consecutive warning
})

test_that("as_tmat returns an mstate transition matrix", {
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
  expect_false(is.na(tm["A", "B"]))             # first-move-only transition kept
  expect_true(is.na(tm["B", "A"]))
  if (requireNamespace("mstate", quietly = TRUE))
    expect_identical(tm, mstate::transMat(list(A = 2L, B = 3L, C = integer(0)),
                                          names = st))
})

test_that("sojourn_to_panel warns when a visit rounds to 0 and is dropped", {
  df <- data.frame(id = 1:2, t.nosp = c(3, 2), t.sp = c(0.3, 0.5), t.nimv = c(2, 0),
                   disch.s = c(1, 1), death.s = c(0, 0))
  segs <- c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv")
  expect_warning(panel <- sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s")),
                 "1 visit\\(s\\).*SP: 1")
  ## subject 1 loses SP (0.3 -> 0): NSP -> NIMV directly; subject 2 keeps SP (0.5 -> 1)
  expect_equal(as.character(panel[id == 1, state]), c("NSP", "NSP", "NSP", "NIMV", "NIMV", "Disch"))
  expect_true("SP" %in% panel[id == 2, state])
  ## half-day durations, as in DIVINE, lose nothing
  expect_no_warning(sojourn_to_panel(df[2, ], "id", segs, c(Disch = "disch.s", Death = "death.s")))
})

test_that("from_msdata warns when two transitions round to the same time", {
  msdat <- data.frame(id = c(1, 1, 1), from = c(1, 2, 3), to = c(2, 3, 4),
                      Tstart = c(0, 2.2, 2.4), Tstop = c(2.2, 2.4, 5), status = c(1, 1, 1))
  expect_warning(panel <- from_msdata(msdat), "1 state visit\\(s\\)")
  expect_false(2 %in% panel$state)          # state 2 (2.2 to 2.4) never occupies a day
  expect_no_warning(from_msdata(transform(msdat, Tstop = c(2, 4, 5), Tstart = c(0, 2, 4))))
})
