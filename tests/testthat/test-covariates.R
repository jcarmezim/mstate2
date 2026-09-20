test_that("prep2 carries baseline covariates through to $triples", {
  panel <- toy_panel(n = 500, seed = 3)
  ids   <- unique(panel$id)
  cov   <- data.frame(id = ids, age = sample(40:80, length(ids), replace = TRUE))
  panel <- merge(panel, cov, by = "id")

  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")

  expect_true("age" %in% names(d$triples))
  expect_identical(d$covariates, "age")
  expect_equal(nrow(d$triples), d$ntriples)

  one <- ids[1]
  expect_true(all(d$triples$age[d$triples$id == one] == cov$age[cov$id == one]))
})

test_that("prep2 errors on a covariate column not present in data", {
  panel <- toy_panel(n = 200, seed = 4)
  expect_error(prep2(panel, states = c("A", "B", "C"), covariates = "nope"),
               "Covariate column")
})

test_that("prep2 warns when a covariate is not constant within subject", {
  panel <- toy_panel(n = 200, seed = 4)
  panel$noisy <- panel$time   # varies within subject by construction
  expect_warning(prep2(panel, states = c("A", "B", "C"), covariates = "noisy"),
                 "not constant")
})

test_that("print.msm2data reports covariates", {
  panel <- toy_panel(n = 200, seed = 4)
  panel$age <- 50
  d <- prep2(panel, states = c("A", "B", "C"), covariates = "age")
  expect_output(print(d), "covariates")
  expect_output(print(d), "age")

  d0 <- prep2(toy_panel(n = 200, seed = 4), states = c("A", "B", "C"))
  expect_output(print(d0), "none")
})

test_that("from_msdata(keep = ) replicates baseline covariates across rows", {
  md <- data.frame(
    id     = c(1, 1, 1, 2, 2),
    from   = c(1, 1, 2, 1, 1), to = c(2, 3, 3, 2, 3),
    Tstart = c(0, 0, 3, 0, 0), Tstop = c(3, 3, 5, 2, 2),
    status = c(1, 0, 1, 0, 1),
    sex    = c("M", "M", "M", "F", "F"))
  class(md) <- c("msdata", "data.frame")

  pan <- from_msdata(md, keep = "sex")
  expect_true("sex" %in% names(pan))
  expect_true(all(pan$sex[pan$id == 1] == "M"))
  expect_true(all(pan$sex[pan$id == 2] == "F"))
})

test_that("prep2 auto-dispatches msdata objects with covariates", {
  md <- data.frame(
    id     = c(1, 1, 1, 2, 2),
    from   = c(1, 1, 2, 1, 1), to = c(2, 3, 3, 2, 3),
    Tstart = c(0, 0, 3, 0, 0), Tstop = c(3, 3, 5, 2, 2),
    status = c(1, 0, 1, 0, 1),
    sex    = c("M", "M", "M", "F", "F"))
  class(md) <- c("msdata", "data.frame")

  d <- prep2(md, states = c("1", "2", "3"), covariates = "sex")
  expect_s3_class(d, "msm2data")
  expect_true(all(d$triples$sex[d$triples$id == 1] == "M"))
  expect_true(all(d$triples$sex[d$triples$id == 2] == "F"))
})
