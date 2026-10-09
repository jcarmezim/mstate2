## The same four patients of an illness-death model in the three layouts.
sojourn_data <- function() {
  data.frame(id = 1:4, t_healthy = c(3, 5, 2, 7), t_ill = c(4, 0, 8, 0), dead = c(1, 1, 0, 0))
}
wide_data <- function() {
  data.frame(id = 1:4, ill_time = c(3, 5, 2, 6), ill_status = c(1, 0, 1, 0),
             dead_time = c(7, 5, 9, 6), dead_status = c(1, 1, 0, 0), age = c(60, 72, 55, 49))
}
same_panel <- function(a, b) {
  isTRUE(all.equal(as.data.frame(a$panel[1:3]), as.data.frame(b$panel[1:3]), check.attributes = FALSE))
}

test_that("sojourn data become a daily panel that prep2() accepts", {
  x <- msprep2(sojourn_data(), durations = c(healthy = "t_healthy", ill = "t_ill"), outcome = c(dead = "dead"))
  expect_s3_class(x, "msm2prep")
  expect_equal(x$settings$format, "sojourn")
  expect_equal(x$states, c("healthy", "ill", "dead"))
  expect_equal(x$absorbing, "dead")
  p <- dplyr::filter(x$panel, id == 1)
  expect_equal(as.character(p$state), rep(c("healthy", "ill", "dead"), c(3, 4, 1)))
  expect_equal(p$time, 0:7)
  expect_equal(x$subjects$status, c("absorbed", "absorbed", "censored", "censored"))
  expect_equal(nrow(x$issues), 0)
  d <- prep2(x)
  expect_s3_class(d, "msm2data")
  expect_equal(d$states, x$states)
  expect_output(print(x), "discrete-time panel")
  expect_output(summary(x), "Observed transitions")
})

test_that("`unit` also applies to sojourn data", {
  s <- data.frame(id = 1:2, t_a = c(14, 7), t_b = c(7, 21), dead = c(1, 0))
  x <- msprep2(s, durations = c(a = "t_a", b = "t_b"), outcome = c(dead = "dead"), unit = 7)
  expect_equal(as.character(x$panel$state), c("a", "a", "b", "dead", "a", "b", "b", "b"))
})

test_that("wide data with Surv() give the same panel, keep the covariates and need no library(survival)", {
  p1 <- msprep2(sojourn_data(), durations = c(healthy = "t_healthy", ill = "t_ill"), outcome = c(dead = "dead"))
  x <- msprep2(wide_data(),
               states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
               trans = c("healthy -> ill -> dead", "healthy -> dead"))
  expect_true(same_panel(p1, x))
  expect_equal(names(x$panel), c("id", "time", "state", "age"))
  expect_equal(x$states, c("healthy", "ill", "dead"))
  expect_equal(x$absorbing, "dead")
  expect_true(all(x$transitions$allowed))
  ## without trans: the absorbing states are the states nobody leaves
  y <- msprep2(wide_data(), states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                          dead = Surv(dead_time, dead_status)))
  expect_identical(x$panel, y$panel)
  expect_null(y$trans)
  ## survival::Surv() and expressions; `keep` chooses the covariates
  z <- msprep2(wide_data(), keep = character(0),
               states = list(healthy = NULL, ill = survival::Surv(ill_time, ill_status == 1),
                             dead = Surv(dead_time, dead_status)))
  expect_true(same_panel(x, z))
  expect_equal(names(z$panel), c("id", "time", "state"))
  ## without an id column, the rows are numbered
  expect_true(same_panel(x, msprep2(wide_data()[-1], states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                                                    dead = Surv(dead_time, dead_status)))))
})

test_that("msdata objects of mstate give the same panel, with their covariates", {
  skip_if_not_installed("mstate")
  tmat <- mstate::transMat(x = list(c(2, 3), 3, c()), names = c("healthy", "ill", "dead"))
  ms <- mstate::msprep(time = c(NA, "ill_time", "dead_time"), status = c(NA, "ill_status", "dead_status"),
                       data = wide_data(), trans = tmat, keep = "age")
  x <- msprep2(ms)
  p1 <- msprep2(sojourn_data(), durations = c(healthy = "t_healthy", ill = "t_ill"), outcome = c(dead = "dead"))
  expect_equal(x$settings$format, "msdata")
  expect_true(same_panel(p1, x))
  expect_equal(names(x$panel), c("id", "time", "state", "age"))
  expect_equal(x$absorbing, "dead")
})

test_that("the transitions can be text or a transMat matrix", {
  a <- msprep2(wide_data(), states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
               trans = c("healthy -> ill", "ill -> dead", "healthy -> dead"))
  b <- msprep2(wide_data(), states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
               trans = c("healthy -> ill -> dead", "healthy -> dead"))   # a chain A -> B -> C
  expect_identical(a$trans, b$trans)
  skip_if_not_installed("mstate")
  tmat <- mstate::transMat(x = list(c(2, 3), 3, c()), names = c("healthy", "ill", "dead"))
  c2 <- msprep2(wide_data(), states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                trans = tmat)
  expect_identical(a$trans, c2$trans)
  expect_identical(a$panel, c2$panel)
})

test_that("the mstate::msprep() workflow of ebmt3 becomes one call", {
  skip_if_not_installed("mstate")
  ebmt3 <- get(utils::data("ebmt3", package = "mstate", envir = environment()))
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "PR", "RelDeath"))
  msbmt <- mstate::msprep(time = c(NA, "prtime", "rfstime"), status = c(NA, "prstat", "rfsstat"),
                          data = ebmt3, trans = tmat, keep = c("dissub", "age", "drmatch", "tcd", "prtime"))
  ref <- msprep2(msbmt)
  x <- msprep2(ebmt3, states = list(Tx = NULL, PR = Surv(prtime, prstat), RelDeath = Surv(rfstime, rfsstat)),
               trans = c("Tx -> PR -> RelDeath", "Tx -> RelDeath"))
  expect_true(same_panel(ref, x))
  expect_equal(names(x$panel), c("id", "time", "state", "dissub", "age", "drmatch", "tcd"))
  expect_equal(names(ref$panel), c("id", "time", "state", "dissub", "age", "drmatch", "tcd", "prtime"))
  ## the same counts as mstate::events()
  n <- stats::setNames(x$transitions$n, paste(x$transitions$from, x$transitions$to))
  expect_equal(unname(n[c("Tx PR", "Tx RelDeath", "PR RelDeath")]), c(1169L, 458L, 383L))
  ## months as time unit
  m <- suppressWarnings(msprep2(ebmt3, states = list(Tx = NULL, PR = Surv(prtime, prstat), RelDeath = Surv(rfstime, rfsstat)),
                                unit = 30.4375))
  expect_lt(max(m$panel$time), max(x$panel$time))
})

test_that("the sojourn layout gives the same panel as sojourn_to_panel()", {
  df <- data.frame(id = 1:3, t.a = c(2.5, 1, 0), t.b = c(0, 2, 3), disch.s = c(1, 0, 0), death.s = c(0, 1, 0))
  segs <- c(A = "t.a", B = "t.b")
  ref <- sojourn_to_panel(df, "id", segs, c(Disch = "disch.s", Death = "death.s"))
  x <- msprep2(df, durations = segs, outcome = c(Disch = "disch.s", Death = "death.s"))
  expect_equal(as.data.frame(x$panel[1:3]) |> transform(state = as.character(state)),
               as.data.frame(ref) |> transform(state = as.character(state)), ignore_attr = TRUE)
  expect_equal(msprep2(df, durations = segs, outcome = c(Disch = "disch.s", Death = "death.s"),
                       states = c("B", "A", "Disch", "Death"))$states, c("B", "A", "Disch", "Death"))
})

test_that("every record dropped or changed is logged, with one warning", {
  w <- wide_data()
  ## rounded_to_zero (sojourn)
  s <- sojourn_data()
  s$t_ill[2] <- 0.4
  expect_warning(x <- msprep2(s, durations = c(healthy = "t_healthy", ill = "t_ill"), outcome = c(dead = "dead")),
                 "rounded_to_zero")
  expect_equal(x$issues$id, 2L)
  ## missing: status 1 without a time, and an invalid status code
  bad <- w
  bad$ill_time[1] <- NA
  bad$dead_status[2] <- 9
  x <- suppressWarnings(msprep2(bad, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                                   dead = Surv(dead_time, dead_status))))
  expect_equal(x$issues$issue, c("missing", "missing"))
  expect_equal(x$issues$id, 1:2)
  ## same_unit: in weeks
  expect_warning(x <- msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                              unit = 7), "same_unit")
  expect_equal(x$issues$id, c(1L, 3L))
  ## after_absorbing and before_start
  bad <- w
  bad$ill_time[1] <- 9
  bad$ill_time[2] <- -1
  bad$ill_status[2] <- 1
  x <- suppressWarnings(msprep2(bad, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                                trans = c("healthy -> ill -> dead", "healthy -> dead")))
  expect_setequal(x$issues$issue, c("after_absorbing", "before_start"))
  ## not_allowed: kept and reported
  expect_warning(x <- msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                              trans = "healthy -> ill -> dead"), "not_allowed")
  expect_equal(x$issues$id, 2L)
  expect_false(all(x$transitions$allowed))
  expect_output(summary(x), "not_allowed")
})

test_that("Surv() checks the status as survival does", {
  w <- wide_data()
  expect_warning(x <- msprep2(transform(w, dead_status = dead_status + 1),
                              states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))),
                 "1/2")
  expect_true(same_panel(x, msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                                     dead = Surv(dead_time, dead_status)))))
  wt <- transform(w, dead_status = ifelse(dead_status == 1, "yes", "no"))
  expect_error(msprep2(wt, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))),
               "status == \"yes\"")
  lvl <- "yes"
  x <- msprep2(wt, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status == lvl)))
  expect_equal(nrow(x$panel), 31)
  expect_error(msprep2(transform(w, ill_time = as.character(ill_time)),
                       states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))),
               "numeric")
})

test_that("wrong arguments stop with clear messages", {
  w <- wide_data()
  expect_error(msprep2(w), "states = list")
  expect_error(msprep2(data.frame(id = numeric())), "no rows")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(t_ill, ill_status), dead = Surv(dead_time, dead_status))), "t_ill")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = NULL, dead = Surv(dead_time, dead_status))), "exactly one state as NULL")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = "ill_time", dead = Surv(dead_time, dead_status))), "Surv\\(time, status\\)")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time), dead = Surv(dead_time, dead_status))), "needs a status")
  expect_error(msprep2(w, states = list(NULL, Surv(ill_time, ill_status))), "Name every state")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                       trans = "healty -> ill"), "healty")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                       trans = "healthy ill"), "from -> to")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                       unit = "day"), "positive number")
  expect_error(msprep2(w, states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status)),
                       keep = "nope"), "nope")
  expect_error(msprep2(cbind(w, time = 1), states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                                      dead = Surv(dead_time, dead_status))), "clash")
  expect_error(msprep2(sojourn_data(), durations = c(healthy = "t_healthy")), "outcome")
  expect_error(msprep2(sojourn_data(), durations = c(healthy = "t_nope"), outcome = c(dead = "dead")), "t_nope")
  expect_error(msprep2(rbind(w, w), states = list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))),
               "duplicates")
})
