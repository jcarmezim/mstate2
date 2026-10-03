## Raw records with dates, as collected: one row per state change.
raw_events <- function() {
  data.frame(id    = c(1, 1, 1, 2, 2, 3, 3, 3),
             date  = c("2020-03-01", "2020-03-03", "2020-03-06",
                       "2020-03-02", "2020-03-09",
                       "2020-03-05", "2020-03-06", "2020-03-12"),
             state = c("NSP", "SP", "Disch", "SP", "Death", "NSP", "SP", "SP"))
}

test_that("events with text dates become a daily panel that prep2() accepts", {
  x <- msprep2(raw_events(), time = "date", state = "state", absorbing = c("Disch", "Death"))
  expect_s3_class(x, "msm2prep")
  p <- dplyr::filter(x$panel, id == 1)
  expect_equal(as.character(p$state), c("NSP", "NSP", "SP", "SP", "SP", "Disch"))
  expect_equal(p$time, 0:5)
  ## a repeated record of the same state extends follow-up (subject 3, censored)
  expect_equal(max(dplyr::filter(x$panel, id == 3)$time), 7)
  expect_equal(x$subjects$status, c("absorbed", "absorbed", "censored"))
  expect_equal(x$states, c("NSP", "SP", "Disch", "Death"))
  d <- prep2(x)
  expect_s3_class(d, "msm2data")
  expect_equal(d$states, x$states)
  expect_equal(d$absorbing, c("Disch", "Death"))
  expect_equal(nrow(x$issues), 0)
  expect_output(print(x), "msm2prep")
  expect_output(summary(x), "Observed transitions")
})

test_that("dd/mm/yyyy text dates, POSIXct times and named units are handled", {
  r <- raw_events()
  r$date <- format(as.Date(r$date), "%d/%m/%Y")
  a <- msprep2(r, time = "date", absorbing = c("Disch", "Death"))
  b <- msprep2(raw_events(), time = "date", absorbing = c("Disch", "Death"))
  expect_equal(a$panel, b$panel)
  r2 <- raw_events()
  r2$date <- as.POSIXct(paste(r2$date, "10:00"), tz = "UTC")
  expect_equal(msprep2(r2, time = "date", absorbing = c("Disch", "Death"))$panel, b$panel)
  w <- suppressWarnings(msprep2(raw_events(), time = "date", unit = "week", absorbing = c("Disch", "Death")))
  expect_lte(max(w$panel$time), 2)
  expect_error(msprep2(data.frame(id = 1, time = 1, state = "A"), unit = "week"), "dates")
  expect_error(msprep2(data.frame(id = 1, time = "yesterday", state = "A")), "dates")
})

test_that("every record dropped or changed is logged, with one warning", {
  r <- data.frame(
    id    = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3),
    time  = c(0, 1.5, 1.58, 4, 6, -5, 0, NA, 3, 0, 17),
    state = c("1", "2", "3", "D", "2", "1", "1", "2", "3", "2", "3"),
    fin   = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 7, 7))
  tm <- matrix(FALSE, 4, 4, dimnames = list(c("NSP", "SP", "IMV", "Death"), c("NSP", "SP", "IMV", "Death")))
  tm["NSP", "SP"] <- tm["SP", "IMV"] <- tm["SP", "Death"] <- tm["IMV", "Death"] <- TRUE
  expect_warning(
    x <- msprep2(r, start = 0, end = "fin", trans = tm,
                 recode = c("1" = "NSP", "2" = "SP", "3" = "IMV", "D" = "Death")),
    "7 record")
  expect_setequal(x$issues$issue, c("missing", "before_start", "after_end", "same_unit",
                                    "after_absorbing", "not_allowed"))
  expect_equal(as.character(dplyr::filter(x$panel, id == 1)$state), c("NSP", "NSP", "IMV", "IMV", "Death"))
  expect_equal(max(dplyr::filter(x$panel, id == 3)$time), 7)          # censored at `end`
  expect_equal(x$absorbing, "Death")                                   # from `trans`
  expect_error(suppressWarnings(msprep2(r, start = 0, trans = tm, check = "error",
                                        recode = c("1" = "NSP", "2" = "SP", "3" = "IMV", "D" = "Death"))),
               "not allowed")
  expect_error(msprep2(r, states = c("NSP", "SP")), "recode")
})

test_that("an absorbing state wins over a later record in the same unit", {
  r <- data.frame(id = 1, time = c(0, 4, 6), state = c("A", "Death", "B"))
  x <- suppressWarnings(msprep2(r, unit = 7, absorbing = "Death"))
  expect_equal(as.character(x$panel$state), c("A", "Death"))
  expect_equal(x$issues$issue, "after_absorbing")
})

test_that("ties = 'first' keeps the first state of the unit", {
  r <- data.frame(id = 1, time = c(0, 2.2, 2.4, 5), state = c("A", "B", "C", "D"))
  last <- suppressWarnings(msprep2(r, absorbing = "D"))
  first <- suppressWarnings(msprep2(r, absorbing = "D", ties = "first"))
  expect_equal(as.character(last$panel$state[3]), "C")
  expect_equal(as.character(first$panel$state[3]), "B")
})

test_that("the initial state and baseline covariates can be given", {
  r <- data.frame(id = c(1, 2), adm = c("2020-03-01", "2020-03-02"), ini = c("1", "2"),
                  t = c("2020-03-04", "2020-03-03"), s = c("D", "D"), age = c(70, 55))
  x <- msprep2(r, time = "t", state = "s", start = "adm", initial = "ini", keep = "age",
               recode = c("1" = "NSP", "2" = "SP", "D" = "Death"), absorbing = "Death")
  expect_equal(as.character(x$panel$state), c("NSP", "NSP", "NSP", "Death", "SP", "Death"))
  expect_equal(unique(x$panel$age), c(70, 55))
})

test_that("the sojourn layout gives the same panel as sojourn_to_panel()", {
  df <- data.frame(id = 1:3, t.a = c(3, 0.5, 2), t.b = c(1, 2, 0.3),
                   disch = c(1, 0, 0), death = c(0, 1, 0))
  ref <- suppressWarnings(sojourn_to_panel(df, "id", c(A = "t.a", B = "t.b"), c(Disch = "disch", Death = "death")))
  x <- suppressWarnings(msprep2(df, durations = c(A = "t.a", B = "t.b"), outcome = c(Disch = "disch", Death = "death")))
  expect_equal(as.data.frame(dplyr::mutate(x$panel, state = as.character(state))), as.data.frame(ref),
               ignore_attr = TRUE)
  expect_equal(x$issues$issue, "rounded_to_zero")
  expect_equal(x$subjects$status, c("absorbed", "absorbed", "censored"))
})

test_that("the wide layout and mstate msdata objects give the same panel as mstate", {
  skip_if_not_installed("mstate")
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "PR", "RelDeath"))
  ebmt3 <- get(utils::data("ebmt3", package = "mstate", envir = environment()))
  ms <- mstate::msprep(data = ebmt3, trans = tmat, time = c(NA, "prtime", "rfstime"),
                       status = c(NA, "prstat", "rfsstat"))
  m1 <- msprep2(ms)
  m2 <- msprep2(ebmt3, times = c(PR = "prtime", RelDeath = "rfstime"),
                status = c(PR = "prstat", RelDeath = "rfsstat"), initial = "Tx", trans = tmat)
  expect_equal(as.data.frame(m1$panel), as.data.frame(m2$panel), ignore_attr = TRUE)
  ## daily units: the same transition counts as mstate::events()
  fr <- mstate::events(ms)$Frequencies
  n <- function(a, b) m1$transitions$n[m1$transitions$from == a & m1$transitions$to == b]
  expect_equal(c(n("Tx", "PR"), n("Tx", "RelDeath"), n("PR", "RelDeath")),
               unname(c(fr["Tx", "PR"], fr["Tx", "RelDeath"], fr["PR", "RelDeath"])))
  expect_equal(m1$absorbing, "RelDeath")
})
