## Tests for the 0.2.0 extensions: duration / pairs in prep2(), the LRT version
## of markov_test(), P1est(), P2boot(), probtrans2(), compare_order() with a
## discrete baseline, divergence_table(), the flexible P2reg() baseline,
## P2reg_all(), plot methods and msprep2().

## A first-order chain written as a second-order tensor (identical slices):
## used for the internal-consistency checks of the research plan (3.5.2).
first_order_tensor <- function() {
  st <- c("A", "B", "C")
  P1 <- matrix(c(.70, .25, .05,
                 .20, .60, .20,
                 0,   0,   1), 3, byrow = TRUE, dimnames = list(st, st))
  list(T = array(P1, c(3, 3, 3), dimnames = list(st, st, st)), P1 = P1)
}

test_that("prep2 stores the time spent in the current state and all one-step moves", {
  p <- data.frame(id = c(1, 1, 1, 1, 1, 2, 2, 2), time = c(0:4, 0:2),
                  state = c("A", "A", "B", "B", "C", "B", "B", "C"))
  d <- prep2(p, states = c("A", "B", "C"))
  expect_true(all(c("d") %in% names(d$triples)))
  expect_equal(d$triples$d, c(2, 1, 2, 2))            # days already in j at s - 1
  ## the pairs table includes the first move of each subject
  expect_equal(d$pairs[from == "A" & to == "A", N], 1L)
  expect_equal(sum(d$pairs$N), nrow(p) - 2L)
  expect_error(prep2(transform(p, d = 1), covariates = "d"), "reserved")
})

test_that("markov_test(method = 'lrt') matches the likelihood-ratio statistic by hand", {
  fit <- toy_fit(n = 1500, seed = 3)
  mt  <- markov_test(fit, method = "lrt")
  e   <- fit$estimate
  N <- e$n.trans[e$l == "C"]; Y <- e$at.risk[e$l == "C"]
  pbar <- sum(N) / sum(Y)
  G2 <- 2 * sum(N * log(N / (Y * pbar)) + (Y - N) * log((Y - N) / (Y * (1 - pbar))))
  expect_equal(mt$statistic[mt$l == "C"], G2, tolerance = 1e-10)
  expect_equal(mt$pooled[mt$l == "C"], pbar)
  ## with two destinations the joint (per-state) test is the same statistic
  ms <- markov_test(fit, method = "lrt", by = "state")
  expect_equal(ms$statistic, G2, tolerance = 1e-10)
  expect_equal(ms$df, 1L)
  expect_s3_class(ms, "markov_test")
  expect_error(markov_test(fit, by = "state"), "lrt")
  expect_output(print(ms), "per current state")
})

test_that("the LRT pooled estimate is the first-order P1est estimate", {
  d  <- prep2(toy_panel(1500, seed = 4), states = c("A", "B", "C"))
  mt <- markov_test(P2est(d), method = "lrt")
  p1 <- P1est(d, B = 0)
  expect_equal(mt$pooled[mt$l == "C"], p1$P["B", "C"])
  expect_null(p1$boot)
  expect_equal(unname(rowSums(p1$P)), c(0, 1, 1))     # A only at time 0: no row
  expect_output(print(p1), "P1est")
})

test_that("with first-order data the second-order estimates reproduce the first-order ones", {
  ft <- first_order_tensor()
  set.seed(21)
  pan <- simulate2(3000, ft$T, ft$P1, init = c(A = 1, B = 0, C = 0))
  d   <- prep2(pan, states = c("A", "B", "C"))
  fit <- P2est(d)
  p1  <- P1est(d, B = 0)
  for (jj in c("A", "B")) for (hh in c("A", "B")) {
    expect_equal(fit$P[jj, , hh], p1$P[jj, ], tolerance = 0.05)
    expect_equal(p1$P[jj, ], ft$P1[jj, ], tolerance = 0.05)
  }
  ## and the Markov tests do not reject the (true) first-order hypothesis
  expect_true(all(markov_test(fit, method = "lrt", by = "state")$p.value > 0.001))
  expect_true(all(markov_test(fit)$p.value > 0.001))
})

test_that("P2boot gives reproducible percentile bands inside [0, 1]", {
  d  <- prep2(toy_panel(800, seed = 5), states = c("A", "B", "C"))
  b1 <- P2boot(d, B = 60, seed = 9)
  b2 <- P2boot(d, B = 60, seed = 9)
  expect_s3_class(b1, "P2boot"); expect_s3_class(b1, "P2est")
  expect_identical(b1$boot, b2$boot)
  expect_equal(dim(b1$boot), c(3, 3, 3, 60))
  expect_equal(b1$P, P2est(d)$P)
  ## each replicate is a proper set of transition probabilities
  rs <- apply(b1$boot, c(1, 3, 4), sum)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
  ck <- ckequations(b1, "A", "B", "C", nsteps = 5, bounds = TRUE)
  expect_true(all(ck$lower >= 0 & ck$upper <= 1 & ck$lower <= ck$estimate + 1e-12))
  cmp <- compare2(b1, h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  expect_identical(attr(cmp, "bands"), "bootstrap")
  expect_true("se.boot" %in% names(b1$estimate))
  expect_output(print(b1), "bootstrap")
  expect_error(P2boot(d, B = 1), "B")
})

test_that("probtrans2 has mstate's probtrans structure and plots", {
  d   <- prep2(toy_panel(800, seed = 6), states = c("A", "B", "C"))
  fit <- P2est(d)
  pt  <- probtrans2(fit, h = "A", j = "B", nsteps = 5)
  expect_s3_class(pt, "probtrans")
  expect_identical(names(pt[[2]]), c("time", paste0("pstate", 1:3), paste0("se", 1:3)))
  expect_equal(pt$predt, 1); expect_equal(pt$direction, "forward")
  expect_equal(unname(rowSums(pt[[2]][, 2:4])), rep(1, 6))
  expect_equal(pt[[2]]$pstate3[-1], ckequations(fit, "A", "B", "C", nsteps = 5))
  expect_null(pt[[1]])
  ## first-order predictions from every transient state, with bootstrap SEs
  p1 <- probtrans2(P1est(d, B = 40, seed = 1), nsteps = 4)
  expect_false(is.null(p1[[2]]))
  expect_true(all(p1[[2]]$se3[-1] > 0))
  pdf(NULL); on.exit(dev.off())
  for (ty in c("filled", "stacked", "single", "separate")) expect_silent(plot(pt, type = ty))
  expect_error(plot(pt, from = 1), "No predictions")
  expect_error(probtrans2(fit, j = "B"), "required")
})

test_that("compare_order accepts a P1est baseline; divergence_table covers all transitions", {
  d   <- prep2(toy_panel(800, seed = 7), states = c("A", "B", "C"))
  fit <- P2est(d)
  p1  <- P1est(d, B = 60, seed = 2)
  co  <- compare_order(fit, p1, h = c("A", "B"), j = "B", l = "C", nsteps = 5)
  expect_true("1st order (discrete)" %in% co$h)
  expect_equal(attr(co, "baseline"), "1st order (discrete)")
  expect_equal(nrow(divergence(co)), 2L)
  expect_error(compare_order(fit, P1est(d, B = 0), h = "A", j = "B", l = "C"), "bootstrap")
  dt <- divergence_table(fit, p1, nsteps = 5)
  expect_true(all(c("j", "l", "h", "separated_steps", "max_abs_diff") %in% names(dt)))
  expect_setequal(unique(dt$l), c("B", "C"))           # A is never reached from B
  expect_true(all(diff(dt$separated_steps) <= 0))
  expect_error(divergence_table(fit, p1, min.h = 3), "min.h")
})

test_that("P2reg accepts the time scales s and d, and P2reg_all reproduces P2est", {
  pan <- toy_panel(1500, seed = 8)
  ids <- unique(pan$id); set.seed(1)
  pan$age <- stats::setNames(sample(40:80, length(ids), TRUE), ids)[as.character(pan$id)]
  d <- prep2(pan, states = c("A", "B", "C"), covariates = "age")
  r <- P2reg(d, h = "B", j = "B", formula = ~ age + log(d) + s)
  expect_true(all(c("log(d)", "s") %in% summary(r)$term))
  expect_equal(dim(predict(r, newdata = data.frame(age = 60, d = 2, s = 4))), c(1L, 1L))
  ## splines in the formula are found through the formula's environment
  expect_s3_class(P2reg(d, h = "A", j = "B", formula = ~ splines::ns(age, 2)), "P2reg")

  all1 <- P2reg_all(d, formula = ~1)
  expect_equal(predict(all1, newdata = data.frame(age = 50)), P2est(d)$P)
  all2 <- P2reg_all(d, formula = ~ age)
  tl <- predict(all2, newdata = data.frame(age = c(40, 80)))
  expect_length(tl, 2L)
  for (P in tl) {
    rs <- apply(P, c(1, 3), sum)
    expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
  }
  expect_equal(nrow(summary(all2)), 4L)                # 2 histories x (intercept, age)
  expect_output(print(all2), "P2reg_all")
  pdf(NULL); on.exit(dev.off())
  expect_invisible(plot(r))
  expect_error(plot(P2reg(d, h = "A", j = "B")), "No covariate")
})

test_that("msprep2 builds the panel from msprep-style wide data", {
  wide <- data.frame(id = 1:3, ill.t = c(2, 5, 4), ill.s = c(1, 0, 1),
                     dth.t = c(6, 5, 4.5), dth.s = c(1, 0, 1), age = c(60, 72, 55))
  p <- msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
               data = wide, states = c("H", "I", "D"), absorbing = "D", id = "id", keep = "age")
  expect_equal(p[id == 1, state], c("H", "H", "I", "I", "I", "I", "D"))
  expect_equal(p[id == 2, state], rep("H", 6))          # censored at 5
  expect_equal(p[id == 3, state], c("H", "H", "H", "H", "I", "D"))   # rnd(4.5) = 5
  expect_equal(unique(p[id == 2, age]), 72)
  expect_warning(msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
                         data = wide), "absorbing")
  expect_error(msprep2(time = c(NA, "x"), status = c(NA, "ill.s"), data = wide,
                       absorbing = "2"), "not found")
})

test_that("msprep2 warns when a state is entered on the same rounded day as the next", {
  wide <- data.frame(id = 1, ill.t = 3.2, ill.s = 1, dth.t = 3.4, dth.s = 1)
  expect_warning(p <- msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
                              data = wide, states = c("H", "I", "D"), absorbing = "D"),
                 "1 state visit\\(s\\)")
  expect_equal(p$state, c("H", "H", "H", "D"))          # I (3.2 to 3.4) is lost
})

test_that("msprep2 agrees with mstate::msprep followed by from_msdata", {
  skip_if_not_installed("mstate")
  ebmt3 <- NULL
  utils::data("ebmt3", package = "mstate", envir = environment())
  ebmt3 <- ebmt3[1:300, ]
  tmat <- mstate::trans.illdeath()
  ms <- mstate::msprep(time = c(NA, "prtime", "rfstime"), status = c(NA, "prstat", "rfsstat"),
                       data = ebmt3, trans = tmat)
  a <- from_msdata(ms)
  a$state <- rownames(tmat)[as.integer(a$state)]
  b <- msprep2(time = c(NA, "prtime", "rfstime"), status = c(NA, "prstat", "rfsstat"),
               data = ebmt3, trans = tmat, id = "id")
  expect_equal(as.data.frame(a), as.data.frame(b), ignore_attr = TRUE)
})

test_that("P1est reuses P2boot replicates; probtrans2 and msprep2 edge cases", {
  d  <- prep2(toy_panel(600, seed = 12), states = c("A", "B", "C"))
  bt <- P2boot(d, B = 30, seed = 3)
  p1 <- P1est(bt)
  expect_identical(p1$boot, bt$boot1)
  expect_equal(p1$P, P1est(d, B = 0)$P)
  expect_error(P1est(P2est(d)), "msm2data")
  pt <- probtrans2(p1, j = "B", nsteps = 3)
  expect_null(pt[[1]]); expect_false(is.null(pt[[2]]))
  expect_true(is.na(pt$h))
  expect_true(all(is.na(probtrans2(P2est(d), "A", "B", nsteps = 2)[[2]]$se3[-1])))
  expect_error(probtrans2(p1, j = "Z"), "not found")
  expect_error(probtrans2(list(states = "A"), j = 1), "must be")
  expect_error(probtrans2(p1, nsteps = 0), "nsteps")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(pt, from = "B", type = "stacked", ord = 3:1))

  wide <- data.frame(ill.t = c(2, 5), ill.s = c(1, 0), dth.t = c(6, 5), dth.s = c(1, 0))
  tm <- matrix(c(NA, NA, NA, 1, NA, NA, 2, 3, NA), 3,
               dimnames = list(c("H", "I", "D"), c("H", "I", "D")))
  ## trans gives names and absorbing states; ids default to 1..n; custom start
  p <- msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
               data = wide, trans = tm, start = list(state = "H", time = 1))
  expect_equal(p[id == 1, time][1], 1)
  expect_equal(p[id == 1, state], c("H", "I", "I", "I", "I", "D"))
  ## an observed move that trans forbids triggers a warning
  tm2 <- tm; tm2["H", "I"] <- NA
  expect_warning(msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
                         data = wide, trans = tm2), "not allowed")
  expect_error(msprep2(time = c(NA, "ill.t"), status = c(NA, "ill.s", "dth.s"), data = wide),
               "same length")
  expect_error(msprep2(time = c(NA, "ill.t", "dth.t"), status = c(NA, "ill.s", "dth.s"),
                       data = wide, states = c("a", "b"), absorbing = "b"), "one label")
})

test_that("P2reg_all reports terms that are constant within a history's risk set", {
  pan <- toy_panel(800, seed = 13)
  ids <- unique(pan$id)
  ## a covariate that is 1 only for subjects absorbed straight from (A, B):
  ## it is constant (0) in the (B, B) risk set, so its effect is not estimable there
  set.seed(4)
  first_exit <- pan[, .(ex = state[3] == "C"), by = id]
  first_exit[, z := as.integer(ex & stats::runif(.N) < 0.5)]   # 1 only for some early exits
  pan <- merge(pan, first_exit[, .(id, z)], by = "id")
  d <- prep2(pan, states = c("A", "B", "C"), covariates = "z")
  a <- P2reg_all(d, formula = ~ z)
  expect_true(any(nzchar(a$pairs$aliased)))
  expect_output(print(a), "not estimable")
  expect_warning(P <- predict(a, newdata = data.frame(z = 1)), "not estimable")
  rs <- apply(P, c(1, 3), sum)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-12))
})

test_that("markov_test adds Holm-adjusted p-values", {
  fit <- P2est(prep2(toy_panel(2000, seed = 4), states = c("A", "B", "C")))
  mt  <- markov_test(fit, method = "lrt")
  expect_true("p.adj" %in% names(mt))
  expect_equal(mt$p.adj, stats::p.adjust(mt$p.value, "holm"))
  expect_equal(markov_test(fit, method = "lrt", p.adjust = "none")$p.adj,
               markov_test(fit, method = "lrt")$p.value)
})

test_that("paired bootstrap difference: compare2/compare_order with shared replicates", {
  d  <- prep2(toy_panel(2000, seed = 5), states = c("A", "B", "C"))
  bt <- P2boot(d, B = 100, seed = 1)
  os <- overlap_step(compare2(bt, h = c("A", "B"), j = "B", l = "C", nsteps = 5))
  expect_false(is.na(os$diff_steps))
  expect_equal(nrow(os$diff), 5L)
  ## the true curves differ a lot here: the paired test detects it at least
  ## as long as the overlap criterion (which is conservative)
  expect_gte(os$diff_steps, os$separated_steps)
  ## evolution intervals: no replicates, no paired test
  expect_true(is.na(overlap_step(compare2(P2est(d), c("A", "B"), "B", "C"))$diff_steps))
  ## P1est built from the same P2boot -> paired; from its own bootstrap -> not
  dv1 <- divergence(compare_order(bt, P1est(bt), h = c("A", "B"), j = "B", l = "C", nsteps = 5))
  expect_false(anyNA(dv1$diff_steps))
  dv2 <- divergence(compare_order(bt, P1est(d, B = 50, seed = 2), h = c("A", "B"),
                                  j = "B", l = "C", nsteps = 5))
  expect_true(all(is.na(dv2$diff_steps)))
  expect_true("diff_steps" %in% names(divergence_table(bt, P1est(bt), nsteps = 5)))
})
