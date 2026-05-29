test_that("RPE/CPE recover known 1-step second-order probabilities", {
  set.seed(1)
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4   # p_ABB, p_ABC
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7   # p_BBB, p_BBC
  tens["C", "C", ] <- 1                                     # C absorbing
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1

  panel <- simulate2(4000, tens, first, init = c(A = 1, B = 0, C = 0))
  d <- prep2(panel, states = st)
  e <- P2est(d, "RPE")$estimate
  get <- function(h, j, l) e$p[e$h == h & e$j == j & e$l == l]

  expect_equal(get("A", "B", "B"), 0.6, tolerance = 0.05)
  expect_equal(get("B", "B", "B"), 0.3, tolerance = 0.05)
  expect_equal(get("A", "B", "B") + get("A", "B", "C"), 1, tolerance = 1e-8)
})

test_that("package SEs match the hand-rolled relative/conditional formulas", {
  set.seed(2)
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1
  panel <- simulate2(3000, tens, first, init = c(A = 1, B = 0, C = 0))

  d   <- prep2(panel, states = st)
  rpe <- P2est(d, "RPE"); cpe <- P2est(d, "CPE")

  ## hand-rolled (Illustration.R style) for B -> B -> B
  P <- data.table::as.data.table(panel)
  data.table::setorder(P, id, time)
  P[, c("hh", "jj") := list(data.table::shift(state, 2), data.table::shift(state, 1)), by = id]
  tr <- P[!is.na(hh) & !is.na(jj)]
  M <- max(tr$time); Y <- Nn <- numeric(0)
  for (s in 2:M) {
    Y[s - 1]  <- nrow(tr[time == s & hh == "B" & jj == "B"])
    Nn[s - 1] <- nrow(tr[time == s & hh == "B" & jj == "B" & state == "B"])
  }
  k <- Y > 0; Y <- Y[k]; Nn <- Nn[k]
  rpe_se <- sqrt((sum(Nn)/sum(Y)) * (1 - sum(Nn)/sum(Y)) / sum(Y))
  cpe_p  <- mean(Nn / Y)
  cpe_se <- sqrt(cpe_p * (1 - cpe_p) * sum(1 / Y) / length(Y)^2)

  er <- rpe$estimate; ec <- cpe$estimate
  expect_equal(er$se[er$h=="B" & er$j=="B" & er$l=="B"], rpe_se, tolerance = 1e-10)
  expect_equal(ec$se[ec$h=="B" & ec$j=="B" & ec$l=="B"], cpe_se, tolerance = 1e-10)
})

test_that("RPE is no less efficient than CPE (Corollary 3) under staggered entry", {
  set.seed(3)
  st <- c("A", "B", "C")
  tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
  tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
  tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
  tens["C", "C", ] <- 1
  first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1

  R <- 40; rpe_est <- cpe_est <- numeric(R)
  for (r in seq_len(R)) {
    panel <- simulate2(800, tens, first, entry = c(A = 0.1))
    d <- prep2(panel, states = st)
    er <- P2est(d, "RPE")$estimate; ec <- P2est(d, "CPE")$estimate
    rpe_est[r] <- er$p[er$h=="B" & er$j=="B" & er$l=="B"]
    cpe_est[r] <- ec$p[ec$h=="B" & ec$j=="B" & ec$l=="B"]
  }
  expect_lt(sd(rpe_est), sd(cpe_est))            # RPE less variable
})
