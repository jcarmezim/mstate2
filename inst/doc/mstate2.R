## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## ----setup--------------------------------------------------------------------
library(mstate2)

## -----------------------------------------------------------------------------
st   <- c("A", "B", "C")
tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4   # from (A,B): stay 0.6
tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7   # from (B,B): stay 0.3
tens["C", "C", ]    <- 1                                  # C absorbing
first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1

set.seed(1)
panel <- simulate2(5000, tens, first, init = c(A = 1, B = 0, C = 0))
head(panel)

## -----------------------------------------------------------------------------
d   <- prep2(panel, id = "id", time = "time", state = "state")
fit <- P2est(d, estimator = "RPE")
fit$estimate

## -----------------------------------------------------------------------------
ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6)

## -----------------------------------------------------------------------------
cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 9)
summary(cmp)

## ----fig.alt = "Evolution intervals for the two preceding states"-------------
plot(cmp, type = "interval")

