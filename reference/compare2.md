# Compare n-step transitions across preceding states (evolution intervals)

For a fixed current state \\j\\ and target \\\ell\\, computes the n-step
transition probabilities for several preceding states \\h\\, with their
evolution intervals. The construction behind Figures 4-5: it assesses
whether, and for how long, the preceding state affects the future
trajectory.

## Usage

``` r
compare2(object, h, j, l, nsteps = 9L, bounds = TRUE)
```

## Arguments

- object:

  A "P2est" object.

- h:

  Vector of preceding states to compare (labels or indices).

- j:

  Current state shared across the compared paths.

- l:

  Target state.

- nsteps:

  Number of steps (default 9).

- bounds:

  If TRUE (default) also compute evolution-interval bounds; FALSE
  returns curves only (e.g. Figure-4 overlays) and skips two thirds of
  the work. If `object` is a
  [`P2boot`](https://jcarmezim.github.io/mstate2/reference/P2boot.md)
  fit, the bounds are percentile bootstrap intervals instead of
  evolution intervals.

## Value

An object of class "msm2pred" (a data frame): columns h, n, estimate
and, when bounds = TRUE, lower and upper.

## Examples

``` r
st   <- c("A", "B", "C")                                 # C is absorbing
tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
tens["C", "C", ]    <- 1
first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1

set.seed(1)
panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
fit <- P2est(prep2(panel))
cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
summary(cmp)
#> Trajectory comparison (RPE, 95% evolution intervals)
#>   target 'B' via current state 'B'; preceding: A, B
#>   intervals first overlap at step 3 (time s = 4); significant for the first 2 step(s).
```
