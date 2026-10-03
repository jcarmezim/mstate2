# Plot n-step evolution curves or evolution intervals

Plot n-step evolution curves or evolution intervals

## Usage

``` r
# S3 method for class 'msm2pred'
plot(
  x,
  type = NULL,
  col = NULL,
  lty = 1,
  lwd = 2,
  alpha = 0.2,
  add = FALSE,
  legend = TRUE,
  dualaxis = TRUE,
  mark.overlap = TRUE,
  ylim = NULL,
  main = NULL,
  xlab = "",
  ylab = "probability",
  ...
)
```

## Arguments

- x:

  A "msm2pred" object.

- type:

  "interval" (shaded bands) or "curve" (lines only). Defaults to
  "interval" when bounds exist, else "curve".

- col, lty, lwd, alpha:

  Appearance of lines and interval shading. The colours are used in the
  order of `h` in
  [`compare2`](https://jcarmezim.github.io/mstate2/reference/compare2.md)

- add:

  Overlay on the current plot. Default FALSE.

- legend, dualaxis:

  Draw the legend / the dual step-and-time axis.

- mark.overlap:

  For a two-group interval plot, draw a vertical line at the first
  overlap step. Default TRUE.

- ylim, main, xlab, ylab:

  Standard graphical parameters.

- ...:

  Passed to the initial plot().

## Examples

``` r
st <- c("A", "B", "C") # C is absorbing
tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
tens["B", "B", "A"] <- 0.6
tens["B", "C", "A"] <- 0.4
tens["B", "B", "B"] <- 0.3
tens["B", "C", "B"] <- 0.7
tens["C", "C", ] <- 1
first <- matrix(0, 3, 3, dimnames = list(st, st))
first["A", "B"] <- 1

set.seed(1)
panel <- simulate2(500, tens, first, init = c(A = 1, B = 0, C = 0))
fit <- P2est(prep2(panel))
cmp <- compare2(fit, h = c("A", "B"), j = "B", l = "B", nsteps = 6)
plot(cmp, type = "interval")

plot(cmp, type = "curve")
```
