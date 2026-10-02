# Round half away from zero (the rounding of the paper's DIVINE analysis)

The function used in the code of the methods paper (Najera-Zuloaga,
Besalu and Gomez Melis, 2025) to turn the sojourn times of the DIVINE
cohort into whole days. Unlike
[`round`](https://rdrr.io/r/base/Round.html), which rounds halves to the
even number (`round(0.5) = 0`, `round(2.5) = 2`), `rnd` always rounds
halves away from zero (`rnd(0.5) = 1`, `rnd(2.5) = 3`). Only this
convention reproduces Table 2 of the paper.

## Usage

``` r
rnd(x)
```

## Arguments

- x:

  numeric vector.

## Value

`trunc(x + sign(x) * 0.5)`.

## Examples

``` r
rnd(c(-2.5, -0.4, 0.4, 0.5, 2.5))
#> [1] -3  0  0  1  3
round(c(-2.5, -0.4, 0.4, 0.5, 2.5))   # base R: halves to even
#> [1] -2  0  0  0  2
```
