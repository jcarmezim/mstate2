# Round half away from zero

Function to turn the sojourn times into whole days. Unlike
[`round`](https://rdrr.io/r/base/Round.html), which rounds halves to the
even number (`round(0.5) = 0`, `round(2.5) = 2`), `rnd` always rounds
halves away from zero (`rnd(0.5) = 1`, `rnd(2.5) = 3`).

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
```
