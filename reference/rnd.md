# Round half up (the convention used in the DIVINE analysis)

Round half up (the convention used in the DIVINE analysis)

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
rnd(c(-2.5, -0.4, 0.4, 2.5))
#> [1] -3  0  0  3
```
