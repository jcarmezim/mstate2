# First step at which two evolution intervals overlap

Scans the two curves of a
[`compare2`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
comparison forward from \\n = 1\\ and returns the first step at which
their intervals overlap: before it, the state at the previous time
changes the prediction.

## Usage

``` r
overlap_step(x)
```

## Arguments

- x:

  A two-group "msm2pred" with interval bounds.

## Value

A list: first overlapping step `n` and time `s = n + 1`; the leading
separated-step count `separated_steps`; per-step `overlap` and signed
`separation` (= max lower - min upper, \> 0 when separated); and the two
`groups`.

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
overlap_step(cmp)
#> $n
#> [1] 3
#> 
#> $s
#> [1] 4
#> 
#> $separated_steps
#> [1] 2
#> 
#> $overlap
#>     1     2     3     4     5     6 
#> FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
#> 
#> $separation
#>            1            2            3            4            5            6 
#>  0.224601660  0.030557457 -0.002568936 -0.004577489 -0.002628205 -0.001208199 
#> 
#> $groups
#> [1] "A" "B"
#> 
```
