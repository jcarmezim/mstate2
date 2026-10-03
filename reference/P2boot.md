# Subject-level bootstrap of the second-order estimates

Resamples subjects with replacement and recomputes the RPE of every
1-step second-order transition probability in each replicate. The
resulting object behaves like a
[`P2est`](https://jcarmezim.github.io/mstate2/reference/P2est.md) fit,
but carries the bootstrap tensors, so that
[`ckequations`](https://jcarmezim.github.io/mstate2/reference/ckequations.md)
and
[`compare2`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
report **percentile bootstrap intervals** for the \\n\\-step
probabilities instead of the evolution intervals obtained by propagating
the one-step confidence limits.

## Usage

``` r
P2boot(object, B = 200, conf.level = 0.95, seed = NULL)
```

## Arguments

- object:

  An "msm2data" object from
  [`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md).

- B:

  Number of bootstrap replicates. Default 200.

- conf.level:

  Confidence level of the intervals. Default 0.95.

- seed:

  Optional integer seed, for reproducible replicates.

## Value

An object of class `c("P2boot", "P2est")`: all the components of
`P2est(object, conf.level)` plus `boot`, an \\M \times M \times M \times
B\\ array of replicate tensors (layout `[j, l, h, b]`), and `B`. The
`estimate` table gains a column `se.boot` (bootstrap standard
deviation).

## Details

Resampling whole subjects keeps the dependence between the several
instants a subject contributes to the same or to different \\(h, j)\\
risk sets. The implementation never re-reads the data: the per-subject
counts of every observed \\(h, j, \ell)\\ triple are tabulated once, and
each replicate is a weighted sum of those counts, so thousands of
replicates are cheap even for large cohorts.

In a replicate where an observed \\(h, j)\\ pair happens to have no
subject at risk, its row is taken from the point estimate, so that the
replicate tensor stays a proper set of transition probabilities.

## Why a bootstrap of subjects

Patients are the independent units; the days of one patient are not.
Resampling whole patients keeps that within-patient dependence, which is
the standard non-parametric bootstrap for clustered data (Davison and
Hinkley, 1997; Field and Welsh, 2007). The n-step predictions are
non-linear (polynomial) functions of all the estimated probabilities, so
a delta-method variance would be cumbersome; propagating each replicate
and taking percentiles (Efron and Tibshirani, 1993) gives their
intervals directly.

## References

Davison, A. C. and Hinkley, D. V. (1997). *Bootstrap Methods and their
Application*. Cambridge University Press.

Efron, B. and Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman and Hall.

Field, C. A. and Welsh, A. H. (2007). Bootstrapping clustered data.
*Journal of the Royal Statistical Society: Series B*, 69(3), 369-390.

## See also

[`ckequations`](https://jcarmezim.github.io/mstate2/reference/ckequations.md),
[`compare2`](https://jcarmezim.github.io/mstate2/reference/compare2.md)

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
bt <- P2boot(prep2(panel), B = 100, seed = 1)
bt
#> <P2boot>  RPE estimates with 100 subject-level bootstrap replicates
#>   states: A, B, C
#>   95% percentile intervals for n-step predictions; 500 subjects
#>   4 estimated transition probabilities (h -> j -> l)
ckequations(bt, h = "A", j = "B", l = "C", nsteps = 4, bounds = TRUE)
#> # A tibble: 4 × 4
#>       n estimate lower upper
#>   <int>    <dbl> <dbl> <dbl>
#> 1     1    0.376 0.343 0.416
#> 2     2    0.804 0.777 0.828
#> 3     3    0.938 0.921 0.951
#> 4     4    0.981 0.972 0.986
```
