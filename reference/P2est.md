# Estimate 1-step second-order transition probabilities

Estimates \\P\_{hj\ell} = P(X_s = \ell \mid X\_{s-1} = j, X\_{s-2} =
h)\\ from an
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md) object
with the relative probability estimator (RPE) of the methods paper
(Najera-Zuloaga, Besalu and Gomez Melis, 2025), with its standard error
and confidence interval. The RPE (Eq. 9) pools all transitions and all
subject-instants at risk of each pair \\(h, j)\\: \$\$\tilde P\_{hj\ell}
= \sum_s \tilde N\_{hj\ell}(s) / \sum_s \tilde Y\_{hj}(s-1).\$\$ Its
estimated asymptotic variance (Theorems 4-5, Eqs. 13-14) gives the
standard error \$\$se = \sqrt{\tilde P\_{hj\ell}(1 - \tilde P\_{hj\ell})
/ \sum_s \tilde Y\_{hj}(s-1)},\$\$ and the Wald interval \\\tilde
P\_{hj\ell} \pm z\_{1-\alpha/2}\\ se\\ is that of Corollary 2.

## Usage

``` r
P2est(object, conf.level = 0.95, ci = c("wald", "logit"), clip = TRUE)
```

## Arguments

- object:

  An "msm2data" object from
  [`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md).

- conf.level:

  Confidence level. Default 0.95.

- ci:

  Confidence-interval type: "wald" (default) or "logit". The logit
  interval is built on the log-odds scale (delta method) and stays
  inside (0, 1), behaving better for probabilities near 0 or 1. At the
  boundary (p = 0 or p = 1, e.g. an absorbing-state row) the logit
  transform is undefined and the interval degenerates to the point
  estimate \[p, p\] rather than \[0, 1\].

- clip:

  For ci = "wald", clip the interval to \[0, 1\] (default TRUE).

## Value

An object of class "P2est": the `estimate` tibble (one row per observed
transition \\h \to j \to \ell\\, with columns
`h, j, l, p, se, lower, upper, n.trans, at.risk`), the point/CI/se
tensors `P`, `P.lower`, `P.upper`, `P.se` (layout `P[j, l, h]`, the
layout of the paper's code), and metadata.

## Details

The estimates, standard errors and Wald limits (clipped to \[0, 1\]) are
those of the code of the paper's illustration, except that
\\z\_{1-\alpha/2}\\ is computed exactly (1.959964 for 95%; the paper's
code uses 1.96). Pairs \\(h, j)\\ with nobody at risk get probability 0,
as in the paper. For an absorbing state \\a\\, \\P\_{haa} = 1\\ for
every \\h\\, including pairs that never occur (the paper's code leaves
those at 0; they have probability 0 of being reached, so predictions do
not change).

## Interval choice

Wald intervals for a proportion are known to behave poorly when the
probability is close to 0 or 1 or the number at risk is small (Brown,
Cai and DasGupta, 2001); `ci = "logit"` is then preferable. The
estimator itself is the relative probability estimator of the methods
paper, which is more efficient than the conditional probability
estimator studied there and is therefore the only one implemented.

## References

Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation
for a binomial proportion. *Statistical Science*, 16(2), 101-133.

Najera-Zuloaga, J., Besalu, M. and Gomez Melis, G. (2025). Second-order
Markov multistate models: nonparametric estimation and inference.
Manuscript submitted for publication.

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
fit$estimate
#> # A tibble: 4 × 9
#>   h     j     l         p     se lower upper n.trans at.risk
#>   <fct> <fct> <fct> <dbl>  <dbl> <dbl> <dbl>   <int>   <int>
#> 1 A     B     B     0.624 0.0217 0.582 0.666     312     500
#> 2 A     B     C     0.376 0.0217 0.334 0.418     188     500
#> 3 B     B     B     0.314 0.0218 0.272 0.357     143     455
#> 4 B     B     C     0.686 0.0218 0.643 0.728     312     455
```
