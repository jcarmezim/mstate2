# Estimate 1-step second-order transition probabilities

Estimates \\P\_{hj\ell} = P(X_s = \ell \mid X\_{s-1} = j, X\_{s-2} =
h)\\ from an
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md) object
with the relative probability estimator (RPE, Eq. 9 / Theorem 5), with
standard errors and confidence intervals (Corollaries 1-2). The RPE
pools all transitions and all subject-instants at risk of each pair
\\(h, j)\\: \$\$\hat P\_{hj\ell} = \sum_s \tilde N\_{hj\ell}(s) / \sum_s
\tilde Y\_{hj}(s-1).\$\$

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

An object of class "P2est": the `estimate` data frame, the point/CI/se
tensors `P`, `P.lower`, `P.upper`, `P.se` (layout `P[j, l, h]`), and
metadata.

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
#>   h j l         p         se     lower     upper n.trans at.risk
#> 1 A B B 0.6240000 0.02166213 0.5815430 0.6664570     312     500
#> 2 A B C 0.3760000 0.02166213 0.3335430 0.4184570     188     500
#> 3 B B B 0.3142857 0.02176347 0.2716301 0.3569413     143     455
#> 4 B B C 0.6857143 0.02176347 0.6430587 0.7283699     312     455
```
