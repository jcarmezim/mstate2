# Extended Chapman-Kolmogorov n-step second-order transition probabilities

Computes \\P\_{hj\ell}(1,n) = P(X\_{n+1} = \ell \mid X_1 = j, X_0 = h)\\
for \\n = 1, \dots, nsteps\\ via the extension of the Chapman-Kolmogorov
relation (Eq. 6). The second-order chain is lifted to a first-order
chain on ordered pairs of states, \\Q\_{(a,b)\to(b,c)} = P\_{abc}\\, and
the initial pair distribution is propagated; this is exact, works for
any number of states and any horizon, and is far cheaper than expanding
the path sum.

## Usage

``` r
ckequations(x, h, j, l = NULL, nsteps = 9L, bounds = FALSE)
```

## Arguments

- x:

  A "P2est" object or an \\M \times M \times M\\ tensor in P\[j, l, h\]
  layout.

- h, j:

  Starting states: h at time 0, j at time 1 (label or index).

- l:

  Target state(s). If NULL (default), the full distribution over all
  states is returned for each step.

- nsteps:

  Number of steps n (default 9), giving \\X_2, \ldots, X\_\\n+1\\\\.

- bounds:

  If x is a "P2est" object and a single l is given, also propagate the
  CI tensors to return evolution-interval bounds (clipped to \[0, 1\]).
  If x is a
  [`P2boot`](https://jcarmezim.github.io/mstate2/reference/P2boot.md)
  object, the bounds are instead percentile bootstrap intervals of the
  n-step probability (every replicate tensor is propagated), which have
  the nominal coverage the evolution intervals lack.

## Value

If l is a single state and bounds = FALSE: a numeric vector of length
nsteps. If l is a vector: an \\nsteps \times length(l)\\ matrix. If l is
NULL: an \\nsteps \times M\\ matrix (one column per state). If bounds =
TRUE: a data.frame with columns n, estimate, lower, upper.

## Why the chain on pairs

A second-order chain on \\M\\ states is a first-order chain on the
\\M^2\\ pairs of consecutive states, a standard representation of
higher-order Markov chains (e.g. Benson, Gleich and Lim, 2017).
Multiplying the pair distribution by its transition matrix once per step
gives the extended Chapman-Kolmogorov relation of Najera-Zuloaga, Besalu
and Gomez Melis (2025) exactly, at a cost linear in the horizon, without
enumerating the \\M^{n-1}\\ paths.

## References

Benson, A. R., Gleich, D. F. and Lim, L.-H. (2017). The spacey random
walk: a stochastic process for higher-order data. *SIAM Review*, 59(2),
321-345.

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
ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6)
#> [1] 0.624000000 0.196114286 0.061635918 0.019371289 0.006088119 0.001913409
ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6, bounds = TRUE)
#>   n    estimate        lower       upper
#> 1 1 0.624000000 0.5815429998 0.666457000
#> 2 2 0.196114286 0.1579645768 0.237886054
#> 3 3 0.061635918 0.0429079321 0.084911367
#> 4 4 0.019371289 0.0116550854 0.030308377
#> 5 5 0.006088119 0.0031658719 0.010818313
#> 6 6 0.001913409 0.0008599461 0.003861503
```
