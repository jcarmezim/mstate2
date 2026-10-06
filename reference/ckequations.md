# Extended Chapman-Kolmogorov n-step second-order transition probabilities

Computes \\P\_{hj\ell}(1,n) = P(X\_{n+1} = \ell \mid X_1 = j, X_0 = h)\\
for \\n = 1, \dots, nsteps\\ via the extension of the Chapman-Kolmogorov
relation. The second-order chain is lifted to a first-order chain on
ordered pairs of states, \\Q\_{(a,b)\to(b,c)} = P\_{abc}\\, and the
initial pair distribution is propagated; this is exact, works for any
number of states and any horizon, and is far cheaper than expanding the
path sum.

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

  Number of steps n (default 9), giving \\X_2, \ldots, X\_{n+1}\\.

- bounds:

  If x is a "P2est" object and a single l is given, also the CI tensors
  to return evolution-interval bounds. If x is a
  [`P2boot`](https://jcarmezim.github.io/mstate2/reference/P2boot.md)
  object, the bounds are instead bootstrap intervals of the n-step
  probability (every replicate tensor propagated), which have the
  nominal coverage the evolution intervals lack.

## Value

If l is a single state and bounds = FALSE: a numeric vector of length
nsteps. If l is a vector: an \\nsteps \times length(l)\\ matrix. If l is
NULL: an \\nsteps \times M\\ matrix (one column per state). If bounds =
TRUE: a tibble with columns n, estimate, lower, upper.

## Details

Step \\n\\ is the probability of being in \\\ell\\ at time \\s = n +
1\\; \\n = 1\\ is the one-step probability \\P\_{hj\ell}\\. The
evolution intervals are obtained by propagating the tensors of lower and
upper confidence limits in the same way.

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
st <- c("A", "B", "C")                                 # C is absorbing
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
ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6)
#> [1] 0.624000000 0.196114286 0.061635918 0.019371289 0.006088119 0.001913409
ckequations(fit, h = "A", j = "B", l = "B", nsteps = 6, bounds = TRUE)
#> # A tibble: 6 × 4
#>       n estimate    lower   upper
#>   <int>    <dbl>    <dbl>   <dbl>
#> 1     1  0.624   0.582    0.666  
#> 2     2  0.196   0.158    0.238  
#> 3     3  0.0616  0.0429   0.0849 
#> 4     4  0.0194  0.0117   0.0303 
#> 5     5  0.00609 0.00317  0.0108 
#> 6     6  0.00191 0.000860 0.00386
```
