# Simulate a second-order Markov multistate process

Generates discrete-time panel data from a second-order Markov model
defined by a 1-step second-order transition tensor, a first-step
(first-order) transition matrix, and entry probabilities. This
generalises the simulation in Section 5 of the paper and produces output
ready for
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md).

## Usage

``` r
simulate2(
  n,
  tensor,
  first,
  init = NULL,
  entry = NULL,
  states = NULL,
  maxT = 1000
)
```

## Arguments

- n:

  Number of individuals.

- tensor:

  An \\M\times M\times M\\ tensor in the `P[j,l,h]` \\= P\_{hj\ell}\\
  layout (e.g. the `P` component of a
  [`P2est`](https://jcarmezim.github.io/mstate2/reference/P2est.md)
  object). Absorbing states must satisfy `tensor[a,a,h] = 1`.

- first:

  An \\M\times M\\ first-order matrix with `first[h,l] =` \\P(X_1 = \ell
  \mid X_0 = h)\\.

- init:

  Probability vector over states for the entry state \\X_0\\ (used when
  `entry = NULL`). Defaults to uniform over non-absorbing states.

- entry:

  Optional named probability vector over entry states; if its sum is
  below 1 the remainder is the per-step probability of *not* yet
  entering, giving staggered entry. If `NULL`, all individuals enter at
  time 0 according to `init`.

- states:

  Optional state labels (defaults to the tensor's dimnames or `1:M`).

- maxT:

  Safety cap on the number of global time steps. Default 1000.

## Value

A `data.table` in panel format (`id`, `time`, `state`), suitable for
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md).

## Details

Each individual receives an entry state \\X_0\\ (its personal time
origin), then \\X_1\\ from the first-step matrix, and subsequent states
from the second-order tensor \\X_s \sim P\_{X\_{s-2}, X\_{s-1},
\cdot}\\. Sampling stops when an absorbing state is reached.

With `entry` supplied, individuals enter the process stochastically over
*global* time (the staggered-entry mechanism of the paper's auxiliary
state 0): at each step every not-yet-entered individual enters state
\\h\\ with probability `entry[h]`, otherwise waits, so that the number
of individuals at risk varies across global time. With `entry = NULL`
all individuals start at global time 0.

## Examples

``` r
st   <- c("A", "B", "C")                                 # C is absorbing
tens <- array(0, c(3, 3, 3), dimnames = list(st, st, st))
tens["B", "B", "A"] <- 0.6; tens["B", "C", "A"] <- 0.4
tens["B", "B", "B"] <- 0.3; tens["B", "C", "B"] <- 0.7
tens["C", "C", ]    <- 1
first <- matrix(0, 3, 3, dimnames = list(st, st)); first["A", "B"] <- 1

set.seed(1)
panel <- simulate2(200, tens, first, init = c(A = 1, B = 0, C = 0))
head(panel)
#>       id  time  state
#>    <int> <int> <fctr>
#> 1:     1     0      A
#> 2:     1     1      B
#> 3:     1     2      C
#> 4:     2     0      A
#> 5:     2     1      B
#> 6:     2     2      B
```
