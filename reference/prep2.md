# Build second-order counting processes from panel data

`prep2` takes discrete-time panel data (one row per subject-time with
the occupied state) and constructs the second-order counting processes
needed to estimate 1-step second-order homogeneous transition
probabilities \\P\_{hj\ell} = P(X_s = \ell \mid X\_{s-1} = j, X\_{s-2} =
h)\\.

## Usage

``` r
prep2(
  data,
  id = "id",
  time = "time",
  state = "state",
  states = NULL,
  absorbing = NULL,
  drop.na = FALSE,
  check.consecutive = TRUE
)
```

## Arguments

- data:

  A data frame in long/panel format (one row per subject and time
  point), e.g. the output of
  [`sojourn_to_panel`](https://jcarmezim.github.io/mstate2/reference/sojourn_to_panel.md)
  or
  [`simulate2`](https://jcarmezim.github.io/mstate2/reference/simulate2.md),
  or an `"msm2prep"` object from
  [`msprep2`](https://jcarmezim.github.io/mstate2/reference/msprep2.md)
  (its panel is used, with its states and absorbing states unless
  `states` or `absorbing` are given).

- id, time, state:

  Column names for subject id, discrete time, and state.

- states:

  Optional state space / ordering. Defaults to sorted observed states.

- absorbing:

  Optional absorbing states; inferred if NULL.

- drop.na:

  If TRUE, drop rows with NA in id/time/state (with a warning) instead
  of erroring. Default FALSE.

- check.consecutive:

  If TRUE (default), warn when any subject's times are not consecutive
  integers.

## Value

An object of class "msm2data": a list with the tibbles `N` (the counts
\\\tilde N\_{hj\ell}(s)\\, columns `h, j, l, s, N`), `Y` (the at-risk
counts \\\tilde Y\_{hj}(s-1)\\, columns `h, j, s, Y`) and `triples` (one
row per subject per observed triple, columns `id, h, j, l, s`; used by
[`P2boot`](https://jcarmezim.github.io/mstate2/reference/P2boot.md) to
resample subjects), where `h, j, l` are factors with levels `states`;
and `states`, `absorbing`, `n` (subjects), `ntriples` and `time.range`.

## Details

For each subject the function forms every consecutive triple of
observations \\(X\_{s-2}, X\_{s-1}, X_s) = (h, j, \ell)\\ and tabulates:

- \\\tilde N\_{hj\ell}(s)\\: number of subjects observing the path \\h
  \to j \to \ell\\ ending at time \\s\\;

- \\\tilde Y\_{hj}(s-1)\\: number of subjects at risk, i.e. occupying
  \\h\\ at \\s-2\\ and \\j\\ at \\s-1\\.

These counting processes are indexed by the time \\s \ge 2\\ of the
destination state. The at-risk count is computed as \\\tilde
Y\_{hj}(s-1) = \sum\_\ell \tilde N\_{hj\ell}(s)\\: with complete
follow-up, as in DIVINE, this is exactly the number of subjects in \\h\\
at \\s-2\\ and \\j\\ at \\s-1\\; a subject whose follow-up stops in
\\j\\ at \\s-1\\ (right-censored) is not counted at risk at \\s\\,
because the next state is not observed.

Consecutive *observations* are treated as consecutive time steps.

## Why these counts

The likelihood of a second-order Markov chain depends on the data only
through the transition counts \\\tilde N\_{hj\ell}\\ and the at-risk
counts \\\tilde Y\_{hj}\\ (Anderson and Goodman, 1957), so they are
computed once and every other function reuses them.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. *Annals of Mathematical Statistics*, 28(1), 89-110.

Najera-Zuloaga, J., Besalu, M. and Gomez Melis, G. (2025). Second-order
Markov multistate models: nonparametric estimation and inference.
Manuscript submitted for publication.

## Examples

``` r
st <- c("A", "B", "C")   # C is absorbing
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
d <- prep2(panel, id = "id", time = "time", state = "state")
d
#> <msm2data>  second-order counting processes
#>   subjects        : 500
#>   observed triples: 955
#>   time range      : 0 - 7
#>   states (3)      : A, B, C
#>   absorbing       : C
#>   distinct (h,j)  : 2
```
