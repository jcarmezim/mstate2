# Prepare multistate data for the analysis

Turns multistate data into the discrete-time panel that
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md) needs
(one row per subject and time unit, with the state occupied), and
reports every record it had to drop or change. It plays the role of
[`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html) for
second-order models, with three ways of giving the data:

- Wide (as
  [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html)):

  One row per subject, with the time of entry into each state, or of
  censoring when its status is 0. `states` is a named list with, in
  order, each state and its `Surv(time, status)`; the state given as
  `NULL` is the initial state of everybody. The column names are written
  without quotes, and survival does not need to be attached.

- Sojourn (as the DIVINE data):

  One row per subject, with the time spent in each transient state
  (`durations`, in visiting order) and 0/1 indicators of the absorbing
  state that ended follow-up (`outcome`).

- msdata:

  An object made by
  [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html); its
  transition matrix is used.

## Usage

``` r
msprep2(
  data,
  states = NULL,
  trans = NULL,
  durations = NULL,
  outcome = NULL,
  id = "id",
  keep = NULL,
  unit = 1
)
```

## Arguments

- data:

  A data frame, or an `msdata` object from
  [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html).

- states:

  Wide data: a named list with the states in order, each one
  `Surv(time, status)` or `NULL` (the initial state), e.g.
  `list(healthy = NULL, ill = Surv(ill_time, ill_status), dead = Surv(dead_time, dead_status))`.
  Sojourn data: optionally, the order of the states (a character vector;
  default the order of `durations` and `outcome`).

- trans:

  Optional allowed transitions, written as text,
  `c("healthy -> ill -> dead", "healthy -> dead")` (a chain
  `"A -> B -> C"` gives A -\> B and B -\> C), or a matrix from
  [`mstate::transMat()`](https://rdrr.io/pkg/mstate/man/transMat.html).
  It also gives the absorbing states (no exit). An `msdata` object
  carries its own.

- durations:

  Sojourn data: named character vector mapping each transient state to
  the column with the time spent in it, in visiting order, e.g.
  `c(healthy = "t_healthy", ill = "t_ill")`.

- outcome:

  Sojourn data: named character vector mapping each absorbing state to
  its 0/1 indicator column, e.g. `c(dead = "dead")`.

- id:

  Name of the subject id column (if missing, the rows are numbered).

- keep:

  Optional names of the baseline covariate columns to carry into the
  panel (default: all the columns that are not used).

- unit:

  Length of one time unit in the units of the data, e.g. `30.4375` to go
  from days to months. Default 1.

## Value

An object of class `"msm2prep"`, which
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md)
accepts directly: a list with

- `panel`:

  tibble `(id, time, state, covariates...)`, one row per subject and
  time unit, `state` a factor with levels `states`;

- `states`, `absorbing`, `trans`:

  state space, absorbing states and the allowed-transition matrix
  (logical, or `NULL`);

- `transitions`:

  tibble `(from, to, n, allowed)` with the number of each observed
  transition;

- `subjects`:

  tibble with one row per subject: first and last time, entry and exit
  state, number of rows and `status` (`"absorbed"` or `"censored"`);

- `issues`:

  tibble `(id, issue, state, time, detail)` with every record dropped or
  changed;

- `settings`:

  the layout and the unit.

## Details

`Surv(time, status)` checks its arguments as
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) does:
the time must be numeric and the status 0/1, `TRUE`/`FALSE` or 1/2 (1 =
censored, 2 = event, with a warning); other codes become `NA` with a
warning. A status written as text must say which value is the event,
e.g. `Surv(dead_time, dead_status == "yes")`.

The times are divided by `unit` and rounded with
[`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md). A subject
occupies each state from the unit in which it was entered until the unit
before the next state; follow-up ends at the entry into an absorbing
state or, otherwise, at the last recorded time (entry or censoring). The
absorbing states are the states with no exit in `trans`, the states of
`outcome`, or the states nobody is seen leaving. All the columns that
are not used (`id`, the columns in
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), `durations`,
`outcome`, or the columns of `msdata`) are kept as baseline covariates.

Every record that is dropped or changed is listed in the `issues` table
of the result, and a single warning gives the counts:

- `missing`:

  a state with status 1 but no time, or a time without a valid status
  (wide);

- `rounded_to_zero`:

  a stay with a positive duration rounds to 0 units (sojourn);

- `same_unit`:

  two different states fall in the same time unit: the state entered
  last (or an absorbing state) is kept, so the other one occupies no
  unit and its transitions are lost;

- `after_absorbing`:

  a state entered after the absorbing state;

- `before_start`:

  a negative time;

- `not_allowed`:

  a transition not allowed by `trans` (kept);

- `no_data`:

  a subject with no usable record.

## Differences with mstate

:msprep(): [`msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html)
needs the transition matrix from state numbers, the time and status
columns in two parallel vectors (`NA` for the initial state) and the
covariates listed in `keep`, and returns one row per subject and
possible transition (the counting-process format of the Cox model).
`msprep2()` gives each state its own `Surv(time, status)`, reads the
transitions as text, keeps the other columns as covariates, and returns
the discrete-time panel of second-order models, with a report of the
data.

## References

de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package
for the analysis of competing risks and multi-state models. *Journal of
Statistical Software*, 38(7), 1-30.

## See also

[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md),
[`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md)

## Examples

``` r
## Illness-death model: healthy -> ill -> dead and healthy -> dead

## Sojourn: the structure of the DIVINE data
sojourn <- data.frame(id = 1:4, t_healthy = c(3, 5, 2, 7), t_ill = c(4, 0, 8, 0),
                      dead = c(1, 1, 0, 0))
p1 <- msprep2(sojourn, durations = c(healthy = "t_healthy", ill = "t_ill"),
              outcome = c(dead = "dead"))
p1
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : sojourn
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 31 (time 0 - 9)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 3 types, 4 in total
#>   issues          : none

## Wide: the data that go into mstate::msprep()
wide <- data.frame(id = 1:4, ill_time = c(3, 5, 2, 6), ill_status = c(1, 0, 1, 0),
                   dead_time = c(7, 5, 9, 6), dead_status = c(1, 1, 0, 0),
                   age = c(60, 72, 55, 49))
p2 <- msprep2(wide,
              states = list(healthy = NULL,
                            ill     = Surv(ill_time, ill_status),
                            dead    = Surv(dead_time, dead_status)),
              trans  = c("healthy -> ill -> dead", "healthy -> dead"))
p2$panel                       # `age` is kept as a covariate
#> # A tibble: 31 × 4
#>       id  time state     age
#>    <int> <int> <fct>   <dbl>
#>  1     1     0 healthy    60
#>  2     1     1 healthy    60
#>  3     1     2 healthy    60
#>  4     1     3 ill        60
#>  5     1     4 ill        60
#>  6     1     5 ill        60
#>  7     1     6 ill        60
#>  8     1     7 dead       60
#>  9     2     0 healthy    72
#> 10     2     1 healthy    72
#> # ℹ 21 more rows
summary(p2)
#> <msm2prep summary>
#> 
#> Observed transitions (from rows to columns):
#>          to
#> from      healthy ill dead
#>   healthy       0   2    1
#>   ill           0   0    1
#>   dead          0   0    0
#> 
#> Time units observed per subject (rows of the panel), by status:
#> # A tibble: 2 × 5
#>   status   subjects   min median   max
#>   <chr>       <int> <dbl>  <dbl> <dbl>
#> 1 absorbed        2     6    7       8
#> 2 censored        2     7    8.5    10
#> 
#> Records dropped or changed:
#>   none

## msdata: the output of mstate::msprep()
if (requireNamespace("mstate", quietly = TRUE)) {
  tmat <- mstate::transMat(x = list(c(2, 3), 3, c()), names = c("healthy", "ill", "dead"))
  ms <- mstate::msprep(time = c(NA, "ill_time", "dead_time"),
                       status = c(NA, "ill_status", "dead_status"),
                       data = wide, trans = tmat)
  p3 <- msprep2(ms)
  p3
}
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : msdata
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 31 (time 0 - 9)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 3 types, 4 in total (0 not allowed by `trans`)
#>   issues          : none

## The three give the same panel
identical(p1$panel$state, p2$panel$state)
#> [1] TRUE

## What is reported: a stay of 0.4 days rounds to 0 days ...
bad <- sojourn
bad$t_ill[2] <- 0.4
msprep2(bad, durations = c(healthy = "t_healthy", ill = "t_ill"),
        outcome = c(dead = "dead"))$issues
#> Warning: 1 record(s) were dropped or changed while building the panel (rounded_to_zero: 1); see the `issues` table of the result.
#> # A tibble: 1 × 5
#>      id issue           state time  detail                        
#>   <int> <chr>           <chr> <chr> <chr>                         
#> 1     2 rounded_to_zero ill   NA    duration 0.4 rounds to 0 units
## ... and a transition not allowed by `trans`
msprep2(wide, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                            dead = Surv(dead_time, dead_status)),
        trans = "healthy -> ill -> dead")$issues
#> Warning: 1 record(s) were dropped or changed while building the panel (not_allowed: 1); see the `issues` table of the result.
#> # A tibble: 1 × 5
#>      id issue       state time  detail                   
#>   <int> <chr>       <chr> <chr> <chr>                    
#> 1     2 not_allowed dead  NA    healthy -> dead at unit 5

## The mstate workflow of ebmt3 in one call (times in days, panel in months)
if (requireNamespace("mstate", quietly = TRUE)) {
  data(ebmt3, package = "mstate")
  bmt <- msprep2(ebmt3,
                 states = list(Tx       = NULL,
                               PR       = Surv(prtime, prstat),
                               RelDeath = Surv(rfstime, rfsstat)),
                 trans  = c("Tx -> PR -> RelDeath", "Tx -> RelDeath"),
                 unit   = 30.4375)
  bmt
}
#> Warning: 110 record(s) were dropped or changed while building the panel (same_unit: 110); see the `issues` table of the result.
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : wide
#>   subjects        : 2204 (841 absorbed, 1363 censored)
#>   panel rows      : 69537 (time 0 - 93)
#>   states (3)      : Tx, PR, RelDeath
#>   absorbing       : RelDeath
#>   transitions     : 3 types, 1900 in total (0 not allowed by `trans`)
#>   issues          : 110 record(s) dropped or changed (same_unit: 110); see x$issues

## The result goes straight into prep2()
prep2(p2)
#> <msm2data>  second-order counting processes
#>   subjects        : 4
#>   observed triples: 23
#>   time range      : 0 - 9
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   distinct (h,j)  : 3
```
