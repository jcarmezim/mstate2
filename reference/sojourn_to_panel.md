# Convert sojourn-time data to discrete-time panel format

Turns one row per subject, with the time spent in each transient state
and an indicator of the absorbing state that ended follow-up, into one
row per subject and time unit with the state occupied. This is the
discretisation used in the illustration of the methods paper: each
sojourn is rounded with
[`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md) and the
states are visited in the order given in `segments`.

## Usage

``` r
sojourn_to_panel(data, id, segments, absorbing, round_fun = rnd)
```

## Arguments

- data:

  A data frame, one row per subject.

- id:

  Name of the subject id column.

- segments:

  Named character vector mapping transient state labels to their
  duration columns, in visiting order, e.g.
  `c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")`.

- absorbing:

  Named character vector mapping absorbing state labels to their 0/1
  indicator columns, e.g. `c(Disch = "disch.s", Death = "death.s")`.

- round_fun:

  Discretisation function for durations. Default `rnd`. Each sojourn is
  rounded once and then expanded into that many time units. A visit with
  a positive duration that rounds to 0 (with `rnd`, a duration below
  0.5) is dropped, together with the transitions into and out of it; a
  warning reports how many visits were lost, by state.

## Value

A tibble in panel format with columns `id`, `time` (0, 1, 2, ... within
each subject) and `state` (character). A subject whose follow-up ended
in an absorbing state gets one final row with that state; a subject with
no indicator equal to 1 (still under observation) has no absorbing row.

## Examples

``` r
d <- data.frame(id = 1:2, t.a = c(2, 1), t.b = c(1, 3),
                disch.s = c(1, 0), death.s = c(0, 1))
sojourn_to_panel(d, id = "id",
                 segments = c(A = "t.a", B = "t.b"),
                 absorbing = c(Disch = "disch.s", Death = "death.s"))
#> # A tibble: 9 × 3
#>      id  time state
#>   <int> <int> <chr>
#> 1     1     0 A    
#> 2     1     1 A    
#> 3     1     2 B    
#> 4     1     3 Disch
#> 5     2     0 A    
#> 6     2     1 B    
#> 7     2     2 B    
#> 8     2     3 B    
#> 9     2     4 Death
```
