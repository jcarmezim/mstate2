# Prepare raw multistate data for the analysis

One entry point that turns multistate data as they are usually collected
into the discrete-time panel that
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md) needs
(one row per subject and time unit, with the state occupied), and
reports every record it had to drop or change. It plays the role of
[`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html), but
works from four input layouts, understands dates, discretises time,
cleans the data and never changes them silently.

## Usage

``` r
msprep2(
  data,
  id = "id",
  format = c("auto", "events", "wide", "sojourn", "msdata"),
  time = "time",
  state = "state",
  times = NULL,
  status = NULL,
  durations = NULL,
  outcome = NULL,
  initial = NULL,
  start = NULL,
  end = NULL,
  unit = 1,
  round_fun = rnd,
  recode = NULL,
  states = NULL,
  absorbing = NULL,
  trans = NULL,
  ties = c("last", "first"),
  check = c("warn", "error"),
  keep = NULL
)
```

## Arguments

- data:

  A data frame with the raw data (or an `msdata` object).

- id:

  Name of the subject id column.

- format:

  Input layout: `"auto"` (default), `"events"`, `"wide"`, `"sojourn"` or
  `"msdata"` (see Details).

- time, state:

  Events layout: names of the time and state columns.

- times:

  Wide layout: named character vector mapping each state to its time
  column, e.g.
  `c(SP = "date_sp", IMV = "date_imv", Death = "date_death")`.

- status:

  Wide layout: optional named character vector mapping each state (names
  as in `times`) to its 0/1 status column.

- durations:

  Sojourn layout: named character vector mapping each transient state to
  its duration column, in visiting order.

- outcome:

  Sojourn layout: named character vector mapping each absorbing state to
  its 0/1 indicator column.

- initial:

  Optional initial state, entered at the origin: a column of `data` or a
  single state label. Useful when the records only list the changes
  after the start (e.g. the state at admission is in a separate column).
  With a list of states, a state with no columns (`NULL`) is the initial
  state of everybody, so `initial` is not needed.

- start:

  Origin of time for each subject: a column of `data` (e.g. the
  admission date) or a single value. Default: the subject's first record
  (events, dates in the wide layout) or 0 (numeric wide layout).

- end:

  End of follow-up for subjects not absorbed: a column of `data` or a
  single value (e.g. the date the data were extracted). Default: none
  (events), or the largest recorded time of the subject when `status` is
  given (wide layout, as in `msprep()`).

- unit:

  Length of one time unit: a number in the units of the times, or a name
  (`"hour"`, `"day"`, `"week"`, `"month"`, `"year"`) when the times are
  dates. Default 1 (one day for dates).

- round_fun:

  Discretisation of the elapsed times. Default
  [`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md).

- recode:

  Optional named character vector to relabel the recorded states,
  `c(old = "new")`, e.g. `c("1" = "NSP", "2" = "SP")`. Applied to
  `state` and `initial`.

- states:

  Either the state space and its order (a character vector; default: the
  row names of `trans`, the order of `times` or `durations` and
  `outcome`, or the sorted observed states), or a named list giving, for
  each state, its time and status columns (wide layout; see Details).

- absorbing:

  Optional absorbing states. Default: the states with no allowed exit in
  `trans`, the states of `outcome`, or the states nobody is seen
  leaving.

- trans:

  Optional allowed transitions, written as text,
  `c("Tx -> PR", "Tx -> RelDeath", "PR -> RelDeath")` (a chain
  `"A -> B -> C"` is allowed); as a named list of destinations,
  `list(Tx = c("PR", "RelDeath"), PR = "RelDeath", RelDeath = NULL)`; or
  as a matrix with the states as row and column names
  ([`mstate::transMat()`](https://rdrr.io/pkg/mstate/man/transMat.html),
  `NA` = not allowed, or logical or 0/1). It also gives the absorbing
  states (no exit).

- ties:

  Which state to keep when several fall in the same time unit: `"last"`
  (default, the state at the end of the unit) or `"first"`.

- check:

  What to do with transitions not allowed by `trans`: `"warn"` (default;
  they are kept and listed) or `"error"`.

- keep:

  Optional names of baseline covariate columns to carry into the panel.
  With conventional names, all the columns that are not `id`, `inistat`,
  `<state>_time` or `<state>_status`.

## Value

An object of class `"msm2prep"`, which
[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md)
accepts directly: a list with

- `panel`:

  tibble `(id, time, state, keep...)`, one row per subject and time
  unit, `time` counted from the origin, `state` a factor with levels
  `states`;

- `states`, `absorbing`, `trans`:

  state space, absorbing states and the allowed-transition matrix
  (logical, or `NULL`);

- `transitions`:

  tibble `(from, to, n, allowed)` with the number of each observed
  transition;

- `subjects`:

  tibble with one row per subject: first and last time, entry and exit
  state, `status` (`"absorbed"` or `"censored"`) and number of rows;

- `issues`:

  tibble `(id, issue, state, time, detail)` with every record dropped or
  changed;

- `settings`:

  the layout, unit and tie rule used.

## Details

**Input layouts** (`format`):

- `"events"`:

  One row per subject and recorded state, with the time at which the
  state was entered or observed (`time`, `state`). Repeated records of
  the same state are merged, so the rows can be state changes or
  repeated observations of the same patient.

- `"wide"`:

  One row per subject, with one time column per state (`times`), as in
  [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html): the
  time the state was entered, or `NA` if it was never visited. With
  `status`, a state is visited when its status is 1, and the times with
  status 0 are censoring times, as in `msprep()`.

- `"sojourn"`:

  One row per subject, with the time spent in each transient state
  (`durations`, in visiting order) and 0/1 indicators of the absorbing
  state that ended follow-up (`outcome`), as in the DIVINE data. Each
  sojourn is rounded separately, as in the code of the methods paper.

- `"msdata"`:

  An object of class `"msdata"` made by
  [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html); its
  transition matrix is used as `trans`.

With `format = "auto"` (default) the layout is taken from the arguments
given: an `msdata` object, `durations`, `times`, conventional names
(below), or else events.

**Conventional names: no arguments needed.** If the data are in the wide
layout with the columns named

- `<state>_time`, `<state>_status`:

  for every state, e.g. `death_time` and `death_status` (time of entry,
  or censoring time when the status is 0; the initial states with time
  0),

- `inistat`:

  the initial state of each individual,

- `id`:

  the individual (if missing, the rows are numbered),

then `msprep2(data)` is enough: the states are the `<state>` prefixes,
in column order, and every other column is kept as a covariate. `trans`,
`absorbing`, `unit` and the other arguments can still be given.

**Any column names: a list of states.** With other names, `states` can
be a named list that assigns to each state, in order, its time and
status columns:

    states = list(healthy = c(time = "t_adm", status = "adm"),
                  ill     = c(time = "t_ill", status = "ill"),
                  dead    = c(time = "t_dth", status = "dth"))

The shortest way is
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html), with the column
names written without quotes:

    states = list(healthy = NULL,
                  ill     = Surv(t_ill, ill),
                  dead    = Surv(t_dth, dth == "yes"))

`Surv(time, status)` is evaluated on the columns of `data`, so it can
contain expressions; the times can be numbers or dates (also text), and
the status is checked by
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) (0/1,
`TRUE`/`FALSE`, or 1/2 = censored/event, with a warning; other values
become `NA` with a warning).
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) and the
column-name forms can be mixed in the same list. Each element can also
be `c("t_ill" = "ill")` (time column = status column),
`c("t_ill", "ill")` (two columns: the one with only 0/1 values is the
status), just `"t_ill"` (no status: the state is visited when its time
is recorded) or `NULL` (no columns: a state that is only entered as
initial state). As with conventional names, `initial` gives the initial
state (default `inistat` if present) and every column not named in the
list, `id`, `initial`, `start` or `end` is kept as a covariate.

**Time.** Times can be numbers or dates (`Date`, `POSIXct`, or text such
as `"2020-03-15"` or `"15/03/2020"`). Each time is measured from the
subject's origin (`start`: a column, a common value, or by default the
subject's first record) in units of `unit` (a number in the units of the
times, or `"hour"`, `"day"`, `"week"`, `"month"` or `"year"` for dates)
and discretised with `round_fun`. A subject occupies each state from the
time unit in which it was entered until the unit before the next state;
follow-up ends at the entry into an absorbing state, at the end of
follow-up (`end`), or, without `end`, at the last record (the subject is
then censored in that unit; a repeated record of the same state counts
as a record).

**Cleaning.** Every record that is dropped or changed is listed in the
`issues` table of the result, and a single warning gives the counts:

- `missing`:

  an id, time or state is missing (events), a state has status 1 but no
  time, or a time has no valid status (wide; e.g. a status code that
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) turns into
  `NA`);

- `before_start`:

  the record is before the subject's origin;

- `after_end`:

  the record is after the end of follow-up;

- `same_unit`:

  two different states fall in the same time unit; `ties = "last"` keeps
  the state entered last (an absorbing state is always kept), so the
  other one occupies no unit and its transitions are lost;

- `rounded_to_zero`:

  a visit with a positive duration rounds to 0 time units (sojourn
  layout);

- `after_absorbing`:

  the record is after the entry into an absorbing state;

- `not_allowed`:

  the transition is not allowed by `trans` (kept, unless
  `check = "error"`, which stops);

- `no_data`:

  the subject has no usable record left.

## Differences with mstate

:msprep(): `msprep()` takes only the wide layout, needs the transition
matrix, works with numeric times and returns one row per subject and
possible transition (the counting-process format of the Cox model).
`msprep2()` returns the discrete-time panel of second-order models,
accepts four layouts and dates, infers the states and absorbing states
when no matrix is given, orders the states by the recorded times, and
checks and reports the data instead of failing or changing them
silently.

## References

de Wreede, L. C., Fiocco, M. and Putter, H. (2011). mstate: an R package
for the analysis of competing risks and multi-state models. *Journal of
Statistical Software*, 38(7), 1-30.

## See also

[`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md),
[`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md)

## Examples

``` r
# Events layout: one row per state change, with dates
raw <- data.frame(
  id    = c(1, 1, 1, 2, 2, 3, 3, 3),
  date  = c("2020-03-01", "2020-03-03", "2020-03-06",
            "2020-03-02", "2020-03-09",
            "2020-03-05", "2020-03-06", "2020-03-12"),
  state = c("NSP", "SP", "Disch", "SP", "Death", "NSP", "SP", "SP"))
x <- msprep2(raw, time = "date", state = "state", absorbing = c("Disch", "Death"))
x
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : events
#>   subjects        : 3 (2 absorbed, 1 censored)
#>   panel rows      : 22 (time 0 - 7)
#>   states (4)      : NSP, SP, Disch, Death
#>   absorbing       : Disch, Death
#>   transitions     : 3 types, 4 in total
#>   issues          : none
x$panel
#> # A tibble: 22 × 3
#>       id  time state
#>    <dbl> <int> <fct>
#>  1     1     0 NSP  
#>  2     1     1 NSP  
#>  3     1     2 SP   
#>  4     1     3 SP   
#>  5     1     4 SP   
#>  6     1     5 Disch
#>  7     2     0 SP   
#>  8     2     1 SP   
#>  9     2     2 SP   
#> 10     2     3 SP   
#> # ℹ 12 more rows
prep2(x)
#> <msm2data>  second-order counting processes
#>   subjects        : 3
#>   observed triples: 16
#>   time range      : 0 - 7
#>   states (4)      : NSP, SP, Disch, Death
#>   absorbing       : Disch, Death
#>   distinct (h,j)  : 3

# Conventional names: x_time, x_status and inistat; nothing else to say
conv <- data.frame(id = 1:4, inistat = c("healthy", "healthy", "ill", "healthy"),
                   healthy_time = 0, healthy_status = c(1, 1, 0, 1),
                   ill_time = c(2, 6, 0, 3), ill_status = c(1, 0, 1, 1),
                   dead_time = c(5, 6, 4, 8), dead_status = c(1, 0, 1, 0),
                   age = c(60, 72, 55, 49))
y <- msprep2(conv, absorbing = "dead")
y
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : wide
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 27 (time 0 - 8)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 2 types, 4 in total
#>   issues          : none
y$panel
#> # A tibble: 27 × 4
#>       id  time state     age
#>    <int> <int> <fct>   <dbl>
#>  1     1     0 healthy    60
#>  2     1     1 healthy    60
#>  3     1     2 ill        60
#>  4     1     3 ill        60
#>  5     1     4 ill        60
#>  6     1     5 dead       60
#>  7     2     0 healthy    72
#>  8     2     1 healthy    72
#>  9     2     2 healthy    72
#> 10     2     3 healthy    72
#> # ℹ 17 more rows

# Any column names: a list giving the time and status columns of each state
own <- data.frame(patient = 1:4, start_state = c("healthy", "healthy", "ill", "healthy"),
                  t_ill = c(2, 6, 0, 3), ill = c(1, 0, 1, 1),
                  t_dth = c(5, 6, 4, 8), dth = c(1, 0, 1, 0), age = c(60, 72, 55, 49))
msprep2(own, id = "patient", initial = "start_state", absorbing = "dead",
        states = list(healthy = NULL,
                      ill     = c(time = "t_ill", status = "ill"),
                      dead    = c(time = "t_dth", status = "dth")))$panel
#> # A tibble: 27 × 4
#>       id  time state     age
#>    <int> <int> <fct>   <dbl>
#>  1     1     0 healthy    60
#>  2     1     1 healthy    60
#>  3     1     2 ill        60
#>  4     1     3 ill        60
#>  5     1     4 ill        60
#>  6     1     5 dead       60
#>  7     2     0 healthy    72
#>  8     2     1 healthy    72
#>  9     2     2 healthy    72
#> 10     2     3 healthy    72
#> # ℹ 17 more rows

# The same with Surv(): column names without quotes, expressions allowed
library(survival)
msprep2(own, id = "patient", initial = "start_state", absorbing = "dead",
        states = list(healthy = NULL,
                      ill     = Surv(t_ill, ill),
                      dead    = Surv(t_dth, dth == 1)))$panel
#> # A tibble: 27 × 4
#>       id  time state     age
#>    <int> <int> <fct>   <dbl>
#>  1     1     0 healthy    60
#>  2     1     1 healthy    60
#>  3     1     2 ill        60
#>  4     1     3 ill        60
#>  5     1     4 ill        60
#>  6     1     5 dead       60
#>  7     2     0 healthy    72
#>  8     2     1 healthy    72
#>  9     2     2 healthy    72
#> 10     2     3 healthy    72
#> # ℹ 17 more rows

# The mstate::msprep() workflow of the ebmt3 data in one call
if (requireNamespace("mstate", quietly = TRUE)) {
  data(ebmt3, package = "mstate")
  bmt <- msprep2(ebmt3,
                 states = list(Tx = NULL,
                               PR = Surv(prtime, prstat),
                               RelDeath = Surv(rfstime, rfsstat)),
                 trans = c("Tx -> PR -> RelDeath", "Tx -> RelDeath"),
                 unit = 30.4375)                       # times in days, panel in months
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

# Wide layout, as in mstate::msprep(): entry time and status per state
wide <- data.frame(id = 1:3,
                   ill.t = c(2, 5, 4), ill.s = c(1, 0, 1),
                   dth.t = c(6, 5, 4.5), dth.s = c(1, 0, 1))
msprep2(wide, times = c(ill = "ill.t", dead = "dth.t"),
        status = c(ill = "ill.s", dead = "dth.s"), initial = "healthy",
        states = c("healthy", "ill", "dead"), absorbing = "dead")$panel
#> # A tibble: 19 × 3
#>       id  time state  
#>    <int> <int> <fct>  
#>  1     1     0 healthy
#>  2     1     1 healthy
#>  3     1     2 ill    
#>  4     1     3 ill    
#>  5     1     4 ill    
#>  6     1     5 ill    
#>  7     1     6 dead   
#>  8     2     0 healthy
#>  9     2     1 healthy
#> 10     2     2 healthy
#> 11     2     3 healthy
#> 12     2     4 healthy
#> 13     2     5 healthy
#> 14     3     0 healthy
#> 15     3     1 healthy
#> 16     3     2 healthy
#> 17     3     3 healthy
#> 18     3     4 ill    
#> 19     3     5 dead   
```
