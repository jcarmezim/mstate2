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
library(survival)                      # Surv()
# Same (id, time, state) in two results?
same <- function(a, b) isTRUE(all.equal(as.data.frame(a$panel[1:3]), as.data.frame(b$panel[1:3]),
                                        check.attributes = FALSE))

## ---- 1. The same four patients in every layout -------------------------------
## Illness-death model: healthy -> ill -> dead and healthy -> dead
##   patient 1: healthy until day 3, ill from day 3, dies on day 7
##   patient 2: healthy until day 5, dies on day 5
##   patient 3: ill from day 2, alive at day 9 (censored)
##   patient 4: healthy until day 6 (censored)

# Sojourn (as DIVINE): days spent in each state, in visiting order, and the outcome
sojourn <- data.frame(id = 1:4, t_healthy = c(3, 5, 2, 7), t_ill = c(4, 0, 8, 0),
                      dead = c(1, 1, 0, 0))
p1 <- msprep2(sojourn, durations = c(healthy = "t_healthy", ill = "t_ill"),
              outcome = c(dead = "dead"))
p1                                     # short report
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : sojourn
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 31 (time 0 - 9)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 3 types, 4 in total
#>   issues          : none
p1$panel                               # one row per patient and day
#> # A tibble: 31 × 3
#>       id  time state  
#>    <int> <int> <fct>  
#>  1     1     0 healthy
#>  2     1     1 healthy
#>  3     1     2 healthy
#>  4     1     3 ill    
#>  5     1     4 ill    
#>  6     1     5 ill    
#>  7     1     6 ill    
#>  8     1     7 dead   
#>  9     2     0 healthy
#> 10     2     1 healthy
#> # ℹ 21 more rows
p1$subjects                            # entry, exit and status of each patient
#> # A tibble: 4 × 7
#>      id first  last entry   exit     rows status  
#>   <int> <int> <int> <chr>   <chr>   <int> <chr>   
#> 1     1     0     7 healthy dead        8 absorbed
#> 2     2     0     5 healthy dead        6 absorbed
#> 3     3     0     9 healthy ill        10 censored
#> 4     4     0     6 healthy healthy     7 censored

# Wide (as mstate::msprep()): time of entry, or of censoring when the status is 0
wide <- data.frame(id = 1:4, ill_time = c(3, 5, 2, 6), ill_status = c(1, 0, 1, 0),
                   dead_time = c(7, 5, 9, 6), dead_status = c(1, 1, 0, 0),
                   age = c(60, 72, 55, 49))
p2 <- msprep2(wide, states = list(healthy = NULL,     # no columns: the initial state
                                  ill     = Surv(ill_time, ill_status),
                                  dead    = Surv(dead_time, dead_status)),
              trans = c("healthy -> ill -> dead", "healthy -> dead"))
p2$panel                               # `age` is kept as a covariate
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

# Events: one row per change of state, dates as text and an end of follow-up
events <- data.frame(id    = c(1, 1, 1, 2, 2, 3, 3, 4),
                     date  = c("01/03/2020", "04/03/2020", "08/03/2020", "01/03/2020",
                               "06/03/2020", "01/03/2020", "03/03/2020", "01/03/2020"),
                     state = c("healthy", "ill", "dead", "healthy", "dead",
                               "healthy", "ill", "healthy"),
                     end   = c(NA, NA, NA, NA, NA, "10/03/2020", "10/03/2020", "07/03/2020"))
p3 <- msprep2(events, time = "date", state = "state", end = "end", absorbing = "dead")

# Codes instead of labels: relabel them with `recode`
coded <- transform(events, state = c(healthy = "H", ill = "I", dead = "D")[state])
p3b <- msprep2(coded, time = "date", state = "state", end = "end",
               recode = c(H = "healthy", I = "ill", D = "dead"), absorbing = "dead",
               states = c("healthy", "ill", "dead"))

# msdata: the wide data already prepared by mstate::msprep()
if (requireNamespace("mstate", quietly = TRUE)) {
  tmat <- mstate::transMat(x = list(c(2, 3), 3, c()), names = c("healthy", "ill", "dead"))
  ms <- mstate::msprep(time = c(NA, "ill_time", "dead_time"),
                       status = c(NA, "ill_status", "dead_status"), data = wide, trans = tmat)
  p4 <- msprep2(ms)                    # recognised by its class; `trans` from its matrix
  print(same(p1, p4))
}
#> [1] TRUE

# All the layouts give the same panel
c(wide = same(p1, p2), events = same(p1, p3), recoded = same(p1, p3b))
#>    wide  events recoded 
#>    TRUE    TRUE    TRUE 

## ---- 2. Other ways of writing the wide layout --------------------------------
# Conventional names (<state>_time, <state>_status, inistat): no arguments needed
conv <- data.frame(wide, inistat = "healthy")
msprep2(conv, absorbing = "dead")
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : wide
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 31 (time 0 - 9)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 3 types, 4 in total
#>   issues          : none

# Column names in quotes (e.g. from a Shiny app), in any of these forms: c(time =, status =),
# c("time column" = "status column"), c("ill_status", "ill_time") (the 0/1 column is the
# status) or just the time column (visited when its time is recorded); Surv() can be mixed in
# and contain expressions, e.g. Surv(dead_time, dead_status == 1)
q <- msprep2(wide, states = list(healthy = NULL,
                                 ill     = c(time = "ill_time", status = "ill_status"),
                                 dead    = c("dead_time" = "dead_status")))
same(p2, q)
#> [1] TRUE

# Or with the arguments `times` and `status`
q <- msprep2(wide, times = c(ill = "ill_time", dead = "dead_time"),
             status = c(ill = "ill_status", dead = "dead_status"), initial = "healthy",
             states = c("healthy", "ill", "dead"), absorbing = "dead")
same(p2, q)
#> [1] TRUE

# The transitions as text, as a list of destinations, as a matrix (mstate::transMat(),
# logical or 0/1) or not given (the absorbing states are then the states nobody leaves)
tr <- msprep2(wide, states = list(healthy = NULL, ill = "ill_time", dead = "dead_time"),
              trans = list(healthy = c("ill", "dead"), ill = "dead", dead = NULL))
#> Warning: 2 record(s) were dropped or changed while building the panel (same_unit: 2); see the `issues` table of the result.
tr$trans                               # `dead` has no exit: absorbing
#>         healthy   ill  dead
#> healthy   FALSE  TRUE  TRUE
#> ill       FALSE FALSE  TRUE
#> dead      FALSE FALSE FALSE

## ---- 3. What msprep2() reports ----------------------------------------------
## Every record dropped or changed is listed in `issues`, with one warning.

# rounded_to_zero: patient 2 spent 0.4 days ill, which rounds to 0 days
bad <- sojourn
bad$t_ill[2] <- 0.4
msprep2(bad, durations = c(healthy = "t_healthy", ill = "t_ill"),
        outcome = c(dead = "dead"))$issues
#> Warning: 1 record(s) were dropped or changed while building the panel (rounded_to_zero: 1); see the `issues` table of the result.
#> # A tibble: 1 × 5
#>      id issue           state time  detail                        
#>   <int> <chr>           <chr> <chr> <chr>                         
#> 1     2 rounded_to_zero ill   NA    duration 0.4 rounds to 0 units
rnd(c(0.4, 0.5, 1.5))                  # halves away from zero ...
#> [1] 0 1 2
round(c(0.4, 0.5, 1.5))                # ... not to the even integer
#> [1] 0 0 2

# same_unit: in weeks (unit = 7; with dates, also unit = "week"), patients 1 and 3
# become ill in the week of admission, so `healthy` occupies no week
msprep2(wide, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                            dead = Surv(dead_time, dead_status)), unit = 7)$issues
#> Warning: 2 record(s) were dropped or changed while building the panel (same_unit: 2); see the `issues` table of the result.
#> # A tibble: 2 × 5
#>      id issue     state   time  detail                        
#>   <int> <chr>     <chr>   <chr> <chr>                         
#> 1     1 same_unit healthy 0     same time unit as ill (unit 0)
#> 2     3 same_unit healthy 0     same time unit as ill (unit 0)

# missing: status 1 without a time (patient 1), and a status code 9 that Surv()
# turns into NA (patient 2)
bad <- wide
bad$ill_time[1] <- NA
bad$dead_status[2] <- 9
msprep2(bad, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                           dead = Surv(dead_time, dead_status)))$issues
#> Warning: Invalid status value, converted to NA
#> Warning: 2 record(s) were dropped or changed while building the panel (missing: 2); see the `issues` table of the result.
#> # A tibble: 2 × 5
#>      id issue   state time  detail                                              
#>   <int> <chr>   <chr> <chr> <chr>                                               
#> 1     1 missing ill   NA    visited (status 1) but without a time               
#> 2     2 missing dead  5     time recorded but status missing: not taken as a vi…

# A status coded 1/2 is read as in survival::Surv() (1 = censored, 2 = event)
w12 <- transform(wide, dead_status = dead_status + 1)
same(p2, msprep2(w12, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                    dead = Surv(dead_time, dead_status))))
#> Warning: A status coded 1/2 is read as in survival::Surv(): 1 = censored, 2 = event.
#> [1] TRUE

# A status written as text is not guessed: say which value is the event
wtxt <- transform(wide, dead_status = ifelse(dead_status == 1, "yes", "no"))
try(msprep2(wtxt, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                dead = Surv(dead_time, dead_status))))
#> Error : The status of Surv() must be 0/1, 1/2 or TRUE/FALSE, not text (values: yes, no). Say which value is the event, e.g. Surv(time, status == "yes").
msprep2(wtxt, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                            dead = Surv(dead_time, dead_status == "yes")))
#> <msm2prep>  discrete-time panel ready for prep2()
#>   layout          : wide
#>   subjects        : 4 (2 absorbed, 2 censored)
#>   panel rows      : 31 (time 0 - 9)
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   transitions     : 3 types, 4 in total
#>   issues          : none

# not_allowed: healthy -> dead (patient 2) is not in `trans`: kept and listed ...
msprep2(wide, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                            dead = Surv(dead_time, dead_status)),
        trans = "healthy -> ill -> dead")$issues
#> Warning: 1 record(s) were dropped or changed while building the panel (not_allowed: 1); see the `issues` table of the result.
#> # A tibble: 1 × 5
#>      id issue       state time  detail                   
#>   <int> <chr>       <chr> <chr> <chr>                    
#> 1     2 not_allowed dead  NA    healthy -> dead at unit 5
# ... or an error with check = "error"
try(msprep2(wide, states = list(healthy = NULL, ill = Surv(ill_time, ill_status),
                                dead = Surv(dead_time, dead_status)),
            trans = "healthy -> ill -> dead", check = "error"))
#> Error : 1 transition(s) not allowed by `trans`: healthy -> dead.

# Events: a record without a date (missing), one before admission (before_start),
# one after the end of follow-up (after_end), one after death (after_absorbing),
# two states on the day of admission (same_unit) and a patient with no usable
# record (no_data)
bad <- rbind(events,
             data.frame(id    = c(1, 2, 4, 1, 5),
                        date  = c(NA, "20/02/2020", "09/03/2020", "12/03/2020", NA),
                        state = "ill",
                        end   = c(NA, NA, "07/03/2020", NA, NA)))
bad$date[bad$id == 3 & bad$state == "ill"] <- "01/03/2020"
y <- msprep2(bad, time = "date", state = "state", end = "end", start = "01/03/2020",
             absorbing = "dead")
#> Warning: 7 record(s) were dropped or changed while building the panel (after_absorbing: 1, after_end: 1, before_start: 1, missing: 2, no_data: 1, same_unit: 1); see the `issues` table of the result.
y$issues
#> # A tibble: 7 × 5
#>      id issue           state   time       detail                               
#>   <dbl> <chr>           <chr>   <chr>      <chr>                                
#> 1     1 missing         ill     NA         missing id, time or state            
#> 2     5 missing         ill     NA         missing id, time or state            
#> 3     2 before_start    ill     2020-02-20 before the origin                    
#> 4     4 after_end       ill     2020-03-09 after the end of follow-up           
#> 5     1 after_absorbing ill     2020-03-12 after the entry into an absorbing st…
#> 6     3 same_unit       healthy 2020-03-01 same time unit as ill (unit 0)       
#> 7     5 no_data         NA      NA         no usable record                     
summary(y)                             # transitions, follow-up and issues by type
#> <msm2prep summary>
#> 
#> Observed transitions (from rows to columns):
#>          to
#> from      healthy ill dead
#>   healthy       0   1    1
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
#> # A tibble: 6 × 2
#>   issue               n
#>   <chr>           <int>
#> 1 after_absorbing     1
#> 2 after_end           1
#> 3 before_start        1
#> 4 missing             2
#> 5 no_data             1
#> 6 same_unit           1

## ---- 4. The mstate::msprep() workflow of ebmt3 in one call -----------------
if (requireNamespace("mstate", quietly = TRUE)) {
  data(ebmt3, package = "mstate")
  bmt <- msprep2(ebmt3,
                 states = list(Tx       = NULL,
                               PR       = Surv(prtime, prstat),
                               RelDeath = Surv(rfstime, rfsstat)),
                 trans  = c("Tx -> PR -> RelDeath", "Tx -> RelDeath"),
                 unit   = 30.4375)     # times in days, panel in months
  print(bmt)
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

## ---- 5. The result goes straight into prep2() --------------------------------
prep2(p2)
#> <msm2data>  second-order counting processes
#>   subjects        : 4
#>   observed triples: 23
#>   time range      : 0 - 9
#>   states (3)      : healthy, ill, dead
#>   absorbing       : dead
#>   distinct (h,j)  : 3
```
