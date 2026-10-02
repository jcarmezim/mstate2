# mstate2 0.2.0

New functionality following the PhD research plan (Objectives 1-3):

* **Testing the Markov assumption.** `markov_test()` gains
  `method = "lrt"`, the likelihood-ratio test of Markov order (Anderson and
  Goodman, 1957) computed from the RPE counts, and `by = "state"`, a joint
  test of all destinations from each current state. Unlike the Wald version,
  the LRT does not drop histories with zero events.
* **First-order baseline from the same data.** New `P1est()`: the first-order
  discrete-time RPE obtained by pooling the second-order counts over the
  preceding state (the estimate a second-order analysis reduces to when the
  first-order property holds), with optional subject-level bootstrap.
* **Systematic first- versus second-order comparison.** `compare_order()` now
  also accepts a `P1est()` fit as baseline (no `mstate` fit needed), and the
  new `divergence_table()` applies it to every current state and target,
  reporting the memory horizon and the largest difference for each history.
* **Bootstrap intervals for multi-step predictions.** New `P2boot()`
  resamples subjects and recomputes all estimates (via per-subject count
  tables, so thousands of replicates are cheap). `ckequations()`,
  `compare2()`, `overlap_step()` and `compare_order()` then use percentile
  bootstrap intervals, which have close to nominal coverage, instead of the
  (conservative) evolution intervals.
* **Flexible baseline and semi-Markov extension in `P2reg()`.** `prep2()` now
  stores, for every triple, the time already spent in the current state
  (`d`). `P2reg()` formulas may use `s` (time; non-homogeneous baseline, e.g.
  `~ age + splines::ns(s, 3)`) and `d` (duration; exploratory semi-Markov
  baseline, e.g. `~ age + log(d)`). Formulas are evaluated in their own
  environment, and `P2reg()` no longer requires baseline covariates.
* **Covariate-adjusted multi-step prediction.** New `P2reg_all()` fits the
  regression for every history and `predict()` returns the covariate-adjusted
  tensor for any profile, ready for `ckequations()` or `simulate2()`.
  `P2reg_all(~1)` reproduces `P2est()` exactly.
  Terms that cannot be estimated in a history (a covariate constant in its
  risk set) are reported by `print()` and act as 0 in predictions.
* **mstate-compatible data preparation and output.** New `msprep2()` builds
  the panel from the wide format of `mstate::msprep()` (same arguments); it
  reproduces `msprep()` followed by `from_msdata()` exactly. New
  `probtrans2()` returns n-step predictions with the structure of an `mstate`
  `probtrans` object, with a `plot()` method using mstate's plot types.
* **Visualisation.** `plot()` for `P2reg` fits (forest plot of hazard ratios)
  and for `probtrans2` objects.
* Bug fix: `as_tmat()` missed transitions that only occur as a subject's
  first move (they have no preceding state, hence no triple). `prep2()` now
  records every one-step move in `$pairs`, and `as_tmat()` uses it. It no
  longer needs the `mstate` package (its output is identical to
  `mstate::transMat()`).
* The conditional probability estimator (CPE) has been removed: the
  relative probability estimator (RPE) is more efficient (Najera-Zuloaga,
  Besalú and Gómez Melis, 2025) and is now the only estimator. `P2est()`
  loses its `estimator` argument, the `t.hj` column of its estimate table
  and the `t_hj` column of `summary()` for `msm2data` objects.
* Methodological review of the new functionality, with references in the
  help pages and a table of methodological choices in `vignette("mstate2")`:
  * `P2reg()` models by default only the moves `l != j` (the discrete-time
    cause-specific hazards); staying is their complement and can still be
    requested with `l = j`. The documentation states precisely when `exp(coef)`
    is a hazard ratio.
  * `markov_test()` adds Holm-adjusted p-values (`p.adj`; argument
    `p.adjust`).
  * With bootstrap replicates (`compare2()` on a `P2boot` fit, or
    `compare_order()` with `P1est()` of the same fit), `overlap_step()`,
    `divergence()` and `divergence_table()` also test the paired difference
    between curves directly (`diff_steps`), since non-overlap of two
    intervals is a conservative criterion.
  * Subsetting a `summary()` of a `P2reg` fit keeps its header attributes.
  * `sojourn_to_panel()`, `msprep2()` and `from_msdata()` warn when
    discretisation drops a visit: a positive sojourn that rounds to 0 time
    units, or a state entered on the same rounded time as the next one. The
    transitions into and out of such a visit would otherwise be lost
    silently. (No visit is lost in DIVINE: all durations are whole or half
    days, and `rnd()` rounds 0.5 up to 1.)
* Both vignettes now use the DIVINE cohort, the real-data illustration of the
  methods paper. Since the data are not distributed with the package, the
  vignettes are precomputed from `vignettes/*.Rmd.orig` with
  `vignettes/precompile.R`.
* Package infrastructure: João Carmezim is the new maintainer. The README
  has been rewritten around a worked example on the DIVINE cohort, the citation now refers to the
  installed package version, and the package-level help page (`?mstate2`)
  describes the whole workflow.
* Covariate names that clash with internal columns (`id`, `time`, `state`,
  `h`, `j`, `l`, `s`, `d`) are now rejected by `prep2()`.

# mstate2 0.1.0

* Initial CRAN submission: `prep2()`, `P2est()` (RPE with variances and
  confidence intervals), `ckequations()`, `compare2()`/`overlap_step()`,
  `simulate2()`, `sojourn_to_panel()`, `as_tmat()`, `from_msdata()`.
* `prep2(covariates = )` and `from_msdata(keep = )` carry baseline covariates
  through to a new individual-level `$triples` table.
* `markov_test()`: a formal, per-transition heterogeneity test of the
  first-order Markov assumption, built on `P2est()`'s own asymptotics.
* `P2reg()`: a covariate-adjusted discrete-time cause-specific hazard
  regression (complementary log-log by default, the discrete-time analogue of
  Cox) for the 1-step second-order transition probabilities, with
  `predict()`, `summary()` and `print()` methods.
* `compare_order()`/`divergence()`: compare a second-order prediction against
  a first-order `mstate::probtrans()` fit and quantify how many steps the
  difference persists.
* `tests/DIVINE_reproduction.R`: reproduces Table 2 and the Section 6.3
  evolution intervals of Najera-Zuloaga, Besalu and Gomez Melis (2025) from
  the real, published DIVINE cohort data (not bundled; the script skips
  quietly, with no network access, when the data file isn't found locally).
* Every exported function now ships a runnable `@examples` block.
* Bug fix: `compare_order()` could silently reuse (flat-line) the last
  available row of `pt` when `nsteps` requested more steps than `pt` actually
  covered, instead of erroring; it now checks the full requested horizon
  against `pt`'s time range upfront.
* Bug fix: `simulate2()` no longer burns through the full `maxT` loop for
  individuals whose next state has zero probability under `tensor`/`first`
  (e.g. an underspecified row); their path now ends immediately instead of
  being retried every step until `maxT`.
* Test coverage raised from 90% to 98% by adding tests for previously
  untested error paths (unresolved states, malformed `pt`/tensor input,
  missing columns) and `print`/`summary`/`plot` methods.
* Bug fix: with staggered entry (`simulate2(entry = )`), an individual who
  entered at global step `s` was also advanced in that same step, so its
  entry state and first move were both recorded at time `s`. Newly entered
  individuals now make their first move at `s + 1`, giving one row per
  subject per time step.
* `print.markov_test()` no longer errors when printing a column subset of a
  `markov_test` object.
* The vignette has been rewritten as a complete, function-by-function guide:
  every exported function is explained (arguments, internal mechanics,
  output, interpretation, caveats) and run on a synthetic hospital cohort
  with known true parameters, so each step can be checked against the truth.
  It also shows how to apply the workflow to the DIVINE cohort.
* Bug fix: evolution-interval bounds from `ckequations(bounds = TRUE)` and
  `compare2()` are now clipped to [0, 1]. The rows of the confidence-limit
  tensors do not sum to 1, so an upper bound feeding an absorbing state could
  exceed 1. First-overlap results (`overlap_step()`) are unaffected.
* New vignette `vignette("reference")`: a complete reference with the
  definition, arguments, value, algorithm, errors/warnings and an example for
  every exported function, S3 method, object class and internal helper.
