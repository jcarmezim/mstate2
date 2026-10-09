# mstate2 (development version, branch `extensiones`)

Extensions beyond the methods paper, rebuilt on the current version of the
package:

* `prep2()` gains `covariates` (baseline covariates carried into `triples`;
  by default those of an `msprep2()` result), and stores in `triples` the time
  already spent in the current state (`d`) and, in `pairs`, every one-step
  move, including each subject's first move.
* `markov_test()`: tests of the first-order Markov assumption, transition by
  transition or state by state, with a Wald (Cochran's Q) or a
  likelihood-ratio (Anderson and Goodman) statistic and Holm-adjusted
  p-values.
* `P1est()`: the first-order RPE from the same counts, pooling over the
  previous state, with a subject bootstrap. `P2boot()` also keeps the
  first-order matrices of its resamples (`boot1`), which `P1est()` reuses.
* `probtrans2()` and its `plot()` method: n-step predictions with the
  structure of an `mstate::probtrans()` object.
* `compare_order()`, `divergence()` and `divergence_table()`: for how long the
  second-order predictions diverge from a first-order baseline (`P1est()` or
  an `mstate::probtrans()` object). `as_tmat()` gives the transition matrix of
  the observed moves in `mstate::transMat()` format.
* `compare2()` on a `P2boot` fit keeps the replicate curves, and
  `overlap_step()` then also tests the paired difference between the two
  curves (`diff_steps`).
* `P2reg()` and `P2reg_all()`: discrete-time cause-specific hazard regression
  (complementary log-log link) for each history, with baseline covariates and
  the time scales `s` (non-homogeneous) and `d` (semi-Markov), cluster-robust
  standard errors, and covariate-adjusted n-step predictions through
  `predict()` and `ckequations()`.

# mstate2 0.1.0

First release. Implements the methods of Najera-Zuloaga, Besalú and Gómez
Melis (2025) and reproduces their analysis of the DIVINE cohort
(`vignette("paper")`).

* `msprep2()` and `rnd()`: build the daily panel from the days spent in each
  state (as in DIVINE), from wide data with one `Surv(time, status)` per state
  (as for `mstate::msprep()`) or from an `msdata` object, rounding half days
  up and reporting every record dropped or changed.
* `prep2()`: second-order counting processes $\tilde N_{hj\ell}(s)$ and
  $\tilde Y_{hj}(s-1)$.
* `P2est()`: relative probability estimator (RPE) of the 1-step second-order
  transition probabilities, with variances and confidence intervals
  (Table 2 of the paper).
* `P2boot()`: subject-level bootstrap of the estimates, which gives
  percentile intervals for the n-step predictions of `ckequations()` and
  `compare2()` (an addition to the methods of the paper).
* `ckequations()`: n-step prediction with the extended Chapman-Kolmogorov
  relation and evolution intervals.
* `compare2()` and `overlap_step()`: comparison of the n-step curves across
  the states at the previous time and first overlap of their intervals
  (Section 6.3 and Figures 4-5), with `print()`, `summary()` and `plot()`
  methods.
* `simulate2()`: simulation from a second-order model, with staggered entry
  (Section 5).
