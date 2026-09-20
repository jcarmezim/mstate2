# mstate2 0.1.0

* Initial CRAN submission: `prep2()`, `P2est()` (RPE/CPE with variances and
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
