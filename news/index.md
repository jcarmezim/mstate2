# Changelog

## mstate2 0.1.0

First release. Implements the methods of Najera-Zuloaga, Besalú and
Gómez Melis (2025) and reproduces their analysis of the DIVINE cohort
([`vignette("paper")`](https://jcarmezim.github.io/mstate2/articles/paper.md)).

- [`msprep2()`](https://jcarmezim.github.io/mstate2/reference/msprep2.md)
  and [`rnd()`](https://jcarmezim.github.io/mstate2/reference/rnd.md):
  build the daily panel from the days spent in each state (as in
  DIVINE), from wide data with one `Surv(time, status)` per state (as
  for [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html))
  or from an `msdata` object, rounding half days up and reporting every
  record dropped or changed.
- [`prep2()`](https://jcarmezim.github.io/mstate2/reference/prep2.md):
  second-order counting processes $`\tilde N_{hj\ell}(s)`$ and
  $`\tilde Y_{hj}(s-1)`$.
- [`P2est()`](https://jcarmezim.github.io/mstate2/reference/P2est.md):
  relative probability estimator (RPE) of the 1-step second-order
  transition probabilities, with variances and confidence intervals
  (Table 2 of the paper).
- [`P2boot()`](https://jcarmezim.github.io/mstate2/reference/P2boot.md):
  subject-level bootstrap of the estimates, which gives percentile
  intervals for the n-step predictions of
  [`ckequations()`](https://jcarmezim.github.io/mstate2/reference/ckequations.md)
  and
  [`compare2()`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
  (an addition to the methods of the paper).
- [`ckequations()`](https://jcarmezim.github.io/mstate2/reference/ckequations.md):
  n-step prediction with the extended Chapman-Kolmogorov relation and
  evolution intervals.
- [`compare2()`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
  and
  [`overlap_step()`](https://jcarmezim.github.io/mstate2/reference/overlap_step.md):
  comparison of the n-step curves across the states at the previous time
  and first overlap of their intervals (Section 6.3 and Figures 4-5),
  with [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
- [`simulate2()`](https://jcarmezim.github.io/mstate2/reference/simulate2.md):
  simulation from a second-order model, with staggered entry (Section
  5).
