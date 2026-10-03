# mstate2 0.1.0

First release. Implements the methods of Najera-Zuloaga, Besalú and Gómez
Melis (2025) and reproduces their analysis of the DIVINE cohort
(`vignette("paper")`).

* `sojourn_to_panel()` and `rnd()`: build the daily panel from the days spent
  in each state, rounding half days up.
* `prep2()`: second-order counting processes $\tilde N_{hj\ell}(s)$ and
  $\tilde Y_{hj}(s-1)$.
* `P2est()`: relative probability estimator (RPE) of the 1-step second-order
  transition probabilities, with variances and confidence intervals
  (Table 2 of the paper).
* `ckequations()`: n-step prediction with the extended Chapman-Kolmogorov
  relation and evolution intervals.
* `compare2()` and `overlap_step()`: comparison of the n-step curves across
  the states at the previous time and first overlap of their intervals
  (Section 6.3 and Figures 4-5), with `print()`, `summary()` and `plot()`
  methods.
* `simulate2()`: simulation from a second-order model, with staggered entry
  (Section 5).
