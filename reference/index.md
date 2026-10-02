# Package index

## Package overview

- [`mstate2`](https://jcarmezim.github.io/mstate2/reference/mstate2-package.md)
  [`mstate2-package`](https://jcarmezim.github.io/mstate2/reference/mstate2-package.md)
  : mstate2: Second-Order Markov Multistate Models

## Building the data

Build a daily panel from the days spent in each state and compute the
second-order counting processes.

- [`sojourn_to_panel()`](https://jcarmezim.github.io/mstate2/reference/sojourn_to_panel.md)
  : Convert sojourn-time data to discrete-time panel format
- [`rnd()`](https://jcarmezim.github.io/mstate2/reference/rnd.md) :
  Round half away from zero (the rounding of the paper's DIVINE
  analysis)
- [`prep2()`](https://jcarmezim.github.io/mstate2/reference/prep2.md) :
  Build second-order counting processes from panel data

## Estimation

Relative probability estimator of the 1-step second-order transition
probabilities, and subject-level bootstrap.

- [`P2est()`](https://jcarmezim.github.io/mstate2/reference/P2est.md) :
  Estimate 1-step second-order transition probabilities
- [`P2boot()`](https://jcarmezim.github.io/mstate2/reference/P2boot.md)
  : Subject-level bootstrap of the second-order estimates

## Prediction

Multi-step prediction via the extended Chapman-Kolmogorov relation.

- [`ckequations()`](https://jcarmezim.github.io/mstate2/reference/ckequations.md)
  : Extended Chapman-Kolmogorov n-step second-order transition
  probabilities

## Does the previous time matter?

Compare trajectories across the states occupied at the previous time,
with evolution or bootstrap intervals.

- [`compare2()`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
  : Compare n-step transitions across preceding states (evolution
  intervals)
- [`overlap_step()`](https://jcarmezim.github.io/mstate2/reference/overlap_step.md)
  : First step at which two evolution intervals overlap
- [`plot(`*`<msm2pred>`*`)`](https://jcarmezim.github.io/mstate2/reference/plot.msm2pred.md)
  : Plot n-step evolution curves or evolution intervals

## Simulation

- [`simulate2()`](https://jcarmezim.github.io/mstate2/reference/simulate2.md)
  : Simulate a second-order Markov multistate process
