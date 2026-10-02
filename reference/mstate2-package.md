# mstate2: Second-Order Markov Multistate Models

Nonparametric estimation and inference for second-order Markov
multistate models in discrete time, in which the next state depends on
the state occupied at the current time *and* at the previous time.
Complements the first-order mstate framework and follows its
conventions.

## Details

The package implements the methods of Najera-Zuloaga, Besalu and Gomez
Melis (2025), so that their analysis of the DIVINE cohort can be
reproduced, and adds subject-level bootstrap intervals for the
multi-step predictions.

## Workflow

1.  **Data.**
    [`sojourn_to_panel`](https://jcarmezim.github.io/mstate2/reference/sojourn_to_panel.md)
    builds a daily panel from sojourn times
    ([`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md)
    discretises them);
    [`prep2`](https://jcarmezim.github.io/mstate2/reference/prep2.md)
    computes the second-order counting processes.

2.  **Estimation.**
    [`P2est`](https://jcarmezim.github.io/mstate2/reference/P2est.md)
    (relative probability estimator with variances and confidence
    intervals).

3.  **Prediction.**
    [`ckequations`](https://jcarmezim.github.io/mstate2/reference/ckequations.md)
    (extended Chapman-Kolmogorov relation, with evolution intervals).

4.  **Does the previous time matter, and for how long?**
    [`compare2`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
    and
    [`overlap_step`](https://jcarmezim.github.io/mstate2/reference/overlap_step.md).

5.  **Bootstrap intervals.**
    [`P2boot`](https://jcarmezim.github.io/mstate2/reference/P2boot.md)
    resamples patients;
    [`ckequations()`](https://jcarmezim.github.io/mstate2/reference/ckequations.md)
    and
    [`compare2()`](https://jcarmezim.github.io/mstate2/reference/compare2.md)
    then report percentile bootstrap intervals.

6.  **Simulation.**
    [`simulate2`](https://jcarmezim.github.io/mstate2/reference/simulate2.md)
    (the simulation design of the paper).

See
[`vignette("mstate2")`](https://jcarmezim.github.io/mstate2/articles/mstate2.md)
for a guided analysis and
[`vignette("reference")`](https://jcarmezim.github.io/mstate2/articles/reference.md)
for the complete function reference.

## See also

Useful links:

- <https://jcarmezim.github.io/mstate2/>

- <https://github.com/jcarmezim/mstate2>

- Report bugs at <https://github.com/jcarmezim/mstate2/issues>

## Author

**Maintainer**: João Carmezim <joaocarmezimcorreia@gmail.com>

Authors:

- Josu Najera-Zuloaga <josu.najera@ehu.eus>

- Mireia Besalú

- Cristian Tebé

- Guadalupe Gómez Melis
