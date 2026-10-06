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
    [`msprep2`](https://jcarmezim.github.io/mstate2/reference/msprep2.md)
    turns sojourn data (as DIVINE), wide data with one
    `Surv(time, status)` per state (as for
    [`mstate::msprep()`](https://rdrr.io/pkg/mstate/man/msprep.html)) or
    an mstate `msdata` object into a daily panel and reports every
    record it drops or changes
    ([`rnd`](https://jcarmezim.github.io/mstate2/reference/rnd.md)
    discretises the times);
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

## Implementation

See
[`vignette("mstate2")`](https://jcarmezim.github.io/mstate2/articles/mstate2.md)
for a guided analysis,
[`vignette("paper")`](https://jcarmezim.github.io/mstate2/articles/paper.md)
for the reproduction of the methods paper and
[`vignette("reference")`](https://jcarmezim.github.io/mstate2/articles/reference.md)
for the complete function reference.

## See also

Useful links:

- <https://jcarmezim.github.io/mstate2/>

- <https://github.com/jcarmezim/mstate2>

- Report bugs at <https://github.com/jcarmezim/mstate2/issues>

## Author

**Maintainer**: Jo\<U+00E3\>o Carmezim <joaocarmezimcorreia@gmail.com>

Authors:

- Josu Najera-Zuloaga <josu.najera@ehu.eus> (creator of the package)

- Mireia Besal\<U+00FA\>

- Cristian Teb\<U+00E9\>

- Guadalupe G\<U+00F3\>mez Melis
