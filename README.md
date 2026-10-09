
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# mstate2 <img src="man/figures/logo.png" align="right" width="200"/>

[![R-CMD-check](https://github.com/jcarmezim/mstate2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jcarmezim/mstate2/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/jcarmezim/mstate2/graph/badge.svg?token=ptcXJwFrji)](https://codecov.io/gh/jcarmezim/mstate2)

**mstate2** implements **second-order Markov multistate models** in R.

The package implements the methods of Najera-Zuloaga, Besalú and Gómez
Melis (2025): nonparametric estimation of the second-order transition
probabilities, multi-step prediction with the extended
Chapman–Kolmogorov relation and evolution intervals, and the comparison
of trajectories according to the state at the previous time.

## Installation

mstate2 is not yet on CRAN. Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("jcarmezim/mstate2", build_vignettes = TRUE)
```

## Example: the DIVINE cohort

The example uses the DIVINE cohort of 2076 patients hospitalised with
COVID-19 (the real-data illustration of the methods paper). The data are
not distributed with mstate2; download `MSM_Data.RData` from the [DIVINE
repository](https://github.com/bruigtp/DIVINE/tree/main/msm_data). Each
row is a patient with the days spent in each state and how follow-up
ended.

``` r
library(mstate2)
load("MSM_Data.RData")     # data frame MSM
```

**1. From sojourn times to a daily panel, and the second-order counts.**

``` r
x <- msprep2(MSM,
             durations = c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov"),
             outcome   = c(Disch = "disch.s", Death = "death.s"),
             states    = c("NSP", "SP", "Recov", "NIMV", "IMV", "Disch", "Death"))
d <- prep2(x)
d
#> <msm2data>  second-order counting processes
#>   subjects        : 2076
#>   observed triples: 23584
#>   time range      : 0 - 138
#>   states (7)      : NSP, SP, Recov, NIMV, IMV, Disch, Death
#>   absorbing       : Disch, Death
#>   distinct (h,j)  : 12
#>   covariates      : inistat
```

**2. Estimation.** The probability of moving from severe pneumonia (SP) to non-invasive (NIMV) or invasive (IMV) ventilation depends strongly on the state at the previous time:

``` r
fit <- P2est(d)
fit$estimate |>
  dplyr::filter(j == "SP" & l %in% c("NIMV", "IMV")) |>
  dplyr::select("h", "j", "l", "p", "lower", "upper")
#> # A tibble: 4 × 6
#>   h     j     l          p  lower  upper
#>   <fct> <fct> <fct>  <dbl>  <dbl>  <dbl>
#> 1 NSP   SP    NIMV  0.224  0.184  0.264 
#> 2 NSP   SP    IMV   0.151  0.116  0.185 
#> 3 SP    SP    NIMV  0.0255 0.0195 0.0315
#> 4 SP    SP    IMV   0.0184 0.0133 0.0235
```

**3. For how long does the previous time matter?** `compare2()` predicts the probability of being in NIMV $n$ days later, for patients in SP whose state at the previous time was NSP or SP, with the evolution intervals of the paper. The intervals first overlap on day 5 (dotted line): the state at the previous time matters for the first 4 days.

``` r
cmp <- compare2(fit, h = c("NSP", "SP"), j = "SP", l = "NIMV")
summary(cmp)
#> Trajectory comparison (RPE, 95% evolution intervals)
#>   target 'NIMV' via current state 'SP'; preceding: NSP, SP
#>   intervals first overlap at step 5 (time s = 6); significant for the first 4 step(s).
plot(cmp, dualaxis = FALSE, xlab = "days ahead (n)")
```

<img src="man/figures/README-compare-1.png" width="75%" />

## Main functions

| Task                                             | Functions                                                               |
|--------------------------------------------------|-------------------------------------------------------------------------|
| Prepare the data                                 | `msprep2()`, `rnd()`, `prep2()`                                         |
| Estimate                                         | `P2est()` (RPE, with Wald or logit intervals)                           |
| Predict                                          | `ckequations()` (extended Chapman-Kolmogorov, with evolution intervals) |
| Does the previous time matter, and for how long? | `compare2()`, `overlap_step()`, with `summary()` and `plot()`           |

## Documentation

- `vignette("mstate2")`: a complete guided analysis, function by
  function.
- `vignette("reference")`: the full reference (signature, arguments,
  value, algorithm and error conditions of every function).
- The package website: <https://jcarmezim.github.io/mstate2/>.

## Methods

The estimators, their variances and the extended Chapman-Kolmogorov
relation are described in:

- Najera-Zuloaga J, Besalú M, Gómez Melis G (2025). *Second-Order Markov
  Multistate Models: Nonparametric Estimation and Inference*. Manuscript
  submitted for publication.
- Besalú M, Gómez Melis G (2024). Second order Markov multistate models.
  *SORT*, 48(2), 209–234.

## Citation

``` r
citation("mstate2")
To cite mstate2 in publications, please cite the methods paper. Citing
the package itself (below) is optional.

  Najera-Zuloaga J, Besalú M, Gómez Melis G (2025). "Second-Order
  Markov Multistate Models: Nonparametric Estimation and Inference."
  Manuscript submitted for publication.

  Najera-Zuloaga J, Carmezim J, Besalú M, Tebé C, Gómez Melis G (2026).
  mstate2: Second-Order Markov Multistate Models. R package version
  0.1.0. https://github.com/jcarmezim/mstate2

To see these entries in BibTeX format, use 'print(<citation>,
bibtex=TRUE)', 'toBibtex(.)', or set
'options(citation.bibtex.max=999)'.
```

## Authors

- Josu Najera-Zuloaga
- João Carmezim (maintainer)
- Mireia Besalú
- Cristian Tebé
- Guadalupe Gómez Melis

## Getting help

Please report bugs, with a minimal reproducible example, on [GitHub Issues](https://github.com/jcarmezim/mstate2/issues).

## License

GPL (\>= 3)
