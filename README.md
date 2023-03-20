# survmob

<!-- badges: start -->

<!-- badges: end -->

`survmob` is an R library that makes it easier to benchmark survival models across multiomics datasets.
It builds upon [mlr3proba](https://github.com/mlr-org/mlr3proba/) and other packages from the [mlr3](https://github.com/mlr-org/) ML ecosystem using [R6](https://github.com/r-lib/R6/) classes.

Key classes and functions implemented:
- Preprocessing - [link](https://github.com/bblodfon/survmob/blob/main/R/preprocess.R)
- Survival learners and hyperparameter spaces class - [link](https://github.com/bblodfon/survmob/blob/main/R/SurvLPS.R)
- Survival measures - [link](https://github.com/bblodfon/survmob/blob/main/R/measures.R)
- Ensemble Feature Selection (eFS) class (uses Random Survival Forest learners and RFE algorithm for finding best feature subsets) - [link](https://github.com/bblodfon/survmob/blob/main/R/eFS.R)
- Multi-omics benchmarking - [link](https://github.com/bblodfon/survmob/blob/main/R/MOBenchmark.R)

## Installation

Install the latest development version of `survmob` using:
```r
remotes::install_github('bblodfon/survmob')
```
