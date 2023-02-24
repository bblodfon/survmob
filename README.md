# survmob

<!-- badges: start -->

<!-- badges: end -->

`survmob` is an R library that makes it easier to benchmark survival models across multiomics datasets.
It builds upon [mlr3proba](https://github.com/mlr-org/mlr3proba/) and other packages from the [mlr3](https://github.com/mlr-org/) ML ecosystem.
Key [R6](https://github.com/r-lib/R6/) classes:

- Preprocess classes
- Survival learners and hyperparameter spaces class
- Ensemble Feature Selection (eFS) class (uses Random Survival Forest learners and RFE algorithm for finding best feature subsets)

## Installation

Install the latest development version of `survmob` using:
```r
remotes::install_github('bblodfon/survmob')
```
