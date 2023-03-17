#' @title Measures for Benchmarking Survival Models
#'
#' @description A convenient function that returns a table with measures that
#' can be used for benchmarking survival learners.
#'
#' @details The following measures are provided via [mlr3proba]:
#' - Discrimination measures (input is ranking prediction - `crank`)
#'    - [Harrell's C-index][mlr3proba::MeasureSurvCindex] (doesn't account
#'    for censoring)
#'    - [Uno's C-index][mlr3proba::MeasureSurvCindex] (estimates censoring
#'    distribution using Kaplan-Meier)
#' - Calibration measures (input is survival predictions - `distr`)
#'    - [Integrated Brier Score][mlr3proba::MeasureSurvGraf] (with
#'    *proper = TRUE, method = 2*, estimates the censoring distribution using
#'    Kaplan-Meier, is also a measure of discrimination)
#'    - [Right-Censored Log-Loss][mlr3proba::MeasureSurvRCLL] (proper loss score)
#'    - [D-Calibration][mlr3proba::MeasureSurvDCalibration] (with *B = 10,
#'    chisq = FALSE*, returning the measure and not the p-value)
#'
#'  The Integrated Brier Score and Right-Censored Log-Loss are also offered in
#'  Explained Residual Variation (ERV) versions, where the resulting score can
#'  be interpreted as the percentage increase in performance over the baseline
#'  Kaplan-Meier (KM) model.
#'  So the best possible value for ERV-standardized scores is 1 and negative
#'  values correspond to worse performance than KM.
#'
#'  Any of the above measures **can be used for model evaluation** (but depends
#'  on several other factors if they *should* be used, e.g. censoring
#'  distribution, dataset size, proportion of events in the test set, etc.).
#'  For model **optimization/tuning**, we suggest the Right-Censored Log-Loss
#'  for survival predictions and Uno's C-index for risk predictions
#'  (e.g. linear predictors).
#'
#' @return a [data.table][data.table::data.table] with the suggested survival
#' measures used for benchmarking purposes
#'
#' @examples
#' library(mlr3proba)
#' ms = bench_msrs()
#'
#' # Available measures ids
#' ms$id
#'
#' # Get specific mlr3 measure objects
#' ms[id %in% c('ibrier','harrell_c')]$measure
#'
#' @export
bench_msrs = function() {
  # Harrell's C-index
  harrell_c = msr('surv.cindex', label = 'HarrellC', id = 'harrell_c')

  # Uno's C-index
  uno_c = msr('surv.cindex', weight_meth = 'G2', label = 'UnoC', id = 'uno_c')

  # Integrated Brier Score (with proper = TRUE, method = 2)
  ibrier = msr('surv.graf', proper = TRUE, label = 'IBrier', id = 'ibrier')

  # Integrated Brier Score (ERV version)
  ibrier_erv = msr('surv.graf', proper = TRUE, ERV = TRUE, label = 'IBrier-ERV',
    id = 'ibrier_erv', minimize = FALSE, range = c(-Inf, 1))

  # Right-Censored Log-Loss
  rcll = msr('surv.rcll', id = 'rcll', label = 'RCLL')

  # Right-Censored Log-Loss (ERV version)
  rcll_erv = msr('surv.rcll', ERV = TRUE, id = 'rcll_erv', label = 'RCLL-ERV',
    minimize = FALSE, range = c(-Inf, 1))

  # D-calibration
  dcal = msr('surv.dcalib', B = 10, chisq = FALSE, label = 'Dcalibration',
    id = 'dcal')

  data.table(
    id = c('harrell_c', 'uno_c', 'ibrier', 'ibrier_erv', 'rcll', 'rcll_erv', 'dcal'),
    measure = list(harrell_c, uno_c, ibrier, ibrier_erv, rcll, rcll_erv, dcal)
  )
}
