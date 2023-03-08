#' @title Measures for Benchmarking Survival Models
#'
#' @description A convenient function that returns a table with measures that
#' can be used for benchmarking survival learners.
#'
#' @details The following measures are provided via [mlr3proba]:
#' - Discrimination measures (ranking `crank` prediction type)
#'    - [Harrell's C-index][mlr3proba::MeasureSurvCindex] (doesn't account for censoring)
#'    - [Uno's C-index][mlr3proba::MeasureSurvCindex] (estimates censoring distribution using KM)
#' - Calibration measures (survival `distr` prediction type)
#'    - [Integrated Brier Score][mlr3proba::MeasureSurvGraf] (with *proper = TRUE, method = 2*, estimates KM censoring distribution, is also a measure of discrimination)
#'    - [Right-Censored Log-Loss][mlr3proba::MeasureSurvRCLL] (proper calibration score)
#'    - [D-Calibration][mlr3proba::MeasureSurvDCalibration] (with *B = 10, chisq = FALSE*, returning the measure and not the p-value)
#'
#'  Any of the above measures **can** be used for model *evaluation* (but depends
#'  on several other factors if they **should** be used, e.g. censoring
#'  distribution, dataset size, proportion of events in the test set, etc.).
#'  For model *optimization/tuning*, we suggest the Right-Censored Log-Loss
#'  for survival predictions and Uno's C-index for risk predictions
#'  (e.g. linear predictors).
#'
#' @return a [data.table][data.table::data.table] with the suggested survival
#' measures used for benchmarking purposes
#'
#' @examples
#' library(mlr3proba)
#' ms = bench_msrs()
#' ms[id %in% c('ibrier','harrell_c')]$measure
#'
#' @export
bench_msrs = function() {
  # Harrell's C-index
  harrell_c = msr('surv.cindex')
  harrell_c$label = 'HarrellC'
  harrell_c$id = 'harrell_c'

  # Uno's C-index
  uno_c = msr('surv.cindex', weight_meth = 'G2')
  uno_c$label = 'UnoC'
  uno_c$id = 'uno_c'

  # Integrated Brier Score (with proper = TRUE, method = 2)
  ibrier = msr('surv.graf')
  ibrier$label = 'IBrier'
  ibrier$id = 'ibrier'
  ibrier$param_set$values$proper = TRUE

  # Right-Censored Log-Loss
  rcll = msr('surv.rcll')
  rcll$label = 'RCLL'
  rcll$id = 'rcll'

  # D-calibration
  dcal = msr('surv.dcalib', B = 10, chisq = FALSE)
  dcal$label = 'Dcal'
  dcal$id = 'dcal'

  data.table(
    id      = c('harrell_c', 'uno_c', 'ibrier', 'rcll', 'dcal'),
    measure = list(harrell_c, uno_c, ibrier, rcll, dcal)
  )
}
