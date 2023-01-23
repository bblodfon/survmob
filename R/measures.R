#' @title Measures for Benchmarking Survival Models
#'
#' @description A convenient function that returns a table with measures that
#' can be used for benchmarking survival learners
#'
#' @details The following measures are provided via [mlr3proba]:
#' - Discrimination measures (ranking `crank` prediction type)
#'    - [Harrell's C-index][mlr3proba::MeasureSurvCindex] (doesn't account for censoring)
#'    - [Uno's C-index][mlr3proba::MeasureSurvCindex] (estimates censoring distribution using KM)
#' - Calibration measures (survival `distr` prediction type)
#'    - [Integrated Brier Score][mlr3proba::MeasureSurvGraf] (with *proper = TRUE, method = 2*, estimates KM censoring distribution)
#'    - [Right-Censored Log-Loss][mlr3proba::MeasureSurvRCLL] (strictly proper calibration score)
#'    - [D-Calibration][mlr3proba::MeasureSurvDCalibration] (with *B = 10, chisq = FALSE*, returning the measure and not the p-value)
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

  # Uno's C-index
  uno_c = msr('surv.cindex', weight_meth = 'G2')
  uno_c$label = 'UnoC'

  # Integrated Brier Score (with proper = TRUE, method=2)
  ibrier = msr('surv.graf')
  ibrier$label = 'IBrier'
  ibrier$param_set$values$proper = TRUE

  # Ring-Censored Log-Loss
  rcll = msr('surv.rcll')

  # D-calibration
  dcal = msr('surv.dcalib', B = 10, chisq = FALSE)

  data.table(
    id      = c('harrell_c', 'uno_c', 'ibrier', 'rcll', 'dcal'),
    measure = list(harrell_c, uno_c, ibrier, rcll, dcal)
  )
}

#' @title Out-of-bag (OOB) Survival Error Measure
#'
#' @name mlr_measures_surv.oob_error
#'
#' @description
#' Returns the out-of-bag error of a [LearnerSurv][mlr3proba::LearnerSurv] if it
#' supports it (e.g. a Random Forest Survival learner that has the property
#' `'oob_error'`).
#' The `oob_error` is usually (1 - C-index) but it depends on the learner.
#' Returns `NA` for unsupported learners.
#'
#' @section Dictionary:
#' This [Measure][mlr3::Measure] can be instantiated via the [dictionary][mlr3misc::Dictionary]
#' [mlr_measures][mlr3::mlr_measures] or with the associated sugar function [msr()][mlr3::msr]:
#' ```
#' MeasureSurvOOBError$new()
#' mlr_measures$get('surv.oob_error')
#' msr('surv.oob_error')
#' ```
#'
#' @examples
#' library(mlr3proba)
#' library(mlr3misc)
#' library(mlr3extralearners)
#'
#' l = lrn('surv.ranger')
#' t = tsk('rats')
#' l$train(t)
#' p = l$predict(t)
#' m = msr('surv.oob_error')
#' p$score(m, learner = l)
#'
#' @export
MeasureSurvOOBError = R6::R6Class('MeasureSurvOOBError',
  inherit = Measure,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      super$initialize(
        id = 'surv.oob_error',
        task_type = NA_character_,
        properties = c('na_score', 'requires_learner'),
        predict_type = 'crank', # every survival learner should have this
        range = c(-Inf, Inf),
        minimize = TRUE,
        label = 'OOB-Error'
      )
    }
  ),

  private = list(
    .score = function(prediction, learner, ...) {
      learner = learner$base_learner()
      if ('oob_error' %nin% learner$properties) {
        return(NA_real_)
      }

      # whatever the learner supports (e.g. 1 - Cindex or some kind of loss)
      return(learner$oob_error())
    }
  )
)
