#' @title Early Stopping Callback for Graph XGBoost Learner
#'
#' @description
#' This [mlr3tuning::CallbackTuning] integrates early stopping into the
#' hyperparameter tuning of an XGBoost learner.
#' Early stopping estimates the optimal number of trees (`nrounds`) for a given
#' hyperparameter configuration.
#' Since early stopping is performed in each resampling iteration, there are
#' several optimal `nrounds` values.
#' The callback writes the maximum value to the archive in the `max_nrounds`
#' column.
#' In the best hyperparameter configuration (`instance$result_learner_param_vals`),
#' the value of `nrounds` is replaced by `max_nrounds` and early stopping is
#' deactivated.
#'
#' @details
#' This callback works only with xgboost `GraphLearner`s available from this
#' package (mainly due to using manually chosen xgboost learner ids and
#' hyperparameter names).
#' Inspired by the default early stopping callback (see [code](https://github.com/mlr-org/mlr3tuning/blob/02e9f338fad6c921ef040df8f63a8249ce6b783b/R/mlr_callbacks.R#L48-L88) and
#' [issue](https://github.com/mlr-org/mlr3tuning/issues/376)).
#' The callback is compatible with the [AutoTuner].
#' The final model is fitted with the best hyperparameter configuration and
#' `max_nrounds` as `nrounds` i.e. early stopping is not performed.
#'
#' @examples
#' library(mlr3proba) # loads mlr3 as well
#' library(mlr3tuning)
#' library(mlr3pipelines)
#' library(mlr3extralearners)
#'
#' # task
#' task = tsk('lung')
#' poe = po('encode')
#' task = poe$train(list(task))[[1L]]
#'
#' # xgboost graph learner with distr prediction + early stopping on the test set
#' s = SurvLPS$new(nthreads_xgb = 2, ids = c('xgboost_cox_early'))
#' dt = s$lrn_tbl()
#' xgb_learner = dt$learner[[1L]]
#' xgb_ps = dt$param_set[[1L]]
#'
#' xgb_at = AutoTuner$new(
#'   learner = xgb_learner,
#'   resampling = rsmp('holdout'),
#'   measure = msr('surv.rcll'),
#'   search_space = xgb_ps,
#'   terminator = trm('evals', n_evals = 5),
#'   tuner = tnr('random_search'),
#'   callbacks = xgb_es_callback
#'   # early stopping callback to get the proper nrounds at the end
#' )
#' xgb_at$train(task)
#'
#' as.data.table(xgb_at$archive)
#' xgb_at$learner
#'
#' @export
xgb_es_callback = mlr3tuning::callback_tuning('mlr3tuning.early_stopping_graph_learner',
  on_optimization_begin = function(callback, context) {
    learner = context$instance$objective$learner # learner has an id here but not later

    # store models temporarily
    callback$state$store_models = context$instance$objective$store_models
    context$instance$objective$store_models = TRUE
  },

  on_eval_after_benchmark = function(callback, context) {
    callback$state$max_nrounds = mlr3misc::map_dbl(
      context$benchmark_result$resample_results$resample_result, function(rr) {
      max(mlr3misc::map_dbl(
        mlr3misc::get_private(rr)$
          .data$learner_states(mlr3misc::get_private(rr)$.view), function(state) {
            # Use path to XGBoost model in graph learner
            if (!is.null(state$model$XGBoostCox)) {
              state$model$XGBoostCox$model$best_iteration
            } else {
              state$model$XGBoostAFT$model$best_iteration
            }
      }))
    })
  },

  on_eval_before_archive = function(callback, context) {
    data.table::set(context$aggregated_performance, j = 'max_nrounds', value = callback$state$max_nrounds)
    if (!callback$state$store_models) context$benchmark_result$discard(models = TRUE)
  },

  on_result = function(callback, context) {
    # Prefix parameters with appropriate learner id
    lrn_id = context$instance$objective$learner$id
    # print(lrn_id) # for checking if it works

    context$result$learner_param_vals[[1]][[paste0(lrn_id, '.early_stopping_rounds')]] = NULL
    context$result$learner_param_vals[[1]][[paste0(lrn_id, '.nrounds')]] = context$instance$archive$best()$max_nrounds
    context$result$learner_param_vals[[1]][[paste0(lrn_id, '.early_stopping_set')]] = 'none'
    context$instance$objective$store_models = callback$state$store_models
  }
)
