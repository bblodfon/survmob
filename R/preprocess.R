#' @title PipeOpSurvShuffle
#' @template param_pipelines
#' @description Shuffles a Survival Task's targets (time, status)
#' @section Parameters:
#' - `replace` :: `logical(1)` \cr
#' Whether sampling should be done with replacement or not. Default: FALSE
#' @examples
#' library(mlr3)
#' library(mlr3proba)
#'
#' task = tsk('lung')
#' task$head(3)
#' poss = PipeOpSurvShuffle$new()
#' task_shuffled = poss$train(list(task))[[1L]]
#' task_shuffled$head(3)
#' @export
PipeOpSurvShuffle = R6Class('PipeOpSurvShuffle',
  inherit = PipeOpTaskPreproc,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function(id = 'survshuffle', param_vals = list()) {
      p = ps(replace = p_lgl(tags = 'required'))
      p$values = list(replace = FALSE)
      super$initialize(id = id, param_set = p, param_vals = param_vals,
        can_subset_cols = FALSE, task_type = 'TaskSurv'
      )
    }
  ),
  private = list(
    .train_task = function(task) {
      pvals = self$param_set$get_values()
      surv  = task$data(cols = c('time', 'status'))
      if (nrow(surv) > 1) {  # `sample` 'misbehaves' when 1st argument has length 1!
        surv$time   = sample(surv$time,   replace = pvals$replace)
        surv$status = sample(surv$status, replace = pvals$replace)
      }
      # $cbind() overwrites old task columns
      task$cbind(surv)
    },
    .predict_task = function(task) task # don't shuffle task during prediction!
  )
)

#' @title PipeOpRemoveNAs
#' @template param_pipelines
#' @description Removes features from a Survival Task for which the proportion
#' of `NA` values are above a certain `cutoff` threshold
#' @section Parameters:
#' @param cutoff `logical(1)`\cr
#' Features with more than `cutoff` percentage of NAs will be removed.
#' Default value: 0.2.
#' Higher `cutoff` values remove less features (less strict).
#' @examples
#' library(mlr3proba)
#'
#' task = tsk('lung')
#' task
#' pona = PipeOpRemoveNAs$new()
#' pona$train(list(task))[[1L]] # meal.cal removed
#'
#' pona$param_set$values$cutoff = 0 # remove every feature with at least 1 NA
#' pona$train(list(task))[[1L]]
#'
#' @export
PipeOpRemoveNAs = R6Class('PipeOpRemoveNAs',
  inherit = PipeOpTaskPreproc,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function(id = 'remove_na', param_vals = list()) {
      p = ps(cutoff = p_dbl(lower = 0, upper = 1, tags = 'required'))
      p$values = list(cutoff = 0.2)
      super$initialize(id = id, param_set = p, param_vals = param_vals,
        can_subset_cols = FALSE, task_type = 'TaskSurv'
      )
    },

    #' @description How many features (columns) have NAs (i.e. at least one)?
    #' @param task [Task][mlr3::Task]
    ncolsNA = function(task) {
      sum(task$missings(cols = task$feature_names) > 0)
    }
  ),
  private = list(
    .train_task = function(task) {
      pvals = self$param_set$get_values()
      per_NA = task$missings(cols = task$feature_names)/task$nrow
      task$select(cols = names(per_NA)[!per_NA > pvals$cutoff])
    },
    .predict_task = function(task) task # Do nothing during prediction
  )
)


#' @title Minimize Backend of Task
#'
#' @description Use this function when the backend of a task hasn't changed
#' after a lot of preprocessing, making the it to heavy to use in further
#' modeling steps
#'
#' @param task [TaskSurv][mlr3proba::TaskSurv]
#' @export
minimize_backend = function(task) {
  mlr3proba::as_task_surv(
    x = task$data(), id = task$id,
    time = 'time', event = 'status'
  )
}
