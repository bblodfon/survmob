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
