#' @title PipeOpSurvShuffle
#' @template param_pipelines
#' @description Shuffles a Survival Task's targets (time, status)
#' @section Parameters:
#' - `replace` :: `logical(1)` \cr
#' Whether sampling should be done with replacement or not. Default: FALSE
#' @section Initialization:
#' ```
#' PipeOpSurvShuffle$new()
#' poss = po('survshuffle')
#' ```
#' @examples
#' library(mlr3)
#' library(mlr3proba)
#' library(mlr3pipelines)
#'
#' task = tsk('lung')
#' task$head(3)
#' poss = PipeOpSurvShuffle$new() # poss = po('survshuffle')
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
#' @description Removes features from a Task for which the proportion
#' of `NA` values are above a certain `cutoff` threshold.
#' @section Parameters:
#' - `cutoff` :: `numeric(1)`\cr
#' Features with more than `cutoff` percentage of NAs will be removed.
#' Default value: 0.2.
#' @section Initialization:
#' ```
#' PipeOpRemoveNAs$new()
#' pona = po('removenas')
#' ```
#' @examples
#' library(mlr3proba)
#' library(mlr3pipelines)
#'
#' task = tsk('lung')
#' task
#' pona = PipeOpRemoveNAs$new()
#' pona$train(list(task))[[1L]] # meal.cal removed
#'
#' pona$param_set$values$cutoff = 0 # remove every feature with at least 1 NA
#' pona$train(list(task))[[1L]]
#'
#' # another way to do the same
#' pona = po('removenas', cutoff = 0.1)
#' pona$train(list(task))[[1L]] # meal.cal removed
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
        can_subset_cols = FALSE, task_type = 'Task'
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
      lgl_vec = per_NA > pvals$cutoff
      # keep track of number of columns/features removed
      self$state = list(removed_column_num = sum(lgl_vec))
      # return task
      task$select(cols = names(per_NA)[!lgl_vec])
    },
    .predict_task = function(task) task # Do nothing during prediction
  )
)

#' @title PipeOpRemoveZeros
#' @template param_pipelines
#' @description
#' Removes features from a Task for which the proportion
#' of zero values (`0`) are above a certain `cutoff` threshold.
#' Useful when the task has read count expression data (e.g. mRNA data).
#' @section Parameters:
#' - `cutoff` :: `numeric(1)`\cr
#' Features with more than `cutoff` percentage of zeros will be removed.
#' Default value: 0.2.
#' @section Initialization:
#' ```
#' PipeOpRemoveZeros$new()
#' porz = po('removezeros')
#' ```
#' @examples
#' library(mlr3proba)
#' library(mlr3pipelines)
#'
#' df = data.frame(
#'   time = c(1,2,3,4),
#'   status = c(0,1,0,1),
#'   X1 = c(999,0,0,0),
#'   X2 = c(23,NA,0,0)
#' )
#'
#' task = as_task_surv(x = df, id = 'test', time = 'time', event = 'status')
#'
#' porz = PipeOpRemoveZeros$new(param_vals = list(cutoff = 0.7))
#' porz$train(list(task))[[1L]] # X1 removed, X2 remained
#'
#' # another way to do the same
#' porz = po('removezeros', cutoff = 0.7)
#' porz$train(list(task))[[1L]] # X1 removed, X2 remained
#'
#' @export
PipeOpRemoveZeros = R6Class('PipeOpRemoveZeros',
  inherit = PipeOpTaskPreproc,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function(id = 'remove_zeros', param_vals = list()) {
      p = ps(cutoff = p_dbl(lower = 0, upper = 1, tags = 'required'))
      p$values = list(cutoff = 0.2)
      super$initialize(id = id, param_set = p, param_vals = param_vals,
        can_subset_cols = FALSE, task_type = 'Task'
      )
    },

    #' @description How many features (columns) have zeros (i.e. at least one)?
    #' @param task [Task][mlr3::Task]
    ncolsZero = function(task) {
      nzeros = private$.get_nzeros(task)
      sum(nzeros > 0)
    }
  ),
  private = list(
    # data.table of 1 row: for every column/feature of the provided task,
    # the number of zeros
    .get_nzeros = function(task) {
      dt = task$data(cols = task$feature_names)
      nzeros = dt[, lapply(.SD, function(X) { sum(X == 0, na.rm = T) })]
      nzeros
    },
    .train_task = function(task) {
      pvals = self$param_set$get_values()
      per_zeros = private$.get_nzeros(task)/task$nrow
      lgl_vec = per_zeros > pvals$cutoff
      # keep track of number of columns/features removed
      self$state = list(removed_column_num = sum(lgl_vec))
      task$select(cols = names(per_zeros)[!lgl_vec])
    },
    .predict_task = function(task) task # Do nothing during prediction
  )
)

#' @title PipeOpLogTransform
#' @template param_pipelines
#' @description
#' Applies the [log()][base::log] function to every numeric feature of a Task.
#' NA's are ignored.
#' @section Parameters:
#' - `base` :: `numeric(1)`\cr
#' The base of the logarithm. Default is to use a base of 2.
#' - `offset` :: `numeric(1)`\cr
#' Offset value to add to each (numeric) feature.
#' Since this is primarily intended to be used with count mRNA data where there
#' are many zeros, default `offset` is 1 (to get a log value of 0).
#' @section Initialization:
#' ```
#' PipeOpLogTransform$new()
#' polog = po('logtransform')
#' ```
#' @examples
#' library(mlr3proba)
#' library(mlr3pipelines)
#'
#' df = data.frame(
#'   time = c(1,2,3,4),
#'   status = c(0,1,0,1),
#'   X1 = c(999,0,0,0),
#'   X2 = c(23,NA,0,0),
#'   X3 = c('a','b','c','d')
#' )
#'
#' task = as_task_surv(x = df, id = 'test', time = 'time', event = 'status')
#'
#' polog = po('logtransform')
#' polog$train(list(task))[[1L]] # X3 didn't change, X1 and X2 got log-transformed
#'
#' @export
PipeOpLogTransform = R6Class('PipeOpLogTransform',
  inherit = PipeOpTaskPreproc,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function(id = 'log_transform', param_vals = list()) {
      p = ps(
        base = p_dbl(lower = 0, upper = Inf, tags = 'required'),
        offset = p_dbl(tags = 'required')
      )
      p$values = list(base = 2, offset = 1)
      super$initialize(id = id, param_set = p, param_vals = param_vals,
        can_subset_cols = FALSE, task_type = 'Task'
      )
    }
  ),
  private = list(
    .train_task = function(task) {
      # original features
      features = task$col_roles$feature

      # only operate on numeric columns
      numeric_cols = task$feature_types[type %in% c('integer', 'numeric'), id]
      if (!length(numeric_cols)) return(task)
      dt = task$data(cols = numeric_cols)

      # get the hyperparameters and do the log-transform
      pvals = self$param_set$get_values()
      dt_trans = log(dt + pvals$offset, base = pvals$base)

      # $cbind() overwrites old task columns
      task$select(setdiff(task$feature_names, numeric_cols))$cbind(dt_trans)
      # keep the order the same
      task$col_roles$feature = features

      task
    },

    .predict_task = function(dt, levels) dt # Do nothing during prediction
  )
)
