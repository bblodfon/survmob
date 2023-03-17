#' @noRd

MOBenchRes = R6Class('MOBenchmarkResult',
  public = list(
    tune_rsmp = NULL,
    tune_measure_id = NULL,
    test_measure_ids = NULL,
    keep_models = FALSE,
    result = NULL, # tibble result, can always rewrite it!

    initialize = function(result) {
      self$result = result

      # fill in the rest
    },

    #' reduce to {task_id, lrn_id, score_{msr_id}}
    summarize = function(msr_id, mode = 'median') {

    },

    # combine with another class results, check for same fields
    combine = function() {
      # checks for similar task_ids, lrn_ids

    },

    # for specific {task_id + lrn_id}
    hpc_plot = function(task_id = NULL, lrn_id = NULL) {

    }
  )
)
