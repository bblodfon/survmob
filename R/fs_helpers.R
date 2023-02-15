#' @title Create task powerset
#'
#' @description Given a list of [tasks][mlr3::Task], this function produces the
#' set (list) of all possible combinations of tasks (powerset), combining their
#' respective features via [mlr3pipelines::mlr_pipeops_featureunion].
#' Every task in the given list must have **different feature names and task ids**,
#' but the same target names.
#'
#' @param `tasks` list of [mlr3::Task]s
#'
#' @export
task_powerset = function(tasks) {
  po_featureunion = mlr3pipelines::po('featureunion')

  # check that all task ids are different
  tsk_ids = mlr3misc::map_chr(tasks, `[[`, 'id')
  stopifnot(length(unique(tsk_ids)) == length(tsk_ids))

  task_subsets = lapply(1:length(tasks), combinat::combn, x = tasks,
    simplify = FALSE) %>% unlist(recursive = FALSE)

  powerset = list()
  for (task_subset in task_subsets) {
    task = po_featureunion$train(task_subset)[[1L]]
    task_id = paste0(sapply(task_subset, function(l) l$id), collapse = '-')
    task$id = task_id
    powerset[[task_id]] = task
  }

  # check that the empty subset is excluded
  stopifnot(length(powerset) == (2^length(tasks)-1))

  powerset
}
