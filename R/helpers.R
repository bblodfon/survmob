#' @title Create task powerset
#'
#' @description Given a list of [tasks][mlr3::Task], this function produces the
#' set (list) of all possible combinations of tasks (powerset), combining their
#' respective features via [mlr3pipelines::mlr_pipeops_featureunion].
#' Every task in the given list must have:
#'
#' - **Different feature names and task ids**
#' - Same target names (target columns only from the first task are kept)
#' - Same number of observations/rows
#'
#' @param `tasks` list of [mlr3::Task]s
#'
#' @export
task_powerset = function(tasks) {
  po_featureunion = mlr3pipelines::po('featureunion')

  # check that all task ids are different
  tsk_ids = mlr3misc::map_chr(tasks, `[[`, 'id')
  stopifnot(length(unique(tsk_ids)) == length(tsk_ids))

  # check that tasks have same number of rows
  stopifnot(length(unique(mlr3misc::map_dbl(tasks, `[[`, 'nrow'))) == 1)

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

#' @title Minimize Backend of Survival Task
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

#' @title Powerset intersect counts
#'
#' @description Use on a 0/1 matrix, where e.g. rows are patients and columns
#' are different data modalities (omics).
#' This function will return for all possible combination of omics,
#' the number of patients that have information on ALL these omics (intersection).
#'
#' @param `df` A 0/1 [data.frame][data.frame] or [matrix][matrix].
#' E.g. a (patient, omic) has a value of `1` if the specific patient has that
#' particular data modality, otherwise `0`.
#' Column names should be the name of the different omics in that case.
#'
#' @return a [tibble][tibble] with rows different omic combinations and columns:
#' - `omics` => list of omic names that make up the combo
#' - `n_omics` => number of omics in the combo
#' - `intersect_count` => number of patients who have all omics in the combo
#'
#' @examples
#' library(dplyr)
#'
#' set.seed(42)
#' m = matrix(sample(x = c(0,1), size = 20, replace = TRUE), ncol = 5, nrow = 4)
#' colnames(m) = LETTERS[1:5]
#' pics = powerset_icounts(m)
#' pics %>%
#'   arrange(desc(intersect_count)) %>%
#'   filter(n_omics > 3)
#'
#' @export
powerset_icounts = function(df) {
  omic_names = colnames(df)

  powerset = lapply(1:length(omic_names), combinat::combn, x = omic_names,
    simplify = FALSE) %>% unlist(recursive = FALSE)

  d = list()
  for (index in 1:length(powerset)) {
    as = powerset[[index]]
    n_omics = length(as)

    if (n_omics == 1) {
      intersect_count = sum(df[,as])
    } else {
      intersect_count = sum(rowSums(df[,as]) == length(as))
    }

    d[[index]] = tibble::tibble(
      omic_combo = list(as),
      n_omics = n_omics,
      intersect_count = intersect_count
    )
  }

  dplyr::bind_rows(d)
}
