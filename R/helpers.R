#' @title Create task powerset
#'
#' @description Given a list of [tasks][mlr3::Task], this function produces the
#' set (list) of all possible combinations of tasks (powerset), combining their
#' respective features via
#' [mlr_pipeops_featureunion][mlr3pipelines::mlr_pipeops_featureunion].
#' Every task in the given list must have:
#'
#' - **Different feature names and task ids**
#' - Same target names (target columns only from the first task are kept)
#' - Same number of observations/rows
#'
#' @param tasks list of [Tasks][mlr3::Task]
#' @param check_targets logical (Default: TRUE).
#' Check if target columns have the same values (same target names is by
#' default checked).
#'
#' @examples
#' library(mlr3proba)
#'
#' # toy data
#' time  = 1:5
#' status = c(1,1,0,0,1)
#' d1 = data.frame(status = status, time = time, a1 = LETTERS[1:5], b1 = rep(1,5))
#' d2 = data.frame(status = status, time = time, a2 = LETTERS[1:5], b2 = rep(1,5))
#' d3 = data.frame(status = status, time = time, a3 = LETTERS[1:5], b3 = rep(1,5))
#'
#'  # survival tasks
#' task1 = as_task_surv(x = d1, time = 'time', event = 'status')
#' task2 = as_task_surv(x = d2, time = 'time', event = 'status')
#' task3 = as_task_surv(x = d3, time = 'time', event = 'status')
#' tasks = list(task1, task2, task3)
#' powset = task_powerset(tasks)
#' powset
#'
#' @export
task_powerset = function(tasks, check_targets = TRUE) {
  # at least 2 tasks
  stopifnot(length(tasks) > 1)

  # check that all task ids are different
  tsk_ids = mlr3misc::map_chr(tasks, `[[`, 'id')
  if (!length(unique(tsk_ids)) == length(tsk_ids)) {
    stop('Task ids must be different')
  }

  # check that tasks have same number of rows
  if (!length(unique(mlr3misc::map_dbl(tasks, `[[`, 'nrow'))) == 1) {
    stop('Tasks must have the same number of rows')
  }

  # check that target names are the same
  trg_names = mlr3misc::map_br(tasks, `[[`, 'target_names')
  for (index in 1:ncol(trg_names)) {
    if (!length(unique(trg_names[,index])) == 1) {
      stop('Target names must be the same')
    }
  }

  # check that target values (truth) are the same
  if (check_targets) {
    for (tsk_index in 1:(length(tasks) - 1)) {
      truth1 = tasks[[tsk_index]]$data(
        cols = tasks[[tsk_index]]$target_names
      )
      truth2 = tasks[[tsk_index + 1]]$data(
        cols = tasks[[tsk_index + 1]]$target_names
      )
      # ?all.equal.data.table
      if (!all.equal(truth1, truth2)) {
        stop('Target values must be the same')
      }
    }
  }

  task_subsets = lapply(1:length(tasks), combinat::combn, x = tasks,
    simplify = FALSE) %>% unlist(recursive = FALSE)

  po_featureunion = mlr3pipelines::po('featureunion')
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

#' @title Powerset intersection counts
#'
#' @description Use on a 0/1 matrix, where e.g. rows are patients and columns
#' are different data modalities (omics).
#' This function will return for all possible combination of omics,
#' the number of patients that have information on **ALL** these omics
#' (intersection counts).
#'
#' @param df A 0/1 [data.frame][data.frame] or [matrix][matrix].
#' E.g. a (patient, omic) has a value of `1` if the specific patient has that
#' particular data modality, otherwise `0`.
#' Column names should be the name of the different omics in that case.
#'
#' @return a [tibble][tibble] with rows different omic combinations and columns:
#' - `combo_name` => collapsed omic name (hyphen '-' is used for concatenation)
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
  if (!any(class(df) %in% c('matrix', 'data.frame'))) {
    stop('Provide matrix or data.frame-like object')
  }

  if (is.null(colnames(df))) {
    stop('Provide column names')
  }

  omic_names = colnames(df)

  powerset = lapply(1:length(omic_names), combinat::combn, x = omic_names,
    simplify = FALSE) %>% unlist(recursive = FALSE)

  d = list()
  for (index in 1:length(powerset)) {
    omic_combo = powerset[[index]]
    n_omics = length(omic_combo)
    combo_name = paste0(omic_combo, collapse = '-')

    if (n_omics == 1) {
      intersect_count = sum(df[,omic_combo])
    } else {
      intersect_count = sum(rowSums(df[,omic_combo]) == n_omics)
    }

    d[[index]] = tibble::tibble(
      combo_name = combo_name,
      omic_combo = list(omic_combo),
      n_omics = n_omics,
      intersect_count = intersect_count
    )
  }

  dplyr::bind_rows(d)
}

# for internal use: `part` is a list with train and test indices,
# used to split a dataset to two disjoint sets
assert_part = function(part) {
  stopifnot(class(part) == 'list')
  stopifnot(length(part) == 2)
  stopifnot(names(part) == c('train', 'test'))
  stopifnot(length(intersect(part$train, part$test)) == 0)
}
