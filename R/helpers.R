#' @title Create task powerset
#'
#' @description Given a list of [tasks][mlr3::Task], this function produces the
#' set (list) of all possible combinations of tasks (powerset), combining their
#' respective features via
#' [mlr_pipeops_featureunion][mlr3pipelines::mlr_pipeops_featureunion].
#' Every task in the given list must have:
#'
#' - **Different** feature names and task ids
#' - **Same** target names
#' - **Same** number of observations/rows
#'
#' @details Only the target columns from the first task of a combination are
#' kept. You need to add the `stratum` property to **all input tasks** if you
#' want every possible task combination to have it as well.
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
#' # survival tasks
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
      res = all.equal(truth1, truth2, check.attributes = FALSE)
      if (class(res) == 'character') {
        stop('Target values must be the same! ', res)
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
#' after a lot of preprocessing, making it too heavy to use in further
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

#' @title Reshape MOBenchmark result object to tibble
#'
#' @description
#' Reshape the output `result` object from [survmob::MOBenchmark]
#'
#' @param res a `tibble` with columns `task_id`, `lrn_id` and `boot_res`.
#' Every other column is discarded.
#' @param add_modality_columns whether to generate individual columns for each
#' modality listed in the `task_id` column.
#' Modality names should be separated by a dash ('-') in multi-modal `task_id`s
#' and have different names.
#' The modality columns will have values 0 or 1, indicating the presence or
#' absence of a corresponding modality in the `task_id` column.
#' Default: `TRUE`.
#' @param endfix string to add to the end of each newly generated modality
#' column.
#' Default: `_omic`.
#'
#' @return `tibble` with columns `task_id`, `lrn_id`, `rsmp_id`, `measure`,
#' `value` and possibly one column per omic/modality with values 1 or 0.
#'
#' @examples
#' library(mlr3verse)
#' library(mlr3proba)
#'
#' # Logging
#' lgr::get_logger('bbotk')$set_threshold('warn')
#' lgr::get_logger('mlr3')$set_threshold('warn')
#'
#' # task lung
#' task = tsk('lung')
#' pre = po('encode', method = 'treatment') %>>%
#'       po('imputelearner', lrn('regr.rpart'))
#' task = pre$train(task)[[1]]
#'
#' # partition to train and test sets
#' part = partition(task, ratio = 0.8)
#'
#' mob = MOBenchmark$new(
#'   tasks = list(task), part = part,
#'   lrn_ids = c('coxph', 'coxnet'),
#'   tune_nevals = 2, test_nrsmps = 20, test_workers = 1,
#'   tune_rsmp = rsmp('holdout', ratio = 0.8),
#'   quiet = FALSE, keep_models = TRUE
#' )
#'
#' # execute benchmark
#' mob$run()
#'
#' # reshape result tibble
#' reshape_mob_res(mob$result)
#'
#' @export
reshape_mob_res = function(res, add_modality_columns = TRUE, endfix = '_omic') {
  # check number of test bootstrap resamplings
  nrsmps_vec = mlr3misc::map_dbl(res$boot_res, `[[`, 'test_nrsmps')
  if (length(unique(nrsmps_vec)) != 1) {
    stop('Number of test bootstrap resamplings different?!')
  }
  nrsmps = res$boot_res[[1]]$test_nrsmps

  df = res %>%
    select(task_id, lrn_id, boot_res) %>%
    mutate(scores = purrr::map(boot_res, 'scores')) %>%
    select(!matches('boot_res')) %>%
    tibble::add_column(id = list(tibble(rsmp_id = 1:nrsmps))) %>%
    tidyr::unnest(cols = c(id, scores)) %>%
    dplyr::relocate(rsmp_id, .after = lrn_id) %>%
    tidyr::pivot_longer(cols = !matches('task_id|lrn_id|rsmp_id'),
      names_to = 'measure', values_to = 'value')

  if (add_modality_columns) {
    # Extract modality names
    modalities = unique(unlist(strsplit(df$task_id, '-')))

    # Iterate over modalities and create new columns
    for (modality in modalities) {
      mod_column = paste0(modality, endfix)
      df = df %>%
        mutate(!!mod_column := as.integer(str_detect(task_id, modality)))
    }
  }

  df
}

