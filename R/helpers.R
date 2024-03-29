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

#' @title Fit Bayesian LME models to compare multiple models
#'
#' @description
#' Using a general format of benchmarking data, this function fits the following
#' Bayesian linear-mixed effects (LME) model using [stan_glmer][rstanarm::stan_glmer]:
#'
#' \deqn{value \sim -1 + lrn\_id + (1 | task\_id/rsmp\_id)}
#'
#' where:
#' - `lrn_id` is the model used
#' - `task_id` is the dataset the model was used on
#' - `rsmp_id` is the specific resampling of the dataset
#' - `value` is the output value for a particular performance `measure`
#'
#' One LME model is fit per `measure`, possibly in parallel using
#' [future_lapply][future.apply::future_lapply].
#' The aim is to get one posterior distribution per `model`, for each
#' performance `measure` in the benchmarking data, while accounting for
#' possible inter-dependencies between the different datasets as well as
#' the nested within-resample correlation that might exist.
#'
#' See example in [MOBenchmark].
#'
#' @param df a `data.frame`/`tibble` object with the benchmarking results.
#' Columns names must include: `task_id`, `lrn_id`, `rsmp_id`, `measure` and
#' `value`.
#' For example, see the output of [reshape_mob_res].
#' @param nrsmps_perc Percentage of resampling ids to keep for faster
#' bayesian model fit (number between 0 and 1).
#' Default is 1, so we keep all of the resampling ids.
#' @param n_chains Number of Markon chains for MCMC. Default = 4
#' @param n_iters Number of iterations per chain. Default = 5000
#' @param n_cores Number of cores for MCMC parallelization. Default = 4
#' @param seed seed number for MCMC
#' @param verbose whether to print some informative messages. Default = TRUE.
#' @param ... other arguments passed on to [rstanarm::stan_glmer]
#'
#' @return A list of [stanreg-objects][rstanarm::stanreg-objects], one per
#' performance `measure` in the benchmarking data.
#'
#' @export
fit_blme_model_cmp = function(df, nrsmps_perc = 1, n_chains = 4,
  n_iters = 5000, n_cores = 4, seed = 42, verbose = TRUE, ...) {
  # input checks
  col_names = c('lrn_id', 'task_id', 'rsmp_id', 'measure', 'value')
  stopifnot(col_names %in% colnames(df))
  stopifnot(nrsmps_perc > 0, nrsmps_perc <= 1)

  if (verbose) message('\n### Fit Bayesian LMEs to compare learners')

  # reducing rsmp_ids for faster model fit
  if (nrsmps_perc < 1) {
    rsmp_ids = unique(df$rsmp_id)
    nrsmps = length(rsmp_ids)

    set.seed(seed) #' for `sample()`
    df = df %>%
      filter(rsmp_id %in% sample(x = rsmp_ids, size = ceiling(nrsmps * nrsmps_perc)))

    if (verbose) message('Reducing to ', ceiling(nrsmps * nrsmps_perc), ' rsmp_ids')
  }

  # Measures - one model per measure
  msr_ids = unique(df$measure)
  if (verbose) message('Measures: ', paste0(msr_ids, collapse = ', '))

  # progress bar - parallelize across measures
  pb = progressr::progressor(steps = length(msr_ids))

  model_list = future.apply::future_lapply(1:length(msr_ids), function(i) {
    set.seed(i)

    msr_id = msr_ids[i]
    pb(sprintf('%s', msr_id))

    df_sub = df %>%
      filter(measure == msr_id) %>%
      select(all_of(col_names))

    model = rstanarm::stan_glmer(
      data = df_sub,
      # MAIN FORMULA FOR MODEL COMPARISON
      formula = value ~ -1 + lrn_id + (1 | task_id/rsmp_id),
      chains = n_chains, cores = n_cores, iter = n_iters, seed = seed, ...
    )
  }, future.seed = TRUE)

  names(model_list) = msr_ids

  model_list
}
