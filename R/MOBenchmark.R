#' @title Multi-Omics Benchmark
#'
#' @description Use this [R6][R6::R6Class] class to execute a multi-omics
#' benchmark using several [survival learners][SurvLPS].
#'
#' The use-case here is that we have several datasets representing
#' different omic profiles for a particular patient cohort.
#' For example we have mRNA expression, CNA, methylation and clinical data for
#' a breast cancer patient cohort.
#' We can choose to generate all possible combinations of omic datasets
#' (referred to as tasks) or simply benchmark the given tasks (individual omics
#' or other datasets).
#' We train and tune several survival models on the same proportion of patients
#' (**train cohort**) across all tasks.
#' Performance is assessed using various survival measures on bootstrap
#' resamplings of a **separate test cohort**.
#'
#' The most important thing to keep in mind is that this class was made to
#' **benchmark datasets that share the same observations** (maybe not
#' conceptually in case of separate given tasks and `gen_task_powerset = FALSE`
#' but most definitely in number).
#' See example below.
#'
#' @examples
#' library(mlr3verse)
#' library(mlr3proba)
#' library(progressr)
#'
#' # Logging
#' lgr::get_logger('bbotk')$set_threshold('warn')
#' lgr::get_logger('mlr3')$set_threshold('warn')
#'
#' # Progress
#' options(progressr.enable = TRUE)
#' handlers(global = TRUE)
#' handlers('progress')
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
#'   lrn_ids = c('coxph', 'coxnet', 'aorsf'),
#'   tune_nevals = 2, test_nrsmps = 50, test_workers = 5,
#'   tune_rsmp = rsmp('holdout', ratio = 0.8),
#'   quiet = FALSE, keep_models = TRUE
#' )
#'
#' # execute benchmark
#' mob$run()
#'
#' # result tibble
#' mob$result
#'
#' # reshape benchmarking results
#' df = reshape_mob_res(mob$result)
#' df
#'
#' # Run Bayesian LME model
#' plan(multisession, workers = 2)
#' res = fit_blme_model_cmp(df)
#'
#' @export
MOBenchmark = R6Class('MOBenchmark',
  public = list(
    #' @field tasks Survival datasets to benchmark
    tasks = NULL,
    #' @field gen_task_powerset Generate task powerset?
    gen_task_powerset = NULL,
    #' @field part List partitioning the tasks to train/tuning and test sets
    part = NULL,
    #' @field lrn_ids Internal survival learner ids to benchmark
    lrn_ids = NULL,
    #' @field nthreads_rsf Number of cores for random survival forests
    nthreads_rsf = NULL,
    #' @field nthreads_xgb Number of cores for xgboost survival learners
    nthreads_xgb = NULL,
    #' @field tune_rsmp Tuning resampling
    tune_rsmp = NULL,
    #' @field tune_measure_id Survival measure for tuning
    tune_measure_id = NULL,
    #' @field tune_nevals Number of evaluations during Bayesian Optimization
    #' tuning
    tune_nevals = NULL,
    #' @field test_measure_ids Survival measures for testing
    test_measure_ids = NULL,
    #' @field test_nrsmps Number of bootstrap resamplings of the test set
    test_nrsmps = NULL,
    #' @field test_workers Number of cores for bootstrap parallelization
    test_workers = NULL,
    #' @field keep_models Keep the trained/tuned models?
    keep_models = NULL,
    #' @field quiet Show elapsed times for training and testing?
    quiet = NULL,
    #' @field result Tibble result with the benchmark results
    result = NULL,

    #' @description Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param tasks List of [TaskSurv][mlr3proba::TaskSurv].
    #' @param gen_task_powerset Whether to generate the
    #' [task powerset][task_powerset] from the given `tasks`.
    #' If TRUE (default), then the `tasks` should be different omic datasets,
    #' with the same target variables (`time`, `status`), i.e. corresponding to
    #' the same patient cohort.
    #' If FALSE, we use `tasks` as is, but you need to make sure they have the
    #' same number of observations (rows), otherwise the common partition
    #' to train and test sets (`part`) will not be tasks.
    #' @param part List with a disjoint (common) partition of the tasks to
    #' `train` and `test` indexes.
    #' The training/tuning is performed on the `train` set and the
    #' bootstrap testing on the `test` set for each task.
    #' Use the [partition][mlr3::partition()] function to perform stratified
    #' train/test splitting on the `status` indicator variable.
    #' @param lrn_ids Internal learner ids to use for benchmarking.
    #' Default is NULL, i.e. use all available learners.
    #' See [supported_lrn_ids()][SurvLPS] for which ids can be used.
    #' @param nthreads_rsf Number of cores to use in random survival forest
    #' learners (implicit parallelization).
    #' Default: use all available cores.
    #' @param nthreads_xgb Number of cores to use in xgboost survival learners
    #' (implicit parallelization).
    #' Default: use 2 cores.
    #' Never use only 1 core since during bootstrap resampling xgboost
    #' learners behave erratically in terms of CPU usage.
    #' @param tune_rsmp [Resampling][mlr3::Resampling] to use for tuning the
    #' learners on the `part$train` set.
    #' Default is repeated-CV (5 folds, 5 times).
    #' @param tune_measure_id Internal measure id to use for tuning the learners.
    #' Default is Uno's C-index.
    #' See [bench_msrs()] for available measures.
    #' @param tune_nevals Number of evaluations (hyperparameter configurations)
    #' to try during Bayesian Optimization tuning before termination.
    #' Default: 100.
    #' @param test_measure_ids Internal measure ids to use for assessing the
    #' performance of the learners in the test set (`part$test`).
    #' See [bench_msrs()] for available measures.
    #' @param test_nrsmps Number of bootstrap resamplings of the test set.
    #' Must be >= 1. Default: 1000.
    #' @param test_workers Number of workers for parallelization of the bootstrap
    #' resampling on the test set.
    #' This should be configured high enough based on available CPU cores.
    #' Default is 1.
    #' The value is temporarily overridden when benchmarking xgboost learners
    #' and set to 1 to avoid strange overuse of CPU.
    #' @param keep_models Whether to keep the trained models after tuning.
    #' Default: FALSE.
    #' @param quiet Whether to report elapsed timings for tuning and testing.
    #' Default: TRUE (**don't** report).
    initialize = function(
      tasks, gen_task_powerset = TRUE, part, lrn_ids = NULL,
      nthreads_rsf = unname(parallelly::availableCores()), nthreads_xgb = 2,
      tune_rsmp = rsmp('repeated_cv', repeats = 5, folds = 5),
      tune_measure_id = 'uno_c', tune_nevals = 100,
      test_measure_ids = c('uno_c', 'rcll'),
      test_nrsmps = 1000, test_workers = 1,
      keep_models = FALSE, quiet = TRUE
    ) {
      # some simple sanity checks
      if (test_nrsmps <= 1) {
        stop('`test_nrsmps` needs to be larger than 1')
      }
      if (test_workers <= 0) {
        stop('`test_workers` needs to be larger than 0')
      }
      if (tune_nevals <= 0) {
        stop('`tune_nevals` needs to be larger than 0')
      }
      stopifnot(nthreads_rsf > 0)
      stopifnot(nthreads_xgb > 0)

      # check tasks
      if (length(tasks) < 1) stop('Please provide some survival tasks')
      assert_tasks(tasks, task_type = 'surv')

      # check train/test partition
      assert_part(part)
      part_row_ids = union(part$train, part$test)
      # check that partition is matching row_ids in the tasks given
      for (task in tasks) {
        stopifnot('Partition does not match task\'s row ids' =
            all(part_row_ids %in% task$row_ids))
      }

      # check resampling
      assert_resampling(tune_rsmp)

      # check measures
      stopifnot(tune_measure_id %in% bench_msrs()$id)
      stopifnot(all(test_measure_ids %in% bench_msrs()$id))

      # assign fields
      self$tasks = tasks
      self$gen_task_powerset = gen_task_powerset
      self$part = part
      self$lrn_ids = lrn_ids
      self$nthreads_rsf = nthreads_rsf
      self$nthreads_xgb = nthreads_xgb
      self$tune_rsmp = tune_rsmp
      self$tune_measure_id = tune_measure_id
      self$tune_nevals = tune_nevals
      self$test_measure_ids = test_measure_ids
      self$test_nrsmps = test_nrsmps
      self$test_workers = test_workers
      self$keep_models = keep_models
      self$quiet = quiet
    },

    #' @description Execute benchmark
    #'
    #' @details Benchmark steps:
    #' - Generate task powerset (if `task_powerset = TRUE`)
    #' - Generate grid of learners and tasks
    #' - For each (learner, task) combo:
    #'    - Tune each learner in the train set
    #'    - Test the tuned learner's performance by taking multiple bootstrap
    #'    resamplings of the test set in the corresponding task and applying
    #'    several test measures
    #'
    #' @return a [tibble][tibble] with columns:
    #' - `task_id` => which task
    #' - `lrn_id` => which learner
    #' - `model` => the trained learner after tuning (if `keep_models = TRUE`)
    #' - `boot_res` => a [BootstrapResult] object (with the bootstrap results)
    run = function() {
      # tasks
      tasks = self$tasks

      # generate task powerset if specified
      if (self$gen_task_powerset && length(tasks) > 1) {
        message('Creating task powerset...\n')
        tasks = task_powerset(tasks)
      }

      # make sure the task list names are the task ids
      task_ids = unname(mlr3misc::map_chr(tasks, `[[`, 'id'))
      names(tasks) = task_ids

      # learners
      s = SurvLPS$new(
        ids = self$lrn_ids,
        nthreads_rsf = self$nthreads_rsf,
        nthreads_xgb = self$nthreads_xgb
      )
      lrn_tbl = s$lrn_tbl()

      # create benchmark grid
      bench_grid = data.table::CJ(
        task_id = task_ids,
        lrn_id = lrn_tbl$id,
        sorted = FALSE
      )
      bench_num = nrow(bench_grid)

      message('\nLearners: ', nrow(lrn_tbl))
      message('Tasks: ', length(task_ids))
      message('Benchmarks: ', bench_num)

      # execute benchmark
      result = list()
      for(row_id in 1:bench_num) {
        message('\nBenchmark ', row_id, '/', bench_num)
        message('#################')

        task_id = bench_grid[row_id]$task_id
        lrn_id  = bench_grid[row_id]$lrn_id

        task = tasks[[task_id]]
        learner = lrn_tbl[id == lrn_id]$learner[[1L]]$clone(deep = TRUE)
        search_space = lrn_tbl[id == lrn_id]$param_set[[1L]]

        # tune learner
        model = self$tune_learner(task, learner, search_space)

        # hack `test_workers` (set to 1) for bootstrap parallelization when
        # using xgboost learners to avoid strange issue that uses all cores
        # but gets nothing done
        test_workers = ifelse(
          grepl(pattern = 'xgb', x = lrn_id), yes = 1, no = self$test_workers
        )

        # test learner
        message('Bootstrap Testing')
        br = BootstrapResult$new(
          test_measure_ids = self$test_measure_ids,
          test_workers = test_workers,
          test_nrsmps = self$test_nrsmps
        )
        br$calculate(
          task = task,
          # model => trained learner on the best hyper-parameter configuration
          learner = model,
          part = self$part,
          quiet = self$quiet
        )

        result[[row_id]] = tibble(
          task_id  = task_id,
          lrn_id   = lrn_id,
          model    = base::switch(self$keep_models, list(model)),
          boot_res = list(br)
        )
      }

      res = bind_rows(result)
      self$result = res
      invisible(res)
    },

    #' @description Tune `learner` on the `given` task.
    #' The tuning is done on the train set of the task (parameter `part`),
    #' using a specific performance measure (parameter `tuning_measure_id`) and
    #' resampling (parameter `tune_rsmp`).
    #'
    #' The tuning is performed by constructing an [AutoTuner][mlr3tuning::AutoTuner]
    #' object and the optimization strategy is by default **Bayesian Optimization**.
    #' Optimization stops after a specific number of evaluations, see parameter
    #' `tune_nevals`.
    #' The hyperparameter configuration with the **best** average resampled
    #' performance is chosen to train a final model which is returned by this
    #' function.
    #'
    #' @param task a [TaskSurv][mlr3proba::TaskSurv]
    #' @param learner a [LearnerSurv][mlr3proba::LearnerSurv]
    #' @param search_space a [ParamSet][paradox::ParamSet] indicating the
    #' learner's hyperparameters that need to be tuned.
    #' If `NULL`, we just `train` the learner on the `task`.
    tune_learner = function(task, learner, search_space) {
      # get train set
      train_set = self$part$train

      # get measure
      measure = bench_msrs()[id == self$tune_measure_id]$measure[[1L]]

      # If NULL `search_space` just train, no tuning required
      if (is.null(search_space)) {
        message('Training: (', learner$id, ', ', task$id, ')')

        tic()
        learner$train(task, row_ids = train_set)
        toc(quiet = self$quiet)

        return(learner)
      }

      # Add early stopping callback for specific xgboost learners
      callbacks = list()
      es_set_cox = learner$param_set$values$XGBoostCox.early_stopping_set
      es_set_aft = learner$param_set$values$XGBoostAFT.early_stopping_set
      if ((!is.null(es_set_cox) && es_set_cox == 'test') ||
          (!is.null(es_set_aft) && es_set_aft == 'test')) {
        callbacks[[1]] = xgb_es_callback
      }

      # Initialize AutoTuner
      at = AutoTuner$new(
        learner = learner,
        resampling = self$tune_rsmp,
        measure = measure,
        search_space = search_space,
        terminator = trm('evals', n_evals = self$tune_nevals),
        tuner = tnr('mbo'),
        callbacks = callbacks
      )

      # Train AutoTuner
      message('Tuning: (', learner$id, ', ', task$id, ')')

      tic()
      at$train(task, row_ids = train_set)
      toc(quiet = self$quiet)

      at$learner
    },

    #' @description remove tasks from the class object to reduce size
    drop_tasks = function() {
      self$tasks = NULL
    },

    #' @description remove models from result tibble to reduce size
    drop_models = function() {
      self$result$model = NULL
    }
  )
)
