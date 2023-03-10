#' @title Bootstrap Result Class
#'
#' @description Use this class to measure the performance of a survival learner
#' on a test dataset using bootstrapping.
#' See example.
#'
#' @examples
#' library(mlr3proba)
#' library(mlr3pipelines)
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' # Lung task - impute missing values
#' task = tsk('lung')
#' pre = po('encode', method = 'treatment') %>>%
#'       po('imputelearner', lrn('regr.rpart'))
#' task = pre$train(task)[[1]]
#' # partition to train and test sets
#' part = partition(task)
#'
#' # train model on the train set
#' cox = lrn('surv.coxph', id = 'cox')
#' cox$train(task, row_ids = part$train)
#'
#' # bootstrap the test set 10 times and measure performance
#' # using all available test metrics from `bench_msrs()`
#' # can be parallelized by increasing `test_workers`
#' brs = BootRes$new(test_nrsmps = 10)
#' brs$calculate(task = task, learner = cox, part = part)
#'
#' # performance on the test set
#' brs$score
#'
#' # bootstrapped performance scores
#' brs$scores
#'
#' # Median summary score per metric
#' brs$score_median()
#'
#' # Percentile Confidence Intervals
#' brs$percent_ci()
#' brs$percent_ci(conf = 0.9)
#'
#' @export
BootRes = R6Class('BootstrapResult',
  public = list(
    #' @field task_id Task id
    task_id = NULL,
    #' @field lrn_id Learner id
    lrn_id = NULL,
    #' @field test_measure_ids Measure ids
    test_measure_ids = NULL,
    #' @field test_nrsmps Number of bootstrapped resamplings
    test_nrsmps = NULL,
    #' @field test_workers Number of workers for bootstrap parallelization
    test_workers = NULL,
    #' @field score Performance on the original test set
    score = NULL,
    #' @field scores Performance on bootstrapped test sets
    scores = NULL,

    #' @description Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param test_measure_ids Measure ids.
    #' Default: use all available measures from [bench_msrs()] function.
    #' @param test_nrsmps Number of bootstrapped resamplings. Default: 1000.
    #' @param test_workers Number of workers for bootstrap parallelization.
    #' Default is 1 (sequential).
    #' We advise careful consideration of how many workers to use.
    #' Several factor might influence this, such as the total number of CPUs
    #' and memory available, the size of the dataset, the learner's memory
    #' footprint, the learner itself (e.g. `xgboost` might have problems), etc.
    initialize = function(test_measure_ids = bench_msrs()$id,
      test_nrsmps = 1000, test_workers = 1) {
      if (!all(test_measure_ids %in% bench_msrs()$id)) {
        stop('Use only supported test measure ids from `bench_msrs()`')
      }
      if (test_nrsmps <= 1) {
        stop('`test_nrsmps` needs to be larger than 1')
      }
      if (test_workers <= 0) {
        stop('`test_workers` needs to be larger than 0')
      }

      self$test_measure_ids = test_measure_ids
      self$test_nrsmps = test_nrsmps
      self$test_workers = test_workers
    },

    #' @description Calculate bootstrap performance scores on the test set of
    #' the given survival task.
    #'
    #' @param task A [survival task][mlr3proba::TaskSurv]
    #' @param learner A [learner][mlr3proba::LearnerSurv], **trained** on the
    #' `part$train` set of the given task.
    #' @param part A [partition][mlr3::partition] of the given task to train
    #' and test sets.
    #' @param quiet Default: TRUE (**don't** report time elapsed)
    #'
    #' @details Parallelization: we check if the system allows for
    #' [multicore][future::multicore]
    #' parallelization backend strategy (Linux systems mostly) and use that,
    #' otherwise we use [multisession][future::multisession].
    #'
    #' @return (invisibly) a [tibble][tibble::tibble] (`scores`) with columns the
    #' metrics used and rows the performance scores measured on the different
    #' bootstrapped resamplings.
    #' Also the performance per metric on the original test set is recorded in
    #' the variable `score`.
    calculate = function(task, learner, part, quiet = TRUE) {
      stopifnot('`calculate()` has already been called!' = is.null(self$scores))

      assert_task(task)
      assert_learner(learner)
      assert_part(part)

      if (is.null(learner$model)) {
        stop('Cannot predict, learner has not been trained yet')
      }

      # assign task and learner ids
      self$task_id = task$id
      self$lrn_id  = learner$id

      # get measures
      measures = bench_msrs()[id %in% self$test_measure_ids]$measure

      # get score on original test set
      self$score =
        learner$
        predict(task, row_ids = part$test)$
        score(measures, task = task, train_set = part$train) %>%
        as_tibble_row()

      # decide on parallelization strategy
      strategy = ifelse(parallelly::supportsMulticore(), 'multicore', 'multisession')
      future::plan(strategy, workers = self$test_workers)

      # get scores on bootstrap test sets
      tic()
      scores = future.apply::future_lapply(1:self$test_nrsmps, function(i) {
        # generate bootstrap test dataset
        set.seed(i)
        row_ids = sample(x = part$test, size = length(part$test), replace = TRUE)
        test_data = task$data(rows = row_ids)

        # predictions
        p = learner$predict_newdata(test_data)

        # scores
        p$score(measures, task = task, train_set = part$train)
      }, future.seed = TRUE) %>% bind_rows()
      toc(quiet = quiet)
      future::plan('sequential') # disable parallelization

      self$scores = scores
      invisible(scores)
    },

    #' @description Percentile confidence intervals per measure
    #'
    #' @param conf Confidence level. Default: 0.95 (i.e. the 2.5th and 97.5th
    #' percentile points are used to define the confidence interval range)
    #'
    #' @return a [tibble][tibble::tibble] with confidence intervals for the
    #' performance score (per measure used).
    percent_ci = function(conf = 0.95) {
      stopifnot('Run `calculate()` first'= !is.null(self$scores))
      stopifnot(conf > 0, conf < 1)

      alphas = (1 + c(-conf, conf))/2
      apply(self$scores, 2, quantile, probs = alphas) %>%
        as_tibble(rownames = 'perc')
    },

    #' @description Median summary score for the bootstrapped scores
    #' (one per measure)
    #'
    #' @return a [tibble row][tibble::tibble] with the median value of a score
    #' across all bootstrap resamplings.
    score_median = function() {
      stopifnot('Run `calculate()` first'= !is.null(self$scores))

      apply(self$scores, 2, median, na.rm = TRUE) %>%
        as_tibble_row()
    }
  )
)
