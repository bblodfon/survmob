#' @title Bootstrap Result Class
#'
#' @description Use this class to measure the performance of a survival learner
#' on a test dataset using bootstrapping.
#'
#' @examples
#' library(mlr3proba)
#' library(tibble)
#' library(dplyr)
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
#' # can be parallelized by increaaing `test_workers`
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
#'
#' @export
BootRes = R6Class('BootstrapResult',
  public = list(
    #' @field task_id mpla
    task_id = NULL,
    #' @field lrn_id mpla
    lrn_id = NULL,
    #' @field test_measure_ids mpla
    test_measure_ids = NULL,
    #' @field test_nrsmps mpla
    test_nrsmps = NULL,
    #' @field test_workers mpla
    test_workers = NULL,
    #' @field score mpla
    score = NULL,
    #' @field scores mpla
    scores = NULL,

    #' @description mpla mpla
    #' @param test_measure_ids mpla
    #' @param test_nrsmps mpla
    #' @param test_workers mpla
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

    #' @description Calculate bootstrapping performance scores
    #' @param task mpla
    #' @param learner mpla
    #' @param part mpla
    #' @param quiet Default: TRUE => don't report time elapsed
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

      # get scores on bootstrap test sets
      future::plan('multisession', workers = self$test_workers)
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
    #' @param conf Confidence level. Default: 0.95 (i.e. the 2.5th and 97.5th
    #' percentile points are used to define the confidence interval range)
    percent_ci = function(conf = 0.95) {
      stopifnot('Run `calculate()` first'= !is.null(self$scores))
      stopifnot(conf > 0, conf < 1)

      alphas = (1 + c(-conf, conf))/2
      apply(self$scores, 2, quantile, probs = alphas) %>%
        as_tibble(rownames = 'perc')
    },

    #' @description Median summary score for the bootstrapping scores
    #' (one per measure)
    score_median = function() {
      stopifnot('Run `calculate()` first'= !is.null(self$scores))

      apply(self$scores, 2, median, na.rm = TRUE) %>%
        as_tibble_row()
    }
  )
)
