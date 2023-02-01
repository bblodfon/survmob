#' @title Ensemble Feature Selection (eFS)
#'
#' @description This class stores the configuration and result of the proposed
#' ensemble feature selection approach, as well as the function (`run()`) to
#' execute it.
#'
#' We use a couple of different **random survival forests** (RSFs) as learners
#' and a **recursive feature elimination algorithm** (RFE, wrapper-based feature
#' selection) to find the most predictive feature subset in the given task
#' (dataset) for each learner.
#' This process repeats a number of times for each RSF learner and therefore
#' returns different best feature subsets for further analysis and processing
#' (e.g. finding robust features).
#'
#' @examples
#' library(mlr3proba)
#' library(mlr3fselect)
#' library(mlr3pipelines)
#' library(mlr3extralearners)
#'
#' # less logging
#' lgr::get_logger('bbotk')$set_threshold('warn')
#' lgr::get_logger('mlr3')$set_threshold('warn')
#'
#' # Lung task - impute missing values
#' task = tsk('lung')
#' pre = po('encode', method = 'treatment') %>>%
#'       po('imputelearner', lrn('regr.rpart'))
#' task = pre$train(task)[[1]]
#'
#' # supported lrn ids
#' eFS$new()$supported_lrn_ids()
#'
#' # create eFS object where every RSF is trained on each RFE feature subset
#' # using all features (train set equals the test set) and the out-of-bag
#' # error is used to assess predictive performance (1 - Cindex)
#' efs = eFS$new(lrn_ids = c('rsf_logrank', 'rsf_cindex'), nthreads_rsf = 4,
#'   feature_fraction = 0.9, n_features = 1, mtry_ratio = 0.5,
#'   repeats = 3, msr_id = 'oob_error', resampling = rsmp('insample')
#' )
#'
#' # useful info for RFE (adaptive mtry.ratio, subset sizes)
#' efs$rfe_info(task)
#'
#' # execute ensemble feature selection
#' efs$run(task)
#'
#' # result tibble
#' res = efs$result
#'
#' res
#' @export
eFS = R6Class('EnsembleFeatureSelection',
  public = list(
    #' @field task_id (`character(1)`)\cr
    #' Task id
    task_id = NULL,

    #' @field lrn_ids (`character()`)\cr
    #' Learner ids
    lrn_ids = NULL,

    #' @field msr_id (`character(1)`)\cr
    #' Measure id
    msr_id = NULL,

    #' @field resampling (`Resampling`)\cr
    #' [Resampling][mlr3::Resampling]
    resampling = NULL,

    #' @field repeats (`int(1)`)\cr
    #' Number of times to run RFE on each RSF learner
    repeats = NULL,

    #' @field n_features (`int(1)`)\cr
    #' Number of features that signals the termination of RFE
    n_features = NULL,

    #' @field feature_fraction (`int(1)`)\cr
    #' Fraction of features to retain in each iteration of RFE
    feature_fraction = NULL,

    #' @field nthreads_rsf (`int(1)`)\cr
    #' Number of cores to use in the random forest survival learners
    nthreads_rsf = NULL,

    #' @field num_trees (`int(1)`)\cr
    #' Number of trees to use in the random forest survival learners
    num_trees = NULL,

    #' @field mtry_ratio (`double(1)`)\cr
    #' Percentage of features to try at each node split (should be between 0
    #' and 1)
    mtry_ratio = NULL,

    #' @field adaptive_mr (`logical(1)`)\cr
    #' Whether to have an adaptive `mtry_ratio` or not in the RFE algorithm
    adaptive_mr = NULL,

    #' @field result (`tibble`)\cr
    #' A [tibble] with the results from the eFS (see `run()` method)
    result = NULL,

    #' @description Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param lrn_ids Internal learner ids.
    #' These will let us use the corresponding [learner][mlr3::Learner]s
    #' in the [AutoFSelector][mlr3fselect::AutoFSelector].
    #' We get the learners using the [SurvLPS] class.
    #' All ids correspond to random survival forests and by default we use all
    #' of them.
    #' Check available learner ids with the method `supported_lrn_ids()`.
    #'
    #' @param msr_id Internal measure id.
    #' This will let us get the [measure][mlr3::Measure] to use in the
    #' [AutoFSelector][mlr3fselect::AutoFSelector].
    #' Must be one of the available ids provided via the [bench_msrs] function
    #' or the [oob_error][mlr3::mlr_measures_oob_error].
    #' Default is Harrell's C-index.
    #'
    #' @param resampling [Resampling][mlr3::Resampling] for the
    #' [AutoFSelector][mlr3fselect::AutoFSelector].
    #' Default resampling is 5 times 5-fold CV.
    #'
    #' @param repeats Number of times to run the RFE algorithm on each RSF
    #' learner. Defaults to 100.
    #'
    #' @param n_features Number of features that signals the termination of the
    #' RFE algorithm. Defaults to 2.
    #'
    #' @param feature_fraction Fraction of features to retain in each iteration
    #' of the RFE algorithm. Defaults to 0.8.
    #'
    #' @param nthreads_rsf Number of cores to use in the random forest survival
    #' learners.
    #' By default we use all available cores to take full advantage of
    #' the implicit parallelization of the RSF learners.
    #'
    #' @param num_trees Number of trees to use in the random forest survival
    #' learners. Defaults to 250.
    #'
    #' @param mtry_ratio Percentage of features to try at each tree node split
    #' (should be between 0 and 1).
    #' Default value is 0.05 which means that 500 features will be randomly
    #' selected in a node split, considering a dataset of 10000 features.
    #'
    #' @param adaptive_mr Whether to use an adaptive `mtry_ratio` or not.
    #' Default: TRUE.
    #'
    #' @details
    #' - If `msr_id` is `oob_error`:
    #'    - Make sure to check which is the `oob_error` provided by each RSF
    #' learner (should be the same ideally, usually 1 - Cindex)
    #'    - It's more efficient to have an [insample resampling][mlr_resamplings_insample]
    #'    in that case.
    #' - `adaptive_mr`:
    #'    - During the execution of the RFE algorithm, progressively smaller
    #'    feature subsets are used, as dictated by the `feature_fraction`
    #'    parameter.
    #'    - Setting `adaptive_mr` to FALSE means that the same value of `mtry_ratio`
    #'    is used, irrespective of the number of features in each subset.
    #'    This creates the following issue: when the feature subsets get really
    #'    smaller (e.g. less 40 features) and with an example `mtry_ratio` of
    #'    0.05, we get only 2 features per tree node split, which results in
    #'    extremely randomized trees.
    #'    - Setting `adaptive_mr` to TRUE (default), the value of `mtry_ratio`
    #'    changes in each iteration of the RFE algorithm according to the
    #'    following formula:
    #'    \deqn{mr^{log(s)/log(n)}}, where `mr` is the `mtry_ratio` provided
    #'    at initialization, `subset_size` is the current number of
    #'    features in the RFE subset and `n` is the total number of features
    #'    in the original dataset.
    #'    The formula results in an increase of the `mtry_ratio`s in the RSFs
    #'    as the feature subsets selected by the RFE algorithm get smaller.
    #'    Therefore more bagged trees are being created with smaller feature
    #'    subsets.
    initialize = function(
      lrn_ids = self$supported_lrn_ids(), msr_id = 'harrell_c',
      resampling = mlr3::rsmp('repeated_cv', repeats = 5, folds = 5),
      repeats = 100, n_features = 2, feature_fraction = 0.8,
      nthreads_rsf = parallelly::availableCores(), num_trees = 250,
      mtry_ratio = 0.05, adaptive_mr = TRUE) {

      # learner ids
      supp_lrn_ids = self$supported_lrn_ids()
      if (!all(lrn_ids %in% supp_lrn_ids) || length(lrn_ids) == 0) {
        stop(paste0('You have used some unsupported learner ids or none at all.
          \nAvailable ids: ', paste0(supp_lrn_ids, collapse = ', ')))
      }
      self$lrn_ids = lrn_ids

      # RFE/AutoFSelector parameters
      supp_msr_ids = c(bench_msrs()$id, 'oob_error')
      if (!msr_id %in% supp_msr_ids) {
        stop(paste0('You have used an unsupported measure id. \n
          Available ids: ', paste0(supp_msr_ids, collapse = ',')))
      }
      self$msr_id = msr_id

      assert_resampling(resampling)
      self$resampling = resampling

      # 'oob_error' should be used with 'insample' resampling
      if (self$msr_id == 'oob_error' && self$resampling$id != 'insample')
        warning('Note: \'oob_error\' is used but not \'insample\' resampling')

      if (repeats > 0)
        self$repeats = repeats
      else
        stop('\'repeats\' needs to be larger than 0')

      if (n_features > 0)
        self$n_features = n_features
      else
        stop('\'n_features\' needs be larger than 0')

      if (feature_fraction >= 0 && feature_fraction < 1)
        self$feature_fraction = feature_fraction
      else
        stop('\'feature_fraction\' needs to be between 0 (inclusive) and 1
          (non-inclusive)')

      if (is.logical(adaptive_mr))
        self$adaptive_mr = adaptive_mr
      else
        stop('\'adaptive_mr\' can only be TRUE or FALSE')

      # RSF parameters
      if (nthreads_rsf > 0)
        self$nthreads_rsf = nthreads_rsf
      else
        stop('\'nthreads_rsf\' needs to be larger than 0')

      if (num_trees > 0)
        self$num_trees = num_trees
      else
        stop('\'num_trees\' needs to be larger than 0')

      if (mtry_ratio >= 0 && mtry_ratio <= 1)
        self$mtry_ratio = mtry_ratio
      else
        stop('\'mtry_ratio\' needs to be between 0 and 1')
    },

    #' @description Returns a vector of internal ids corresponding to the
    #' random survival forest learners that can be used in method `run()`
    supported_lrn_ids = function() {
      s = SurvLPS$new()
      ids = s$lrn_ids()
      ids[grepl(pattern = 'rsf_', ids)]
    },

    #' @description Returns the feature subset sizes that the RFE algorithm
    #' will use and the corresponding `mtry_ratio` and `mtry` values for the
    #' RSF learners.
    #' These depend on the total number of features of the given task,
    #' as well as the `n_features` and `feature_fraction` parameters
    #' (initialized upon class construction).
    #'
    #' @details
    #' The subset sizes are calculated from the respective
    #' [code](https://github.com/mlr-org/mlr3fselect/blob/HEAD/R/FSelectorRFE.R#L123).
    #'
    #' @param task [mlr3::Task]
    #' @return a [tibble] with columns `subset_size`, `mtry_ratio` and `mtry`
    rfe_info = function(task) {
      n = length(task$feature_names) # number of total features
      n_features = self$n_features
      feature_fraction = self$feature_fraction
      mr = self$mtry_ratio

      subset_sizes = unique(floor(cumprod(
        c(n, rep(feature_fraction, log(n_features / n) / log(feature_fraction))))
      ))

      dt = tibble(subset_size = subset_sizes)
      if (self$adaptive_mr) {
        dt = dt %>% mutate(mtry_ratio = mr^(log(subset_size)/log(n)))
      } else {
        dt = dt %>% mutate(mtry_ratio = mr)
      }

      dt = dt %>% mutate(mtry = ceiling(mtry_ratio * subset_size))

      dt
    },

    #' @description Runs the ensemble feature selection on the given task.
    #'
    #' @details For every RSF learner, we run `repeats` times the RFE algorithm
    #' using a properly constructed [AutoFSelector][mlr3fselect::AutoFSelector]).
    #' From each `repeat`ition we get the best feature subset.
    #' The aggregated result is invisibly returned as a tibble and is available
    #' also in the `result` field of this object.
    #'
    #' @param task [TaskSurv][mlr3proba::TaskSurv]
    #' @param verbose Write log messages or not? Default: TRUE.
    #' @param store_archive Whether to also store the [ArchiveFSelect][mlr3fselect::ArchiveFSelect]
    #' archive object created by `AutoFSelector`, for debugging purposes.
    #' Default: FALSE.
    #'
    #' @return a tibble with columns:
    #'    - `lrn_id` => which learner was used
    #'    - `iter` => in this particular execution of RFE (out of a total `repeats`)
    #'    - `selected_features` => the best feature subset selected by RFE
    #'    - `nfeatures` => how many features were selected
    #'    - `score` => the performance score of the chosen best feature subset
    #'     (depends on the `msr_id`)
    #'    - `ArchiveFSelect` => additional object for validation and checking
    #'    (**not included** by default due to large size)
    run = function(task, verbose = TRUE, store_archive = FALSE) {
      # task
      assert_task(task)
      self$task_id = task$id

      # RSF learners
      learners = private$.get_lrns()

      # mtry_ratio callback
      mr = self$mtry_ratio
      n = length(task$feature_names) # total #features

      if (self$adaptive_mr) {
        mr_clbk = callback_fselect(id = 'mtry',
          on_eval_after_design = function(callback, context) {
            # #features in RFE subset
            nfeats = length(context$design$task[[1]]$feature_names)

            # adaptive mtry.ratio formula
            mtry.ratio = mr^(log(nfeats)/log(n))

            context$design$learner[[1]]$param_set$set_values(mtry.ratio = mtry.ratio)
          })
      } else {
        mr_clbk = callback_fselect(id = 'empty')
      }

      # performance measure
      if (self$msr_id == 'oob_error') {
        measure = msr('oob_error')
      } else {
        measure = bench_msrs()[id == self$msr_id]$measure[[1L]]
      }

      # eFS
      # make grid of RFE runs
      rfe_grid = full_join(
        tibble(lrn_id = self$lrn_ids),
        tibble(iter   = 1:self$repeats),
        by = character()
      )
      if (verbose) {
        message('A total of ', nrow(rfe_grid), ' RFE runs')
      }

      result = list()
      for (index in 1:nrow(rfe_grid)) {
        lrn_id = rfe_grid[index,]$lrn_id
        learner = learners[[lrn_id]]$reset() # un-train
        iter = rfe_grid[index,]$iter

        if (verbose) {
          message('### Learner: ', learner$label, ' (', iter, '/',
            self$repeats, ') ###')
        }

        at = AutoFSelector$new(
          learner = learner,
          resampling = self$resampling,
          measure = measure,
          terminator = trm('none'),
          fselector = fs('rfe', n_features       = self$n_features,
                                feature_fraction = self$feature_fraction),
          store_models = TRUE, # necessary for RFE
          callbacks = mr_clbk
        )

        at$train(task)

        selected_features = at$fselect_instance$result_feature_set
        nfeatures = length(selected_features)
        score = at$archive$best()[[measure$id]]

        result[[index]] = tibble(
          lrn_id = lrn_id,
          iter   = iter,
          selected_features = list(selected_features),
          nfeatures = nfeatures,
          score = score
        )

        if (store_archive) {
          result[[index]] = tibble(
            result[[index]],
            archive = list(at$archive)
          )
        }
        # store/save intermediate result[[index]] somehow?
      }

      res = bind_rows(result)
      self$result = res
      invisible(res)
    }
  ),
  private = list(
    # Returns a list of survival random forest learners with the appropriate
    # hyperparameters
    .get_lrns = function() {
      rsf_lrns =
        SurvLPS$
        new(ids = self$lrn_ids, nthreads_rsf = self$nthreads_rsf)$
        lrns()

      # Apply the following parameter list to each RSF learner
      param_list = list(
        num.trees = self$num_trees,
        mtry.ratio = self$mtry_ratio,
        min.node.size = 3, # default for survival (RSF)
        importance = 'permutation' # don't expose it to the user (yet)
      )

      for (learner in rsf_lrns) {
        learner$param_set$values = param_list
      }

      rsf_lrns
    }
  )
)
