#' @title Ensemble Feature Selection (eFS)
#'
#' @description This class stores the configuration and result of the proposed
#' hybrid ensemble feature selection approach, as well as the function (`run()`)
#' to execute it.
#'
#' We use a couple of different **random survival forests** (RSFs) as learners
#' and a **recursive feature elimination algorithm** (RFE, wrapper-based feature
#' selection) to find the most predictive feature subset in the given task
#' (dataset) for each learner.
#' The dataset can be subsampled before the execution of each RFE run.
#' This process repeats a number of times for each RSF learner and therefore
#' returns different best feature subsets for further analysis and processing
#' (e.g. finding robust features).
#'
#' @examples
#' library(mlr3proba)
#' library(mlr3fselect)
#' library(mlr3pipelines)
#' library(mlr3extralearners)
#' library(ggplot2)
#' set.seed(42)
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
#' # create eFS object
#' efs = eFS$new(lrn_ids = c('rsf_logrank', 'rsf_cindex'), nthreads_rsf = 2,
#'   feature_fraction = 0.6, n_features = 1, repeats = 3, subsample_ratio = 0.8
#' )
#'
#' # useful info for RFE (feature subset sizes, mtry used)
#' efs$rfe_info(task)
#'
#' # Execute ensemble feature selection
#' # In each RFE iteration, every RSF will be trained on all the features
#' # available and the out-of-bag error will be used to assess predictive
#' # performance (1 - Harrell's Cindex)
#' efs$run(task, verbose = TRUE)
#'
#' # Get result in a tibble format
#' res = efs$result
#' res
#'
#' # Get consensus performance score and number of features
#' median(res$score) # 1 - C-index
#' median(res$nfeatures)
#'
#' # get frequency selection stats (per learner and consensus)
#' fss = efs$fs_stats()
#' fss$consensus # ranked consensus features
#'
#' # Barplots: Feature Selection Frequency
#' efs$ffs_plot(lrn_id = 'consensus', title = 'RSF consensus')
#' efs$ffs_plot(lrn_id = 'rsf_logrank', title = 'RSF logrank')
#'
#' # Performance plot (per RSF learner)
#' efs$res_plot(msr_label = 'OOB (1 - C-index)')
#'
#' # Number of selected features plot (per RSF learner)
#' efs$res_plot(type = 'nfeat', ylimits = c(3,7))
#'
#' # Stability plot
#' efs$res_plot(type = 'stab', task = task)
#'
#' @export
eFS = R6Class('eFS',
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

    #' @field subsample_ratio (`double(1)`)\cr
    #' Ratio of observations to be selected in each RFE repeat
    subsample_ratio = NULL,

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
    #' Default is 'oob_error'.
    #'
    #' @param resampling [Resampling][mlr3::Resampling] for the
    #' [AutoFSelector][mlr3fselect::AutoFSelector].
    #' Default resampling is [insample][mlr3::ResamplingInsample].
    #'
    #' @param repeats Number of times to run the RFE algorithm on each RSF
    #' learner. Defaults to 100.
    #'
    #' @param subsample_ratio Ratio of observations to use in each RFE run.
    #' The selection is done randomly each time and is stratified on the `status`
    #' target variable (so that the censoring distribution is the same as the
    #' original dataset).
    #' This option enables the ensemble feature selection algorithm to be used
    #' on different subsets of a given dataset.
    #' Defaults to 0.9.
    #' A value of 1 means that we always use all observations and no subsampling
    #' is performed.
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
    #' @details
    #' - By default, `msr_id` is `oob_error`:
    #'    - Each RSF provides (1 - Harrell's C-index) as Out-Of-Bag error
    #'    - An [insample resampling][mlr3::ResamplingInsample] is used in that
    #'    case in order to use all training data for bootstrapping in the RSFs.
    #' - `mtry` used in RSFs:
    #'    - During the execution of the RFE algorithm, progressively smaller
    #'    feature subsets are used, as dictated by the `feature_fraction`
    #'    parameter.
    #'    - For each feature subset, the RSF learner tries by default
    #'    `ceiling(sqrt(#features-in-subset))` features as candidate variables
    #'    for splitting (`mtry`).
    #'    So `mtry` values also get progressively smaller.
    #'    See `rfe_info()` for more details.
    initialize = function(
      lrn_ids = self$supported_lrn_ids(),
      msr_id = 'oob_error',
      resampling = mlr3::rsmp('insample'),
      repeats = 100, subsample_ratio = 0.9,
      n_features = 2, feature_fraction = 0.8,
      nthreads_rsf = unname(parallelly::availableCores()), num_trees = 250
    ) {
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
          Available ids: ', paste0(supp_msr_ids, collapse = ', ')))
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

      if (subsample_ratio >= 0 && subsample_ratio <= 1)
        self$subsample_ratio = subsample_ratio
      else
        stop('\'subsample_ratio\' needs to be between 0 and 1 (inclusive)')

      if (n_features > 0)
        self$n_features = n_features
      else
        stop('\'n_features\' needs be larger than 0')

      if (feature_fraction >= 0 && feature_fraction < 1)
        self$feature_fraction = feature_fraction
      else
        stop('\'feature_fraction\' needs to be between 0 (inclusive) and 1
          (non-inclusive)')

      # RSF parameters
      if (nthreads_rsf > 0)
        self$nthreads_rsf = nthreads_rsf
      else
        stop('\'nthreads_rsf\' needs to be larger than 0')

      if (num_trees > 0)
        self$num_trees = num_trees
      else
        stop('\'num_trees\' needs to be larger than 0')
    },

    #' @description Returns a vector of internal ids corresponding to the
    #' random survival forest learners that can be used in method `run()`
    supported_lrn_ids = function() {
      s = SurvLPS$new()
      ids = s$lrn_ids()
      ids[grepl(pattern = 'rsf', ids)]
    },

    #' @description Returns the feature subset sizes that the RFE algorithm
    #' will use and the corresponding `mtry` values for the
    #' RSF learners.
    #' These depend on the total number of features of the given task,
    #' as well as the `n_features` and `feature_fraction`
    #' parameters (initialized upon class construction).
    #'
    #' @details
    #' The subset sizes are calculated from the respective
    #' [code](https://github.com/mlr-org/mlr3fselect/blob/HEAD/R/FSelectorRFE.R)
    #' (see `rfe_subsets()` function, `ml3fselect` v0.10.0).
    #'
    #' @param task [mlr3::Task]
    #' @return a [tibble] with columns `subset_size` and `mtry`
    rfe_info = function(task) {
      n = length(task$feature_names) # number of total features
      n_features = self$n_features
      feature_fraction = self$feature_fraction

      subset_sizes = unique(floor(cumprod(
        c(n, rep(feature_fraction, log(n_features / n) / log(feature_fraction))))
      ))

      tibble(subset_size = subset_sizes) %>%
        mutate(mtry = ceiling(sqrt(subset_size)))
    },

    #' @description Runs the ensemble feature selection on the given task.
    #'
    #' @details For every RSF learner, we run `repeats` times the RFE algorithm
    #' using a properly constructed [AutoFSelector][mlr3fselect::AutoFSelector]).
    #' From each `repeat`-ition we get the best feature subset.
    #' The aggregated result is invisibly returned as a tibble and is available
    #' also in the `result` field of this object.
    #'
    #' @param task [TaskSurv][mlr3proba::TaskSurv]
    #' @param verbose Write log messages or not? Default: TRUE.
    #' @param store_archive Whether to also store the
    #' [ArchiveFSelect][mlr3fselect::ArchiveFSelect]
    #' archive object created by `AutoFSelector`, for debugging purposes.
    #' Default: FALSE.
    #'
    #' @return a tibble with columns:
    #'    - `lrn_id` => which learner was used
    #'    - `iter` => which iteration of RFE
    #'    - `selected_features` => the best feature subset selected by RFE
    #'    - `nfeatures` => how many features were selected
    #'    - `score` => the performance score of the best chosen feature subset
    #'     (depends on the `msr_id`)
    #'    - `ArchiveFSelect` => additional object for validation and checking
    #'    (**not included** by default due to large size)
    run = function(task, verbose = TRUE, store_archive = FALSE) {
      # task
      assert_task(task)
      self$task_id = task$id

      # RSF learners
      learners = private$.get_lrns()

      # performance measure
      if (self$msr_id == 'oob_error') {
        measure = msr('oob_error')
      } else {
        measure = bench_msrs()[id == self$msr_id]$measure[[1L]]
      }

      # eFS
      # make grid of RFE runs
      rfe_grid = dplyr::cross_join(
        tibble(lrn_id = self$lrn_ids),
        tibble(iter   = 1:self$repeats)
      )

      # initial message
      if (verbose) {
        message(length(self$lrn_ids), ' RSF learner(s) x ',
                self$repeats, ' repeats')
        message('A total of ', nrow(rfe_grid), ' RFE runs')
      }

      result = list()
      for (index in 1:nrow(rfe_grid)) {
        lrn_id = rfe_grid[index,]$lrn_id
        learner = learners[[lrn_id]]$reset() # un-train
        iter = rfe_grid[index,]$iter

        if (verbose) {
          message('### Learner: ', learner$id, ' (', iter, '/',
            self$repeats, '), Iter: ', index, '/', nrow(rfe_grid),
            ' (', round(100*index/nrow(rfe_grid), digits = 1), '%)')
        }

        at = AutoFSelector$new(
          learner = learner,
          resampling = self$resampling,
          measure = measure,
          terminator = trm('none'),
          fselector = fs('rfe', n_features = self$n_features,
            feature_fraction = self$feature_fraction),
          store_models = store_archive # hacked :)
        )

        part = partition(task, ratio = self$subsample_ratio, stratify = TRUE)

        # envelop in try/catch in case something goes wrong
        tryCatch(
          {
            at$train(task, row_ids = part$train)

            selected_features = at$fselect_instance$result_feature_set
            nfeatures = length(selected_features)
            score = at$archive$best()[[measure$id]]

            result[[index]] = tibble(
              lrn_id = lrn_id,
              iter   = iter,
              selected_features = list(selected_features),
              nfeatures = nfeatures,
              score = score,
              archive = base::switch(store_archive, list(at$archive))
            )
          }, error = function(e) {
            message('### Error: ', e$message)
          }
        )
      }
      res = bind_rows(result)

      self$result = res
      invisible(res)
    },

    #' @description Runs the ensemble feature selection on the given task in
    #' parallel
    #'
    #' @details Same as `run()` but with the possibility to parallelize the RFE
    #' runs.
    #' Should be used with `progressr` and `future` libraries and it's a bit
    #' *experimental* as this parallelization can conflict with the
    #' parallelization of some RSF learners.
    #' See examples.
    #'
    #' @param task [TaskSurv][mlr3proba::TaskSurv]
    #' @param verbose Write log messages or not? Default: TRUE.
    #' @param store_archive Whether to also store the
    #' [ArchiveFSelect][mlr3fselect::ArchiveFSelect]
    #' archive object created by `AutoFSelector`, for debugging purposes.
    #' Default: FALSE.
    #'
    #' @return a tibble with columns, same as in `run()`.
    run_parallel = function(task, verbose = TRUE, store_archive = FALSE) {
      # task
      assert_task(task)
      self$task_id = task$id

      # RSF learners
      learners = private$.get_lrns()

      # performance measure
      if (self$msr_id == 'oob_error') {
        measure = msr('oob_error')
      } else {
        measure = bench_msrs()[id == self$msr_id]$measure[[1L]]
      }

      # eFS
      # make grid of RFE runs
      rfe_grid = dplyr::cross_join(
        tibble(lrn_id = self$lrn_ids),
        tibble(iter   = 1:self$repeats)
      )

      # initial message
      if (verbose) {
        message(length(self$lrn_ids), ' RSF learner(s) x ',
          self$repeats, ' repeats')
        message('A total of ', nrow(rfe_grid), ' RFE runs')
      }

      # future and progressr prerequisites
      is_sequential = inherits(plan(), 'sequential')
      stdout = if (is_sequential) NA else TRUE
      lgr_thres = mlr3misc::map_int(mlr_reflections$loggers, 'threshold')
      pb = progressr::progressor(steps = nrow(rfe_grid))

      res = future.apply::future_lapply(1:nrow(rfe_grid), function(index) {
        # restore logger thresholds to avoid unnecessary logging
        for (package in names(lgr_thres)) {
          logger = lgr::get_logger(package)
          threshold = lgr_thres[package]
          logger$set_threshold(threshold)
        }

        lrn_id = rfe_grid[index,]$lrn_id
        learner = learners[[lrn_id]]$reset() # un-train
        iter = rfe_grid[index,]$iter

        if (verbose) {
          pb(sprintf('%s (%i/%i)', learner$id, iter, self$repeats))
        }

        at = AutoFSelector$new(
          learner = learner,
          resampling = self$resampling,
          measure = measure,
          terminator = trm('none'),
          fselector = fs('rfe', n_features = self$n_features,
            feature_fraction = self$feature_fraction),
          store_models = store_archive # hacked :)
        )

        part = partition(task, ratio = self$subsample_ratio, stratify = TRUE)

        # envelop in try/catch in case something goes wrong
        tryCatch(
          {
            at$train(task, row_ids = part$train)

            selected_features = at$fselect_instance$result_feature_set
            nfeatures = length(selected_features)
            score = at$archive$best()[[measure$id]]

            res_tbl = tibble(
              lrn_id = lrn_id,
              iter   = iter,
              selected_features = list(selected_features),
              nfeatures = nfeatures,
              score = score,
              archive = base::switch(store_archive, list(at$archive)),
            )
          }, error = function(e) {
            message('### Error: ', e$message)
          }
        )
        res_tbl
      }, future.seed = TRUE, future.conditions = 'message', future.stdout = stdout
      ) %>% bind_rows()

      self$result = res
      invisible(res)
    },

    #' @description Frequency selection statistics (feature ranking).
    #' This function uses the best feature sets found by RFE
    #' and creates one table per RSF learner with the features in descending
    #' order of selection frequency.
    #' A consensus frequency results table across all feature sets generated
    #' by all RSFs during the RFE is also returned.
    #'
    #' @return A list of `tibble`s with columns:
    #'    - `feat_name` => feature name
    #'    - `times` => how many times a feature was chosen in the best feature
    #' subsets across all RFE runs?
    #'    - `freq` => selection frequency
    fs_stats = function() {
      result = self$result
      if (is.null(result)) stop('Need to execute run() first')

      # 1 frequency table per RSF learner
      freq_list = list()
      lrn_ids = unique(result$lrn_id)
      for (id in lrn_ids) {
        res_subset = result %>% filter(lrn_id == id)

        sf_list = res_subset$selected_features

        n = length(sf_list)
        res = sort(table(unlist(sf_list)), decreasing = TRUE)
        times = as.vector(unname(res))
        freq_list[[id]] =
          tibble(feat_name = names(res), times = times, freq = times/n)
      }

      # 1 frequency table for all RSFs together (consensus)
      if (length(lrn_ids) > 1) {
        sf_list = result$selected_features

        n = length(sf_list) # total number of best subsets
        res = sort(table(unlist(sf_list)), decreasing = TRUE)
        times = as.vector(unname(res))
        freq_list[['consensus']] =
          tibble(feat_name = names(res), times = times, freq = times/n)
      }

      freq_list
    },

    #' @description Stability assessment of the ensemble feature selection.
    #' Currently two stability metrics are supported via the [stabm] R package,
    #' namely Jaccard and Nogueira's measure.
    #' Stability is assessed on the feature sets produced per inidividual RSF
    #' learner used, as well as on all of them combined (consensus stability).
    #'
    #' @param task [TaskSurv][mlr3proba::TaskSurv] that was used in `run()`
    #' @param stab_metrics vector of stability metrics
    #'
    #' @return a tibble with a `lrn_id` column and as many columns as the
    #' stability measures chosen.
    #' The `lrn_id` column includes the `consensus` stability when the feature
    #' set from all RSFs are merged (and when more than one RSF learner was
    #' used).
    stab = function(task, stab_metrics = NULL) {
      result = self$result
      if (is.null(result)) stop('Need to execute run() first')
      if (self$repeats == 1) stop('Need at least 2 RFE `repeats`, i.e. 2 feature
        sets, to measure stability')

      # which metrics to use (default => all)
      avail_metrics = private$.get_stab_metrics()
      if (is.null(stab_metrics)) {
        stab_metrics = avail_metrics
      } else {
        stopifnot(all(stab_metrics %in% avail_metrics))
      }

      # Nogueira measure needs the task to get the number of features (p)
      if ('nogueira' %in% stab_metrics)
        assert_task(task)

      lrn_ids = unique(result$lrn_id)
      stab_list = list()
      for (id in lrn_ids) {
        res_subset = result %>% filter(lrn_id == id)
        sf_list    = res_subset$selected_features

        jaccard = NULL
        if ('jaccard' %in% stab_metrics)
          jaccard = stabm::stabilityJaccard(features = sf_list,
            correction.for.chance = 'none')

        nogueira = NULL
        if ('nogueira' %in% stab_metrics)
          nogueira = stabm::stabilityNogueira(features = sf_list,
            p = length(task$feature_names))

        stab_list[[id]] = tibble(lrn_id = id, jaccard = jaccard, nogueira = nogueira)
      }

      # consensus stability from all feature sets generated by all RSFs
      if (length(lrn_ids) > 1) {
        sf_list = result$selected_features

        jaccard = NULL
        if ('jaccard' %in% stab_metrics)
          jaccard = stabm::stabilityJaccard(features = sf_list,
            correction.for.chance = 'none')

        nogueira = NULL
        if ('nogueira' %in% stab_metrics)
          nogueira = stabm::stabilityNogueira(features = sf_list,
            p = length(task$feature_names))

        stab_list[['consensus']] = tibble(lrn_id = 'consensus',
          jaccard = jaccard, nogueira = nogueira)
      }

      bind_rows(stab_list)
    },

    #' @description Feature selection frequency barplot
    #' @param lrn_id One of RSF learner ids that was used during initialization
    #' or 'consensus' (default)
    #' @param top_n plot only the `n` features with the higher selection
    #' frequency
    #' @param title title of barplot (if NULL, the `lrn_id` is used)
    ffs_plot = function(lrn_id = 'consensus', top_n = 10, title = NULL) {
      fss = self$fs_stats()
      lrn_ids = names(fss)

      if (!lrn_id %in% lrn_ids)
        stop(lrn_id, ' not included in ', paste0(lrn_ids, collapse = ', '))

      if (is.null(title)) {
        title = lrn_id
      }

      p = fss[[lrn_id]] %>%
        slice(1:top_n) %>%
        mutate(feat_name = forcats::fct_reorder(feat_name, times, .desc = FALSE)) %>%
        ggplot(aes(x = feat_name, y = freq)) +
        geom_bar(stat = 'identity', fill = '#377EB8', show.legend = FALSE) +
        theme_bw(base_size = 14) +
        labs(x = 'Feature name', y = 'Selection Frequency', title = title) +
        scale_y_continuous(labels = scales::label_percent()) +
        coord_flip()

      p
    },

    #' @description Produce various plots of the eFS results
    #'
    #' @param type Type of plot to produce. Can be `perf`, `nfeat` or `stab`
    #' @param msr_label Measure name as the y-axis label (default: `self$msr_id`)
    #' @param ylimits Passing on to `ylim()` (y-axis limits)
    #' @param title Title for plot (default is `self$task_id`, the task that was
    #' used during `run()`)
    #' @param include_legend By default no legend is included
    #' @param task this is needed for `type` = `stab`. Should be the same used
    #' when executing `run()`
    #' @param stab_metrics this is needed for `type` = `stab`. By default all
    #' stability metrics are used.
    #'
    #' @details 3 types of plots can be produced:
    #' 1. **Performance boxplot**
    #' 2. **Number of features boxplot**
    #' 3. **Stability metric barplots (1 per metric)**, where also the consensus
    #' feature set is included as a comparison category. The results are ordered
    #' according to the values of the first metric used in `stab()` ('jaccard'
    #' by default) decided by the value of `stab_metrics`.
    res_plot = function(type = 'perf', msr_label = NULL, ylimits = NULL,
      title = NULL, include_legend = FALSE, task = NULL, stab_metrics = NULL) {
      # checks
      result = self$result
      if (is.null(result)) stop('Need to execute run() first')

      if (is.null(msr_label)) {
        msr_label = self$msr_id
      }

      if (is.null(title)) {
        title = self$task_id
      }

      if (type == 'stab' && is.null(task)) {
        stop('Provide task when type = \'stab\'')
      }

      if (type == 'perf') {
        p = result %>%
          ggplot(aes(x = lrn_id, y = score, fill = lrn_id)) +
          geom_boxplot() +
          labs(y = msr_label, x = NULL, title = title)
      } else if (type == 'nfeat') {
        p = result %>%
          ggplot(aes(x = lrn_id, y = nfeatures, fill = lrn_id)) +
          geom_boxplot() +
          labs(y = 'Number of selected features', x = NULL, title = title)
      } else if (type == 'stab') {
        stab_tbl = self$stab(task, stab_metrics)

        # pick the first stability measure
        sm = colnames(stab_tbl)[2]
        # check it's part of the available ones
        stopifnot(sm %in% private$.get_stab_metrics())

        stab_tbl = stab_tbl %>%
          mutate(lrn_id = fct_reorder(lrn_id, .data[[sm]], .desc = TRUE)) %>%
          pivot_longer(cols = !c('lrn_id'), names_to = 'stab_measure')

        p = stab_tbl %>%
          ggplot(aes(x = lrn_id, y = value, fill = lrn_id)) +
          geom_bar(position = 'dodge', stat = 'identity') +
          facet_wrap(~stab_measure) +
          labs(x = NULL, y = 'Similarity score', title = title)
      } else {
        stop('Wrong type, can be only: perf, nfeat or stab')
      }

      p = p +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5))

      if (!is.null(ylimits)) {
        p = p + ylim(ylimits)
      }

      if (!include_legend) {
        p = p + theme(legend.position = 'none')
      }

      p
    },

    #' @description
    #' Opens the help page for this object.
    help = function() {
      mlr3misc::open_help('survmob::eFS')
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

      # Apply the following parameter list to each RSF learner (ranger)
      param_list = list(
        num.trees = self$num_trees,
        min.node.size = 3, # default for survival (RSF)
        importance = 'permutation' # don't expose it to the user (yet)
      )

      # Apply the following parameter list to each Oblique RSF learner
      param_list_oblique = list(
        n_tree = self$num_trees,
        leaf_min_obs = 3, # default for survival (RSF)
        importance = 'anova',
        attach_data = TRUE # needed for importance
      )

      for (learner in rsf_lrns) {
        if (startsWith(x = learner$id, prefix = 'Oblique')) {
          learner$param_set$values = param_list_oblique
        } else {
          learner$param_set$values = param_list
        }
      }

      rsf_lrns
    },

    # convenience function to get ids for the stability measures
    .get_stab_metrics = function() {
      c('jaccard', 'nogueira')
    }
  )
)
