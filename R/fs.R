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
#'   feature_fraction = 0.6, n_features = 1, repeats = 3
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
#' fss$consensus
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
      repeats = 100, n_features = 2, feature_fraction = 0.8,
      nthreads_rsf = parallelly::availableCores(), num_trees = 250
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
      rfe_grid = full_join(
        tibble(lrn_id = self$lrn_ids),
        tibble(iter   = 1:self$repeats),
        by = character()
      )
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
            self$repeats, '), Iter: ', index, '/', nrow(rfe_grid))
        }

        at = AutoFSelector$new(
          learner = learner,
          resampling = self$resampling,
          measure = measure,
          terminator = trm('none'),
          fselector = fs('rfe', n_features       = self$n_features,
                                feature_fraction = self$feature_fraction),
          store_models = store_archive # hacked :)
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
        geom_bar(stat = "identity", fill = '#377EB8', show.legend = FALSE) +
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

#' @title Subset tasks class
#'
#' @description Subset (survival) tasks to the most stable features identified
#' by the ensemble feature approach ([eFS]).
#' For a summary table of the number of features that can be selected via
#' various methods for cutoffs, etc. see `feat_stats()` and for actually
#' subsetting the tasks, see `subset_tasks()`.
#'
#' @export
TskSub = R6Class('FSTaskSubsettor',
  public = list(
    #' @field nfeats (`numeric(1)`)\cr
    #' Top N features
    nfeats = NULL,
    #' @field perc (`numeric(1)`)\cr
    #' Top N% (percentage cutoff)
    perc = NULL,
    #' @field cutoff (`numeric(1)`)\cr
    #' Consensus Selection Frequency cutoff
    cutoff = NULL,
    #' @field max_nfeats (`numeric(1)`)\cr
    #' Max #features
    max_nfeats = NULL,

    #' @description Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param nfeats Top N features
    #' @param perc Top N% (percentage cutoff)
    #' @param cutoff Selection Frequency cutoff
    #' @param max_nfeats Max #features
    initialize = function(nfeats = 10, perc = 0.05, cutoff = 0.6,
      max_nfeats = 15) {
      # basic checks
      stopifnot(nfeats > 0)
      stopifnot(perc >= 0, perc <= 1)
      stopifnot(cutoff >= 0, cutoff <= 1)
      stopifnot(max_nfeats > 0)

      self$nfeats = nfeats
      self$perc = perc
      self$cutoff = cutoff
      self$max_nfeats = max_nfeats
    },

    #' @description Summary table with the number of features that can be
    #' selected based on available methods. See `subset_tasks()` for more
    #' details.
    #'
    #' @param eFSlist list of [eFS] objects
    #' @param tasks list of [TaskSurv][mlr3proba::TaskSurv] (1-1 correspondence with
    #' each eFS object)
    #'
    #' @details Main sanity checks:
    #' - There has to be 1-1 correspondence between the tasks and the [eFS]
    #' classes
    #' - The tasks have to be different (different id at least :)
    #' - The [eFS] objects must have executed the method `$run()` so that there
    #' is a `result` tibble
    #' - No task should have more than `max_nfeats` (see `initialize()`) features
    #'
    #' @return A [tibble][dplyr::tibble], with columns:
    #' - `task_id` => different task id per row
    #' - `features` => features names to select from (different vector per row,
    #' ordered by consensus selection frequency)
    #' - `top_n` => #features to select based on `top_n` method
    #' - `top_nperc` => #features to select based on `top_nperc` method
    #' - `gt_freq` => #features to select based on `gt_freq` method
    #' - `score_median` => median performance score out of all feature subsets
    #' used in eFS
    #' - `nfeats_median` => median #features selected out of all feature subsets
    #' used in eFS
    #' - In case of a C-index measure, a combination measure and its normalized
    #' version are included (`comb_score`, `comb_score_norm`) as well as a
    #' column named `optim_nfeats`, see `subset_tasks()` for more info.
    #'
    summary_stats = function(eFSlist, tasks) {
      private$.check(eFSlist, tasks)

      fs_tbl = lapply(eFSlist, function(efs) {
        fss = efs$fs_stats()
        if (!'consensus' %in% names(fss)) {
          stop('Consensus `fs_stats` not found, eFS was run on single RSF?')
        }
        cons_fss = fss$consensus
        res = efs$result
        task_id = efs$task_id
        msr_id = efs$msr_id

        # top N
        n1 = self$nfeats

        # top N%
        total_nfeats = nrow(cons_fss)
        n2 = floor(self$perc * total_nfeats)

        # greater than `cutoff` selection frequency
        n3 = cons_fss %>% filter(freq > self$cutoff) %>% nrow()

        # median score
        score_median = median(res$score)

        # median #features selected
        nfeats_median = floor(median(res$nfeatures))

        tibb = tibble::tibble(
          task_id = task_id,
          features = list(cons_fss$feat_name),
          top_n = n1,
          top_nperc = n2,
          gt_freq = n3,
          score_median = score_median,
          nfeats_median = nfeats_median
        )

        # add extra score if measure is C-index-like
        if (msr_id %in% c('harrell_c', 'uno_c', 'oob_error')) {
          cindex = ifelse(msr_id == 'oob_error', 1 - score_median, score_median)

          # higher C-index with more features gets a higher score
          comb_score = cindex * nfeats_median
          tibb = tibb %>% mutate(comb_score = comb_score)
        }

        tibb
      }) %>% dplyr::bind_rows()

      if ('comb_score' %in% colnames(fs_tbl)) {
        fs_tbl = fs_tbl %>%
          # normalize score
          mutate(comb_score_norm = comb_score/max(comb_score)) %>%
          # calculate #features based on the normalized score and `max_nfeats`
          mutate(optim_nfeats = floor(self$max_nfeats * comb_score_norm))
      }

      fs_tbl
    },

    #' @description Subset tasks (i.e. reduce the number of features) based on
    #' the corresponding consensus feature selection results (eFS)
    #'
    #' @param eFSlist list of [eFS] objects
    #' @param tasks list of [TaskSurv][mlr3proba::TaskSurv] (1-1 correspondence with
    #' each eFS object)
    #' @param method method that decides which features to select from each task.
    #' Available options are `top_n`, `top_nperc`, `gt_freq` or `optim_nfeats`.
    #'
    #' @details To decide which features to select for a particular task, we use
    #' the consensus frequency selection statistics from the corresponding [eFS]
    #' object (features are ranked from most frequently selected to least).
    #' This means that at least two RSFs were used during the execution of
    #' the eFS algorithm.
    #'
    #' Available options for `method`:
    #' - `top_n` => Select the top N (`nfeats`) features.
    #' After subsetting, all tasks will have the same number of features.
    #' - `top_nperc` => Select the top N% (`perc`) features.
    #' This depends on the total number of features selected from each task by
    #' the respective eFS, so the number of the top 5% features for one task
    #' might be different than the number of the top 5% features in another task.
    #' - `gt_freq` => Select features whose selection frequency in the eFS
    #' was above X% (X is the `cutoff`).
    #' - `optim_nfeats` => The number of features for each task is decided based
    #' on the median performance and median number of selected features in the
    #' eFS (**the larger** the both of these are, the more features will be
    #' assigned to a task).
    #' This method can be used only when all the eFS were executed with an
    #' optimization measure that's a **C-index** or **OOB = 1 - C-index**.
    #' So, one task will get a normalized `comb_score_norm` equal to 1 (larger
    #' product of C-index * #features) and this task will get `max_nfeats` (from
    #' class initialization) - the rest of the tasks will get less #features.
    #'
    #' @return list of [TaskSurv][mlr3proba::TaskSurv]s with only the most
    #' frequently selected features based on the consensus eFS results
    #' (how many features is decided by the `method` and the result of
    #' `summary_stats()`)
    #'
    subset_tasks = function(eFSlist, tasks, method = 'top_n') {
      # check 'method'
      avail_methods = c('top_n', 'top_nperc', 'gt_freq', 'optim_nfeats')
      if (!method %in% avail_methods) {
        stop('Method not in available methods: ',
          paste0(avail_methods, collapse = ', '))
      }

      # more checks are done in `.check()`
      fs_tbl = self$summary_stats(eFSlist, tasks)

      if (method == 'optim_nfeats' && (!'optim_nfeats' %in% colnames(fs_tbl))) {
        stop('No column `optim_nfeats` in summary fs stats table')
      }

      lapply(tasks, function(task) {
        # how many task features based on the method chosen should we select?
        nfeats =
          fs_tbl %>%
          filter(task_id == task$id) %>%
          pull(method)

        # cannot subset to 0 features!
        if (nfeats < 1) return(task)

        # get all the available feature names (from consensus eFS)
        features =
          fs_tbl %>%
          filter(task_id == task$id) %>%
          pull(features) %>%
          `[[`(1)

        # subset the features
        features = features[1:nfeats]

        minimize_backend(task$clone()$select(cols = features))
      })
    },

    #' @description
    #' Opens the help page for this object.
    help = function() {
      mlr3misc::open_help('survmob::TskSub')
    }
  ),
  private = list(
    # @param eFSlist list of [eFS] objects
    # @param tasks list of [TaskSurv][mlr3proba::TaskSurv]
    .check = function(eFSlist, tasks) {
      stopifnot(length(eFSlist) == length(tasks))

      # check that task ids are different (from tasks)
      task_ids = mlr3misc::map_chr(tasks, `[[`, 'id')
      if (length(task_ids) > 1 && length(unique(task_ids)) != length(task_ids)) {
        stop('Found same task id among given tasks.
              Did you run eFS on the same task more than once? You shouldn\'t!')
      }

      # check that task ids are not NULL in eFS objects
      for (efs in eFSlist) {
        if (is.null(efs$task_id)) {
          stop('Empty task id in eFS. Probably `eFS$run()` has not been executed')
        }
      }

      # check that task ids match (tasks <=> eFSlist)
      stopifnot(mlr3misc::map_chr(eFSlist, `[[`, 'task_id') == task_ids)

      # check that all measure ids used are the same in eFS objects
      msr_ids = mlr3misc::map_chr(eFSlist, `[[`, 'msr_id')
      stopifnot(length(unique(msr_ids)) == 1)

      # check that we don't have any NULL results (`$run()` not executed)
      res_list = mlr3misc::map(eFSlist, `[[`, 'result')
      lgl_vec = sapply(res_list, is.null)
      if (any(lgl_vec)) {
        stop('Found NULL $result')
      }

      # check that `max_nfeats` does not exceed any of the tasks #nfeats
      if (!any(self$max_nfeats <
          sapply(mlr3misc::map(tasks, `[[`, 'feature_names'), length))) {
        stop('There is a task that has more features than `max_nfeats`!')
      }

      # check that `nfeats` (top N) does not exceed any of the tasks #nfeats
      if (!any(self$nfeats <
          sapply(mlr3misc::map(tasks, `[[`, 'feature_names'), length))) {
        stop('There is a task that has more features than `nfeats`!')
      }
    }
  )
)
