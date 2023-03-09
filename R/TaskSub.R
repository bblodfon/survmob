#' @title Subset tasks class
#'
#' @description Subset (survival) tasks to the most stable features identified
#' by the ensemble feature approach ([eFS]).
#' For a summary table of the number of features that can be selected via
#' various methods for cutoffs, etc. see `feat_stats()` and for actually
#' subsetting the tasks, see `subset_tasks()`.
#'
#' @export
TaskSub = R6Class('FSTaskSubsettor',
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
      mlr3misc::open_help('survmob::TaskSub')
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
