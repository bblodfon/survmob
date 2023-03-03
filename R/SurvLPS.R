#' @title Survival Learners and Parameter Spaces
#'
#' @description Convenient [R6][R6::R6Class] class that offers methods to retrieve
#' survival learners (to be used for benchmarking) as well as suggested
#' hyperparameter spaces for those learners.
#'
#' @examples
#' library(paradox)
#' library(mlr3proba)
#' library(mlr3extralearners)
#'
#' s = SurvLPS$new(nthreads_rsf = 4)
#'
#' # All available learner ids
#' s$lrn_ids()
#'
#' # Get all parameter spaces
#' pss = s$pss()
#'
#' # CoxPH doesn't have hyperparameters to tune
#' pss$coxph
#' # random survival forest hyperparameter space
#' pss$rsf_logrank
#' # XGBoost hyperparameter space with early stopping
#' pss$xgboost_cox_early$trafo(x = list(
#'   XGBoostCox.nrounds = 150,
#'   XGBoostCox.eta = -5,
#'   XGBoostCox.max_depth = 5))
#' # XGBoost AFT hyperparameter space with additional
#' # regularization parameters
#' pss$xgboost_aft_reg
#'
#' # Get only 2 learner objects
#' s$initialize(ids = c('coxnet', 'coxboost'))
#' s$lrns()
#'
#' # Get a convenient data table format
#' dt = s$lrn_tbl()
#' dt
#'
#' @export
SurvLPS = R6Class('SurvLPS',
  public = list(
    #' @field nthreads_rsf (`int(1)`)\cr
    #' Number of cores to use in the random forest survival learners
    nthreads_rsf = NULL,

    #' @field nthreads_xgb (`int(1)`)\cr
    #' Number of cores to use in the xgboost survival learners
    nthreads_xgb = NULL,

    #' @description Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param ids (`character()`)\cr
    #' Survival learner ids for which to request the corresponding
    #' [LearnerSurv][mlr3proba::LearnerSurv] objects and suggested
    #' hyperparameter spaces.
    #' If not provided, the other methods of this class will return objects
    #' which will include all learners currently supported.
    #' See method `lrn_ids()` for available ids to request.
    #' @param nthreads_rsf (`int(1)`)\cr
    #' Number of cores to use in random survival forest learners (implicit
    #' parallelization). Default value: use all available cores.
    #'
    #' @param nthreads_xgb (`int(1)`)\cr
    #' Number of cores to use in xgboost survival learners (implicit
    #' parallelization). Default value: use 1 core.
    #'
    #' @details If implicit parallelization is desired, the more `nthreads_rsf`
    #' the better, but that is not the case with the xgboost learners.
    #' With xgboost, adding more cores (especially when training small datasets)
    #' might create performance issues, so test before you use too many `nthreads_xgb`!
    initialize = function(ids          = NULL,
                          nthreads_rsf = unname(parallelly::availableCores()),
                          nthreads_xgb = 1) {
      lrn_ids = self$supported_lrn_ids()

      if (!is.null(ids)) {
        lrn_ids = lrn_ids[lrn_ids %in% ids]
      }

      private$.ids = lrn_ids
      self$nthreads_rsf = nthreads_rsf
      self$nthreads_xgb = nthreads_xgb
    },

    #' @description Supported Survival Learners IDs
    #' @return (`character()`)\cr
    #' A vector of ALL available survival learner ids
    supported_lrn_ids = function() {
      private$.lrn_ids
    },

    #' @description Survival Learners IDs
    #' @return (`character()`)\cr
    #' A subset of the available survival learner ids that will be used in
    #' functions `lrns()` and `lrn_tbl()` (filtered by the user upon
    #' initialization)
    lrn_ids = function() {
      private$.ids
    },

    #' @description Survival learners table
    #'
    #' @details This function uses `lrns()` and `pss()` methods to return a
    #' data table of available survival learners and tuning parameters for each.
    #'
    #' @return a [data.table::data.table] with 3 columns:
    #' 1. `id`: learner id
    #' 2. `learner`: an [mlr3proba::LearnerSurv] object
    #' 3. `param_set`: a [paradox::ParamSet] object to be used for tuning each
    #' respective survival learner
    lrn_tbl = function() {
      # get survival learners
      learners = self$lrns()

      # get parameter spaces
      pss = self$pss()

      # check that lrn_ids of lists are the same
      stopifnot(names(learners) == names(pss))

      # return table
      data.table(
        id = names(learners),
        learner = learners,
        param_set = pss
      )
    },

    #' @description Survival Learners
    #' @return (`list()`)\cr
    #' A list of mlr3 survival learners able to predict both `crank`/
    #' `lp` (ranks) and `distr` predictions (survival probabilities) on test
    #' data.
    #'
    #' @details
    #' - All learners return the prediction types `crank` and `distr` (so that
    #' both discrimination and calibration performance measures can be used)
    #' - The CoxPH and CoxBoost learners return `distr` prediction by applying a
    #' transformation on the `lp` prediction using the Breslow estimator (as is
    #' done by default from the respective packages - we don't change anything)
    #' - The `SurvivalTree`, `CoxNet`, `XGBoostCox` and `XGBoostAFT` learners
    #' return a `distr` prediction by transforming it from the `crank`
    #' prediction (the latter being equal to `lp` if the learner returns it)
    #' using the [mlr3proba::distrcompositor].
    lrns = function() {
      lrn_ids = private$.ids # user-specified, filtered at initialization
      nthreads_rsf = self$nthreads_rsf
      nthreads_xgb = self$nthreads_xgb

      # suppress warning when the fallback learner (kaplan-meier) and
      # a base learner have different main predict types
      suppressWarnings({
        # returns list
        sapply(lrn_ids, function(lrn_id) {
          if (lrn_id == 'coxph') { # CoxPH
            lrn('surv.coxph', id = 'CoxPH', fallback = lrn('surv.kaplan'))
          } else if (lrn_id == 'coxnet') { # CoxNet
            coxnet = mlr3pipelines::ppl('distrcompositor',
              learner = lrn('surv.glmnet', id = 'CoxNet',
                fallback = lrn('surv.kaplan'),
                standardize = FALSE, maxit = 10^3),
              estimator = 'kaplan',
              form = 'ph',
              overwrite = TRUE,
              graph_learner = TRUE
            )
            coxnet$id = 'CoxNet'
            coxnet$label = 'Elastic Net CoxPH'
            coxnet
          } else if (lrn_id == 'surv_tree') { # Survival Tree
            surv_tree = mlr3pipelines::ppl('distrcompositor',
              learner = lrn('surv.rpart', id = 'SurvivalTree',
                fallback = lrn('surv.kaplan')),
              estimator = 'kaplan',
              form = 'ph',
              overwrite = TRUE,
              graph_learner = TRUE
            )
            surv_tree$id = 'SurvivalTree'
            surv_tree$label = 'Survival Tree'
            surv_tree
          } else if (lrn_id == 'rsf_cindex') { # Random Survival Forests
            lrn('surv.ranger', verbose = FALSE,
              id = 'SurvivalForestCIndex',
              label = 'Random Forest (C-index)',
              fallback = lrn('surv.kaplan'),
              num.threads = nthreads_rsf,
              splitrule = 'C'  # Harrell's C-index
            )
          } else if (lrn_id == 'rsf_logrank') {
            lrn('surv.ranger', verbose = FALSE,
              id = 'SurvivalForestLogRank',
              label = 'Random Forest (Logrank)',
              fallback = lrn('surv.kaplan'),
              num.threads = nthreads_rsf,
              splitrule = 'logrank'
            )
          } else if (lrn_id == 'rsf_maxstat') {
            lrn('surv.ranger', verbose = FALSE,
              id = 'SurvivalForestMaxStat',
              label = 'Random Forest (Maximally selected rank statistics)',
              fallback = lrn('surv.kaplan'),
              num.threads = nthreads_rsf,
              splitrule = 'maxstat'
            )
          } else if (lrn_id == 'rsf_extratrees') {
            lrn('surv.ranger', verbose = FALSE,
              id = 'SurvivalForestExtraTrees',
              label = 'Random Forest (Extremely Randomized Trees)',
              fallback = lrn('surv.kaplan'),
              num.threads = nthreads_rsf,
              splitrule = 'extratrees',
              # default (but keep it here, to show we are not tuning it)
              num.random.splits = 1
            )
          } else if (lrn_id == 'aorsf') {
            lrn('surv.aorsf',
              id = 'ObliqueSurvivalForestFast',
              label = 'Accelerated Oblique Random Forest',
              fallback = lrn('surv.kaplan'),
              control_type = 'fast',
              oobag_pred_type = 'surv',
              importance = 'anova', # very fast, leave it as it is for now
              attach_data = TRUE # this is needed for prediction and importance
            )
          } else if (lrn_id == 'coxboost') { # CoxBoost
            lrn('surv.coxboost', id = 'CoxBoost',
              fallback = lrn('surv.kaplan'),
              standardize = FALSE, # data already standardized
              return.score = FALSE # don't need this in the output
            )
          } else if (grepl(pattern = 'xgboost_cox', x = lrn_id)) { # XGBoost Cox
            xgboost_cox = mlr3pipelines::ppl('distrcompositor',
              learner = lrn('surv.xgboost', nthread = nthreads_xgb,
                booster = 'gbtree', fallback = lrn('surv.kaplan'),
                objective = 'survival:cox', id = 'XGBoostCox'),
              estimator = 'kaplan',
              form = 'ph',
              overwrite = TRUE,
              graph_learner = TRUE
            )
            xgboost_cox$id = 'XGBoostCox'
            xgboost_cox$label = 'Extreme Gradient Boosting (Cox)'

            # force the use of a 'test' set for early validation
            if (grepl(pattern = 'early', x = lrn_id)) {
              xgboost_cox$param_set$values$XGBoostCox.early_stopping_set = 'test'
            }
            xgboost_cox
          } else if (grepl(pattern = 'xgboost_aft', x = lrn_id)) {
            xgboost_aft = mlr3pipelines::ppl('distrcompositor',
              learner = lrn('surv.xgboost', nthread = nthreads_xgb,
                booster = 'gbtree', fallback = lrn('surv.kaplan'),
                objective = 'survival:aft', id = 'XGBoostAFT'),
              estimator = 'kaplan',
              form = 'aft',
              overwrite = TRUE,
              graph_learner = TRUE
            )
            xgboost_aft$id = 'XGBoostAFT'
            xgboost_aft$label = 'Extreme Gradient Boosting (AFT)'

            # force the use of a 'test' set for early validation
            if (grepl(pattern = 'early', x = lrn_id)) {
              xgboost_aft$param_set$values$XGBoostAFT.early_stopping_set = 'test'
            }
            xgboost_aft
          }
        }, simplify = FALSE, USE.NAMES = TRUE)
      })
    },

    #' @description Suggested Parameter Spaces
    #'
    #' @details For the XGBoostCox learner, we provide:
    #'   1) A basic parameter space with 4 hyperparameters
    #'   2) The above 4 hyperparameters + support of early stopping (10% of the
    #'   `nrounds`)
    #'   3) The basic 4 hyperparameters + 3 regularization hyperparameters.
    #'
    #' The above 3 parameter spaces are also offered for the XGBoostAFT learner with
    #' the addition of 2 more hyperparameter to tune. See examples.
    #'
    #' @return (`list()`)\cr
    #' A list with suggested [paradox::ps()] objects (to be used for tuning)
    pss = function() {
      lrn_ids = private$.ids # user-specified, filtered at initialization

      # returns list
      sapply(lrn_ids, function(lrn_id) {
        # Note
        # l$coxph = NULL (by default in R) so no need to specify explicitly
        if (lrn_id == 'coxnet') { # CoxNet
          paradox::ps(
            CoxNet.lambda = p_dbl(1e-03, 10, logscale = TRUE),
            CoxNet.alpha  = p_dbl(0, 1) # from Ridge to Lasso penalty
          )
        } else if (lrn_id == 'surv_tree') { # Survival Tree
          paradox::ps(
            SurvivalTree.minsplit = p_int(1, 64, logscale = TRUE),
            SurvivalTree.cp = p_dbl(1e-04, 1, logscale = TRUE)
          )
        } else if (lrn_id == 'rsf_cindex') { # Random Survival Forests
          paradox::ps(
            num.trees = p_int(100, 1500),
            mtry.ratio = p_dbl(0.1, 0.9),
            min.node.size = p_int(3, 20)
          )
        } else if (lrn_id == 'rsf_logrank') {
          paradox::ps(
            num.trees = p_int(100, 1500),
            mtry.ratio = p_dbl(0.1, 0.9),
            min.node.size = p_int(3, 20)
          )
        } else if (lrn_id == 'rsf_maxstat') {
          paradox::ps(
            num.trees = p_int(100, 1500),
            mtry.ratio = p_dbl(0.1, 0.9),
            min.node.size = p_int(3, 20)
          )
        } else if (lrn_id == 'rsf_extratrees') {
          paradox::ps(
            num.trees = p_int(100, 1500),
            mtry.ratio = p_dbl(0.1, 0.9),
            min.node.size = p_int(3, 20)
          )
        } else if (lrn_id == 'aorsf') {
          paradox::ps(
            n_tree = p_int(100, 1500),
            mtry_ratio = p_dbl(0.1, 0.9),
            leaf_min_obs = p_int(3, 20)
          )
        } else if (lrn_id == 'coxboost') { # CoxBoost
          paradox::ps(
            stepno = p_int(50, 500),
            # leave at default => 9 * sum(status == 1)?
            penalty = p_int(10, 1000, logscale = TRUE),
            # leave at default => 1?
            stepsize.factor = p_dbl(1e-01, 10, logscale = TRUE),
            # use the penalized scores - `hpscore` is faster to calculate
            criterion = p_fct(c('pscore', 'hpscore'))
          )
        } else if (lrn_id == 'xgboost_cox') { # XGBoost Cox
          paradox::ps(
            XGBoostCox.nrounds = p_int(100, 1500),
            XGBoostCox.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostCox.max_depth = p_int(2, 8), # shallow trees
            XGBoostCox.min_child_weight = p_dbl(1, 128, logscale = TRUE)
          )
        } else if (lrn_id == 'xgboost_cox_early') {
          # add support for early stopping
          # should be used with `XGBoostCox.early_stopping_set = 'test'
          paradox::ps(
            XGBoostCox.nrounds = p_int(100, 1500),
            XGBoostCox.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostCox.max_depth = p_int(2, 8), # shallow trees
            XGBoostCox.min_child_weight = p_dbl(1, 128, logscale = TRUE),
            .extra_trafo = function(x, param_set) {
              x$XGBoostCox.early_stopping_rounds =
                as.integer(ceiling(0.1 * x$XGBoostCox.nrounds))
              x
            }
          )
        } else if (lrn_id == 'xgboost_cox_reg') {
          # more explicit control of regularization
          paradox::ps(
            XGBoostCox.nrounds = p_int(100, 1500),
            XGBoostCox.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostCox.max_depth = p_int(2, 8),
            XGBoostCox.min_child_weight = p_dbl(1, 128, logscale = TRUE),
            XGBoostCox.alpha  = p_dbl(1e-03, 10, logscale = TRUE),
            XGBoostCox.lambda = p_dbl(1e-03, 10, logscale = TRUE),
            XGBoostCox.gamma  = p_dbl(0, 5)
          )
        } else if (lrn_id == 'xgboost_aft') { # and the AFT counterparts
          paradox::ps(
            XGBoostAFT.nrounds = p_int(100, 1500),
            XGBoostAFT.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostAFT.max_depth = p_int(2, 8), # shallow trees
            XGBoostAFT.min_child_weight = p_dbl(1, 128, logscale = TRUE),
            XGBoostAFT.aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
            XGBoostAFT.aft_loss_distribution_scale = p_dbl(0.5, 2.0)
          )
        } else if (lrn_id == 'xgboost_aft_early') {
          paradox::ps(
            XGBoostAFT.nrounds = p_int(100, 1500),
            XGBoostAFT.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostAFT.max_depth = p_int(2, 8), # shallow trees
            XGBoostAFT.min_child_weight = p_dbl(1, 128, logscale = TRUE),
            XGBoostAFT.aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
            XGBoostAFT.aft_loss_distribution_scale = p_dbl(0.5, 2.0),
            .extra_trafo = function(x, param_set) {
              x$XGBoostAFT.early_stopping_rounds =
                as.integer(ceiling(0.1 * x$XGBoostAFT.nrounds))
              x
            }
          )
        } else if (lrn_id == 'xgboost_aft_reg') {
          paradox::ps(
            XGBoostAFT.nrounds = p_int(100, 1500),
            XGBoostAFT.eta = p_dbl(1e-04, 0.3, logscale = TRUE),
            XGBoostAFT.max_depth = p_int(2, 8),
            XGBoostAFT.min_child_weight = p_dbl(1, 128, logscale = TRUE),
            XGBoostAFT.alpha  = p_dbl(1e-03, 10, logscale = TRUE),
            XGBoostAFT.lambda = p_dbl(1e-03, 10, logscale = TRUE),
            XGBoostAFT.gamma  = p_dbl(0, 5),
            XGBoostAFT.aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
            XGBoostAFT.aft_loss_distribution_scale = p_dbl(0.5, 2.0)
          )
        }
      }, simplify = FALSE, USE.NAMES = TRUE)
    },

    #' @description
    #' Opens the help page for this object.
    help = function() {
      mlr3misc::open_help('survmob::SurvLPS')
    }
  ), private = list(
    # Learner ids
    # UPDATE here if any more are added
    .lrn_ids = c(
      'coxph',
      'coxnet',
      'surv_tree',
      'rsf_cindex',
      'rsf_logrank',
      'rsf_maxstat',
      'rsf_extratrees',
      'aorsf',
      'coxboost',
      'xgboost_cox',
      'xgboost_cox_early',
      'xgboost_cox_reg',
      'xgboost_aft',
      'xgboost_aft_early',
      'xgboost_aft_reg'
    ),
    # filtered, user-provided `lrn_ids`
    .ids = NULL
  )
)
