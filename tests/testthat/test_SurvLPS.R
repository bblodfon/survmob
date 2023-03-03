test_that('SurvLPS works as expected', {
  # SurvLPS
  s = SurvLPS$new()
  expect_class(s, 'SurvLPS')
  lrn_ids = s$lrn_ids()
  all_pss = s$pss()

  # pss()
  expect_true(all(names(all_pss) %in% lrn_ids))
  expect_equal(length(all_pss), length(lrn_ids))

  s2 = SurvLPS$new(ids = c('notProperId', 'coxnet', 'xgboost_cox_early'))
  pss = s2$pss()
  expect_class(pss, 'list')
  expect_equal(length(pss), 2)
  expect_equal(names(pss), c('coxnet', 'xgboost_cox_early'))

  expect_equal(pss$coxnet$ids(), c('CoxNet.lambda', 'CoxNet.alpha'))
  res_list = pss$xgboost_cox_early$trafo(x = list(XGBoostCox.nrounds = 150,
    XGBoostCox.max_depth = 5))
  expect_equal(names(res_list), c('XGBoostCox.nrounds', 'XGBoostCox.max_depth',
    'XGBoostCox.early_stopping_rounds'))
  expect_equal(res_list$XGBoostCox.nrounds, 150)
  expect_equal(res_list$XGBoostCox.early_stopping_rounds, 15) # 10%

  # lrns()
  s = SurvLPS$new(nthreads_rsf = 5, nthreads_xgb = 4,
    ids = c('coxnet', 'notProperId', 'xgboost_cox'))
  learners = s$lrns()
  expect_equal(length(learners), 2)
  expect_equal(s$nthreads_rsf, 5)
  expect_equal(s$nthreads_xgb, 4)
  expect_equal(names(learners), c('coxnet', 'xgboost_cox'))
  expect_equal(learners$xgboost_cox$param_set$values$XGBoostCox.nthread, 4)

  # lrn_tbl()
  dt = s$lrn_tbl()
  expect_equal(dim(dt), c(2, 3))

  # test some RSFs
  rsf_lrn_ids = c('rsf_maxstat', 'rsf_logrank', 'rsf_cindex', 'rsf_extratrees',
    'aorsf')
  s = SurvLPS$new(nthreads_rsf = 5, ids = rsf_lrn_ids)
  expect_equal(sort(s$lrn_ids()), sort(rsf_lrn_ids))
  dt = s$lrn_tbl()
  expect_equal(dim(dt), c(5, 3))
  rsf_lrn = dt[id == 'rsf_extratrees']$learner[[1]]
  expect_class(rsf_lrn, c('LearnerSurv', 'LearnerSurvRanger'))
  expect_equal(rsf_lrn$param_set$values$num.threads, 5)
  expect_false(rsf_lrn$param_set$values$verbose)
  expect_equal(rsf_lrn$param_set$values$splitrule, 'extratrees')
  expect_equal(rsf_lrn$param_set$values$num.random.splits, 1)

  rsf_lrn = dt[id == 'aorsf']$learner[[1]]
  expect_class(rsf_lrn, c('LearnerSurv', 'LearnerSurvAorsf'))
  expect_equal(rsf_lrn$param_set$values$control_type, 'fast')
  expect_equal(rsf_lrn$param_set$values$oobag_pred_type, 'surv')
  expect_equal(rsf_lrn$param_set$values$importance, 'anova')
  expect_true(rsf_lrn$param_set$values$attach_data)
})

test_that('Tune CoxNet using Uno\'s C-index', {
  s = SurvLPS$new(ids = 'coxnet')
  dt = s$lrn_tbl()

  coxnet_at = AutoTuner$new(
    learner = dt$learner[[1L]],
    resampling = rsmp('holdout'),
    measure = msr('surv.cindex', weight_meth = 'G2'),
    search_space = dt$param_set[[1L]],
    terminator = trm('evals', n_evals = 5),
    tuner = tnr('random_search')
  )
  coxnet_at$train(taskv)
  expect_true('surv.cindex' %in% colnames(coxnet_at$archive$data))
  expect_numeric(coxnet_at$archive$data$surv.cindex)

  p = coxnet_at$predict(taskv)
  expect_class(p$distr, c('Distribution', 'Matdist')) # distr prediction exists!

  # All following scores should return something
  expect_number(p$score(msr('surv.cindex')))
  expect_number(p$score(msr('surv.rcll')))
  expect_number(p$score(msr('surv.brier')))
})

test_that('Tune XGBoost (Cox, AFT) using RCLL measure + early stopping', {
  s = SurvLPS$new(nthreads_xgb = 1, ids = c('xgboost_cox_early', 'xgboost_aft_early'))
  dt = s$lrn_tbl()
  xgb_learner = dt$learner[[1L]]
  xgb_learner_aft = dt$learner[[2L]]
  # IMPORTANT that this is set if using an *_early xgboost learner
  expect_equal(xgb_learner$param_set$values$XGBoostCox.early_stopping_set, 'test')
  expect_equal(xgb_learner_aft$param_set$values$XGBoostAFT.early_stopping_set, 'test')

  xgb_ps = dt$param_set[[1L]]
  expect_equal(xgb_ps$trafo(x =
    list(XGBoostCox.nrounds = 100))$XGBoostCox.early_stopping_rounds,
    10
  )
  xgb_aft_ps = dt$param_set[[2L]]
  expect_true('XGBoostAFT.aft_loss_distribution' %in% xgb_aft_ps$ids())
  expect_true('XGBoostAFT.aft_loss_distribution_scale' %in% xgb_aft_ps$ids())

  # rewrite pss for speed (use less `nrounds`)
  xgb_ps = paradox::ps(
    XGBoostCox.nrounds = p_int(3, 15),
    XGBoostCox.eta = p_dbl(0.1, 0.3, logscale = TRUE),
    XGBoostCox.max_depth = p_int(2, 8), # shallow trees
    XGBoostCox.min_child_weight = p_dbl(1, 128, logscale = TRUE),
    .extra_trafo = function(x, param_set) {
      x$XGBoostCox.early_stopping_rounds =
        as.integer(ceiling(0.1 * x$XGBoostCox.nrounds))
      x
    }
  )

  xgb_aft_ps = paradox::ps(
    XGBoostAFT.nrounds = p_int(3, 15),
    XGBoostAFT.eta = p_dbl(0.1, 0.3, logscale = TRUE),
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

  # XGBoostCox tuning with early stopping
  xgb_at = AutoTuner$new(
    learner = xgb_learner,
    resampling = rsmp('holdout'),
    measure = msr('surv.rcll'),
    search_space = xgb_ps,
    terminator = trm('evals', n_evals = 2),
    tuner = tnr('random_search'),
    callbacks = xgb_es_callback
  )
  xgb_at$train(taskv)

  # check tuning results
  expect_true('surv.rcll' %in% colnames(xgb_at$archive$data))
  # max_nrounds => max nrounds used in early stopping (best iteration)
  expect_true('max_nrounds' %in% colnames(xgb_at$archive$data))
  expect_numeric(xgb_at$archive$data$surv.rcll)
  # check that the `max_rounds` of the best archive config was used as the
  # normal `nrounds` of the final learner
  expect_equal(xgb_at$learner$param_set$values$XGBoostCox.nrounds, xgb_at$archive$best()$max_nrounds)
  expect_equal(xgb_at$learner$param_set$values$XGBoostCox.early_stopping_set, 'none')
  expect_null(xgb_at$learner$param_set$values$XGBoostCox.early_stopping_rounds)

  # check predictions
  p = xgb_at$predict(taskv)
  expect_class(p, 'PredictionSurv')
  expect_numeric(p$crank)
  expect_numeric(p$lp)
  expect_class(p$distr, c('Distribution', 'Matdist'))

  # tricky to get the importance scores!
  expect_numeric(xgb_at$learner$graph_model$pipeops$
      XGBoostCox$learner_model$importance())

  expect_number(p$score(msr('surv.cindex')))
  expect_number(p$score(msr('surv.rcll')))
  expect_number(p$score(msr('surv.brier')))

  # XGBoostAFT tuning with early stopping
  xgb_aft_at = AutoTuner$new(
    learner = xgb_learner_aft,
    resampling = rsmp('holdout'),
    measure = msr('surv.rcll'),
    search_space = xgb_aft_ps,
    terminator = trm('evals', n_evals = 2),
    tuner = tnr('random_search'),
    callbacks = xgb_es_callback
  )
  xgb_aft_at$train(taskv)

  # check tuning results
  expect_true('surv.rcll' %in% colnames(xgb_aft_at$archive$data))
  # max_nrounds => max nrounds used in early stopping (best iteration)
  expect_true('max_nrounds' %in% colnames(xgb_aft_at$archive$data))
  expect_numeric(xgb_aft_at$archive$data$surv.rcll)
  # check that the `max_rounds` of the best archive config was used as the
  # normal `nrounds` of the final learner
  expect_equal(xgb_aft_at$learner$param_set$values$XGBoostAFT.nrounds, xgb_aft_at$archive$best()$max_nrounds)
  expect_equal(xgb_aft_at$learner$param_set$values$XGBoostAFT.early_stopping_set, 'none')
  expect_null(xgb_aft_at$learner$param_set$values$XGBoostAFT.early_stopping_rounds)

  # check predictions
  p = xgb_aft_at$predict(taskv)
  expect_class(p, 'PredictionSurv')
  expect_numeric(p$crank)
  expect_numeric(p$lp)
  expect_class(p$distr, c('Distribution', 'Matdist'))

  expect_number(p$score(msr('surv.cindex')))
  expect_number(p$score(msr('surv.rcll')))
  expect_number(p$score(msr('surv.brier')))
})
