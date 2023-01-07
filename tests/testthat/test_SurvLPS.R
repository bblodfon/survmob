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
  s = SurvLPS$new(nthreads = 5, ids = c('coxnet', 'notProperId', 'xgboost_cox'))
  learners = s$lrns()
  expect_equal(length(learners), 2)
  expect_equal(s$nthreads, 5)
  expect_equal(names(learners), c('coxnet', 'xgboost_cox'))
  expect_equal(learners$xgboost_cox$param_set$values$XGBoostCox.nthread, 5)

  # lrn_tbl()
  tbl = s$lrn_tbl()
  expect_equal(dim(tbl), c(2, 3))
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

test_that('Tune XGBoost using RCLL', {
  s = SurvLPS$new(ids = c('xgboost_cox_early', 'xgboost_aft_early'))
  dt = s$lrn_tbl()
  xgb_learner = dt$learner[[1L]]
  xgb_learner_aft = dt$learner[[2L]]
  # IMPORTANT that this is set if using an *_early xgboost learner
  expect_equal(xgb_learner$param_set$values$XGBoostCox.early_stopping_set,
              'test')
  expect_equal(xgb_learner_aft$param_set$values$XGBoostAFT.early_stopping_set,
              'test')

  xgb_ps = dt$param_set[[1L]]
  expect_equal(xgb_ps$trafo(x =
    list(XGBoostCox.nrounds = 100))$XGBoostCox.early_stopping_rounds,
    10
  )
  xgb_aft_ps = dt$param_set[[2L]]
  expect_true('XGBoostAFT.aft_loss_distribution' %in% xgb_aft_ps$ids())
  expect_true('XGBoostAFT.aft_loss_distribution_scale' %in% xgb_aft_ps$ids())

  skip('XGBoost tuning can be very slow')
  # set part of the dataset for validation (IMPORTANT - otherwise there are
  # internal errors in the xgboost learners)
  taskv2 = taskv$clone(deep = TRUE)
  split = partition(taskv2, ratio = 0.8)
  taskv2$set_row_roles(split$test, 'test')

  xgb_at = AutoTuner$new(
    learner = xgb_learner,
    resampling = rsmp('holdout'),
    measure = msr('surv.brier'),
    search_space = xgb_ps,
    terminator = trm('evals', n_evals = 5),
    tuner = tnr('random_search')
  )
  xgb_at$train(taskv2)
  expect_true('surv.rcll' %in% colnames(xgb_at$archive$data))
  expect_numeric(xgb_at$archive$data$surv.rcll)

  p = xgb_at$predict(taskv2)
  expect_class(p, 'PredictionSurv')
  expect_numeric(p$crank)
  expect_numeric(p$lp)
  expect_class(p$distr, c('Distribution', 'Matdist'))

  expect_numeric(xgb_at$learner$graph_model$pipeops$
      XGBoostCox$learner_model$importance())

  expect_number(p$score(msr('surv.cindex')))
  expect_number(p$score(msr('surv.rcll')))
  expect_number(p$score(msr('surv.brier')))
})
