test_that('bench_msrs returns the expected measures', {
  ms = bench_msrs()
  expect_equal(dim(ms), c(7, 2))

  # check external measure ids
  expect_equal(ms$id, c('harrell_c', 'uno_c', 'ibrier', 'ibrier_erv', 'rcll',
    'rcll_erv', 'dcal'))

  # Uno's C-index
  uno_c = ms[id == 'uno_c']$measure[[1L]]
  expect_class(uno_c, c('R6', 'MeasureSurv', 'MeasureSurvCindex'))
  expect_equal(uno_c$param_set$values$weight_meth, 'G2')
  expect_equal(uno_c$param_set$values$eps, 1e-3)
  expect_equal(uno_c$id, 'uno_c')
  expect_equal(uno_c$label, 'UnoC')
  expect_false(uno_c$minimize)

  # Integrated Brier Score
  ibrier = ms[id == 'ibrier']$measure[[1L]]
  expect_class(ibrier, c('R6', 'MeasureSurv', 'MeasureSurvGraf'))
  expect_true(ibrier$param_set$values$integrated)
  expect_equal(ibrier$param_set$values$method, 2)
  expect_true(ibrier$param_set$values$proper)
  expect_false(ibrier$param_set$values$ERV)
  expect_equal(ibrier$id, 'ibrier')
  expect_equal(ibrier$label, 'IBrier')
  expect_true(ibrier$minimize)

  # Integrated Brier Score (ERV version)
  ibrier_erv = ms[id == 'ibrier_erv']$measure[[1L]]
  expect_class(ibrier_erv, c('R6', 'MeasureSurv', 'MeasureSurvGraf'))
  expect_true(ibrier_erv$param_set$values$integrated)
  expect_equal(ibrier_erv$param_set$values$method, 2)
  expect_true(ibrier_erv$param_set$values$proper)
  expect_true(ibrier_erv$param_set$values$ERV)
  expect_equal(ibrier_erv$id, 'ibrier_erv')
  expect_equal(ibrier_erv$label, 'IBrier-ERV')
  expect_false(ibrier_erv$minimize)

  # RCLL
  rcll = ms[id == 'rcll']$measure[[1L]]
  expect_class(rcll, c('R6', 'MeasureSurv', 'MeasureSurvRCLL'))
  expect_equal(rcll$param_set$values$eps, 1e-15)
  expect_false(rcll$param_set$values$ERV)
  expect_equal(rcll$id, 'rcll')
  expect_equal(rcll$label, 'RCLL')
  expect_true(rcll$minimize)

  # RCLL (ERV version)
  rcll_erv = ms[id == 'rcll_erv']$measure[[1L]]
  expect_class(rcll_erv, c('R6', 'MeasureSurv', 'MeasureSurvRCLL'))
  expect_equal(rcll_erv$param_set$values$eps, 1e-15)
  expect_true(rcll_erv$param_set$values$ERV)
  expect_equal(rcll_erv$id, 'rcll_erv')
  expect_equal(rcll_erv$label, 'RCLL-ERV')
  expect_false(rcll_erv$minimize)

  # check that the two C-indexes' internal ids are different
  msr_ids = sapply(ms$measure, function(measure) measure$id)
  expect_equal(msr_ids[1:2], c('harrell_c', 'uno_c'))
})

test_that('oob_error with ranger() works', {
  learner = lrn('surv.ranger')
  learner$train(taskv)
  pred = learner$predict(taskv)
  measure = msr('oob_error')
  expect_equal(unname(pred$score(measure, learner = learner)),
               learner$oob_error())
})
