test_that('bench_msrs returns the expected measures', {
  ms = bench_msrs()
  expect_equal(dim(ms), c(5, 2))

  # check external measure ids
  expect_equal(ms$id, c('harrell_c', 'uno_c', 'ibrier', 'rcll', 'dcal'))

  # Uno's C-index
  uno_c = ms[id == 'uno_c']$measure[[1L]]
  expect_class(uno_c, c('R6', 'MeasureSurv', 'MeasureSurvCindex'))
  expect_equal(uno_c$param_set$values$weight_meth, 'G2')
  expect_equal(uno_c$param_set$values$eps, 1e-3)
  expect_equal(uno_c$id, 'uno_c')

  # Integrated Brier Score
  ibrier = ms[id == 'ibrier']$measure[[1L]]
  expect_class(ibrier, c('R6', 'MeasureSurv', 'MeasureSurvGraf'))
  expect_equal(ibrier$param_set$values$integrated, TRUE)
  expect_equal(ibrier$param_set$values$method, 2)
  expect_equal(ibrier$param_set$values$proper, TRUE)
  expect_equal(ibrier$id, 'ibrier')

  # RCLL
  rcll = ms[id == 'rcll']$measure[[1L]]
  expect_class(rcll, c('R6', 'MeasureSurv', 'MeasureSurvRCLL'))
  expect_equal(rcll$param_set$values$eps, 1e-15)
  expect_equal(rcll$id, 'rcll')

  # check labels
  msr_labels = sapply(ms$measure, function(measure) measure$label)
  expect_equal(msr_labels, c('HarrellC', 'UnoC', 'IBrier', 'RCLL', 'Dcal'))

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
