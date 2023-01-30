test_that('bench_msrs returns the expected measures', {
  ms = bench_msrs()
  expect_equal(dim(ms), c(5, 2))

  expect_equal(ms$id, c('harrell_c', 'uno_c', 'ibrier', 'rcll', 'dcal'))

  ibrier = ms[id == 'ibrier']$measure[[1L]]
  expect_class(ibrier, c('R6', 'MeasureSurv', 'MeasureSurvGraf'))
  expect_equal(ibrier$param_set$values$integrated, TRUE)
  expect_equal(ibrier$param_set$values$method, 2)
  expect_equal(ibrier$param_set$values$proper, TRUE)
})

test_that('oob_error with ranger() works', {
  learner = lrn('surv.ranger')
  learner$train(taskv)
  pred = learner$predict(taskv)
  measure = msr('oob_error')
  expect_equal(unname(pred$score(measure, learner = learner)),
               learner$oob_error())
})
