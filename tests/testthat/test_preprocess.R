test_that('time-status suffling works', {
  task = tsk('lung')
  poss = PipeOpSurvShuffle$new()

  expect_equal(poss$param_set$values$replace, FALSE)

  task_shuffled = poss$train(list(task))[[1L]]
  dt_shuffled = task_shuffled$data(cols = c('time', 'status'))
  dt_original = task$data(cols = c('time', 'status'))

  expect_false(all(dt_original == dt_shuffled))
})

test_that('PipeOpRemoveNAs works', {
  task = tsk('lung')
  expect_equal(length(task$feature_names), 8)

  pona = PipeOpRemoveNAs$new()
  expect_equal(pona$ncolsNA(task), 6) # 6 features with at least 1 NA
  expect_equal(pona$param_set$values$cutoff, 0.2) # default cutoff value

  task2 = pona$train(list(task))[[1L]]
  expect_equal(length(task2$feature_names), 7) # meal.cal removed

  pona$param_set$values$cutoff = 0.05
  task3 = pona$train(list(task))[[1L]]
  expect_equal(length(task3$feature_names), 6) # meal_cal + wt_loss removed

  # remove all features with at least 1 NA
  pona$param_set$values$cutoff = 0
  task4 = pona$train(list(task))[[1]]
  expect_equal(length(task4$feature_names), 2)

  # create another task to test what happens when all features are removed
  df = data.frame(
    time = c(1,2,3,4),
    status = c(0,1,0,1),
    X1 = c(NA,0,0,0),
    X2 = c(NA,NA,0,0),
    X3 = c(NA,NA,NA,0),
    X4 = c(NA,NA,NA,NA)
  )
  ttask = as_task_surv(x = df, id = 'test', time = 'time', event = 'status')
  expect_equal(pona$ncolsNA(ttask), 4)

  pona$param_set$values$cutoff = 0.75
  ttask1 = pona$train(list(ttask))[[1L]]
  expect_equal(length(ttask1$feature_names), 3)

  pona$param_set$values$cutoff = 0.5
  ttask2 = pona$train(list(ttask))[[1L]]
  expect_equal(length(ttask2$feature_names), 2)

  pona$param_set$values$cutoff = 0.25
  ttask3 = pona$train(list(ttask))[[1L]]
  expect_equal(length(ttask3$feature_names), 1)

  pona$param_set$values$cutoff = 0.2
  ttask4 = pona$train(list(ttask))[[1L]]
  expect_equal(length(ttask4$feature_names), 0) # no features remain
})

test_that('PipeOpRemoveZeros works', {
  # create survival task from read count test data
  df = data.frame(
    time = c(1,2,3,4),
    status = c(0,1,0,1),
    X1 = c(999,0,0,0),
    X2 = c(23,NA,0,0),
    X3 = c(NA,3,NA,0),
    X4 = c(NA,1231,32,134),
    X5 = c(0,0,0,0),
    X6 = c(12,13,14,15)
  )
  task = as_task_surv(x = df, id = 'test', time = 'time', event = 'status')
  expect_equal(length(task$feature_names), 6)

  poz = PipeOpRemoveZeros$new()
  expect_equal(poz$ncolsZero(task), 4) # 4 features with at least 1 zero
  expect_equal(poz$param_set$values$cutoff, 0.2) # default cutoff value
  nzeros = poz$.__enclos_env__$private$.get_nzeros(task) # zero counts
  expect_equal(nzeros, data.table(X1 = 3, X2 = 2, X3 = 1, X4 = 0, X5 = 4, X6 = 0))
  # nzeros/task$nrow # percentage of zeros per column

  poz$param_set$values$cutoff = 1
  task2 = poz$train(list(task))[[1L]]
  expect_equal(length(task2$feature_names), 6) # no features removed

  poz$param_set$values$cutoff = 0.75
  task3 = poz$train(list(task))[[1L]]
  expect_equal(length(task3$feature_names), 5) # 1 feature removed
  expect_equal(task3$feature_names, c('X1', 'X2', 'X3', 'X4', 'X6'))

  poz$param_set$values$cutoff = 0.5
  task4 = poz$train(list(task))[[1L]]
  expect_equal(length(task4$feature_names), 4) # 2 features removed
  expect_equal(task4$feature_names, c('X2', 'X3', 'X4', 'X6'))

  poz$param_set$values$cutoff = 0.1
  task5 = poz$train(list(task))[[1L]]
  expect_equal(length(task5$feature_names), 2) # 4 features removed
  expect_equal(task5$feature_names, c('X4', 'X6'))

  # remove all features with at least 1 zero
  poz$param_set$values$cutoff = 0
  task6 = poz$train(list(task))[[1]]
  expect_equal(length(task6$feature_names), 2)
  expect_equal(task6$feature_names, c('X4', 'X6'))
})

test_that('minimize_backend works', {
  task = tsk('lung')

  pona = PipeOpRemoveNAs$new(param_vals = list(cutoff = 0))
  task1 = pona$train(list(task))[[1L]]
  expect_equal(length(task1$feature_names), 2)

  # backend hasn't changed!
  expect_equal(task1$backend$ncol, 11) # 8 feats, 2 targets (time+status) + row_id
  expect_equal(task1$backend$nrow, 228)

  # now it is
  task2 = minimize_backend(task1)
  expect_equal(task2$backend$ncol, 5) # 2 feats, 2 targets (time+status) + row_id
  expect_equal(task2$backend$nrow, 228)

  # original task didn't change
  expect_equal(task1$backend$ncol, 11)
})
