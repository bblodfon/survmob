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
  expect_equal(length(task2$feature_names), 7)
  expect_equal(pona$state$removed_column_num, 1) # meal.cal removed

  pona$param_set$values$cutoff = 0.05
  task3 = pona$train(list(task))[[1L]]
  expect_equal(length(task3$feature_names), 6)
  expect_equal(pona$state$removed_column_num, 2) # meal_cal + wt_loss removed

  # remove all features with at least 1 NA
  pona$param_set$values$cutoff = 0
  task4 = pona$train(list(task))[[1]]
  expect_equal(length(task4$feature_names), 2)
  expect_equal(pona$state$removed_column_num, 6)

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
  expect_equal(pona$state$removed_column_num, 3)

  pona$param_set$values$cutoff = 0.2
  ttask4 = pona$train(list(ttask))[[1L]]
  expect_equal(length(ttask4$feature_names), 0) # no features remain
  expect_equal(pona$state$removed_column_num, 4) # all features removed
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
  expect_equal(length(task2$feature_names), 6)
  expect_equal(poz$state$removed_column_num, 0) # no features removed

  poz$param_set$values$cutoff = 0.75
  task3 = poz$train(list(task))[[1L]]
  expect_equal(length(task3$feature_names), 5)
  expect_equal(poz$state$removed_column_num, 1) # 1 feature removed
  expect_equal(task3$feature_names, c('X1', 'X2', 'X3', 'X4', 'X6'))

  poz$param_set$values$cutoff = 0.5
  task4 = poz$train(list(task))[[1L]]
  expect_equal(length(task4$feature_names), 4)
  expect_equal(poz$state$removed_column_num, 2) # 2 features removed
  expect_equal(task4$feature_names, c('X2', 'X3', 'X4', 'X6'))

  poz$param_set$values$cutoff = 0.1
  task5 = poz$train(list(task))[[1L]]
  expect_equal(length(task5$feature_names), 2)
  expect_equal(poz$state$removed_column_num, 4) # 4 features removed
  expect_equal(task5$feature_names, c('X4', 'X6'))

  # remove all features with at least 1 zero
  poz$param_set$values$cutoff = 0
  task6 = poz$train(list(task))[[1]]
  expect_equal(length(task6$feature_names), 2)
  expect_equal(poz$state$removed_column_num, 4)
  expect_equal(task6$feature_names, c('X4', 'X6'))
})

test_that('PipeOpLogTransform works', {
  # create survival task from read count test data
  df = data.frame(
    time = c(1,2,3,4),
    status = c(0,1,0,1),
    X1 = c(999,0,0,0),
    X2 = c('i',NA,'r','t'), # char column
    X3 = c(NA,3,NA,0),
    X4 = c(NA,1231,32,134),
    X5 = c(0,0,0,0),
    X6 = LETTERS[1:4] # char column
  )
  task = as_task_surv(x = df, id = 'test', time = 'time', event = 'status')
  expect_equal(length(task$feature_names), 6)

  polog = PipeOpLogTransform$new()
  # default hyperparameters are set correctly
  expect_equal(polog$param_set$values$base, 2)
  expect_equal(polog$param_set$values$offset, 1)

  task2 = polog$train(list(task))[[1L]]
  # non-numeric columns didn't change
  expect_equal(
    task$data (cols = c('time', 'status', 'X2', 'X6')),
    task2$data(cols = c('time', 'status', 'X2', 'X6'))
  )

  # order of features didn't change
  expect_equal(task2$feature_names, task$feature_names)

  # log(x+1, 2) was applied correctly
  expect_equal(
    task2$data(cols = 'X1', rows = 1)$X1, # 999
    log(task$data(cols = 'X1', rows = 1)$X1 + 1, 2)
  )
  expect_equal(
    task2$data(cols = 'X4', rows = 3)$X4, # 32
    log(task$data(cols = 'X4', rows = 3)$X4 + 1, 2)
  )
  expect_equal(
    task2$data(cols = 'X5')$X5, # all zeros
    log(task$data(cols = 'X5')$X5 + 1, 2)
  )

  # different base and offset => log(x+2, exp(1))
  polog2 = PipeOpLogTransform$new(param_vals = list(base = exp(1), offset = 2))
  task3 = polog2$train(list(task))[[1L]]

  expect_equal(
    task3$data(cols = 'X1', rows = 1)$X1, # 999
    log(task$data(cols = 'X1', rows = 1)$X1 + 2, exp(1))
  )
  expect_equal(
    task3$data(cols = 'X4', rows = 3)$X4, # 32
    log(task$data(cols = 'X4', rows = 3)$X4 + 2, exp(1))
  )
  expect_equal(
    task3$data(cols = 'X5')$X5, # all zeros
    log(task$data(cols = 'X5')$X5 + 2, exp(1))
  )

  # order the same
  expect_equal(task3$col_roles$feature, task$col_roles$feature)
})

test_that('initialization with po (PipeOp) shorthand constructor works', {
  poss = po('survshuffle')
  expect_equal(poss$param_set$values$replace, FALSE)

  pona = po('removenas')
  expect_equal(pona$param_set$values$cutoff, 0.2)

  porz = po('removezeros')
  expect_equal(porz$param_set$values$cutoff, 0.2)

  polog = po('logtransform', base = exp(1), offset = 1)
  expect_equal(polog$param_set$values$base, exp(1))
  expect_equal(polog$param_set$values$offset, 1)
})
