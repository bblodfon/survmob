test_that('time-status suffling works', {
  task = tsk('lung')
  poss = PipeOpSurvShuffle$new()

  expect_equal(poss$param_set$values$replace, FALSE)

  task_shuffled = poss$train(list(task))[[1L]]
  dt_shuffled = task_shuffled$data(cols = c('time', 'status'))
  dt_original = task$data(cols = c('time', 'status'))

  expect_false(all(dt_original == dt_shuffled))
})
