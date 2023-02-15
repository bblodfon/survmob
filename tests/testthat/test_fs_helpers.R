test_that('task_powerset() works', {
  # generate 3 survival tasks
  gen = mlr3::mlr_task_generators$get('simsurv')
  task1 = gen$generate(5)
  task2 = gen$generate(5)
  task3 = gen$generate(5)

  # change feature names (need to be different)
  task2$rename(old = task2$feature_names, new = paste0(task2$feature_names, '2'))
  task3$rename(old = task3$feature_names, new = paste0(task3$feature_names, '3'))

  tasks = list(task1, task2, task3)
  # task ids are all the same
  expect_true(task1$id == task2$id)
  expect_true(task2$id == task3$id)
  expect_error(task_powerset(tasks))

  # `tasks` list is modified by reference
  task1$id = 'task1'
  task2$id = 'task2'
  task3$id = 'task3'

  # check 3-element list
  powset = task_powerset(tasks)

  # empty subset is excluded
  expect_equal(length(powset), 2^3 - 1)
  # all features are included
  expect_equal(powset[['task1-task3']]$feature_names,
               c(task1$feature_names, task3$feature_names))

  # check 2-element list
  powset2 = task_powerset(tasks[1:2])
  expect_equal(length(powset2), 2^2 - 1)
  expect_equal(names(powset2), c(task1$id, task2$id, paste0(task1$id, '-', task2$id)))

  # check 1-element list
  powset1 = task_powerset(tasks[1])
  expect_equal(length(powset1), 1)
  expect_equal(names(powset1), task1$id)
})
