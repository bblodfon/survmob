test_that('task_powerset works', {
  # generate 3 survival tasks
  gen = mlr3::mlr_task_generators$get('simsurv')
  task1 = gen$generate(5)
  task2 = gen$generate(5)
  task3 = gen$generate(5)
  task4 = gen$generate(6)

  # change feature names (need to be different)
  task2$rename(old = task2$feature_names, new = paste0(task2$feature_names, '2'))
  task3$rename(old = task3$feature_names, new = paste0(task3$feature_names, '3'))

  tasks = list(task1, task2, task3)
  # task ids are all the same (error)
  expect_true(task1$id == task2$id)
  expect_true(task2$id == task3$id)
  expect_error(task_powerset(tasks))

  # `tasks` list is modified by reference
  task1$id = 'task1'
  task2$id = 'task2'
  task3$id = 'task3'

  # number of rows must be the same
  expect_error(task_powerset(list(task1, task4)))

  # check 3-element list
  powset = task_powerset(tasks, check_targets = FALSE)

  # empty subset is excluded
  expect_equal(length(powset), 2^3 - 1)
  # all features are included
  expect_equal(powset[['task1-task3']]$feature_names,
    c(task1$feature_names, task3$feature_names))

  # check 2-element list
  powset2 = task_powerset(tasks[1:2], check_targets = FALSE)
  expect_equal(length(powset2), 2^2 - 1)
  expect_equal(names(powset2), c(task1$id, task2$id, paste0(task1$id, '-', task2$id)))

  # check 1-element list
  expect_error(task_powerset(tasks[1]))

  # test tasks with exactly the same targets
  time  = 1:5
  status = c(1,1,0,0,1)
  d1 = data.frame(status = status, time = time, a1 = LETTERS[1:5], b1 = rep(1,5))
  d2 = data.frame(status = status, time = time, a2 = LETTERS[1:5], b2 = rep(1,5))
  d3 = data.frame(status = status, time = time, a3 = LETTERS[1:5], b3 = rep(1,5))
  # survival tasks
  task1 = as_task_surv(x = d1, time = 'time', event = 'status')
  task2 = as_task_surv(x = d2, time = 'time', event = 'status')
  task3 = as_task_surv(x = d3, time = 'time', event = 'status')
  # classification tasks (2 features per task)
  task4 = as_task_classif(x = d1 %>% select(-time), id = 'd1', target = 'status')
  task5 = as_task_classif(x = d2 %>% select(-time), id = 'd2', target = 'status')
  task6 = as_task_classif(x = d3 %>% select(-time), id = 'd3', target = 'status')

  # check survival tasks
  powset3 = task_powerset(list(task1, task2, task3)) # check_targets = TRUE
  expect_equal(length(powset3), 2^3 - 1)
  expect_equal(powset3[['d1-d2']]$feature_names,
    c(task1$feature_names, task2$feature_names))
  expect_equal(powset3[['d1-d3']]$feature_names,
    c(task1$feature_names, task3$feature_names))
  expect_equal(powset3[['d1-d2-d3']]$feature_names,
    c(task1$feature_names, task2$feature_names, task3$feature_names))

  # check classification tasks
  powset4 = task_powerset(list(task4, task5, task6))
  expect_equal(length(powset4$d1$feature_names), 2)
  expect_equal(length(powset4$d3$feature_names), 2)
  expect_equal(length(powset4$`d1-d2`$feature_names), 4)
  expect_equal(length(powset4$`d1-d2-d3`$feature_names), 6)
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

test_that('powerset_icounts works', {
  # need matrix or data.frame
  expect_error(powerset_icounts(3))
  # test matrix
  m = matrix(data = c(1,0,1, 1,1,1, 1,0,0), nrow = 3, byrow = TRUE)
  # no column names
  expect_error(powerset_icounts(m))

  colnames(m) = LETTERS[1:3]
  pics = powerset_icounts(m)

  expect_equal(nrow(pics), 2^ncol(m)-1)
  expect_equal(pics %>% filter(combo_name == 'A') %>% pull(intersect_count), 3)
  expect_equal(pics %>% filter(combo_name == 'B') %>% pull(intersect_count), 1)
  expect_equal(pics %>% filter(combo_name == 'C') %>% pull(intersect_count), 2)
  expect_equal(pics %>% filter(combo_name == 'A-B') %>% pull(intersect_count), 1)
  expect_equal(pics %>% filter(combo_name == 'B-C') %>% pull(intersect_count), 1)
  expect_equal(pics %>% filter(combo_name == 'A-C') %>% pull(intersect_count), 2)
  expect_equal(pics %>% filter(combo_name == 'A-B-C') %>% pull(intersect_count), 1)
})

test_that('assert_part works', {
  expect_error(assert_part(list(a = 1)))
  expect_error(assert_part(list(a = 1, b = 2)))
  expect_error(assert_part(list(train = 1)))
  expect_error(assert_part(list(train = 1, test = 1)))
})
