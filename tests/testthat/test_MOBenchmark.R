test_that('MOBenchmark works', {
  part = partition(taskv, ratio = 0.8)
  lrn_ids = c('coxph', 'aorsf')
  mob = MOBenchmark$new(
    lrn_ids = lrn_ids, use_callr = FALSE, tasks = list(taskv), part = part,
    nthreads_rsf = 2, nthreads_xgb = 2,
    tune_nevals = 2, test_nrsmps = 50, test_workers = 4,
    tune_rsmp = rsmp('holdout', ratio = 0.8),
    keep_models = TRUE, quiet = TRUE
  )

  # check initialized parameters
  expect_list(mob$tasks)
  expect_equal(length(mob$tasks), 1)
  expect_equal(mob$lrn_ids, lrn_ids)
  expect_false(mob$use_callr)
  expect_equal(mob$part, part)
  expect_equal(mob$nthreads_rsf, 2)
  expect_equal(mob$nthreads_xgb, 2)
  expect_true(mob$gen_task_powerset) # but not used since we use only 1 task
  expect_true(mob$quiet)
  expect_true(mob$keep_models)
  expect_equal(mob$tune_rsmp$param_set$values$ratio, 0.8)
  expect_equal(mob$tune_measure_id, 'uno_c')
  expect_equal(mob$tune_nevals, 2)
  expect_equal(mob$test_measure_ids, c('uno_c', 'rcll'))
  expect_equal(mob$test_nrsmps, 50)
  expect_equal(mob$test_workers, 4)
  expect_null(mob$result)

  # execute benchmark
  mob$run()

  # check benchmark result
  res = mob$result
  expect_class(res, 'tbl_df')
  expect_equal(dim(res), c(2,4))
  expect_equal(colnames(res), c('task_id', 'lrn_id', 'model', 'boot_res'))
  expect_equal(res$task_id, rep('veteran', 2)) # same task
  expect_true(all(res$lrn_id %in% lrn_ids)) # learners
  # 2 `test_measure_ids`
  expect_equal(colnames(res$boot_res[[1]]$score), c('uno_c', 'rcll'))
  expect_equal(colnames(res$boot_res[[1]]$scores), c('uno_c', 'rcll'))
  # bootstrap checks
  expect_equal(nrow(res$boot_res[[1]]$scores), mob$test_nrsmps)

  # drop tasks
  mob$drop_tasks()
  expect_null(mob$tasks)

  # drop models
  mob$drop_models()
  expect_equal(dim(mob$result), c(2,3))

  # reshape results
  tbl_res = reshape_mob_res(res = mob$result, add_modality_columns = F)
  expect_equal(colnames(tbl_res), c('task_id', 'lrn_id', 'rsmp_id', 'measure',
    'value'))
  expect_equal(unique(tbl_res$task_id), 'veteran')
  expect_equal(sort(unique(tbl_res$lrn_id)), sort(lrn_ids))
  expect_equal(length(unique(tbl_res$rsmp_id)), 50) #' `test_nrsmps = 50`
  expect_equal(sort(unique(tbl_res$measure)), c('rcll', 'uno_c'))
})
