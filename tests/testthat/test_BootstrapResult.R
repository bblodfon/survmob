test_that('BootRes works', {
  expect_error(BootRes$new(test_nrsmps = 1)) # too low
  expect_error(BootRes$new(test_measure_ids = 'NonId'))
  expect_error(BootRes$new(test_workers = 0)) # too low

  boot_res = BootRes$new(test_nrsmps = 2)

  # Initialization ----
  expect_equal(length(boot_res$test_measure_ids), 5)
  expect_equal(boot_res$test_nrsmps, 2)
  expect_equal(boot_res$test_workers, 1) # sequential
  expect_null(boot_res$task_id)
  expect_null(boot_res$lrn_id)
  expect_null(boot_res$score)
  expect_null(boot_res$scores)

  # calculate bootstrapping scores ----
  cox = lrn('surv.coxph', id = 'cox')
  part = partition(taskv)
  # learner is not trained yet
  expect_error(boot_res$calculate(task = taskv, learner = cox, part = part))
  # bootstrapped scores have not been calculated yet
  expect_error(boot_res$percent_ci())
  expect_error(boot_res$score_median())

  cox$train(taskv, row_ids = part$train)
  boot_res$calculate(task = taskv, learner = cox, part = part)
  expect_equal(boot_res$lrn_id, 'cox')
  expect_equal(boot_res$task_id, 'veteran')

  # score ----
  score = boot_res$score
  expect_class(score, 'tbl_df')
  expect_equal(colnames(score), bench_msrs()$id) # all measures
  expect_numeric(score$ibrier)

  # scores ----
  scores = boot_res$scores
  expect_class(scores, 'tbl_df')
  expect_equal(colnames(scores), bench_msrs()$id) # all measures
  expect_numeric(scores$rcll)
  expect_equal(dim(scores), c(2,5))

  # median scores ----
  score_median = boot_res$score_median()
  expect_class(score_median, 'tbl_df')
  expect_equal(colnames(score_median), bench_msrs()$id) # all measures
  expect_numeric(score_median$dcal)

  # percentile confidence intervals ----
  per_ci = boot_res$percent_ci()
  expect_class(per_ci, 'tbl_df')
  expect_equal(per_ci$perc, c('2.5%', '97.5%'))
  expect_numeric(per_ci$uno_c)
  expect_equal(dim(boot_res$percent_ci()), c(2,6))

  # 1 measure ----
  boot_res2 = BootRes$new(test_nrsmps = 2, test_measure_ids = 'rcll')
  boot_res2$calculate(task = taskv, learner = cox, part = part)
  expect_class(boot_res2$scores, 'tbl_df')
  expect_numeric(boot_res2$scores$rcll)
  expect_numeric(boot_res2$score_median()$rcll)
  # tidy way to get low and upper bounds for the conf. interval
  expect_true(boot_res2$percent_ci()$rcll[1] < boot_res2$score_median()$rcll)
  expect_true(boot_res2$percent_ci()$rcll[2] > boot_res2$score_median()$rcll)
})
