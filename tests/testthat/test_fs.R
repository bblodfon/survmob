test_that('eFS is initialized properly', {
  efs = eFS$new(nthreads_rsf = 7)

  # NULL task at initialization
  expect_null(efs$task_id)

  # RFE/AutoFSelector parameters
  expect_equal(efs$msr_id, 'harrell_c')
  expect_equal(efs$resampling$id, 'repeated_cv')
  expect_equal(efs$resampling$param_set$values$repeats, 5)
  expect_equal(efs$resampling$param_set$values$folds, 5)
  expect_equal(efs$repeats, 100)
  expect_equal(efs$feature_fraction, 0.8)
  expect_equal(efs$n_features, 2)

  # RSF parameters
  expect_equal(efs$nthreads_rsf, 7)
  expect_equal(efs$num_trees, 250)
  expect_equal(efs$mtry_ratio, 0.1)

  # some errors to check
  expect_error(efs$new(msr_id = 'NotProperMeasure'))
  expect_error(eFS$new(resampling = 'NotProperResamplingObj'))
  expect_error(eFS$new(repeats = 0))
  expect_error(eFS$new(n_features = 0))
  expect_error(eFS$new(feature_fraction = 1)) # must be less than 1 strictly
  expect_error(eFS$new(num_trees = 0))
  expect_error(eFS$new(nthreads_rsf = 0))
  expect_error(eFS$new(mtry_ratio = 1.2))
})

test_that('subsetsRFE works properly', {
  efs = eFS$new(feature_fraction = 0, n_features = 1)
  expect_equal(efs$subsetsRFE(taskv), 9) # taskv has 9 features

  efs1 = eFS$new(feature_fraction = 0.3, n_features = 1)
  expect_equal(efs1$subsetsRFE(taskv), c(9,2))

  efs2 = eFS$new(feature_fraction = 0.5, n_features = 1)
  expect_equal(efs2$subsetsRFE(taskv), c(9,4,2,1))

  efs3 = eFS$new(feature_fraction = 0.9, n_features = 1)
  expect_equal(efs3$subsetsRFE(taskv), rev(1:9))
})

test_that('RSF learners and ids are properly initialized', {
  efs = eFS$new(nthreads_rsf = 3)
  supp_lrn_ids = c('rsf_cindex', 'rsf_logrank', 'rsf_maxstat')

  # learner ids are assigned properly (when correctly initialized)
  expect_equal(efs$lrn_ids, supp_lrn_ids)
  expect_equal(efs$supported_lrn_ids(), supp_lrn_ids)
  expect_equal(eFS$new(lrn_ids = 'rsf_logrank')$lrn_ids, 'rsf_logrank')
  expect_equal(eFS$new(lrn_ids = supp_lrn_ids[1:2])$lrn_ids, supp_lrn_ids[1:2])
  expect_error(efs$new(lrn_ids = c('rsf_logrank', 'UnknownID')))

  # check that learners have proper default hyperparameters
  efs = eFS$new(nthreads_rsf = 3)
  rsf_lrns = efs$.__enclos_env__$private$.get_lrns()

  expect_equal(names(rsf_lrns), supp_lrn_ids)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$num.trees, 250)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$mtry.ratio, 0.1)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$min.node.size, )
  rsf_lrns$rsf_cindex$param_set$values$importance

})


