test_that('eFS is initialized properly', {
  efs = eFS$new(nthreads_rsf = 7)

  # NULL task at initialization
  expect_null(efs$task_id)

  # RFE/AutoFSelector parameters
  expect_equal(efs$msr_id, 'oob_error')
  expect_equal(efs$resampling$id, 'insample')
  expect_equal(efs$repeats, 100)
  expect_equal(efs$feature_fraction, 0.8)
  expect_equal(efs$n_features, 2)
  expect_true(efs$adaptive_mr)

  # RSF parameters
  expect_equal(efs$nthreads_rsf, 7)
  expect_equal(efs$num_trees, 250)
  expect_equal(efs$mtry_ratio, 0.05)

  # errors
  expect_error(efs$new(msr_id = 'NotProperMeasure'))
  expect_error(eFS$new(resampling = 'NotProperResamplingObj'))
  expect_error(eFS$new(repeats = 0))
  expect_error(eFS$new(n_features = 0))
  expect_error(eFS$new(feature_fraction = 1)) # must be less than 1 strictly
  expect_error(eFS$new(num_trees = 0))
  expect_error(eFS$new(nthreads_rsf = 0))
  expect_error(eFS$new(mtry_ratio = 1.2))

  # warnings
  expect_warning(eFS$new(msr_id = 'oob_error', resampling = rsmp('cv')))
})

test_that('rfe_info() works properly', {
  efs = eFS$new(feature_fraction = 0, n_features = 1)
  tbl = efs$rfe_info(taskv)
  expect_class(tbl, c('tbl', 'data.frame'))
  expect_equal(tbl$subset_size, 9) # taskv has 9 features

  efs1 = eFS$new(feature_fraction = 0.3, n_features = 1)
  tbl1 = efs1$rfe_info(taskv)
  expect_equal(tbl1$subset_size, c(9,2))

  efs2 = eFS$new(feature_fraction = 0.5, n_features = 1)
  tbl2 = efs2$rfe_info(taskv)
  expect_equal(tbl2$subset_size, c(9,4,2,1))

  efs3 = eFS$new(feature_fraction = 0.9, n_features = 1, mtry_ratio = 0.5)
  tbl3 = efs3$rfe_info(taskv)
  expect_equal(tbl3$subset_size, rev(1:9))
  expect_equal(tbl3$mtry_ratio[1], 0.5)
  expect_equal(tbl3$mtry_ratio[9], 1) # `mtry_ratio` increases

  efs3$adaptive_mr = FALSE
  tbl4 = efs3$rfe_info(taskv)
  expect_equal(tbl4$subset_size, rev(1:9))
  expect_equal(tbl4$mtry_ratio, rep(0.5, 9)) # `mtry_ratio` doesn't change
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
  expect_error(efs$new(lrn_ids = c()))
  expect_error(efs$new(lrn_ids = NULL))

  # check that learners have proper hyperparameters
  efs = eFS$new(nthreads_rsf = 3)
  rsf_lrns = efs$.__enclos_env__$private$.get_lrns()

  expect_equal(names(rsf_lrns), supp_lrn_ids)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$num.trees, 250)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$mtry.ratio, 0.05)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$min.node.size, 3)
  expect_equal(rsf_lrns$rsf_cindex$param_set$values$importance, 'permutation')

  for(rsf_lrn in rsf_lrns) {
    expect_true(all(c('importance', 'oob_error') %in% rsf_lrn$properties))
  }
})

test_that('run() works', {
  efs = eFS$new(lrn_ids = 'rsf_logrank', nthreads_rsf = 1,
    feature_fraction = 0.9, n_features = 1, mtry_ratio = 0.5,
    repeats = 1
  )
  expect_null(efs$result)
  expect_null(efs$task_id)

  # run eFS
  result = efs$run(task = taskv, verbose = FALSE, store_archive = TRUE)
  expect_equal(result, efs$result)

  # check that parameters were updated
  expect_equal(efs$task_id, 'veteran')
  expect_class(result, 'tbl')
  # check output tibble has correct properties
  expect_equal(dim(result), c(1,6))
  expect_equal(colnames(result), c('lrn_id', 'iter', 'selected_features',
                                   'nfeatures', 'score', 'archive'))

  # check subset_sizes and adaptive mtry.ratio is correct
  arch = efs$result$archive[[1]]
  rfe_tbl = efs$rfe_info(taskv)
  subset_sizes = unlist(lapply(as.data.table(arch)$importance, length))
  mtry_ratios = unlist(mlr3misc::map(as.data.table(arch)$resample_result,
    function(rr) rr$learner$param_set$values$mtry.ratio))
  expect_equal(rfe_tbl$subset_size, subset_sizes)
  expect_equal(rfe_tbl$mtry_ratio, mtry_ratios)

  # check that the smallest oob_error is actually selected
  dt = as.data.table(arch)
  oob_error_min = min(dt$oob_error)
  index = which.min(dt$oob_error)
  expect_equal(efs$result$score, oob_error_min)

  # check that the proper features are selected
  sel_feat_mat = as.matrix(dt[index, 1:length(taskv$feature_names)])
  sel_feats = colnames(sel_feat_mat)[which(sel_feat_mat)]
  expect_equal(efs$result$selected_features[[1L]], sel_feats)
  expect_equal(efs$result$nfeatures, length(sel_feats))

  # RCLL may not work with just one observation in the test as censored
  # so we do a stratified CV split to 4 folds where each fold will have
  # (2,2,2,3) censored observations as `taskv` has a total of 9
  efs2 = eFS$new(lrn_ids = c('rsf_logrank'), nthreads_rsf = 1,
    msr_id = 'rcll', resampling = rsmp('cv', folds = 4),
    repeats = 1, mtry_ratio = 0.8)

  # check that msr_id and resampling are different
  expect_equal(efs2$msr_id, 'rcll')
  expect_equal(efs2$resampling$id, 'cv')
  expect_equal(efs2$resampling$param_set$values$folds, 4)

  taskv2 = taskv$clone(deep = TRUE)
  taskv2$col_roles$stratum = 'status'

  result2 = efs2$run(task = taskv2, verbose = FALSE)
  expect_equal(result2, efs2$result)
})

test_that('fs_stats() works', {
  efs = eFS$new(nthreads_rsf = 1)
  expect_error(efs$fs_stats())

  # hacky result object (1 learner)
  efs$result = tibble(
    lrn_id = rep('lrn', 4),
    iter   = 1:4,
    selected_features = list(LETTERS[1:4], LETTERS[2:3], LETTERS[2:6], LETTERS[3]),
    score  = rep(0.98, 4)
  )

  freq_list = efs$fs_stats()
  expect_class(freq_list, 'list')
  expect_equal(length(freq_list), 1)
  expect_equal(names(freq_list), 'lrn')
  expect_equal(freq_list$lrn$feat_name[1:3], c('C', 'B', 'D'))

  # hacky result object (2 learners)
  efs$result = tibble(
    lrn_id = c('lrn1', 'lrn1', 'lrn2', 'lrn2'),
    iter   = 1:4,
    selected_features = list(LETTERS[1:3], LETTERS[2:3], LETTERS[2:6], LETTERS[3]),
    score  = rep(0.98, 4)
  )

  freq_list2 = efs$fs_stats()
  expect_class(freq_list2, 'list')
  expect_equal(length(freq_list2), 3)
  expect_equal(names(freq_list2), c('lrn1', 'lrn2', 'consensus'))
  expect_equal(freq_list2$lrn2$feat_name[1], 'C')
  expect_equal(freq_list2$consensus$feat_name[1:2], c('C', 'B'))
})

test_that('stab() works', {
  efs = eFS$new(nthreads_rsf = 1)
  expect_error(efs$stab()) # need to have a `$result`

  # hacky result object (1 learner)
  efs$result = tibble(
    lrn_id = rep('lrn', 4),
    iter   = 1:4,
    selected_features = list(LETTERS[1:4], LETTERS[2:3], LETTERS[2:6], LETTERS[3]),
    score  = rep(0.98, 4)
  )

  # all measures asked by default so need to provide task (Nogueira measure)
  expect_error(efs$stab())
  expect_error(efs$stab(stab_metrics = 'nogueira'))
  stab_ms = efs$.__enclos_env__$private$.get_stab_metrics()
  expect_equal(stab_ms, c('jaccard', 'nogueira'))
  # not in the available list
  expect_error(efs$stab(stab_metrics = 'NotValidStabilityMeasure'))

  # don't need task for Jaccard
  t1 = efs$stab(stab_metrics = 'jaccard')
  expect_equal(colnames(t1), c('lrn_id', 'jaccard'))
  expect_equal(t1$lrn_id, 'lrn')

  # hacky task
  m = as.data.frame(matrix(data = rep(0,10), nrow = 1))
  colnames(m) = c('time', 'status', LETTERS[1:8])
  task = as_task_surv(x = m, id = 'empty', time = 'time', event = 'status')

  t2 = efs$stab(task)
  expect_equal(colnames(t2), c('lrn_id', 'jaccard', 'nogueira'))
  expect_equal(t2$lrn_id, 'lrn')
})
