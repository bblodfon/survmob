test_that('TaskSub works', {
  ts = TaskSub$new(nfeats = 2, cutoff = 0.7, perc = 0.2, max_nfeats = 4)

  # check initialization
  expect_equal(ts$nfeats, 2)
  expect_equal(ts$cutoff, 0.7)
  expect_equal(ts$perc, 0.2)

  # 2 hacky tasks
  m = as.data.frame(matrix(data = rep(0,10), nrow = 1))
  t1fs = LETTERS[1:8] # task1 feature names
  colnames(m) = c('time', 'status', t1fs)
  task1 = as_task_surv(x = m, id = 'tsk1', time = 'time', event = 'status')
  t2fs = paste0(LETTERS[1:8],'1') # task2 feature names
  colnames(m) = c('time', 'status', t2fs)
  task2 = as_task_surv(x = m, id = 'tsk2', time = 'time', event = 'status')

  # 2 hacky eFS
  efs1 = eFS$new(msr_id = 'harrell_c')
  efs2 = eFS$new(msr_id = 'harrell_c')

  el = list(efs1, efs2)
  tsks = list(task1, task2)

  # Summary stats ----
  # task ids (from tasks) must be different
  task2$id = 'tsk1' # same as task1$id
  expect_error(ts$summary_stats(el, tsks))
  task2$id = 'tsk2' # return to proper, different id

  # no task ids (NULL) in eFS objects
  expect_error(ts$summary_stats(el, tsks))

  # same task ids in eFS objects
  # no 1-1 correspondence between task ids and eFS task ids
  efs1$task_id = task1$id
  efs2$task_id = task1$id
  expect_error(ts$summary_stats(el, tsks))

  # different task ids in eFS objects, but NULL results
  efs2$task_id = task2$id
  expect_error(ts$summary_stats(el, tsks))

  # hacky result object with 2 learners (so that there is a consensus fs table)
  result1 = tibble(
    lrn_id = c('lrn1', 'lrn1', 'lrn2', 'lrn2'),
    iter   = 1:4,
    selected_features = list(t1fs[1:3], t1fs[2:3], t1fs[2:6], t1fs[3]),
    nfeatures = c(3,2,5,1),
    score  = c(0.8, 0.6, 0.7, 0.8) # C-index-like
  )
  efs1$result = result1

  # still one NULL result so we get an error
  expect_error(ts$summary_stats(el, tsks))

  # no error if it's just 1 eFS + 1 task
  expect_class(ts$summary_stats(el[1], tsks[1]), c('tbl_df', 'data.frame'))

  # another hacky result
  result2 = tibble(
    lrn_id = c('lrn1', 'lrn1', 'lrn2', 'lrn2'),
    iter   = 1:4,
    selected_features = list(t2fs[1:2], t2fs[2:4], t2fs[3:5], t2fs[5]),
    nfeatures = c(2,2,4,1),
    score  = c(0.5, 0.5, 0.5, 0.5) # random C-index performance
  )
  efs2$result = result2

  # now it works
  fs_tbl = ts$summary_stats(el, tsks)
  expect_class(fs_tbl, c('tbl_df', 'data.frame'))
  expect_equal(nrow(fs_tbl), 2) # 2 tasks
  # C-index (or oob = 1 - cindex) measure so I expect to have extra columns
  extra_cols = c('comb_score', 'comb_score_norm', 'optim_nfeats')
  expect_true(all(extra_cols %in% colnames(fs_tbl)))

  # top N features
  expect_equal(fs_tbl$top_n, c(ts$nfeats, ts$nfeats))
  # top N%
  fss1 = efs1$fs_stats()$consensus
  fss2 = efs2$fs_stats()$consensus
  top_nperc1 = floor(ts$perc * length(fss1$feat_name))
  top_nperc2 = floor(ts$perc * length(fss2$feat_name))
  expect_equal(fs_tbl$top_nperc, c(top_nperc1, top_nperc2))
  # #feats selected more than (frequency cutoff) times
  gt_freq1 = fss1 %>% filter(freq > ts$cutoff) %>% nrow()
  gt_freq2 = fss2 %>% filter(freq > ts$cutoff) %>% nrow()
  expect_equal(fs_tbl$gt_freq, c(gt_freq1, gt_freq2))
  # median score
  s1 = median(result1$score)
  s2 = median(result2$score)
  expect_equal(fs_tbl$score_median, c(s1, s2))
  # median #features selected
  mfs1 = floor(median(result1$nfeatures))
  mfs2 = floor(median(result2$nfeatures))
  expect_equal(fs_tbl$nfeats_median, c(mfs1, mfs2))
  # check features available for selection (and ordering)
  expect_equal(fs_tbl$features, list(fss1$feat_name, fss2$feat_name))
  # check optim_nfeats
  expect_equal(fs_tbl$optim_nfeats,
    floor(fs_tbl$comb_score_norm * ts$max_nfeats))

  # All msr_ids should be the same!
  efs1$msr_id = 'rcll'
  expect_error(ts$summary_stats(el, tsks))

  # RCLL test (hacked - results are like C-indexes)
  efs1$msr_id = 'rcll'
  efs2$msr_id = 'rcll'
  fs_tbl_rcll = ts$summary_stats(el, tsks) # now it works
  expect_class(fs_tbl_rcll, c('tbl_df', 'data.frame'))
  # no extra columns
  expect_false(any(extra_cols %in% colnames(fs_tbl_rcll)))

  # no consensus in an eFS!
  efs2$result = result2[1:2,]
  expect_error(ts$summary_stats(el, tsks))
  efs2$result = result2 # restore consensus result

  # max_nfeats > task's #features
  ts$max_nfeats = 30
  expect_error(ts$summary_stats(el, tsks))
  ts$max_nfeats = 4 # restore

  # nfeats > task's #features
  ts$nfeats = 20
  expect_error(ts$summary_stats(el, tsks))
  ts$nfeats = 2

  # Subset tasks ----
  # not a proper method
  expect_error(ts$subset_tasks(el, tsks, method = 'de'))
  # no `optim_nfeats` method supported for RCLL score!
  expect_error(ts$subset_tasks(el, tsks, method = 'optim_nfeats'))

  # restore C-index measure
  efs1$msr_id = 'harrell_c'
  efs2$msr_id = 'harrell_c'

  # check that the different methods return correctly subsetted tasks
  list1 = ts$subset_tasks(el, tsks) # top_n by default ($nfeats)
  expect_true(length(list1) == length(tsks))
  f1 = mlr3misc::map(list1, `[[`, 'feature_names') # which features were used
  nf1 = sapply(f1, length) # how many features
  expect_equal(nf1, fs_tbl$top_n)
  # correct features were used
  expect_equal(sort(f1[[1]]),
    fss1 %>% slice(1:nf1[1]) %>% pull(feat_name) %>% sort())
  expect_equal(sort(f1[[2]]),
    fss2 %>% slice(1:nf1[2]) %>% pull(feat_name) %>% sort())
  # backends are minimized (less features, less size)
  expect_true(list1[[1]]$backend$ncol < task1$backend$ncol)
  expect_true(list1[[2]]$backend$ncol < task2$backend$ncol)

  list2 = ts$subset_tasks(el, tsks, method = 'top_nperc')
  expect_true(length(list2) == length(tsks))
  f2 = mlr3misc::map(list2, `[[`, 'feature_names')
  nf2 = sapply(f2, length)
  expect_equal(nf2, fs_tbl$top_nperc)
  # check if the correct features were used
  expect_equal(sort(f2[[1]]),
    fss1 %>% slice(1:nf2[1]) %>% pull(feat_name) %>% sort())
  expect_equal(sort(f2[[2]]),
    fss2 %>% slice(1:nf2[2]) %>% pull(feat_name) %>% sort())
  # backends are minimized (less features, less size)
  expect_true(list2[[1]]$backend$ncol < task1$backend$ncol)
  expect_true(list2[[2]]$backend$ncol < task2$backend$ncol)

  list3 = ts$subset_tasks(el, tsks, method = 'gt_freq')
  expect_true(length(list3) == length(tsks))
  f3 = mlr3misc::map(list3, `[[`, 'feature_names')
  nf3 = sapply(f3, length)
  # replace 0 selected features (`fs_tbl$gt_freq[2]`) with all task's features
  expect_equal(nf3, c(fs_tbl$gt_freq[1], length(t2fs)))
  # correct features were used
  expect_equal(sort(f3[[1]]),
    fss1 %>% slice(1:nf3[1]) %>% pull(feat_name) %>% sort())
  expect_equal(sort(f3[[2]]), sort(t2fs)) # all task2 features
  # backends are minimized (less features, less size)
  expect_true(list3[[1]]$backend$ncol < task1$backend$ncol)
  expect_true(list3[[2]]$backend$ncol == task2$backend$ncol) # task didn't change

  list4 = ts$subset_tasks(el, tsks, method = 'optim_nfeats')
  expect_true(length(list4) == length(tsks))
  f4 = mlr3misc::map(list4, `[[`, 'feature_names')
  nf4 = sapply(f4, length)
  expect_equal(nf4, fs_tbl$optim_nfeats)
  # correct features were used
  expect_equal(sort(f4[[1]]),
    fss1 %>% slice(1:nf4[1]) %>% pull(feat_name) %>% sort())
  expect_equal(sort(f1[[2]]),
    fss2 %>% slice(1:nf4[2]) %>% pull(feat_name) %>% sort())
  # backends are minimized (less features, less size)
  expect_true(list4[[1]]$backend$ncol < task1$backend$ncol)
  expect_true(list4[[2]]$backend$ncol < task2$backend$ncol)

  # works if it's just 1 eFS + 1 task
  list5 = ts$subset_tasks(el[1], tsks[1])
  expect_class(list5, 'list')
  expect_true(length(list5) == 1)
})
