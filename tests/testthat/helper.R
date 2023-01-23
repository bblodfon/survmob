library(mlr3proba) # loads mlr3 as well
library(mlr3tuning)
library(mlr3pipelines)
library(mlr3extralearners)
library(mlr3misc)
library(checkmate) # for more expect_*() functions

# less logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# code/objects available to all tests
## veteran task
taskv = as_task_surv(x = survival::veteran,
  time = 'time', event = 'status')
poe = po('encode')
taskv = poe$train(list(taskv))[[1L]]
