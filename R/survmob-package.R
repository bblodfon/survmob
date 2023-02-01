#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table data.table
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom dplyr mutate
#' @importFrom lgr get_logger
#' @importFrom mlr3 lrn
#' @importFrom mlr3 msr
#' @importFrom mlr3 tsk
#' @importFrom mlr3fselect AutoFSelector
#' @importFrom mlr3fselect callback_fselect
#' @importFrom mlr3fselect fs
#' @importFrom mlr3misc open_help
#' @importFrom mlr3pipelines %>>%
#' @importFrom mlr3pipelines PipeOpTaskPreproc
#' @importFrom mlr3pipelines po
#' @importFrom mlr3proba as_task_surv
#' @importFrom paradox p_dbl
#' @importFrom paradox p_fct
#' @importFrom paradox p_int
#' @importFrom paradox p_lgl
#' @importFrom paradox ps
#' @importFrom parallelly availableCores
#' @importFrom R6 R6Class
#' @importFrom tibble tibble
## usethis namespace: end
NULL

register_mlr3 = function() {
  x = utils::getFromNamespace('mlr_pipeops', ns = 'mlr3pipelines')
  x$add('survshuffle', PipeOpSurvShuffle)
  x$add('removenas', PipeOpRemoveNAs)
  x$add('removezeros', PipeOpRemoveZeros)
  x$add('logtransform', PipeOpLogTransform)
}

.onLoad = function(libname, pkgname) {
  register_mlr3()
}
