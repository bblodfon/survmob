#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table data.table
#' @importFrom dplyr filter
#' @importFrom lgr get_logger
#' @importFrom mlr3 lrn
#' @importFrom mlr3 msr
#' @importFrom mlr3 tsk
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
## usethis namespace: end
NULL

register_mlr3 = function() {
  x = utils::getFromNamespace('mlr_measures', ns = 'mlr3')
  x$add('surv.oob_error', MeasureSurvOOBError)

  x = utils::getFromNamespace('mlr_pipeops', ns = 'mlr3pipelines')
  x$add('survshuffle', PipeOpSurvShuffle)
  x$add('removenas', PipeOpRemoveNAs)
  x$add('removezeros', PipeOpRemoveZeros)
}

.onLoad = function(libname, pkgname) {
  register_mlr3()
}
