#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom combinat combn
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr slice
#' @importFrom forcats fct_reorder
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom lgr get_logger
#' @importFrom mlr3 assert_learner
#' @importFrom mlr3 assert_task
#' @importFrom mlr3 lrn
#' @importFrom mlr3 msr
#' @importFrom mlr3 partition
#' @importFrom mlr3 tsk
#' @importFrom mlr3fselect AutoFSelector
#' @importFrom mlr3fselect callback_fselect
#' @importFrom mlr3fselect fs
#' @importFrom mlr3misc map_chr
#' @importFrom mlr3misc map_dbl
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
#' @importFrom scales label_percent
#' @importFrom stabm stabilityJaccard
#' @importFrom stabm stabilityNogueira
#' @importFrom tibble as_tibble
#' @importFrom tibble as_tibble_row
#' @importFrom tibble tibble
#' @importFrom tictoc tic
#' @importFrom tictoc toc
#' @importFrom tidyr pivot_longer
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
