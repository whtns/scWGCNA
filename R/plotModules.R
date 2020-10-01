#' @title plotModules
#'
#' @description Use a bubbleplot to display expression levels and percentage for
#' genes identified as functioning within a network module.
#'
#' @param object scRNA-seq data object
#' @param modules Gene expression modules to plot, as determined by
#' seuratClusterWGCNA
#' @param ... Additional arguments to pass to bubbleplot
#'
#' @import patchwork
#' @importFrom ggplot2 labs xlab ylab theme aes
#' @importFrom SeuratBubblePlot bubbleplot
#' @importFrom gridExtra grid.arrange
#' @importFrom purrr reduce
#'
#' @return gg
#' @export
#'
#' @examples
plotModules <- function(object, modules, combine = TRUE, ...) {
  modplots <- lapply(
    (1:length(modules)),
    function(x) {
      bp <- bubbleplot(
        object = object,
        features_plot = modules[[x]],
        do.return = TRUE,
        ...
      ) +
        labs(title = names(modules)[[x]],
             x = NULL,
             y = NULL) +
        theme(legend.position = "none")
    }
  )

  if (combine){
    patchwork::wrap_plots(modplots)
  } else {
    modplots
  }
}
