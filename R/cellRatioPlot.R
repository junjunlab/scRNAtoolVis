#' @name cellRatioPlot
#' @author Junjun Lao
#' @title using cellRatioPlot to visualize cell percent ratio across clusters
#'
#' @param object object seurat object.
#' @param sample.name sample name in meta info, default NULL.
#' @param celltype.name celltype name in meta info, default NULL.
#' @param col.width column width, default 0.7.
#' @param flow.alpha flow color alpha, default 0.25.
#' @param flow.curve flow curve rate, default 0.
#' @param fill.col supply your own fill colors, default NULL.
#' @param sample.order the order for samples, default NULL.
#'
#' @return a ggplot object.
#' @export

globalVariables(c("n", "num"))

cellRatioPlot <- function(
    object = NULL,
    sample.name = NULL,
    celltype.name = NULL,
    sample.order = NULL,
    col.width = 0.7,
    flow.alpha = 0.25,
    flow.curve = 0,
    fill.col = NULL) {
  # get metainfo
  meta <- object@meta.data

  # order
  if(!is.null(sample.order)){
    meta$sample.name <- factor(meta$sample.name,levels = sample.order)
  }

  # calculate percent ratio
  ratio.info <- meta %>%
    dplyr::group_by(.data[[sample.name]], .data[[celltype.name]]) %>%
    dplyr::summarise(num = n()) %>%
    dplyr::mutate(rel_num = num / sum(num))

  # color
  if (is.null(fill.col)) {
    fill.col <- jjAnno::useMyCol("paired", n = length(unique(meta[, celltype.name])))
  } else {
    fill.col <- fill.col
  }

  # plot
  p <-
    ggplot2::ggplot(
      ratio.info,
      ggplot2::aes_string(x = sample.name, y = "rel_num")
    ) +
    ggplot2::geom_col(
      ggplot2::aes_string(fill = celltype.name),
      width = col.width
    ) +
    ggalluvial::geom_flow(
      ggplot2::aes_string(
        stratum = celltype.name,
        alluvium = celltype.name,
        fill = celltype.name
      ),
      width = 0.5,
      alpha = flow.alpha,
      knot.pos = flow.curve
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(
      values = fill.col,
      name = "Cell Type"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      legend.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("Cell percent ratio")

  return(p)
}
