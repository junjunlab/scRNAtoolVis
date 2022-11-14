#' @name drawLegend
#' @author Junjun Lao
#' @title using drawLegend to add cluster id in legend
#'
#' @param object object seurat object, default NULL.
#' @param plot ggplot object, default NULL.
#' @param cellType cell type column name, default NULL.
#' @param clusters cluster id column name, default NULL.
#' @param ncol ncols to draw legend(1/2), default 1.
#' @param col point color, default hue_pal().
#' @param pt.size point size, default 8.
#' @param text.size text size, default 4.
#'
#' @return combine plot
#' @export

globalVariables(c("x","y"))

drawLegend <- function(object = NULL,
                       plot = NULL,
                       cellType = NULL,
                       clusters = NULL,
                       ncol = 1,
                       col = NULL,
                       pt.size = 8,
                       text.size = 4){
  # prepare data
  leg.data <- object@meta.data %>%
    dplyr::select(.data[[cellType]],.data[[clusters]]) %>%
    unique()

  colnames(leg.data) <- c("cellType","clusters")

  # reorder
  leg.data <- leg.data[match(levels(leg.data$cellType),leg.data$cellType),]

  # add xy position
  if(ncol > 1){
    leg.data$x <- rep(1:ncol,c(ceiling(nrow(leg.data)/ncol),
                               nrow(leg.data) - ceiling(nrow(leg.data)/ncol)))

    leg.data$y <- c(1:ceiling(nrow(leg.data)/ncol),
                    1:(nrow(leg.data) - ceiling(nrow(leg.data)/ncol)))
  }else{
    leg.data$x <- 1
    leg.data$y <- 1:nrow(leg.data)
  }

  # order
  leg.data$cellType <- factor(leg.data$cellType,levels = rev(levels(leg.data$cellType)))

  # plot
  if(is.null(col)){
    color <- rev(scales::hue_pal()(nrow(leg.data)))
  }else{
    color <- rev(col)
  }

  pleg <-
    ggplot2::ggplot(leg.data,ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = cellType),
                        show.legend = F,
                        size = pt.size) +
    ggplot2::geom_text(ggplot2::aes(label = clusters)) +
    ggplot2::geom_text(ggplot2::aes(label = cellType),
                       hjust = 0,
                       nudge_x = 0.2,
                       size = text.size) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::scale_y_reverse() +
    ggplot2::xlim(0,ncol + 1) +
    ggplot2::theme_void()

  # COMBINE
  cowplot::plot_grid(plot,pleg)
}
