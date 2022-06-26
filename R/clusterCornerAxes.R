#' @name clusterCornerAxes
#' @title Add corner axes on seurat UMAP/tSNE cluster figures
#' @param object seurat object.
#' @param reduction "string", reduction type (umap/tsne).
#' @param groupFacet "string", give the column name in seurat metadata to facet plot.
#' @param clusterCol "string", the point color to group by,cluster name, defaults "seurat_clusters".
#' @param pSize "num", point size.
#' @param noSplit 'logic', whether to split/facet the plot, defaults "TRUE".
#' @param nrow "num", rows to plot when noSplit = FALSE.
#' @param relLength 'num', the corner axis line relative length to plot axis(0-1).
#' @param relDist "num" ,the relative distance of corner axis label to axis.
#' @param axes "string", show multiple corner axis or only one (mul/one), defaults "mul".
#' @param legendPos "string", legend position same as ggplot theme function, defaults "right".
#' @param lineTextcol "string", corner line and label color, defaults "black".
#' @param stripCol "string", facet balckground color, defaults "white".
#' @param arrowType "string", arrow type (open/closed), defaults "closed".
#' @param cornerTextSize "num", the corner label text size, defaults is 5.
#' @param base_size "num", theme base size, defaults is 14.
#' @param themebg Another theme style, defaults is 'default', or 'bwCorner'.
#' @return
#' @export
#' @examples
#' test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
#'
#' tmp <- readRDS(test)
#'
#' # umap
#' clusterCornerAxes(object = tmp,reduction = 'umap',
#'                   noSplit = TRUE)
#'
#' # arrowType
#' clusterCornerAxes(object = tmp,reduction = 'umap',
#'                   noSplit = TRUE,arrowType = 'open')
#'
#' # facet by metadata column "orig.ident"
#' clusterCornerAxes(object = tmp,reduction = 'umap',
#'                   noSplit = FALSE,groupFacet = 'orig.ident',
#'                   relLength = 0.5)
#'
#' # retain only one axes
#' clusterCornerAxes(object = tmp,reduction = 'umap',
#'                   noSplit = FALSE,groupFacet = 'orig.ident',
#'                   relLength = 0.5,
#'                   axes = 'one')
#'
#' # line color
#' clusterCornerAxes(object = tmp,reduction = 'umap',
#'                   noSplit = FALSE,groupFacet = 'orig.ident',
#'                   relLength = 0.5,
#'                   lineTextcol = 'grey50')
#'
#' # tsne
#' clusterCornerAxes(object = tmp,reduction = 'tsne',
#'                   noSplit = FALSE,groupFacet = 'orig.ident',
#'                   relLength = 0.5)

# define viriables
globalVariables(c("x1", "y1", "linegrou","angle","lab"))

# define function
clusterCornerAxes <- function(object,
                              reduction = 'umap',
                              groupFacet = groupFacet,
                              clusterCol = 'seurat_clusters',
                              pSize = 1,
                              noSplit = TRUE,
                              nrow = 1,
                              relLength = 0.25,
                              relDist = 0.1,
                              axes = 'mul',
                              legendPos = 'right',
                              lineTextcol = 'black',
                              stripCol = 'white',
                              arrowType = 'closed',
                              cornerTextSize = 5,
                              base_size = 14,
                              themebg = 'default') {
  # make PC data
  reduc <-
    data.frame(Seurat::Embeddings(object, reduction = reduction))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)

  # data range
  range <- floor(min(min(pc12[, 1]), min(pc12[, 2])))

  # get botomn-left coord
  lower <- range - relDist * abs(range)

  # label reldist to axes
  labelRel <- relDist * abs(lower)

  # get relative line length
  linelen <- abs(relLength * lower) + lower

  # mid point
  mid <- abs(relLength * lower) / 2 + lower

  # give reduction type
  if (reduction == 'umap') {
    axs_label <- paste('UMAP', 2:1, sep = '')
  } else if (reduction == 'tsne') {
    axs_label <- paste('t-SNE', 2:1, sep = '')
  } else{
    print('Please give correct type(umap or tsne)!')
  }

  if (axes == 'mul') {
    # axies data
    axes <- data.frame(
      'x1' = c(lower, lower, lower, linelen),
      'y1' = c(lower, linelen, lower, lower),
      'linegrou' = c(1, 1, 2, 2)
    )
    # axies label
    label <- data.frame(
      'lab' = c(axs_label),
      'angle' = c(90, 0),
      'x1' = c(lower - labelRel, mid),
      'y1' = c(mid, lower - labelRel)
    )
  } else if (axes == 'one') {
    firstFacet <- unique(pc12[, groupFacet])[1]
    # axies data
    axes <- data.frame(
      'x1' = c(lower, lower, lower, linelen),
      'y1' = c(lower, linelen, lower, lower),
      'linegrou' = c(1, 1, 2, 2),
      'group' = rep(firstFacet, 2)
    )
    # axies label
    label <- data.frame(
      'lab' = c(axs_label),
      'angle' = c(90, 0),
      'x1' = c(lower - labelRel, mid),
      'y1' = c(mid, lower - labelRel),
      'group' = rep(firstFacet, 2)
    )

    # rename group name
    colnames(axes)[4] <- groupFacet
    colnames(label)[5] <- groupFacet
  } else{
    print('Please give correct args(mul or one)!')
  }

  # plot
  p <- ggplot2::ggplot(pc12,
                       ggplot2::aes(x = pc12[, 1], y = pc12[, 2])) +
    ggplot2::geom_point(ggplot2::aes_string(color = clusterCol),
                        size = pSize) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = '', y = '') +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(colour = NA, fill = stripCol),
      aspect.ratio = 1,
      legend.position = legendPos,
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    ) +
    ggplot2::geom_line(
      data = axes,
      ggplot2::aes(x = x1, y = y1, group = linegrou),
      color = lineTextcol,
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.1, "inches"),
        ends = "last",
        type = arrowType
      )
    ) +
    ggplot2::geom_text(
      data = label,
      color = lineTextcol,
      ggplot2::aes(
        x = x1,
        y = y1,
        angle = angle,
        label = lab
      ),
      fontface = 'italic',
      size = cornerTextSize,
    )

  # facet plot
  if (noSplit == TRUE) {
    p1 <- p
  } else{
    p1 <- p + ggplot2::facet_wrap( ~ get(groupFacet), nrow = nrow)
  }

  # theme style
  if(themebg == 'bwCorner'){
    p2 <- p1 +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     aspect.ratio = 1,
                     strip.background = ggplot2::element_rect(colour = NA,fill = stripCol))
    return(p2)
  }else if(themebg == 'default'){
    return(p1)
  }
}
