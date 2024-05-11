#' @name clusterCornerAxes
#' @author Junjun Lao
#' @title Add corner axes on seurat UMAP/tSNE cluster figures
#' @param object seurat object.
#' @param reduction "string", reduction type (umap/tsne).
#' @param groupFacet "string", give the column name in seurat metadata to facet plot.
#' @param clusterCol "string", the point color to group by,cluster name, default "seurat_clusters".
#' @param pSize "num", point size.
#' @param aspect.ratio "num", plot width and height ratio, default NULL.
#' @param noSplit 'logic', whether to split/facet the plot, default "TRUE".
#' @param nrow "num", rows to plot when noSplit = FALSE.
#' @param relLength 'num', the corner axis line relative length to plot axis(0-1).
#' @param relDist "num" ,the relative distance of corner axis label to axis.
#' @param axes "string", show multiple corner axis or only one (mul/one), default "mul".
#' @param legendPos "string", legend position same as ggplot theme function, default "right".
#' @param keySize The legned point size, default is 5.
#' @param lineTextcol "string", corner line and label color, default "black".
#' @param stripCol "string", facet balckground color, default "white".
#' @param arrowType "string", arrow type (open/closed), default "closed".
#' @param cornerTextSize "num", the corner label text size, default is 3.
#' @param base_size "num", theme base size, default is 14.
#' @param themebg Another theme style, default is 'default', or 'bwCorner'.
#' @param addCircle Logic, whether add circle on clusters, default is 'FALSE'.
#' @param cicAlpha "num", circle fill color alpha, default is 0.1,
#' @param cicLineSize "num", circle line size, default is 1.
#' @param cicLineColor "num", circle line color, default is 'grey50'.
#' @param cicLineLty "num", circle line type, default is 'dashed'.
#' @param nbin "num", number of points used to shape the hull, default 100.
#' @param nsm "num", number of points used to perform convolution, should less than nbin, default 10.
#' @param addsm "num", number of additional times of convolution performed, default 1.
#' @param sfac "num", quantile of each sector, used to determine the edge of the hull, should less than 1, default 1.
#' @param qval "num", expansion size factor, larger value means bigger hull, default 1.5.
#'
#' @param cellLabel Whether to label cell type on plot, default is FALSE.
#' @param cellLabelSize Cell type label size, default is 6.
#' @param cellLabelColor Cell type label color, default is "balck".
#' @param show.legend Wheher show legend, default is TRUE.
#'
#' @return Return a ggplot object.
#' @export
#' @examples
#' test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
#'
#' tmp <- readRDS(test)
#'
#' # umap
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = TRUE
#' )
#'
#' # arrowType
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = TRUE, arrowType = "open"
#' )
#'
#' # facet by metadata column "orig.ident"
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5
#' )
#'
#' # retain only one axes
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5,
#'   axes = "one"
#' )
#'
#' # line color
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5,
#'   lineTextcol = "grey50"
#' )
#'
#' # tsne
#' clusterCornerAxes(
#'   object = tmp, reduction = "tsne",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5
#' )
#'
# define viriables
globalVariables(c("x1", "y1", "linegrou", "angle", "lab", ".data"))

# define function
clusterCornerAxes <- function(
    object = NULL,
    reduction = "umap",
    groupFacet = groupFacet,
    clusterCol = "seurat_clusters",
    pSize = 1,
    aspect.ratio = NULL,
    noSplit = TRUE,
    nrow = 1,
    relLength = 0.25,
    relDist = 0.1,
    axes = "mul",
    show.legend = TRUE,
    legendPos = "right",
    keySize = 5,
    cellLabel = FALSE,
    cellLabelSize = 6,
    cellLabelColor = "black",
    lineTextcol = "black",
    stripCol = "white",
    arrowType = "closed",
    cornerTextSize = 3,
    base_size = 14,
    themebg = "default",
    addCircle = FALSE,
    cicAlpha = 0.1,
    cicLineSize = 1,
    cicLineColor = "grey50",
    cicLineLty = "dashed",
    nbin = 100,
    nsm = 10,
    addsm = 1,
    qval = 1,
    sfac = 1.5) {
  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)

  #######################################
  # text data
  namePos <- pc12 %>%
    dplyr::group_by(.data[[clusterCol]]) %>%
    dplyr::summarise(
      posMedia1 = stats::median(get(colnames(pc12)[1])),
      posMedia2 = stats::median(get(colnames(pc12)[2]))
    )

  #######################################

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
  if (startsWith(reduction, "umap")) {
    axs_label <- paste("UMAP", 2:1, sep = "")
  } else if (startsWith(reduction, "tsne")) {
    axs_label <- paste("t-SNE", 2:1, sep = "")
  } else {
    print("Please give correct type(umap or tsne)!")
  }

  if (axes == "mul") {
    # axies data
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2)
    )
    # axies label
    label <- data.frame(
      "lab" = c(axs_label),
      "angle" = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel)
    )
  } else if (axes == "one") {
    firstFacet <- unique(pc12[, groupFacet])[1]
    # axies data
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2),
      "group" = rep(firstFacet, 2)
    )
    # axies label
    label <- data.frame(
      "lab" = c(axs_label),
      "angle" = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel),
      "group" = rep(firstFacet, 2)
    )

    # rename group name
    colnames(axes)[4] <- groupFacet
    colnames(label)[5] <- groupFacet
  } else {
    print("Please give correct args(mul or one)!")
  }

  ######################################################
  # plot
  p <- ggplot2::ggplot(
    pc12,
    ggplot2::aes_string(x = colnames(pc12)[1], y = colnames(pc12)[2])
  ) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = clusterCol),
      size = pSize,
      show.legend = show.legend
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(colour = NA, fill = stripCol),
      aspect.ratio = aspect.ratio,
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
      ggplot2::aes(x = x1, y = y1, angle = angle, label = lab),
      color = lineTextcol,
      fontface = "italic",
      size = cornerTextSize
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = keySize)))

  ######################################################
  # add text label
  if (cellLabel == FALSE) {
    plabel <- p
  } else {
    plabel <- p +
      ggrepel::geom_text_repel(
        data = namePos,
        ggplot2::aes_string(x = "posMedia1", y = "posMedia2", label = clusterCol),
        show.legend = FALSE,
        size = cellLabelSize,
        color = cellLabelColor
      )
  }

  ######################################################
  # add circle line
  if (addCircle == FALSE) {
    p0 <- plabel
  } else {
    p0 <- plabel +
      ggunchull::stat_unchull0(
        ggplot2::aes_string(fill = clusterCol),
        alpha = cicAlpha,
        size = cicLineSize,
        color = cicLineColor,
        lty = cicLineLty,
        show.legend = FALSE,
        nbin = nbin,
        nsm = nsm,
        addsm = addsm,
        sfac = sfac,
        qval = qval
      )
  }

  ######################################################
  # facet plot
  if (noSplit == TRUE) {
    p1 <- p0
  } else {
    p1 <- p0 + ggplot2::facet_wrap(facets = groupFacet, nrow = nrow)
  }

  ######################################################
  # theme style
  if (themebg == "bwCorner") {
    p2 <- p1 +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        aspect.ratio = 1,
        strip.background = ggplot2::element_rect(colour = NA, fill = stripCol)
      )
  } else if (themebg == "default") {
    p2 <- p1
  }

  return(p2)
}
