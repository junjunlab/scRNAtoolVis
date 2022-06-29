#' @name FeatureCornerAxes
#' @title Add corner axes on seurat UMAP/tSNE gene FeaturePlot function figures
#' @param object object seurat object.
#' @param reduction "string",reduction type (umap/tsne).
#' @param features "string",the gene you want to plot.
#' @param groupFacet "string",give the column name in seurat metadata to facet plot.
#' @param relLength "num",the corner axis line relative length to plot axis(0-1).
#' @param relDist "num",the relative distance of corner axis label to axis.
#' @param low "string",point color with low expression.
#' @param high "string",point color with high expression.
#' @param axes "string",show multiple corner axis or only one (mul/one),defaults "mul".
#' @param legendPos "string",legend position same as ggplot theme function,defaults "right".
#' @param RowSample "num",rows for groupFacet.
#' @param ColGene "num",cols for features.
#' @param stripCol "string",facet balckground color,defaults "white".
#' @param pSize "num",point size.
#' @param arrowType "string",arrow type (open/closed),defaults "closed".
#' @param lineTextcol "string",facet balckground color,defaults "white".
#' @param cornerTextSize "num", the corner label text size, defaults is 5.
#' @param base_size "num", theme base size, defaults is 14.
#' @param themebg Another theme style, defaults is 'default', or 'bwCorner'.
#' @return Return a ggplot.
#' @export
#' @examples
#'
#' test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
#'
#' tmp <- readRDS(test)
#'
#' # umap
#' FeatureCornerAxes(
#'   object = tmp, reduction = "umap",
#'   groupFacet = "orig.ident",
#'   relLength = 0.5, relDist = 0.2,
#'   features = c("Actb", "Ythdc1", "Ythdf2")
#' )
#'
#' # one axes
#' FeatureCornerAxes(
#'   object = tmp, reduction = "umap",
#'   groupFacet = "orig.ident",
#'   features = c("Actb", "Ythdc1", "Ythdf2"),
#'   relLength = 0.5, relDist = 0.2,
#'   axes = "one",
#'   lineTextcol = "grey50"
#' )
#'
#' # tsne
#' FeatureCornerAxes(
#'   object = tmp, reduction = "tsne",
#'   groupFacet = "orig.ident",
#'   relLength = 0.5, relDist = 0.2,
#'   features = c("Actb", "Ythdc1", "Ythdf2")
#' )
#'
#'
# define viriables
globalVariables(c("x1", "y1", "linegrou", "angle", "lab", "gene_name", "scaledValue"))

# define function
FeatureCornerAxes <- function(object = NULL,
                              reduction = "umap",
                              features = NULL,
                              groupFacet = "orig.ident",
                              relLength = 0.25,
                              relDist = 0.1,
                              low = "lightgrey",
                              high = "red",
                              axes = "mul",
                              legendPos = "right",
                              RowSample = 1,
                              ColGene = 1,
                              stripCol = "white",
                              pSize = 1,
                              arrowType = "closed",
                              lineTextcol = "black",
                              cornerTextSize = 5,
                              base_size = 14,
                              themebg = "default") {
  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)

  # get gene expression
  geneExp <- Seurat::FetchData(object = object, vars = features)

  # cbind
  mer <- cbind(pc12, geneExp)

  # merge data
  megredf <- reshape2::melt(
    mer,
    id.vars = colnames(pc12),
    variable.name = "gene_name",
    value.name = "scaledValue"
  )

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
  if (reduction == "umap") {
    axs_label <- paste("UMAP", 2:1, sep = "")
  } else if (reduction == "tsne") {
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
      angle = c(90, 0),
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

  # sample
  genes <- unique(megredf$gene_name)

  # bacth plot
  lapply(genes, function(x) {
    tmpdf <- dplyr::filter(megredf, gene_name == x)

    # plot
    p1 <- ggplot2::ggplot(
      tmpdf,
      ggplot2::aes(x = tmpdf[, 1], y = tmpdf[, 2])
    ) +
      ggplot2::geom_point(ggplot2::aes(color = scaledValue), size = pSize) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::scale_color_gradient(
        name = "",
        low = low,
        high = high
      ) +
      ggplot2::labs(x = "", y = x) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(colour = NA, fill = stripCol),
        aspect.ratio = 1,
        legend.position = legendPos,
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
      ) +
      ggplot2::facet_wrap(~get(groupFacet), nrow = RowSample) +
      ggplot2::geom_line(
        data = axes,
        color = lineTextcol,
        ggplot2::aes(x = x1, y = y1, group = linegrou),
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
        fontface = "italic",
        size = cornerTextSize
      )

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
      return(p2)
    } else if (themebg == "default") {
      return(p1)
    }
  }) -> plst

  # combine
  p3 <- patchwork::wrap_plots(plst, ncol = ColGene, guides = "collect")
  return(p3)
}
