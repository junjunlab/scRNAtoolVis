#' @name averageHeatmap
#' @author Junjun Lao
#' @title Plot averaged gene expression cross cluster cells
#'
#' @param object object seurat object.
#' @param markerGene Your marker genes.
#' @param group.by Categories for grouping (e.g, ident, replicate, celltype). "ident" by default.
#' @param assays Which assays to use. Default is "RNA" assays.
#' @param slot Slot(s) to use. Default is "data".
#' @param htCol Heatmap colors. Default is c("#0099CC", "white", "#CC0033").
#' @param colseed Cluster annotation colors seed, these colors are produced randomly, so you can give a seed to assure produce same colors.  Default is 666.
#' @param htRange Heatmap values range. Default is c(-2, 0, 2).
#' @param annoCol Whether use your own annotation clusters colors. Default is "FALSE".
#' @param myanCol You can specify your own annotation clusters colors vectors. Default is "null".
#' @param annoColType Cluster annotation colors type (bright, light, dark and random). Default is light.
#' @param annoColTypeAlpha Cluster annotation colors transparency. Default is 0.
#' @param row_title Heatmap row title. Default is "Cluster top Marker genes".
#' @param row_names_side Heatmap gene name side. Default is "left".
#' @param border Whether to shOw heatmap border. Default is "FALSE".
#' @param fontsize Heatmap gene name fontsize. Default is 10.
#' @param column_names_rot Cluster name rotation. Default is 45.
#' @param showRowNames whether to show rownames. Default is "TRUE".
#' @param markGenes Provide your tartget genes to mark on the plot. Default is "NULL".
#' @param clusterAnnoName Whether to add clsuetr column annotation name. Default is "TRUE".
#' @param width The heatmap body width. Default is "NULL".
#' @param height The heatmap body height. Default is "NULL".
#' @param cluster.order The cell clusters order. Default is "NULL".
#'
#' @param cluster_columns Whether cluster columns. Default is "FALSE".
#' @param cluster_rows Whether cluster rows. Default is "FALSE".
#' @param gene.order the gene orders for heatmap. Default is "NULL".
#'
#' @param ... Other arguments passed with ComplexHeatmap::rowAnnotation and ComplexHeatmap::Heatmap.
#' @return Return a plot.
#' @export
#'
#' @examples
#' httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
#' pbmc <- readRDS(httest)
#'
#' # load markergene
#' markergene <- system.file("extdata", "top5pbmc.markers.csv", package = "scRNAtoolVis")
#' markers <- read.table(markergene, sep = ",", header = TRUE)
#'
#' # plot
#' averageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene
#' )
#'
#' # change color
#' averageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene,
#'   htCol = c("#339933", "#FFCC00", "#FF0033")
#' )
#'
# define function
averageHeatmap <- function(
    object = NULL,
    markerGene = NULL,
    group.by = "ident",
    assays = "RNA",
    slot = "data",
    htCol = c("#0099CC", "white", "#CC0033"),
    colseed = 666,
    htRange = c(-2, 0, 2),
    annoCol = FALSE,
    myanCol = NULL,
    annoColType = "light",
    annoColTypeAlpha = 0,
    row_title = "Cluster top Marker genes",
    clusterAnnoName = TRUE,
    showRowNames = TRUE,
    row_names_side = "left",
    markGenes = NULL,
    border = FALSE,
    fontsize = 10,
    column_names_rot = 45,
    width = NULL,
    height = NULL,
    cluster.order = NULL,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    gene.order = NULL,
    ...) {
  # get cells mean gene expression
  # check Seurat version first
  vr <- utils::compareVersion(as.character(utils::packageVersion("Seurat")),"5")
  if(vr == 1){
    mean_gene_exp <- as.matrix(
      data.frame(
        Seurat::AverageExpression(object,
                                  features = markerGene,
                                  group.by = group.by,
                                  assays = assays,
                                  layer = slot
        )
      )
    )
  }else{
    mean_gene_exp <- as.matrix(
      data.frame(
        Seurat::AverageExpression(object,
                                  features = markerGene,
                                  group.by = group.by,
                                  assays = assays,
                                  slot = slot
        )
      )
    )
  }

  # add colnames
  # name1 <- gsub(
  #   pattern = paste0(assays, ".", sep = ""),
  #   replacement = "",
  #   colnames(mean_gene_exp)
  # )
  #
  # colnames(mean_gene_exp) <- gsub(
  #   pattern = "\\.",
  #   replacement = " ", name1
  # )

  colnames(mean_gene_exp) <- levels(Seurat::Idents(object))

  # Z-score
  htdf <- t(scale(t(mean_gene_exp), scale = TRUE, center = TRUE))

  # cluster order
  if (!is.null(cluster.order)) {
    htdf <- htdf[, cluster.order]
  }

  # gene order
  if (!is.null(gene.order)) {
    htdf <- htdf[gene.order, ]
  }

  # color
  col_fun <- circlize::colorRamp2(htRange, htCol)

  # anno color
  if (annoCol == FALSE) {
    set.seed(colseed)
    anno_col <- circlize::rand_color(
      ncol(htdf),
      luminosity = annoColType,
      transparency = annoColTypeAlpha
    )
    print(c("Your cluster annotation color is:", anno_col))
  } else if (annoCol == TRUE) {
    # give your own color vectors
    anno_col <- myanCol
  } else {
    print("Give TRUE or FALSE paramters!")
  }
  names(anno_col) <- colnames(htdf)

  # top annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    cluster = colnames(htdf),
    show_legend = FALSE,
    show_annotation_name = clusterAnnoName,
    col = list(cluster = anno_col)
  )

  # whether mark your genes on plot
  if (!is.null(markGenes)) {
    # all genes
    rowGene <- rownames(htdf)

    # tartget gene
    annoGene <- markGenes

    # get target gene index
    index <- match(annoGene, rowGene)

    # some genes annotation
    geneMark <- ComplexHeatmap::rowAnnotation(
      gene = ComplexHeatmap::anno_mark(
        at = index,
        labels = annoGene,
        labels_gp = grid::gpar(
          fontface = "italic",
          fontsize = fontsize
        )
      ),
      ...
    )

    right_annotation <- geneMark
  } else {
    right_annotation <- NULL
  }

  # control heatmap width and height
  if (is.null(width) || is.null(height)) {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      name = "Z-score",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows,
      row_title = row_title,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      ...
    )
  } else {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      name = "Z-score",
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      row_title = row_title,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      width = ggplot2::unit(width, "cm"),
      height = ggplot2::unit(height, "cm"),
      ...
    )
  }
}
