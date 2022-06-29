#' @name AverageHeatmap
#' @title Plot averaged gene expression cross cluster cells
#' @param object object seurat object.
#' @param markerGene Your marker genes.
#' @param group.by Categories for grouping (e.g, ident, replicate, celltype). 'ident' by default.
#' @param assays Which assays to use. Default is "RNA" assays.
#' @param slot Slot(s) to use. Default is "data"
#' @param htCol Heatmap colors. Default is c("#0099CC", "white", "#CC0033").
#' @param colseed Cluster annotaion colors seed, these colors are produed randomly, so you can give a seed to assure produce same colors.  Default is 666.
#' @param htRange Heatmap values range. Default is c(-2, 0, 2).
#' @param annoCol Whether use your own annotation clusters colors. Default is 'FALSE'.
#' @param myanCol You can specify your own annotation clusters colors vectors. Default is 'null'.
#' @param annoColType Cluster annotaion colors type (bright, light, dark and random). Default is light.
#' @param annoColTypeAlpha Cluster annotaion colors transparency. Default is 0.
#' @param row_title Heatmap row title. Default is "Cluster top Marker genes".
#' @param row_names_side Heatmap gene name side. Default is "left".
#' @param border Whether to shOw heatmap border. Default is "FALSE".
#' @param fontsize Heatmap gene name fontsize. Default is 10.
#' @param column_names_rot Cluster name rotation. Default is 45.
#'
#' @return
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
#' AverageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene
#' )
#'
#' # change color
#' AverageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene,
#'   htCol = c("#339933", "#FFCC00", "#FF0033")
#' )
#'
# define viriables
AverageHeatmap <- function(object = NULL,
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
                           row_names_side = "left",
                           border = FALSE,
                           fontsize = 10,
                           column_names_rot = 45) {
  # get cells mean gene expression
  mean_gene_exp <- as.matrix(data.frame(Seurat::AverageExpression(object,
    features = markerGene,
    group.by = group.by,
    assays = assays,
    slot = slot
  )))

  # add colnames
  name1 <- gsub(
    pattern = paste0(assays, ".", sep = ""),
    replacement = "",
    colnames(mean_gene_exp)
  )

  colnames(mean_gene_exp) <- gsub(
    pattern = "\\.",
    replacement = " ", name1
  )

  # Z-score
  htdf <- t(scale(t(mean_gene_exp), scale = T, center = T))

  # color
  col_fun <- circlize::colorRamp2(htRange, htCol)

  # anno color
  if (annoCol == FALSE) {
    set.seed(colseed)
    anno_col <- circlize::rand_color(ncol(htdf),
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
    show_legend = F,
    col = list(cluster = anno_col)
  )

  # plot
  ComplexHeatmap::Heatmap(htdf,
    name = "Z-score",
    cluster_columns = F,
    cluster_rows = F,
    row_title = row_title,
    # column_title = "Clusters",
    row_names_gp = grid::gpar(fontface = "italic", fontsize = fontsize),
    row_names_side = row_names_side,
    border = border,
    column_names_side = "top",
    column_names_rot = column_names_rot,
    top_annotation = column_ha,
    col = col_fun
  )
}
