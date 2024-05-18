#' @name tracksPlot
#' @author Jun Zhang
#' @title This function generates a track or heatmap plot based on the provided data.
#'
#' @param object An optional object containing the data.
#' @param plot.type The type of plot to generate, either "track" or "heatmap".
#' @param genes A vector or data frame specifying the genes to plot.
#' @param vmin The minimum value for the color scale (only applicable for the heatmap plot).
#' @param vmax The maximum value for the color scale (only applicable for the heatmap plot).
#' @param cell.order An optional vector specifying the order of cells in the plot.
#' @param gene.order An optional vector specifying the order of genes in the plot.
#' @param facet_nested_params A list of additional parameters to customize the facet_nested plot.
#' @param theme_params A list of additional parameters to customize the plot's theme.
#' @param strip_nested_params A list of additional parameters to customize the strip_nested plot.
#'
#' @return A ggplot object representing the track or heatmap plot.
#'
#' @export

globalVariables(c("barcode", "distinct"))

tracksPlot <- function(
    object = NULL,
    plot.type = c("track", "heatmap"),
    genes = NULL,
    vmin = -2, vmax = 2,
    cell.order = NULL,
    gene.order = NULL,
    facet_nested_params = list(),
    theme_params = list(),
    strip_nested_params = list()) {
  plot.type <- match.arg(plot.type, c("track", "heatmap"))

  # check markers gene
  if (is.data.frame(genes)) {
    markers_tmp <- genes |> dplyr::mutate(gene = make.unique(gene))

    markers <- markers_tmp$gene
    names(markers) <- markers_tmp$cluster
  } else {
    markers <- genes
  }

  # get barcode info
  barcode_info <- data.frame(Seurat::Idents(object))
  barcode_info$barcode <- rownames(barcode_info)
  colnames(barcode_info)[1] <- "cell"

  # get normalized matrix
  df <- data.frame(t(as.matrix(object@assays$RNA@data)), check.names = FALSE)[, markers]

  # do zscore
  if (plot.type == "heatmap") {
    df <- scale(df, center = TRUE) |>
      data.frame(check.names = FALSE)
    df <- apply(df, c(1, 2), function(x) {
      if (x > vmax) {
        x <- vmax
      } else if (x < vmin) {
        x <- vmin
      } else {
        x
      }
    }) |>
      data.frame(check.names = FALSE)
  }

  df$barcode <- rownames(df)

  # add cell type
  df <- df |>
    dplyr::left_join(y = barcode_info, by = "barcode")

  # wide to long
  df_long <- reshape2::melt(df, id.vars = c("cell", "barcode"))
  colnames(df_long)[3:4] <- c("gene", "exp")

  # order
  if (!is.null(cell.order)) {
    df_long$cell <- factor(df_long$cell, levels = cell.order)
  }

  if (!is.null(gene.order)) {
    df_long$gene <- factor(df_long$gene, levels = gene.order)
  }

  # whether add cluster
  if (is.data.frame(genes)) {
    plyr::ldply(seq_along(markers), function(x) {
      df_tmp <- df_long |>
        dplyr::filter(gene == markers[x]) |>
        dplyr::mutate(cluster = names(markers[x]))

      return(df_tmp)
    }) -> df_long
  }

  # ============================================================================
  # plot
  # ============================================================================
  # strip color
  strip <- do.call(ggh4x::strip_nested, modifyList(
    list(),
    strip_nested_params
  ))

  # facet layer
  if (is.data.frame(genes)) {
    facet_nested <-
      do.call(
        ggh4x::facet_nested, modifyList(
          list(cluster + gene ~ cell,
            scales = "free",
            space = "fixed",
            switch = "y",
            nest_line = ggplot2::element_line(),
            strip = strip
          ),
          facet_nested_params
        )
      )
  } else {
    facet_nested <-
      do.call(
        ggh4x::facet_nested, modifyList(
          list(gene ~ cell,
            scales = "free",
            space = "fixed",
            switch = "y",
            nest_line = ggplot2::element_line(),
            strip = strip
          ),
          facet_nested_params
        )
      )
  }

  # main layer
  pmain <- ggplot2::ggplot(df_long) +
    ggplot2::theme_bw(base_size = 12) +
    facet_nested +
    # do.call(
    #  facet_grid, modifyList(
    #   list(
    #     gene~cell,
    #     scales = "free",
    #     space = "fixed",switch = "y"
    #   ),
    #     facet_grid_params
    #   )
    # ) +
    do.call(
      ggplot2::theme, modifyList(
        list(
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          strip.placement = "outside",
          strip.background.y = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold.italic"),
          strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)
        ),
        theme_params
      )
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("")

  # add layers
  if (plot.type == "heatmap") {
    p <- pmain +
      ggplot2::geom_tile(ggplot2::aes(x = barcode, y = gene, fill = exp)) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::scale_fill_gradient2(
        low = "#313695", mid = "white", high = "#A50026",
        midpoint = 0, na.value = "white"
      )
  } else {
    p <- pmain +
      ggplot2::geom_col(ggplot2::aes(x = barcode, y = exp, fill = gene), width = 1)
  }

  return(p)
}
