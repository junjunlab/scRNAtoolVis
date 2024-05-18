#' @name jjVolcano
#' @author Junjun Lao
#' @title using jjVolcano to visualize marker genes
#'
#' @param diffData diff results with data.frame format, default NULL.
#' @param myMarkers whether supply your own gene labels, default NULL.
#' @param log2FC.cutoff log2FoldChange cutoff, default 0.25.
#' @param pvalue.cutoff pvalue cutoff to filter, default 0.05.
#' @param adjustP.cutoff adjusted pvalue cutoff to be colored in plot, default 0.01.
#' @param topGeneN top genes to be labeled in plot, default 5.
#' @param col.type point color type("updown/adjustP"), default "updown".
#' @param back.col background color, default "grey93".
#' @param pSize point size, default 0.75.
#' @param aesCol point mapping color, default c("#0099CC","#CC3333").
#' @param legend.position legend position in plot, default c(0.7,0.9).
#' @param base_size theme base size, default 14.
#' @param tile.col cluster tile fill color, default jjAnno::useMyCol("paired",n = 9).
#' @param ... other arguments passed by "geom_text_repel".
#' @param cluster.order whether given your cluster orders, default NULL.
#' @param polar whether make the plot to br polar, default FALSE.
#' @param expand the y axis expand, default c(-1,1).
#' @param flip whether flip the plot, default FALSE.
#'
#' @param order.by top marker gene selection method, how the order is, default c("avg_log2FC").
#'
#' @return a ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' jjVolcano(diffData = pbmc.markers)
#' }
globalVariables(c("p_val", "p_val_adj", "type", "type2"))
jjVolcano <- function(
    diffData = NULL,
    myMarkers = NULL,
    order.by = c("avg_log2FC"), # c("avg_log2FC","p_val")
    log2FC.cutoff = 0.25,
    pvalue.cutoff = 0.05,
    adjustP.cutoff = 0.01,
    topGeneN = 5,
    col.type = "updown",
    back.col = "grey93",
    pSize = 0.75,
    aesCol = c("#0099CC", "#CC3333"),
    legend.position = c(0.7, 0.9),
    base_size = 14,
    tile.col = jjAnno::useMyCol("paired", n = 9),
    cluster.order = NULL,
    polar = FALSE,
    expand = c(-1, 1),
    flip = FALSE,
    ...) {
  # filter data
  diff.marker <- diffData %>%
    dplyr::filter(abs(avg_log2FC) >= log2FC.cutoff & p_val < pvalue.cutoff)

  # assign type
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(avg_log2FC >= log2FC.cutoff, "sigUp", "sigDown")) %>%
    dplyr::mutate(type2 = ifelse(p_val_adj < adjustP.cutoff,
      paste("adjust Pvalue < ", adjustP.cutoff, sep = ""),
      paste("adjust Pvalue >= ", adjustP.cutoff, sep = "")
    ))

  # cluster orders
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster,
      levels = cluster.order
    )
  }

  # get background cols
  purrr::map_df(unique(diff.marker$cluster), function(x) {
    tmp <- diff.marker %>%
      dplyr::filter(cluster == x)

    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC) - 0.2,
      max = max(tmp$avg_log2FC) + 0.2
    )
    return(new.tmp)
  }) -> back.data

  # get top gene
  top.marker.tmp <- diff.marker %>%
    dplyr::group_by(cluster)

  # order
  # if(length(order.by) == 1){
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::group_by(cluster) %>%
  #     dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  #
  # }else{
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_head(n = topGeneN)
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_tail(n = topGeneN)
  # }

  top.marker.max <- top.marker.tmp %>%
    dplyr::slice_max(n = topGeneN, order_by = get(order.by))

  top.marker.min <- top.marker.tmp %>%
    dplyr::slice_min(n = topGeneN, order_by = get(order.by))

  # combine
  top.marker <- rbind(top.marker.max, top.marker.min)

  # whether supply own genes
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>%
      dplyr::filter(gene %in% myMarkers)
  } else {
    top.marker <- top.marker
  }

  # ====================================================================
  # plot
  p1 <- ggplot2::ggplot(
    diff.marker,
    ggplot2::aes(x = cluster, y = avg_log2FC)
  ) +
    # add back cols
    ggplot2::geom_col(
      data = back.data,
      ggplot2::aes(x = cluster, y = min), fill = back.col
    ) +
    ggplot2::geom_col(
      data = back.data,
      ggplot2::aes(x = cluster, y = max), fill = back.col
    )

  # ap1 <- paste("adjust Pvalue >= ", adjustP.cutoff, sep = '')
  # ap2 <- paste("adjust Pvalue < ", adjustP.cutoff, sep = '')

  # color type
  if (col.type == "updown") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize) +
      ggplot2::scale_color_manual(values = c("sigDown" = aesCol[1], "sigUp" = aesCol[2]))
  } else if (col.type == "adjustP") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type2), size = pSize) +
      ggplot2::scale_color_manual(values = c(aesCol[2], aesCol[1]))
  }

  # theme details
  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

  # add tile
  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, fill = cluster),
      color = "black",
      height = log2FC.cutoff * 2,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = tile.col) +
    # add gene label
    ggrepel::geom_text_repel(
      data = top.marker,
      ggplot2::aes(x = cluster, y = avg_log2FC, label = gene),
      max.overlaps = 50,
      ...
    )

  # whether coord_plolar
  if (polar == TRUE) {
    p5 <- p4 +
      geomtextpath::geom_textpath(ggplot2::aes(x = cluster, y = 0, label = cluster)) +
      ggplot2::scale_y_continuous(
        n.breaks = 6,
        expand = ggplot2::expansion(mult = expand)
      ) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        legend.position = legend.position,
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::coord_polar(clip = "off", theta = "x")
  } else {
    # whether flip plot
    if (flip == TRUE) {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = cluster, y = 0, label = cluster)) +
        ggplot2::theme(
          axis.line.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_flip()
    } else {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(x = cluster, y = 0, label = cluster)) +
        ggplot2::theme(
          axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }
  return(p5)
}

###############################
#' This is a test data for this package
#' test data description
#'
#' @name pbmc.markers
#' @docType data
#' @author Junjun Lao
"pbmc.markers"
