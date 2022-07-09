#' @name markerVocalno
#' @author Junjun Lao
#' @title Marker genes vocalno plot
#' @param markers Dataframe marker genes from findAllmarkers function from seurat.
#' @param topn  Numbers top genes to label, defaults is 5.
#' @param log2FC The threshold of log2FC, defaults is 0.25.
#' @param hlineSize Hline size, defaults is 1.
#' @param hlineColor Hline color, defaults is 'grey50'.
#' @param pforce  Positive gene force parameters to avoid overlap gene labels, defaults is 5.
#' @param nforce Negtive gene force parameters to avoid overlap gene labels, defaults is 2.5.
#' @param nudge_x Ajustments on the horizotal of the gene label, defaults is 0.8.
#' @param pnudge_y Ajustments on the horizotal of the positive gene label, defaults is 0.25.
#' @param nnudge_y Ajustments on the horizotal of the negtive gene label, defaults is 0.
#' @param base_size Theme base size, defaults is 14.
#' @param facetColor Facet border color, defaults is NA.
#' @param facetFill Facet fill color, defaults is 'white'.
#' @param ylab Plot y label, defaults is 'Log2-Fold Change'.
#' @param nrow Numbers rows to plot, defaults is 1.
#'
#' @return Return a ggplot.
#' @export
#'
#' @examples
#' test <- system.file("extdata", "pbmc.markers.csv", package = "scRNAtoolVis")
#' markers <- read.csv(test)
#'
#' markerVocalno(markers = markers,
#'               topn = 5,
#'               labelCol = ggsci::pal_npg()(9))

# define viriables
globalVariables(c("avg_log2FC", "cluster", "gene", "pct.1", "pct.2"))

# define function
markerVocalno <- function(markers = NULL,
                          topn = 5,
                          log2FC = 0.25,
                          labelCol = NULL,
                          hlineSize = 1,
                          hlineColor = 'grey50',
                          pforce = 5,
                          nforce = 2.5,
                          nudge_x = 0.8,
                          pnudge_y = 0.25,
                          nnudge_y = 0,
                          base_size = 14,
                          facetColor = NA,
                          facetFill = 'white',
                          ylab = 'Log2-Fold Change',
                          nrow = 1) {
  # top genes
  toppos <- markers %>% dplyr::group_by(cluster) %>%
    dplyr::top_n(n = topn, wt = avg_log2FC)
  topnegtive <- markers %>% dplyr::group_by(cluster) %>%
    dplyr::top_n(n = -topn, wt = avg_log2FC)

  # merge
  topgene <- rbind(toppos, topnegtive)

  # plot
  ggplot2::ggplot(markers,
                  ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC)) +
    ggplot2::geom_point(color = 'grey80') +
    ggplot2::geom_hline(
      yintercept = c(-log2FC, log2FC),
      lty = 'dashed',
      size = hlineSize,
      color = hlineColor
    ) +
    ggrepel::geom_text_repel(
      data = toppos,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        label = gene,
        color = cluster
      ),
      show.legend = F,
      direction = 'y',
      hjust = 1,
      nudge_y = pnudge_y,
      force = pforce,
      nudge_x = -nudge_x - (toppos$pct.1 - toppos$pct.2)
    ) +
    ggrepel::geom_text_repel(
      data = topnegtive,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        label = gene,
        color = cluster
      ),
      show.legend = F,
      direction = 'y',
      hjust = 0,
      nudge_y = nnudge_y,
      force = nforce,
      nudge_x = nudge_x - (topnegtive$pct.1 - topnegtive$pct.2)
    ) +
    ggplot2::geom_point(
      data = topgene,
      show.legend = F,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        color = cluster
      )
    ) +
    ggplot2::scale_color_manual(name = '', values = labelCol) +
    # x y breaks label
    # scale_y_continuous(limits = c(-6,10),breaks = seq(-6,10,2)) +
    # scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(color = facetColor, fill = facetFill)
    ) +
    ggplot2::xlab(expression(Delta ~ 'Percentage Difference')) +
    ggplot2::ylab(ylab) +
    ggplot2::facet_wrap(~ cluster, nrow = nrow, scales = 'fixed')
}
