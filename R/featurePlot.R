#' featurePlot Function
#'
#' This function creates a scatter plot for multiple genes or features from a
#' Seurat object.
#'
#' @author Jun Zhang
#'
#' @param object A Seurat object containing the data.
#' @param dim The dimension to use for plotting, default is "umap".
#' @param genes A character vector of gene names or feature names to be plotted.
#' @param nrow Number of rows in the plot grid.
#' @param ncol Number of columns in the plot grid.
#' @param quantile.val The quantile value to determine the color cutoff for each gene.
#' @param color A vector of colors to be used for plotting, defaults to a
#' predefined set of colors.
#' @param rm.axis Logical value indicating whether to remove axis labels and ticks,
#' defaults to FALSE.
#' @param rm.legend Logical value indicating whether to remove the color legend,
#' defaults to FALSE.
#' @param add.rect Logical value indicating whether to add a rectangle around
#' each plot panel, defaults to FALSE.
#' @param add.corArrow Logical value indicating whether to add arrows indicating
#' the correlation direction, defaults to FALSE.
#' @param add.strip Logical value indicating whether to add a strip at the top of
#' each plot panel, defaults to FALSE.
#' @param corLabel.dist Distance between the corner arrows and the axis labels.
#' @param arrow.len Length of the corner arrows.
#' @param arrow.label.size Font size of the corner arrow labels.
#' @param plot.size Size of each individual scatter plot.
#' @param keep.oneCor Logical value indicating whether to keep only one set of
#' corner arrows across the entire plot grid, defaults to FALSE.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param respect Logical value indicating whether to respect the specified number
#' of rows and columns in the plot grid, defaults to TRUE.
#' @param point.size the point size, default 1.
#'
#' @examples
#' \dontrun{
#' # Assuming "seurat_obj" is a Seurat object
#' featurePlot(object = seurat_obj, genes = c("gene1", "gene2", "gene3"), nrow = 2, ncol = 2)
#' }
#'
#' @import dplyr
#' @import grDevices
#' @import ggplot2
#' @import Seurat
#' @importFrom grid unit viewport pushViewport popViewport grid.rect grid.text
#' grid.points arrow gpar grid.layout
#' @importFrom stats quantile
#'
#' @export

globalVariables(c("col_rg", "tmp_col"))

featurePlot <- function(
    object = NULL,
    dim = "umap",
    genes = NULL,
    nrow = NULL,
    ncol = NULL,
    quantile.val = 1,
    color = NULL,
    rm.axis = FALSE,
    rm.legend = FALSE,
    add.rect = FALSE,
    add.corArrow = FALSE,
    add.strip = FALSE,
    corLabel.dist = 0.08,
    arrow.len = 0.2,
    arrow.label.size = 6,
    plot.size = 0.6,
    keep.oneCor = FALSE,
    xlab = NULL,
    ylab = NULL,
    respect = TRUE,
    point.size = 1) {
  # ============================================================================
  # 1_extract data
  # ============================================================================
  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = dim))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)
  pc12$idents <- Idents(object)

  # get gene expression
  geneExp <- Seurat::FetchData(object = object, vars = genes)

  # cbind
  mer <- cbind(pc12, geneExp)

  # get nrow and ncol
  if (is.null(nrow)) {
    nrow <- ifelse(is.null(ncol), 1, ceiling(length(genes) / ncol))
  }

  if (is.null(ncol)) {
    ncol <- ifelse(is.null(nrow), length(genes), ceiling(length(genes) / nrow))
  }

  gene_mtx <- suppressWarnings(matrix(genes, nrow = nrow, ncol = ncol))

  # assign colors
  if (is.null(color)) {
    cols <- c("grey90", "#57C5B6", "#159895", "#1A5F7A", "#002B5B")
  } else {
    cols <- color
  }

  # ============================================================================
  # 2_draw plot
  # ============================================================================
  if (rm.axis == FALSE) {
    lab.shift <- unit(-2.5, "lines")
  } else {
    lab.shift <- unit(-1, "lines")
  }

  # CANVAS FOR PLOT
  grid.newpage()
  pushViewport(
    viewport(
      x = 0.5, y = 0.5,
      width = 0.9, height = 0.9,
      xscale = range(mer[, 1]), yscale = range(mer[, 2]),
      layout = grid.layout(nrow = nrow, ncol = ncol, respect = respect)
    )
  )

  # loop
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      # check genes numbers
      if (i * j > length(genes)) {
        break
      }

      # ===========================================================
      # 1_panel grid
      pushViewport(
        viewport(layout.pos.row = i, layout.pos.col = j)
      )

      if (add.rect == TRUE) {
        grid.rect()
      }

      # process data
      quantile_val <- quantile(mer[, gene_mtx[i, j]], probs = quantile.val)
      mer <- mer |>
        dplyr::mutate(tmp_col = if_else(.data[[gene_mtx[i, j]]] > quantile_val,
          quantile_val,
          .data[[gene_mtx[i, j]]]
        ))

      tmp_data <- mer |>
        dplyr::arrange(tmp_col)
      col_p <- colorRampPalette(cols)(100)
      cut_range <- cut(tmp_data[, "tmp_col"], 100)

      labs <- levels(cut_range)
      names(labs) <- col_p

      tmp_data <- tmp_data |>
        dplyr::mutate(col_rg = as.character(cut_range)) |>
        dplyr::mutate(col_f = ifelse(col_rg %in% labs, names(labs)[match(col_rg, labs)], "black"))

      # ===========================================================
      # 2_scatter plot
      pushViewport(
        viewport(
          x = 0.5, y = 0.5, width = plot.size, height = plot.size,
          xscale = extendrange(range(tmp_data[, 1]), f = 0.05),
          yscale = extendrange(range(tmp_data[, 2]), f = 0.05)
        )
      )

      # whether add corner arrows
      if (keep.oneCor == TRUE) {
        if (j == 1) {
          if (add.corArrow == TRUE) {
            grid.segments(
              x0 = 0, x1 = arrow.len, y0 = 0, y1 = 0,
              arrow = arrow(length = unit(2, "mm"), type = "closed"),
              gp = gpar(fill = "black")
            )
            grid.text(
              label = paste0(toupper(dim), " 1"), x = arrow.len / 2, y = -corLabel.dist,
              gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
            )
            grid.segments(
              x0 = 0, x1 = 0, y0 = 0, y1 = arrow.len,
              arrow = arrow(length = unit(2, "mm"), type = "closed"),
              gp = gpar(fill = "black")
            )
            grid.text(
              label = paste0(toupper(dim), " 2"),
              x = -corLabel.dist, y = arrow.len / 2, rot = 90,
              gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
            )
          } else {
            grid.rect()
          }
        }
      } else {
        if (add.corArrow == TRUE) {
          grid.segments(
            x0 = 0, x1 = arrow.len, y0 = 0, y1 = 0,
            arrow = arrow(length = unit(2, "mm"), type = "closed"),
            gp = gpar(fill = "black")
          )
          grid.text(
            label = paste0(toupper(dim), " 1"), x = arrow.len / 2, y = -corLabel.dist,
            gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
          )
          grid.segments(
            x0 = 0, x1 = 0, y0 = 0, y1 = arrow.len,
            arrow = arrow(length = unit(2, "mm"), type = "closed"),
            gp = gpar(fill = "black")
          )
          grid.text(
            label = paste0(toupper(dim), " 2"),
            x = -corLabel.dist, y = arrow.len / 2, rot = 90,
            gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
          )
        } else {
          grid.rect()
        }
      }

      grid.points(
        x = tmp_data[, 1], y = tmp_data[, 2], pch = 19, size = unit(point.size, "pt"),
        gp = gpar(col = tmp_data$col_f)
      )

      # whether draw axis
      if (add.corArrow == FALSE) {
        if (rm.axis == FALSE) {
          # grid.xaxis()
          # grid.yaxis()
          jjPlot::grid.xaxis2(label.space = 0.5)
          jjPlot::grid.yaxis2(label.space = 0.25)
        }
      }

      # add strip
      if (add.strip == TRUE) {
        grid.rect(
          x = 0.5, y = 1, width = 1,
          height = 0.15, gp = gpar(fill = "grey85"),
          just = "bottom"
        )
      }

      grid.text(
        label = gene_mtx[i, j], x = 0.5, y = unit(1 + 0.15 / 2, "npc"),
        gp = gpar(fontface = "bold.italic")
      )
      if (add.corArrow == FALSE) {
        # axis labels
        if (!is.null(xlab) || !is.null(ylab)) {
          axis.label.x <- xlab
          axis.label.y <- ylab
        } else {
          axis.label.x <- paste0(toupper(dim), " dimension 1")
          axis.label.y <- paste0(toupper(dim), " dimension 2")
        }

        grid.text(label = axis.label.x, x = 0.5, y = lab.shift)
        grid.text(label = axis.label.y, x = lab.shift, y = 0.5, rot = 90)
      }

      popViewport()

      # ===========================================================
      # 3_draw legend
      if (rm.legend == FALSE) {
        pushViewport(
          viewport(
            x = 0.5 + plot.size / 2 + 0.01, y = 0.5,
            width = 0.025, height = unit(plot.size, "npc"),
            just = "left",
            yscale = range(tmp_data[, gene_mtx[i, j]])
          )
        )
        # grid.rect(x = 0.5, y = unit(seq(0.25,0.75, length = 100), "npc"),
        #           width = unit(1, "npc"), height = unit(0.5, "npc"),
        #           just = "centre",default.units = "npc",
        #           gp = gpar(col = NA, fill = col_p))
        # grid.rect(gp = gpar(fill = NA))
        # # grid.yaxis(main = FALSE)
        # jjPlot::grid.yaxis2(side = "right",tick.len = 0.25)

        jjPlot::grid.colorkey(
          x = tmp_data[, gene_mtx[i, j]],
          color = cols,
          pos = "v",
          ticks.side = "right"
        )

        popViewport()
      }
      popViewport()
    }
  }
}
