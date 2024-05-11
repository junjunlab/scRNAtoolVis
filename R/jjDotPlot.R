#' @name jjDotPlot
#' @author Junjun Lao
#' @title using dotplot to visualize gene expression
#'
#' @param object seurat object, default NULL.
#' @param assay the assay to be selected, default NULL.
#' @param slot gene expression data to be selected, default "data".
#' @param id the cell clusters id in the metadata info, default "seurat_clusters".
#' @param split.by the group name to split, default NULL.
#' @param split.by.aesGroup whether the dot color filled by group, default FALSE.
#' @param gene the genes to be drawed in plot, default NULL.
#' @param markerGene the marker genes with celltype info to be drawed, default NULL.
#' @param point.geom the ggplot "point" geom layer to be shown, default TRUE.
#' @param point.shape the point shape,default 21.
#' @param tile.geom the ggplot "tile" geom layer to be shown, default FALSE.
#' @param dot.col dot colors, default "c("white","#990000")".
#' @param midpoint the midpoint value when using rescale parameters, default 0.5.
#' @param scale whether scale the gene expressions, default TRUE.
#' @param col.min minimum scaled average expression threshold (everything smaller will be set to this), default -2.5.
#' @param col.max maximum scaled average expression threshold (everything larger will be set to this), default 2.5.
#' @param rescale whether rescale the average expression to specific range, default FALSE.
#' @param rescale.min minimum rescaled average expression threshold, default 0.
#' @param rescale.max maximum rescaled average expression threshold, default 1.
#' @param dot.min the minimum size of dot point, default 1.
#' @param dot.max the maximum size of dot point, default 6.
#' @param base_size the theme base size, default 14.
#' @param x.text.angle the x text angle, default 90.
#' @param x.text.hjust the x text hjust, default 1
#' @param x.text.vjust the x text vjust, default 0.
#' @param ytree whether add the y axis dendrogram, default TRUE.
#' @param xtree whether add the x axis dendrogram, default FALSE.
#' @param tree.pos the dendrogram position to be placed, default NULL.
#' @param same.pos.label whether place the x/y axis label on same side, default FALSE.
#' @param size.width the point legend width, default 0.3.
#' @param bar.width the colorbar legend width, default 4.5.
#' @param plot.margin the plot margin, default c(1,1,1,1) which relative to c(t,r,b,l).
#' @param anno whether anno celltype, default FALSE.
#' @param aesGroName the markerGene data.frame column name which refers to celltype, default "cluster".
#' @param segWidth annoSegment width, default 0.8.
#' @param lwd annoSegment line size, default 3.
#' @param textRot annoSegment text angle, default 90.
#' @param textSize annoSegment text size, default 14.
#' @param hjust annoSegment text hjust, default 0.
#' @param legend.position ggplot legend position, default "right".
#' @param bar.legendTitle colorbar legend title, default "Mean expression in group".
#' @param point.lengdTitle point size legend title, default "Fraction of cells in group (%)".
#' @param ... other parameters passed to annoSegment function.
#'
#' @param gene.order supply your own gene orders, default NULL.
#' @param cluster.order supply your own cluster number orders, default NULL.
#'
#' @import tidyverse
#' @import ggdendro
#' @import ggsci
#' @import patchwork
#' @import utils
#'
#' @return Return a ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
#' pbmc <- readRDS(httest)
#'
#' # add groups
#' pbmc$groups <- rep(c("stim", "control"), each = 1319)
#' # add celltype
#' pbmc$celltype <- Seurat::Idents(pbmc)
#'
#' # load markergene
#' markergene <- data("top3pbmc.markers")
#'
#' # ====================================
#' jjDotPlot(
#'   object = pbmc,
#'   gene = top3pbmc.markers$gene
#' )
#'
#' jjDotPlot(
#'   object = pbmc,
#'   gene = top3pbmc.markers$gene,
#'   id = "celltype"
#' )
#'
#' jjDotPlot(
#'   object = pbmc,
#'   markerGene = top3pbmc.markers
#' )
#'
#' jjDotPlot(
#'   object = pbmc,
#'   markerGene = top3pbmc.markers,
#'   xtree = TRUE
#' )
#'
#' jjDotPlot(
#'   object = pbmc,
#'   markerGene = top3pbmc.markers,
#'   anno = TRUE,
#'   plot.margin = c(3, 1, 1, 1)
#' )
#' }
globalVariables(c("%||%", ".", "avg.exp", "avg.exp.scaled", "celltype", "group", "pct.exp", "unit"))
PercentAbove <- utils::getFromNamespace("PercentAbove", "Seurat")

# define function
jjDotPlot <- function(
    object = NULL,
    assay = NULL,
    slot = "data",
    id = "seurat_clusters",
    split.by = NULL,
    split.by.aesGroup = FALSE,
    gene = NULL,
    markerGene = NULL,
    gene.order = NULL,
    cluster.order = NULL,
    point.geom = TRUE,
    point.shape = 21,
    tile.geom = FALSE,
    dot.col = c("white", "#990000"),
    midpoint = 0.5,
    scale = TRUE,
    col.min = -2.5,
    col.max = 2.5,
    rescale = FALSE,
    rescale.min = 0,
    rescale.max = 1,
    dot.min = 1,
    dot.max = 6,
    base_size = 14,
    x.text.angle = 90,
    x.text.hjust = 1,
    x.text.vjust = 0,
    ytree = TRUE,
    xtree = FALSE,
    tree.pos = NULL,
    same.pos.label = FALSE,
    size.width = 0.3,
    bar.width = 4.5,
    plot.margin = c(1, 1, 1, 1),
    anno = FALSE,
    aesGroName = "cluster",
    segWidth = 0.8,
    lwd = 3,
    textRot = 90,
    textSize = 14,
    hjust = 0,
    legend.position = "right",
    bar.legendTitle = "Mean expression \n in group",
    point.lengdTitle = "Fraction of cells \n in group (%)",
    ...) {
  # set assays
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  Seurat::DefaultAssay(object = object) <- assay

  # get gene expression
  if (is.null(gene) && !is.null(markerGene)) {
    geneExp <- Seurat::FetchData(
      object = object,
      vars = unique(markerGene$gene),
      slot = slot
    )
  } else if (!is.null(gene) && is.null(markerGene)) {
    geneExp <- Seurat::FetchData(
      object = object,
      vars = gene,
      slot = slot
    )
  } else {
    print("Please supply  one option!")
  }

  # get cluster number or celtype
  # whether split by groups
  if (is.null(split.by)) {
    if (id %in% colnames(object@meta.data)) {
      # using seurat_clusters or Idents
      geneExp$id <- object@meta.data[[id]]
    } else {
      geneExp$id <- Seurat::Idents(object)
    }
  } else {
    # using seurat_clusters or Idents
    if (id %in% colnames(object@meta.data)) {
      geneExp$id <- paste(object@meta.data[[id]],
        " (", object@meta.data[[split.by]], ")",
        sep = ""
      )
    } else {
      geneExp$id <- paste(Seurat::Idents(object),
        " (", object@meta.data[[split.by]], ")",
        sep = ""
      )
    }
  }

  # calculate mean expressions and percent expression
  data.plot <- lapply(
    X = unique(x = geneExp$id),
    FUN = function(ident) {
      data.use <- geneExp[geneExp$id == ident, 1:(ncol(x = geneExp) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)

      res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
      res$gene <- rownames(res)
      return(res)
    }
  ) %>%
    do.call("rbind", .) %>%
    data.frame()

  # scale mean expressions
  if (scale == TRUE) {
    purrr::map_df(unique(data.plot$gene), function(x) {
      tmp <- data.plot %>%
        dplyr::filter(gene == x)

      # scale value corss groups
      if (rescale == TRUE) {
        avg.exp.scale <- tmp %>%
          dplyr::select(avg.exp) %>%
          scale(.) %>%
          scales::rescale(., to = c(rescale.min, rescale.max))
      } else {
        avg.exp.scale <- tmp %>%
          dplyr::select(avg.exp) %>%
          scale(.) %>%
          Seurat::MinMax(., min = col.min, max = col.max)
      }

      # add scaled exp
      tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
      return(tmp)
    }) -> data.plot
  } else {
    # no scale or rescale using log1p
    data.plot$avg.exp.scaled <- log1p(data.plot$avg.exp)
  }

  if (!is.null(markerGene)) {
    celltype_info <- markerGene

    # add celltype
    purrr::map_df(seq_len(nrow(data.plot)), function(x) {
      tmp <- data.plot[x, ]
      tmp$celltype <- celltype_info[which(celltype_info$gene == tmp$gene), aesGroName][[1]]
      return(tmp)
    }) -> data.plot.res

    # arrange cell type
    data.plot.res <- data.plot.res %>%
      dplyr::arrange(celltype)
  } else {
    data.plot.res <- data.plot
  }

  # gene order
  if (is.null(gene.order)) {
    data.plot.res$gene <- factor(data.plot.res$gene, levels = unique(data.plot.res$gene))
  } else {
    data.plot.res$gene <- factor(data.plot.res$gene, levels = unique(gene.order))
  }

  # cluster order
  if (!is.null(cluster.order)) {
    data.plot.res$id <- factor(data.plot.res$id, levels = unique(cluster.order))
  }

  # add group info
  if (!is.null(split.by)) {
    data.plot.res$group <- sapply(strsplit(as.character(data.plot.res$id), split = "\\(|\\)"), "[", 2)
  }

  # ============================================================================
  # plot
  pmain <- ggplot2::ggplot(
    data.plot.res,
    ggplot2::aes(x = gene, y = id)
  ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        t = plot.margin[1], r = plot.margin[2],
        b = plot.margin[3], l = plot.margin[4],
        unit = "cm"
      ),
      axis.text = ggplot2::element_text(color = "black"),
      legend.direction = "horizontal",
      axis.text.x = ggplot2::element_text(
        angle = x.text.angle,
        hjust = x.text.hjust,
        vjust = x.text.vjust
      ),
      legend.position = legend.position
    )

  # colorbar layer
  colorbar.layer <- ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title = bar.legendTitle,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(bar.width, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.8,
      ticks.colour = "black"
    )
  )
  # point size layer
  point.layer <- ggplot2::guides(
    size = ggplot2::guide_legend(
      title = point.lengdTitle,
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      override.aes = list(
        color = "black",
        fill = "grey50"
      ),
      keywidth = ggplot2::unit(size.width, "cm")
    )
  )

  # chenge colors
  if (is.null(split.by)) {
    if (length(dot.col) == 2) {
      pmain <- pmain +
        ggplot2::scale_fill_gradient(low = dot.col[1], high = dot.col[2])
    } else {
      pmain <- pmain +
        ggplot2::scale_fill_gradient2(
          low = dot.col[1],
          mid = dot.col[2],
          high = dot.col[3],
          midpoint = midpoint
        )
    }
  } else {
    if (split.by.aesGroup == FALSE) {
      pmain <- pmain +
        ggplot2::scale_fill_manual(values = dot.col)
    } else {
      if (length(dot.col) == 2) {
        pmain <- pmain +
          ggplot2::scale_fill_gradient(low = dot.col[1], high = dot.col[2])
      } else {
        pmain <- pmain +
          ggplot2::scale_fill_gradient2(
            low = dot.col[1],
            mid = dot.col[2],
            high = dot.col[3],
            midpoint = midpoint
          )
      }
    }
  }

  # add point or tile layer
  if (point.geom == TRUE) {
    if (is.null(split.by)) {
      pmain <- pmain +
        ggplot2::geom_point(
          ggplot2::aes(fill = avg.exp.scaled, size = pct.exp),
          color = "black", shape = point.shape
        ) +
        colorbar.layer +
        point.layer +
        ggplot2::scale_size(range = c(dot.min, dot.max))
    } else {
      if (split.by.aesGroup == FALSE) {
        pmain <- pmain +
          ggplot2::geom_point(
            ggplot2::aes(fill = group, size = pct.exp),
            color = "black", shape = point.shape
          ) +
          point.layer +
          ggplot2::scale_size(range = c(dot.min, dot.max))
      } else {
        pmain <- pmain +
          ggplot2::geom_point(
            ggplot2::aes(fill = avg.exp.scaled, size = pct.exp),
            color = "black", shape = point.shape
          ) +
          colorbar.layer +
          point.layer +
          ggplot2::scale_size(range = c(dot.min, dot.max))
      }
    }
  } else if (tile.geom == TRUE) {
    if (is.null(split.by)) {
      pmain <- pmain +
        ggplot2::geom_tile(
          ggplot2::aes(fill = avg.exp.scaled),
          color = "black"
        ) +
        colorbar.layer
    } else {
      if (split.by.aesGroup == FALSE) {
        pmain <- pmain +
          ggplot2::geom_tile(
            ggplot2::aes(fill = group),
            color = "black"
          ) +
          colorbar.layer
      } else {
        pmain <- pmain +
          ggplot2::geom_tile(
            ggplot2::aes(fill = avg.exp.scaled),
            color = "black"
          ) +
          colorbar.layer
      }
    }
  }

  # long to wide matrix
  plot.matrix <- data.plot.res %>%
    reshape2::dcast(id ~ gene, value.var = "avg.exp.scaled")

  rownames(plot.matrix) <- plot.matrix$id
  plot.matrix <- plot.matrix %>% dplyr::select(-id)

  # =======================================
  # define label position
  if (is.null(tree.pos)) {
    if (ytree == TRUE && xtree == FALSE) {
      tree.pos <- c("right", "right")
    } else if (xtree == TRUE && ytree == FALSE) {
      tree.pos <- c("top", "top")
    } else {
      tree.pos <- c("right", "top")
    }
  } else {
    tree.pos <- tree.pos
  }

  # add ytree
  if (ytree == TRUE) {
    # whether put label same with tree
    if (same.pos.label == FALSE) {
      pytree <- pmain +
        ggh4x::scale_y_dendrogram(
          hclust = stats::hclust(stats::dist(plot.matrix)),
          position = tree.pos[1],
          labels = NULL
        ) +
        ggplot2::guides(y.sec = ggh4x::guide_axis_manual(labels = function(x) {
          x
        }))
    } else {
      pytree <- pmain +
        ggh4x::scale_y_dendrogram(
          hclust = stats::hclust(stats::dist(plot.matrix)),
          position = tree.pos[1]
        )
    }
  } else {
    pytree <- pmain
  }

  # add xtree
  if (xtree == TRUE) {
    # whether put label same with tree
    if (same.pos.label == FALSE) {
      pxtree <- pytree +
        ggh4x::scale_x_dendrogram(
          hclust = stats::hclust(stats::dist(t(plot.matrix))),
          position = tree.pos[2],
          labels = NULL
        ) +
        ggplot2::guides(x.sec = ggh4x::guide_axis_manual(labels = function(x) {
          x
        }))
    } else {
      pxtree <- pytree +
        ggh4x::scale_x_dendrogram(
          hclust = stats::hclust(stats::dist(t(plot.matrix))),
          position = tree.pos[2]
        )
    }
  } else {
    pxtree <- pytree
  }

  # ================================
  # add celltype annotation
  if (!is.null(markerGene) && anno == TRUE) {
    panno <- jjAnno::annoSegment(
      object = pxtree,
      annoPos = "top",
      aesGroup = TRUE,
      aesGroName = "celltype",
      segWidth = 0.8,
      lwd = 3, addBranch = TRUE, branDirection = -1,
      pCol = rep("black", length(unique(celltype_info$cluster))),
      addText = TRUE, textRot = 90,
      textCol = rep("black", length(unique(celltype_info$cluster))),
      textSize = 14,
      hjust = 0,
      ...
    )
  } else {
    panno <- pxtree
  }
  return(panno)
}


###############################
#' This is a test data for this package
#' test data describtion
#'
#' @name top3pbmc.markers
#' @docType data
#' @author Junjun Lao
"top3pbmc.markers"
