# library(scRNAtoolVis)
# library(Seurat)
#
# test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
# tmp <- readRDS(test)
#
# # add celltype
# tmp$cellType <- Idents(tmp)
#
# ################################clusterCornerAxes
#
# # legend key size
# # umap
# clusterCornerAxes(object = tmp,
#                   reduction = 'umap',
#                   clusterCol = 'cellType',
#                   noSplit = T,
#                   keySize = 8)
#
# # add cell type
# clusterCornerAxes(object = tmp,
#                   reduction = 'umap',
#                   clusterCol = "cellType",
#                   noSplit = T,
#                   cellLabel = T,
#                   cellLabelSize = 5)
#
# # remove legend
# clusterCornerAxes(object = tmp,
#                   reduction = 'umap',
#                   clusterCol = "cellType",
#                   noSplit = T,
#                   cellLabel = T,
#                   cellLabelSize = 5,
#                   show.legend = F)
#
# # split
# clusterCornerAxes(object = tmp,
#                   reduction = 'umap',
#                   clusterCol = "cellType",
#                   groupFacet = 'orig.ident',
#                   noSplit = F,
#                   cellLabel = T,
#                   cellLabelSize = 3,
#                   show.legend = F,
#                   aspect.ratio = 1,
#                   themebg = 'bwCorner')
#
# ################################FeatureCornerAxes
#
# # default
# FeatureCornerAxes(object = tmp,reduction = 'umap',
#                   groupFacet = 'orig.ident',
#                   relLength = 0.5,
#                   relDist = 0.2,
#                   features = c("Actb","Ythdc1", "Ythdf2"))
#
# # remove legend
# FeatureCornerAxes(object = tmp,reduction = 'umap',
#                   groupFacet = 'orig.ident',
#                   relLength = 0.5,
#                   relDist = 0.2,
#                   features = c("Actb","Ythdc1", "Ythdf2"),
#                   show.legend = F)
#
# # no facet group
# FeatureCornerAxes(object = tmp,reduction = 'umap',
#                   groupFacet = NULL,
#                   relLength = 0.5,
#                   relDist = 0.2,
#                   features = c("Actb","Ythdc1", "Ythdf2"),
#                   aspect.ratio = 1)
#
# # specify corner position
# p1 <- FeatureCornerAxes(object = tmp,reduction = 'umap',
#                         groupFacet = 'orig.ident',
#                         relLength = 0.5,
#                         relDist = 0.2,
#                         aspect.ratio = 1,
#                         features = c("Actb","Ythdc1", "Ythdf2"),
#                         axes = 'one')
#
# p2 <- FeatureCornerAxes(object = tmp,reduction = 'umap',
#                         groupFacet = 'orig.ident',
#                         relLength = 0.5,
#                         relDist = 0.2,
#                         aspect.ratio = 1,
#                         features = c("Actb","Ythdc1", "Ythdf2"),
#                         axes = 'one',
#                         cornerVariable = 'ST4')
#
# # combine
# cowplot::plot_grid(p1,p2,ncol = 2,align = 'hv')
#
# # given a range to plot
# p1 <- FeatureCornerAxes(object = tmp,reduction = 'umap',
#                         groupFacet = NULL,
#                         relLength = 0.5,
#                         relDist = 0.2,
#                         features = c("Actb","Ythdc1", "Ythdf2"),
#                         aspect.ratio = 1,
#                         themebg = 'bwCorner')
#
# p2 <- FeatureCornerAxes(object = tmp,reduction = 'umap',
#                         groupFacet = NULL,
#                         relLength = 0.5,
#                         relDist = 0.2,
#                         features = c("Actb","Ythdc1", "Ythdf2"),
#                         aspect.ratio = 1,
#                         themebg = 'bwCorner',
#                         minExp = 0,maxExp = 2)
#
# # combine
# cowplot::plot_grid(p1,p2,ncol = 1,align = 'hv')
#
# ################################AverageHeatmap
# httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
# pbmc <- readRDS(httest)
#
# # load markergene
# markergene <- system.file("extdata", "top5pbmc.markers.csv", package = "scRNAtoolVis")
# markers <- read.table(markergene, sep = ',', header = TRUE)
#
# # remove rownames
# AverageHeatmap(object = pbmc,
#                markerGene = markers$gene,
#                showRowNames = F)
#
# # remove cluster anno name
# AverageHeatmap(object = pbmc,
#                markerGene = markers$gene,
#                clusterAnnoName = F)
#
# # mark some genes
# # tartget gene
# annoGene <- c("LDHB","CCR7","LEF1","NKG7","CST7",
#               "GZMK","HLA-DQA1","HLA-DRB1","HLA-DPA1")
#
# AverageHeatmap(object = pbmc,
#                markerGene = markers$gene,
#                clusterAnnoName = F,
#                markGenes = annoGene)
