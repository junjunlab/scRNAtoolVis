# scRNAtoolVis
 Some useful function to make your scRNA-seq plot more beautiful.
 
 ## install

```R
install.packages('devtools')
devtools::install_github('junjunlab/scRNAtoolVis')

library(scRNAtoolVis)
```

## introduction

> This package mainly is used to replot **seurat** default plot and with other interesting functions. I will add other functions in this package in the future.

## clusterCornerAxes

**clusterCornerAxes** is used to add corner axis on the **left-bottom** UMAP/tSNE principle component plot.

We load test data in **scRNAtoolVis** package:

```R
test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
tmp <- readRDS(test)
```

default plot:

```R
# umap
clusterCornerAxes(object = tmp,reduction = 'umap',
                  noSplit = T)
```

![](https://files.mdnice.com/user/15573/f78954bf-11f6-4a9d-9f45-ce33a58d3e62.png)

We can change arrow type:

```R
# arrowType
clusterCornerAxes(object = tmp,reduction = 'umap',
                  noSplit = T,arrowType = 'open')
```

![](https://files.mdnice.com/user/15573/724a0c5b-742c-468d-be95-2408e4ad5ddf.png)

We can facet by seurat metadata column catogary variable:

```R
# facet by metadata column "orig.ident"
clusterCornerAxes(object = tmp,reduction = 'umap',
                  noSplit = F,groupFacet = 'orig.ident',
                  relLength = 0.5)
```

![](https://files.mdnice.com/user/15573/9289639b-d13d-4036-93ce-0041561dd117.png)

If multiple corner axises will confuse you, you can also set **axes = 'one'** to retain only one axis on the left:

```R
# retain only one axes
clusterCornerAxes(object = tmp,reduction = 'umap',
                  noSplit = F,groupFacet = 'orig.ident',
                  relLength = 0.5,
                  axes = 'one')
```

![](https://files.mdnice.com/user/15573/87b6cdd2-5372-47f1-826b-ca5b752fc7f0.png)

Change the axis and label color:

```R
# line color
clusterCornerAxes(object = tmp,reduction = 'umap',
                  noSplit = F,groupFacet = 'orig.ident',
                  relLength = 0.5,
                  lineTextcol = 'grey50')
```

![](https://files.mdnice.com/user/15573/c2708702-23c0-4590-97c1-cf008ba09afd.png)

Use tSNE reduction data:

```R
# tsne
clusterCornerAxes(object = tmp,reduction = 'tsne',
                  noSplit = F,groupFacet = 'orig.ident',
                  relLength = 0.5)
```

![](https://files.mdnice.com/user/15573/bd1ceafb-ac51-4f8f-b27a-5c544e0952ab.png)

## FeatureCornerAxes

**FeatureCornerAxes** is used to add corner axises on the gene expression reduction map:

```R
# umap
FeatureCornerAxes(object = tmp,reduction = 'umap',
                  groupFacet = 'orig.ident',
                  relLength = 0.5,relDist = 0.2,
                  features = c("Actb","Ythdc1", "Ythdf2"))
```

![](https://files.mdnice.com/user/15573/0a33ee17-8fdc-4ddc-93cb-638c29aa53b7.png)

Keep one axis:

```R
# one axes
FeatureCornerAxes(object = tmp,reduction = 'umap',
                  groupFacet = 'orig.ident',
                  features = c("Actb","Ythdc1", "Ythdf2"),
                  relLength = 0.5,relDist = 0.2,
                  axes = 'one',
                  lineTextcol = 'grey50')
```

![](https://files.mdnice.com/user/15573/c384b247-9f31-4210-acfb-c11b0bade32e.png)

tSNE reduction:

```R
# tsne
FeatureCornerAxes(object = tmp,reduction = 'tsne',
                  groupFacet = 'orig.ident',
                  relLength = 0.5,relDist = 0.2,
                  features = c("Actb","Ythdc1", "Ythdf2"))
```

![](https://files.mdnice.com/user/15573/ec1c34c8-26b2-412b-b7b2-9c222e122865.png)

## AverageHeatmap

**AverageHeatmap** is used to plot averaged expression cross cluster cells.

load data:

```R
httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
pbmc <- readRDS(httest)

# load markergene
markergene <- system.file("extdata", "top5pbmc.markers.csv", package = "scRNAtoolVis")
markers <- read.table(markergene, sep = ',', header = TRUE)
```
plot:

```R
# plot
AverageHeatmap(object = pbmc,
               markerGene = markers$gene)
```

![image](https://user-images.githubusercontent.com/64965509/175778192-0d898fa3-c72e-44e0-8b4f-c47c9d71e3ef.png)

change color:

```R
# change color
AverageHeatmap(object = pbmc,
               markerGene = markers$gene,
               htCol = c("#339933", "#FFCC00", "#FF0033"))
```

![image](https://user-images.githubusercontent.com/64965509/175778256-87ce45b9-45f7-4a25-ba28-8c8335107bc1.png)

## help

More parameters refer to:

```R
?clusterCornerAxes
?FeatureCornerAxes
?AverageHeatmap
```

