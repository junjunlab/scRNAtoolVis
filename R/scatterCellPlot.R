globalVariables(c("idents"))
#' Scatter Cell Plot
#'
#' This function creates a scatter cell plot using the grid package in R.
#'
#' @author Jun Zhang
#'
#' @param object A Seurat object containing the data.
#' @param color A vector of colors for each cell type. If NULL, random colors will be assigned.
#' @param dim The dimension used for plotting (default is "umap").
#' @param rm.axis Logical value indicating whether to remove the x and y-axis (default is FALSE).
#' @param cell.id Name of the column in the metadata that represents cell identity (default is NULL).
#' @param bar.width Width of the barplot (default is 0.2).
#' @param point.size Size of the points in the scatter plot (default is 1).
#'
#' @return None
#'
#' @examples
#' \dontrun{scatterCellPlot(object = seurat_object)}
#'
#' @importFrom grid grid.newpage pushViewport popViewport viewport grid.rect grid.xaxis grid.yaxis grid.points grid.segments grid.text arrow gpar
#'
#' @export
scatterCellPlot <- function(object = NULL,
                            color = NULL,
                            dim = "umap",
                            rm.axis = FALSE,
                            cell.id = NULL,
                            bar.width = 0.2,
                            point.size = 1){
  # ============================================================================
  # 1_extract data
  # ============================================================================
  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = dim))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)
  pc12$idents <- Seurat::Idents(object)

  # summary celltype numbers
  if(is.null(cell.id)){
    cell_num <- pc12 |>
      dplyr::group_by(idents) |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::arrange(n)
  }else{
    cell_num <- pc12 |>
      dplyr::group_by(idents,.data[[cell.id]]) |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::arrange(n)
  }

  # ============================================================================
  # 2_draw plot
  # ============================================================================
  if(rm.axis == FALSE){
    lab.shift = unit(-2.5,"lines")
  }else{
    lab.shift = unit(-1,"lines")
  }

  grid.newpage()
  pushViewport(viewport(x = unit(0.1, "npc"), y = unit(0.5, "npc"),
                        width = unit(0.5, "npc"),
                        height = unit(0.7, "npc"),
                        just = "left",
                        xscale = grDevices::extendrange(range(pc12[,1]),f = 0.05),
                        yscale = grDevices::extendrange(range(pc12[,2]),f = 0.05),
  ))
  grid.rect()
  if(rm.axis == FALSE){
    # grid.xaxis()
    # grid.yaxis()
    jjPlot::grid.xaxis2(label.space = 0.5)
    jjPlot::grid.yaxis2(label.space = 0.25)
  }

  celltype <- cell_num$idents

  if(is.null(color)){
    # create colors
    cols <- circlize::rand_color(n = length(celltype))
  }else{
    cols <- color
  }

  # draw points
  # i = 1
  for (i in seq_along(celltype)) {
    tmp <- pc12 |> dplyr::filter(idents == celltype[i])

    grid.points(x = tmp[,1],y = tmp[,2],pch = 19,size = unit(point.size,"pt"),
                gp = gpar(col = cols[i]))
  }

  # arrow
  if(rm.axis == TRUE){
    grid.segments(x0 = 0.025,x1 = 0.2,y0 = 0.05,y1 = 0.05,
                  arrow = arrow(length = unit(2,"mm"),type = "closed"),
                  gp = gpar(fill = "black"))
    grid.text(label = paste0(toupper(dim)," 1"),x = (0.2+0.025)/2,y = 0.025,
              gp = gpar(fontsize = 6,fontface = "bold.italic"))
    grid.segments(x0 = 0.05,x1 = 0.05,y0 = 0.025,y1 = 0.2,
                  arrow = arrow(length = unit(2,"mm"),type = "closed"),
                  gp = gpar(fill = "black"))
    grid.text(label = paste0(toupper(dim)," 2"),x = 0.025,y = (0.2+0.025)/2,rot = 90,gp =
                gpar(fontsize = 6,fontface = "bold.italic"))
  }

  # labs
  grid.text(label = paste0(toupper(dim)," dimension 1"),x = 0.5,y = lab.shift)
  grid.text(label = paste0(toupper(dim)," dimension 2"),x = lab.shift,y = 0.5,rot = 90)

  popViewport()

  # ============================================================================
  # barplot
  pushViewport(viewport(x = unit(0.61, "npc"), y = unit(0.5, "npc"),
                        width = unit(bar.width, "npc"),
                        height = unit(0.7, "npc"),
                        just = "left",
                        yscale = c(0,nrow(cell_num) + 0.75),
                        xscale = c(0,max(cell_num$n) + 0.1*max(cell_num$n))))

  if(rm.axis == FALSE){
    # grid.xaxis()
    jjPlot::grid.xaxis2(label.space = 0.5,
                        at = c(0,max(cell_num$n)),
                        labels = as.character(c(0,max(cell_num$n))))
  }
  grid.rect(x = rep(0,nrow(cell_num)),y = 1:nrow(cell_num),
            width = cell_num$n,height = unit(0.08,"npc"),
            just = "left",
            gp = gpar(fill = cols,col = NA),
            default.units = "native")
  grid.rect(gp = gpar(fill = "transparent"))
  grid.text(label = "Number of cells",x = 0.5,y = lab.shift)
  popViewport()

  # ============================================================================
  # legend
  pushViewport(viewport(x = unit(0.61 + bar.width, "npc"), y = unit(0.5, "npc"),
                        width = unit(0.2, "npc"),
                        height = unit(0.7, "npc"),
                        just = "left",
                        yscale = c(0,nrow(cell_num) + 0.75)))

  grid.points(x = rep(0.1,nrow(cell_num)),y = 1:nrow(cell_num),pch = 19,
              gp = gpar(col = cols),size = unit(1.5, "char"))
  if(!is.null(cell.id)){
    grid.text(label = as.character(unlist(cell_num[,cell.id])),x = 0.1,y = 1:nrow(cell_num),
              default.units = "native")
  }
  grid.text(label = cell_num$idents,x = 0.2,y = 1:nrow(cell_num),
            just = "left",
            gp = gpar(fontsize = 10),
            default.units = "native")
  popViewport()
}
