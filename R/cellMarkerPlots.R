##' Creates Feature Plots based on average gene expression of a list of genes
##' @param seuratObj Seurat Object
##' @param geneList A list of lists containing gene names
##' @param pt.size Size of each point/cell on the UMAP plot. Default is 0.3
##' @param label Whether to have Ident labels on the UMAP plot. Default is FALSE.
##' @param label.size Size of Indent labels on UMAP plot. Default is 3.
##' @param repel Whether to repel Indent labels on UMAP plot. Default is FALSE.
##' @param combine Whether to stitch/combine all the UMAP plots into 1 ggplot object. Default is TRUE.
##' @author Kimberle Shen
##' @return Feature plot with gene set scores
##' @import Seurat
##' @export
##'

cellMarkerPlots <- function(seuratObj, geneList = gsl, pt.size = 0.3, label = FALSE, label.size = 3, repel = FALSE, combine = TRUE){
  for (name in names(geneList)) {
    log2_norm_umi <- GetAssayData(seuratObj)
    log2_norm_umi_selected <- log2_norm_umi[rownames(log2_norm_umi) %in% geneList[[name]], ]
    averageScore <- colMeans(as.array(log2_norm_umi_selected))
    seuratObj<- AddMetaData(seuratObj, col.name = name, metadata = averageScore)
    print(paste0("Calculating average expression of genes in ", name, " list."))
  }
  plottingInputs <- c(pt.size, label, label.size, repel, combine)
  fPlot <- FeaturePlot(seuratObj, names(geneList), pt.size = plottingInputs[1], label = plottingInputs[2], label.size = plottingInputs[3], repel = plottingInputs[4], combine = plottingInputs[5])
  return(fPlot)
}
