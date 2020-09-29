##' Creates Feature Plots based on average gene expression of a list of genes
##' @param seuratObj Seurat Object
##' @param geneList A list of lists containing gene names
##' @author Kimberle Shen
##' @return Feature plot with gene set scores
##' @import Seurat
##' @export
##'

cellMarkerPlots <- function(seuratObj, geneList = gsl){
  for (name in names(geneList)) {
    log2_norm_umi <- GetAssayData(seuratObj)
    log2_norm_umi_selected <- log2_norm_umi[rownames(log2_norm_umi) %in% geneList[[name]], ]
    averageScore <- colMeans(as.array(log2_norm_umi_selected))
    seuratObj<- AddMetaData(seuratObj, col.name = name, metadata = averageScore)
    print(paste0("Calculating average expression of genes in ", name, " list."))
  }
  fPlot <- FeaturePlot(seuratObj, names(geneList))
  return(fPlot)
}
