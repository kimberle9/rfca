##' Takes a Seurat Object and returns the Seurat Object with predicted CNS cell types as Idents
##'
##' @param seuratObj Seurat Object.
##' @param rfModel Random Forest Object from randomForest package. Default: A Random Forest Object trained with CNS cell types.
##' @author Kimberle Shen
##' @return Seurat Object with predicted cell types as Idents and stored under "predictedCellType" metadata
##' @import Seurat
##' @import randomForest
##' @export
##'

predictCells <- function(seuratObj, rfModel = rfmodel){
  matrix <- as.matrix(GetAssayData(seuratObj, slot = "counts"))
  genelist <- row.names(importance(rfModel))
  m <- subset(matrix, rownames(matrix) %in% genelist)
  m <- t(m)
  pred <- randomForest:::predict.randomForest(rfModel, m)
  seuratObj <- AddMetaData(seuratObj, col.name = "predictedCellType", metadata = pred)
  Idents(seuratObj) <- "predictedCellType"
  print(table(Idents(seuratObj)))
  return (seuratObj)
}
