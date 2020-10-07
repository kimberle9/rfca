##' Creates Random Forest model from Seurat Object's Idents
##'
##' @param labelled.seurat.obj Seurat Object where Idents are the labels for training. Typically cell type/cell state.
##' @param nfeatures Number of most variable features to be used for training the Random Forest. Default 200.
##' @author Kimberle Shen
##' @return A list of two items, a Random Forest model and features used, to be used in predictCells()
##' @import randomForest
##' @import Seurat
##' @export
##'


createRFModel<- function(labelled.seurat.obj, nfeatures = 200){
  features <- FindVariableFeatures(labelled.seurat.obj, selection.method = "vst", nfeatures = nfeatures)
  topfeatures <- VariableFeatures(features)
  topfeaturesedited <- make.names(topfeatures)
  topfeaturescleaned <- intersect(topfeatures, topfeaturesedited)
  matrix <- as.matrix(GetAssayData(labelled.seurat.obj, slot = "counts"))
  m <- subset(matrix, rownames(matrix) %in% topfeaturescleaned)
  m <- t(m)
  table <- as.data.frame(as.matrix(Idents(labelled.seurat.obj)))
  merged <- merge(m, table, by = "row.names", all = TRUE)
  merged$V1 <- as.factor(merged$V1)
  merged$"Row.names" <- merged$"row.names" <- NULL
  rfmodel <- randomForest(V1 ~ ., data=merged, importance=TRUE, proximity=TRUE, na.action=na.roughfix)
  print (rfmodel)
  return (rfmodel)
}

