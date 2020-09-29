##' Creates gene lists frrom Random Forest Model
##' @param rfModel Random Forest Model
##' @param numberPerList Number of genes per classification to return. Default 10.
##' @author Kimberle Shen
##' @return A list of named gene lists
##' @import randomForest
##' @export
##'

createGeneLists <- function(rfModel = rfmodel, numberPerList = 10){
  geneList <- list()
  vec <- colnames(importance(rfModel))
  l <- length(vec)-2
  vec <- vec[1:l]
  df <- importance(rfModel)
  for (cellType in vec){
    cellTypeList <- df[, which(colnames(df) == cellType)]
    cellTypeList <- names(sort(cellTypeList, decreasing = TRUE))[1:numberPerList]
    geneList[[ cellType ]] <- cellTypeList
    print(paste0("Adding gene list for ", cellType))
  }
  return (geneList)
}
