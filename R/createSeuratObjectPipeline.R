##' Creates Seurat object from 10X cell ranger output
##'
##' @param data.dir Directory where data is stored. Default "~/filtered_feature_bc_matrix".
##' @param nFeature_RNA_lower Lower bound of nFeature_RNA. Default 500.
##' @param nFeature_RNA_upper Upper bound of nFeature_RNA. Default 5000.
##' @param percent.mt Upper bound of percentage of mitochondria genes. Default 5.
##' @param nfeatures Number of top most variable genes. Default 2000.
##' @param dims Number of dimensions for finding neighbors and UMAP. Default 20.
##' @param clusterResolution Cluster resolution for finding clusters. Default 0.8.
##' @author Kimberle Shen
##' @return Seurat Object with pca and UMAP calculated, ready for use with predictCells().
##' @import Seurat
##' @export
##'

createSeuratObjectPipeline <- function(data.dir = "~/filtered_feature_bc_matrix", nFeature_RNA_lower = 500, nFeature_RNA_upper = 5000, percent.mt = 5,
                               nfeatures = 2000, dims = 20, clusterResolution = 0.8){
  data <- Read10X(data.dir = data.dir)
  seurat.obj <- CreateSeuratObject(counts = data, project = "mySeuratProj")
  print(seurat.obj)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  p <- c(nFeature_RNA_lower, nFeature_RNA_upper, percent.mt)
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > p[1] & nFeature_RNA < p[2] & percent.mt < p[3])
  print(paste0(ncol(subset(seurat.obj, subset = nFeature_RNA > p[1] & nFeature_RNA < p[2] & percent.mt < p[3])), " cells passed QC and are included"))
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nfeatures)
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:dims) # first 20 PCs based on elbowplot
  seurat.obj <- FindClusters(seurat.obj, resolution = clusterResolution)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:dims)
  matrix <- as.matrix(GetAssayData(seurat.obj, slot = "counts"))
  print(table(Idents(seurat.obj)))
  return (seurat.obj)
}
