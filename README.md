rfca
----

This package contains functions that trains a Random Forest Model with a
labelled Seurat Object, for predicting cell types/states in unlabelled
datasets. It also contains a pre-trained Random Forest model, as well as
example datasets.

Introduction
------------

Manual cell annotation of scRNAseq datasets, typically based on marker
genes, can be time-consuming and biased. Being able to automatically
predict cell types/states in a cell-by-cell and cluster-unbiased way is
useful for fast and accurate phenotyping.

In addition, despite the increasing amounts of scRNAseq datasets being
generated, thorough analyis of these datasets is lagging, and/or done in
silos. This package comes with a preloaded Random Forest model based on
different datasets and cell types/states, that will be constantly
updated.

Examples
--------

``` r
# Load the rfca and Seurat libraries
library(rfca)
library(Seurat)

# Create Seurat Object with PCA and UMAP calculated
mySeuratObject <- createSeuratObjectPipeline(data.dir = "~/filtered_feature_bc_matrix", nFeature_RNA_lower = 500, nFeature_RNA_upper = 5000, percent.mt = 5, nfeatures = 2000, dims = 20, clusterResolution = 0.8)
DimPlot(mySeuratObject)
# Assign cell type/state Idents to mySeuratObject manually, if you want to use it as a training dataset.

# Create Random Forest Model with your labelled Seurat Object
myRandomForestModel <- createRFModel(mySeuratObjectLabelled)

# Create Random Forest Model with example labelled Seurat Object
data("exampleSeuratObjectLabelled")
exampleRandomForestModel <- createRFModel(exampleSeuratObjectLabelled)

# Predict cells based on your own created Random Forest Model created above
automaticallyLabelledSeuratObject <- predictCells(mySeuratObject, myRandomForestModel)
DimPlot(automaticallyLabelledSeuratObject)

# Predict cells based on pre-loaded Random Forest Model
automaticallyLabelledSeuratObject <- predictCells(mySeuratObject)
DimPlot(automaticallyLabelledSeuratObject)

# Predict cells based on example unlabelled Seurat Object
data("exampleSeuratObjectUnlabelled")
automaticallyLabelledSeuratObject <- predictCells(exampleSeuratObjectUnlabelled)
DimPlot(automaticallyLabelledSeuratObject)
```

Installation
------------

You can install this package from GitHub with
`devtools::install_github("kimberle9/rfca")`
