# command line for the clustering,
# next step after integration_reduction_dimension.R commands
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day2/clustering/
# Tue Jul  5 14:39:21 2022

# load the libraries
library(Seurat)
library(clustree)

# load the seurat object previously created by previous R script (see above in comments)
seu_int <- readRDS("seu_int_day2_part1.rds")

# check if we used the integrated data 
DefaultAssay(seu_int)
# [1] "integrated" --> correct

# find the neighbors
seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:25)

# find the clusters, test with different resolution values
seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))

# check the clustering information in the metadata table
head(seu_int@meta.data)

# display the clustering result
clustree::clustree(seu_int@meta.data[,grep("integrated_snn_res", colnames(seu_int@meta.data))],
                   prefix = "integrated_snn_res.")

# display the clustering result on the UMAP for a resolution of 0.1
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.1")

# display the clustering result on the UMAP for a resolution of 0.2
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.2")

# display the clustering result on the UMAP for a resolution of 0.3
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")

# display the clustering result on the UMAP for a resolution of 0.5
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.5")
# number of stable clusters --> 0.3 is good

# answer for the second question:
#
# When do the number of neighbors need to be changed? How does changing the 
# method of clustering in FindClusters affect the output? Which parameter should 
# be changed?
#
# As FindClusters is an unsupervised clustering method on the PCA data and UMAP 
# is a good summary of the PCA dimension selected, clusters and UMAP plot should 
# go along. If one has reasons to change the number of neighbors in the UMAP 
# function, here the same parameter should be adapted.
# 
# The method can be changed with algorithm = 2,3 or 4

