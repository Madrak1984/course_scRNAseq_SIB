# command line for the integration and dimension reduction,
# next step after normalization_scaling_metadata.R commands
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day2/dimensionality_reduction/
# Tue Jul  5 09:15:06 2022


# load the library
library(Seurat)

# load the R object created on previous day
seu <- readRDS("seu_day1.rds")


# dimension reduction -----------------------------------------------------

# run the PCA with Seurat
seu <- Seurat::RunPCA(seu)

# display the PCA
Seurat::DimPlot(seu, reduction = "pca")

# change the color according to new observation that we want to do
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.globin")

# change the color according to the HBA1 expression
Seurat::FeaturePlot(seu, reduction = "pca", features = "HBA1")

# change the color according to the IGKC expression
Seurat::FeaturePlot(seu, reduction = "pca", features = "IGKC")

# change the color according to the IGKC and HBA1 expression
Seurat::FeaturePlot(seu, reduction = "pca", features = c("HBA1", "IGKC"))

# change the color according to the IGKC and HBA1 expression
Seurat::FeaturePlot(seu, reduction = "pca", features = c("HBA1", "IGKC", "percent.globin"))


# display the heatmaps based on the principal component scores calculated on the 
# rotation matrix
Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)
# yellow colour, genes more expressed
# purple colour, genes low expressed

# select the best informative PC based on elbow point finding
Seurat::ElbowPlot(seu, ndims = 40)

# to run a t-sSNE with Seurat
?RunTSNE

# run a UMAP with Seurat
seu <- Seurat::RunUMAP(seu, dims = 1:25)

# display the UMAP
Seurat::DimPlot(seu, reduction = "umap")

# color the UMAP according to the percentage of globin
Seurat::FeaturePlot(seu, reduction = "umap", features = "percent.globin")

# color the UMAP according to the HBA1 expression
Seurat::FeaturePlot(seu, reduction = "umap", features = "HBA1")

# color the UMAP according to the percentage of ribosomal gene
Seurat::FeaturePlot(seu, reduction = "umap", features = "percent.ribo")
# erythrocytes are on the left based on the UMAP1 axis

# change the number of neighbors parameter for UMAP
?Seurat::RunUMAP
# n.neighbors parameter, 30L by default, concretely 30
# to increase or decrease to have a better representation of the dataset

# test UMAP by changing the number of dimension
# 1:5 dimension
seu <- Seurat::RunUMAP(seu, dims = 1:5)
# display the UMAP for 1:5 dimension
Seurat::DimPlot(seu, reduction = "umap")

# 1:50 dimension
seu <- Seurat::RunUMAP(seu, dims = 1:50)
# display the UMAP for 1:50 dimension
Seurat::DimPlot(seu, reduction = "umap")
# it affects the gap between the dots

# it is better to much to be more precise

# test with dim 1:100
seu <- Seurat::RunUMAP(seu, dims = 1:100)
# we have only 50 reduction calculated

# we need to be more precise to detect batch effect.



# integration -------------------------------------------------------------


# go back to the previous UMAP
seu <- Seurat::RunUMAP(seu, dims = 1:25)


# plot again before the integration
Seurat::DimPlot(seu, reduction = "umap")


# preprocess the integration
# split the object per samples
seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")

# perform log-normalization and find variable feature individually for each dataset
# based on vst method
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- Seurat::NormalizeData(seu_list[[i]])
  seu_list[[i]] <- Seurat::FindVariableFeatures(seu_list[[i]], 
                    selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# selection of integration anchor
seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list,
                                              dims = 1:25)
# see ?Seurat::FindIntegrationAnchors for all parameters:
# k.filter must have a number lower than the minimal number of cell in batch


# perform the integration
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors, dims = 1:25)

# switch the default assay by the integrated assay
# important to check the default assay because the integrated assay must not be 
# used for the differential expression step
Seurat::DefaultAssay(seu_int) <- "integrated"


# exercice
# scale the integrated assay
seu_int <- Seurat::ScaleData(seu_int)

# run the PCA for the integrated assay
seu_int <- Seurat::RunPCA(seu_int)

# run the UMAP for the integrated assay
seu_int <- Seurat::RunUMAP(seu_int, dims = 1:25)

# display the UMAP for 1:25 dimension
Seurat::DimPlot(seu_int, reduction = "umap")


# Save the dataset and clear environment -------------------------------------------------
saveRDS(seu_int, "seu_int_day2_part1.rds")
rm(list = ls())
gc()
.rs.restartR()
