# command line for the normalization and scaling,
# next step after analysis_tools_QC_day1_exercice.R commands
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day1/normalization_scaling/
#Mon Jul  4 15:56:25 2022


# load the library


# normalization for cells -------------------------------------------------
# see for the assay before normalization
Seurat::GetAssayData(seu)[1:10,1:10]  

# perform the normalization
seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

# see for the assay after normalization
Seurat::GetAssayData(seu)[1:10,1:10]  
# yes


# scaling for genes (features) -------------------------------------------------

# selection of genes 
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

# display the top 10 selected genes
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

# scaling the genes (for PCA), only the variable genes are scaled by default
seu <- Seurat::ScaleData(seu,
                         features = rownames(seu))


# Bonus exercice with SCTransform
seu <- Seurat::SCTransform(seu, variable.features.n = 2000)
# save in seu@assays$SCT
# update the name
DefaultAssay(seu) <- "RNA"

# Save the dataset and clear environment -------------------------------------------------
saveRDS(seu, "seu_day1.rds")
rm(list = ls())
gc()
.rs.restartR()
