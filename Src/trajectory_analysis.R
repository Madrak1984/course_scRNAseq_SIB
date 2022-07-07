# command line for the trajectory analysis,
# next step after enrichment_analysis.R
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day3/trajectory_analysis/
# Wed Jul  6 16:01:40 2022

# loading step ------------------------------------------------------------

# load the libraries
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)
library(monocle3)

# load the data
deng_SCE <- readRDS("course_data/deng-reads.rds")

# reorder the cell type information
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy",
                                         "early2cell",
                                         "mid2cell",
                                         "late2cell",
                                         "4cell",
                                         "8cell",
                                         "16cell",
                                         "earlyblast",
                                         "midblast",
                                         "lateblast"))

# do the trajectory with slingshot -------------------------------------------------


# run the PCA
deng_SCE <- scater::runPCA(deng_SCE, ncomponents = 50)

# save the result of the PCA
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")
head(pca)

# store the PCA result in the sce object
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

# plot the PCA
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# rank cells by their PC1 score
deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  

# do the jitter plot of the PC1 according to the cell type
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")


# calculate the trajectory with slingshot
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')
# we omit to indicate the starting cluster
# correction: we did not indicate the clusters in the clusterLabels parameters.

# create a custon function to plot the PCA on the slingshot
PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  # set palette for factorial variables
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  # set palette for numeric variables
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  # extract pca from SingleCellExperiment object
  pca <- SingleCellExperiment::reducedDims(sce)$PCA
  
  if(is.null(variable)){
    col <- "black"
  }
  if(is.character(variable)){
    variable <- as.factor(variable)
  }
  if(is.factor(variable)){
    colpal <- palf(length(levels(variable)))
    colors <- colpal[variable]
  }
  if(is.numeric(variable)){
    colpal <- paln(50)
    colors <- colpal[cut(variable,breaks=50)]
  }
  
  # draw the plot
  plot(pca, bg = colors, pch = 21)
  # draw lines
  if(draw_lines){
    lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  # add legend
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)
    
  }
}


# display the plot
PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

# display the pseudotime vs the celltype
ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
                                             y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# improvement
gcdata <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(deng_SCE),
                                     project = "slingshot")

gcdata <- Seurat::NormalizeData(object = gcdata,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

gcdata <- Seurat::FindVariableFeatures(object = gcdata,
                                       mean.function = ExpMean,
                                       dispersion.function = LogVMR)

gcdata <- Seurat::ScaleData(object = gcdata,
                            do.center = T,
                            do.scale = F)

gcdata <- Seurat::RunPCA(object = gcdata,
                         pc.genes = gcdata@var.genes)

gcdata <- Seurat::FindNeighbors(gcdata,
                                reduction = "pca",
                                dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
                               resolution = 0.6)

# display the slingshot with the clusters
deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

sce <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            start.clus = "2")

# check the slingshot evolution
SlingshotDataSet(sce)

# Plot PC1 versus PC2 colored by slingshot pseudotime:
PCAplot_slingshot(sce, variable = sce$slingPseudotime_2)

# Plot Slingshot pseudotime vs cell stage.
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_1 = sce$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_2 = sce$slingPseudotime_2),
       aes(x = slingPseudotime_2, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# plot the slingshot according to the cluster labelling
PCAplot_slingshot(sce,
                  variable = deng_SCE$Seurat_clusters,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

# plot the slingshot according to the cell type
PCAplot_slingshot(sce,
                  variable = deng_SCE$cell_type2,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

# select the ending cluster
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_1 = sce$slingPseudotime_1,
                  Seurat_clusters = deng_SCE$Seurat_clusters),
       aes(x = slingPseudotime_1, y = cell_type2,
           colour = Seurat_clusters)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
# cluster 0

# run the slingshot
sce_end <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            end.clus = "0")
# plot the slingshot according to the cell type
PCAplot_slingshot(sce_end,
                  variable = deng_SCE$cell_type2,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)
# it is similar
# answer:
sce <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            end.clus = c("0", "3", "5")) ## check which would be the best according to bio


# clean the environment -------------------------------------------------

rm(list = ls())
gc()
.rs.restartR()


# do the trajectory with Monocle3 -------------------------------------------------

# load the data
seu_int <- readRDS("seu_int_day2_part2.rds")

# create the monocle3 object
# get matrix and filter for minimum number of cells and features (the latter is a fix for backward compatibility)
mat_tmp <- seu_int@assays$RNA@counts
seu_tmp <- Seurat::CreateSeuratObject(mat_tmp, min.cells = 3,
                                      min.features = 100)

feature_names <- as.data.frame(rownames(seu_tmp))
rownames(feature_names) <- rownames(seu_tmp)
colnames(feature_names) <- "gene_short_name"

seu_int_monocl <- monocle3::new_cell_data_set(seu_tmp@assays$RNA@counts,
                                              cell_metadata = seu_int@meta.data,
                                              gene_metadata = feature_names)

# preprocess the data
?preprocess_cds
seu_int_monocl <- monocle3::preprocess_cds(seu_int_monocl)

# check the elbow point
monocle3::plot_pc_variance_explained(seu_int_monocl)

# perform the UMAP
seu_int_monocl <- monocle3::reduce_dimension(seu_int_monocl, reduction_method = "UMAP")

# plot the monocle3 object
# according to the clustering
monocle3::plot_cells(seu_int_monocl, 
                     color_cells_by = "integrated_snn_res.0.3", 
                     cell_size = 1, 
                     show_trajectory_graph = FALSE)

# according to the CD79A gene expression
monocle3::plot_cells(seu_int_monocl, genes = "CD79A", 
                     show_trajectory_graph = FALSE, 
                     cell_size = 1)

# Cluster cells using monocle3‘s clustering function
seu_int_monocl <- monocle3::cluster_cells(seu_int_monocl, resolution=0.00025)
monocle3::plot_cells(seu_int_monocl, label_cell_groups = F)


# perform the trajectory
seu_int_monocl <- monocle3::learn_graph(seu_int_monocl)
monocle3::plot_cells(seu_int_monocl)

monocle3::plot_cells(seu_int_monocl, genes = c("CD79A", "CD34"), 
                     show_trajectory_graph = FALSE, 
                     cell_size = 1)
# cluster 7

# Select the “initial” cells in the B-cell cluster to calculate pseudotime
seu_int_monocl<-monocle3::order_cells(seu_int_monocl)# issue with the shinyapp

monocle3::plot_cells(seu_int_monocl,
                     color_cells_by = "pseudotime",
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=FALSE,
                     graph_label_size=1.5, cell_size = 1)
