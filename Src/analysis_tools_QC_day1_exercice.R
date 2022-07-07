# command line for the analysis tools and QC
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day1/analysis_tools_qc/
# Mon Jul  4 14:36:40 2022

# load the library
library(Seurat)
library(ggplot2)
library(Matrix)
library(Seurat)

# create functions

# function to display the most expressed genes
most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
  
  # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
  
  # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}


# load the samples info
sample_info <- read.csv("course_data/sample_info_course.csv")

#  save all the files with their corresponding path in a vector
datadirs <- file.path("course_data", "count_matrices", sample_info$SampleName,
                      "outs", "filtered_feature_bc_matrix")

#  names the elements of the vector according to the sample name, Seurat dislikes _ symbol
# generate a warning message if presence of _, replace by - (dashed)
names(datadirs) <- gsub("_", "-", sample_info$SampleName)

#  select the PBMMC files
datadirs <- datadirs[1:3]

# import the raw count from the result of cell ranger in one Seurat Object
sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
# read the 3 files located in the directory filtered_featute_bc_matrix directory:
# - matrix.mtx.gz: contains the count table
# - barcodes.tsv.gz: contains the cell names, corresponds to the columns names in the count table
# - features.tsv.gz: contains the genes names, corresponds to the row names in the count table (ensembl_Id and gene name)

# display the first 30 cells for the specified genes
sparse_matrix[c("PECAM1", "CD8A", "TSPAN1"), 1:30]

# create the Seurat object by importing genes which are expressed at minimum in 
# 3 cells and cells which express at minimum 100 genes
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "pbmmc",
                                  min.cells = 3,
                                  min.features = 100)

# see the count of imported cells and genes
dim(seu@assays$RNA@counts)
# [1] 18673  6946
# 18673 genes (features), 6946 cells
# or simply
seu

# to have more information about the functions available on Seurat
# (https://satijalab.org/seurat/articles/essential_commands.html)

# exercice
View(seu)
# allow to display the content of the Seurat object
# @active.ident contains the identity of selected cells

# @meta.data
seu@meta.data
# provide count of genes and UMI per cell

# histogram
hist(seu@meta.data$nCount_RNA, breaks = 1000)

# functions available on Seurat package to display plot based on SeuratObject
# function which used ggplot2 package and object
# source: https://satijalab.org/seurat/articles/essential_commands.html
# example
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# allow to see heterogeneity among samples

# Visualizing QC per cell and gene
Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))

# calculate relative abundance of mitochondrial, ribosomal and hemoglobin genes
# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")
# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")
# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")
# results are stored in 3 additional columns in meta.data
seu@meta.data

# plot these new results
Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))
# PBMMC-2 is different to the 2 other based on % of ribosomal and hemoglobin genes
# checked if these % are negatively correlated
Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")
# erythrocytes

# check the relative expression of genes
# display the most expressed genes
most_expressed_boxplot(seu)

# cell filtering based on # expressed genes between 200 and 5000, and percentage
# of mitochondrial gene < 8, based on publication, but
# we could determine threshold based on the observed values in previous plot for
# nFeature_RNA and percent.mito
seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)

# visualize the result
Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))