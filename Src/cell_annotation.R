# command line for the cell annotation,
# next step after clustering.R commands
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day2/cell_annotation/
# Tue Jul  5 16:14:57 2022

# load the libraries
library(celldex)
library(SingleR)

# select the clustering with the high sensibility
seu_int <- Seurat::SetIdent(seu_int, value = seu_int$integrated_snn_res.0.3)

# must use the original assay (not the integrated)
DefaultAssay(seu_int) <- "RNA"

# example of expression for a gene in each cluster
Seurat::FeaturePlot(seu_int, "HBA1")


# cell annotation ---------------------------------------------------------

# prepare known marker gene for annotation
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

# see the expression of these cells
Seurat::FeaturePlot(seu_int, tcell_genes, ncol=2)

# check in which clusters are T-cells
Seurat::VlnPlot(seu_int,
                features = tcell_genes,
                ncol = 2)
# clusters 8 and 0

# check in which clusters are monocytes
Seurat::FeaturePlot(seu_int, monocyte_genes, ncol=2)
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")
Seurat::VlnPlot(seu_int,
                features = monocyte_genes,
                ncol = 2)
# cluster 2
                    
# automatisation per cell with addModuleScore
seu_int <- Seurat::AddModuleScore(seu_int,
                                  features = list(tcell_genes),
                                  name = "tcell_genes")
str(seu_int@meta.data)
# the name of the added column is tcell_genes1. It contains numeric values.

Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")
# display the UMAP for this column
Seurat::FeaturePlot(seu_int, reduction = "umap", features = "tcell_genes1")
# display the violinplot for this column
Seurat::VlnPlot(seu_int, features = "tcell_genes1")
# yes, it is similar to the previous result.


# phase annotation --------------------------------------------------------

# use the CellCycleScore function provided by Seurat package

# extract gene markers for cell cycling
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# run the CellCycleScore function
seu_int <- Seurat::CellCycleScoring(seu_int,
                                    s.features = s.genes,
                                    g2m.features = g2m.genes)

# visualization of result
Seurat::DimPlot(seu_int, group.by = "Phase")

#### note: it is possible to remove the phase effect in the scale step by using
# the parameter vars.to.regress of Seurat::ScaleData function (e.g vars.to.regress = "Phase)


# cell annotation with SingleR --------------------------------------------


# take reference RNAseq from celldex (https://bioconductor.org/packages/3.14/data/experiment/vignettes/celldex/inst/doc/userguide.html)
ref <- celldex::NovershternHematopoieticData()
class(ref)
table(ref$label.main)

# run the singleR analysis
seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, slot = "data"),
                                    ref = ref,
                                    labels = ref$label.main)

# display the first rows
head(seu_int_SingleR)

# display singleR quality score
SingleR::plotScoreHeatmap(seu_int_SingleR)
SingleR::plotDeltaDistribution(seu_int_SingleR)

# remove annotation with low cells
singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

# results are added in the meta.data table
seu_int$SingleR_annot <- singleR_labels

# display the singleR annotation on the UMAP
dittoSeq::dittoDimPlot(seu_int, "SingleR_annot", size = 0.7)

# display the count of cell per annotation
dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")

# comparison with manual annotation
dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "integrated_snn_res.0.3")
# it seems based on the T-cells and monocytes

# Save the dataset and clear environment -------------------------------------------------
saveRDS(seu_int, "seu_int_day2_part2.rds")
rm(list = ls())
gc()
.rs.restartR()
