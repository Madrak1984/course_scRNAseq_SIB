# command line for the differential gene expression,
# next step after cell_annotation commands, but must work on raw count
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day3/differential_gene_expression/
# Wed Jul  6 09:13:14 2022


# loading step ------------------------------------------------------------

# load the SeuratObject from previous work
seu_int <- readRDS("seu_int_day2_part2.rds")

# load the libraries
library(Seurat)
library(edgeR) # BiocManager::install("edgeR")
library(limma)
library(dplyr)
library(scuttle)


# differential expression analysis ----------------------------------------

# check if we are using the raw count
DefaultAssay(seu_int)
# yes

# perform the differential expression analysis for each cluster against all other cells
# --> find the marker genes in each cluster
de_genes <- Seurat::FindAllMarkers(seu_int,  min.pct = 0.25,
                                   only.pos = TRUE)

# subset the table to keep only DEG (adj-pval <= 0.05)
de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)
head(de_genes)
# pct.1 indicates the percentage of cells in the first group, and pct.2 in the second group

# write the result in a csv file
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)

# select the top 3 specific markers
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# display the selection
dittoSeq::dittoDotPlot(seu_int, vars = unique(top_specific_markers$gene), 
                       group.by = "integrated_snn_res.0.3")
# significant marker in cluster 0: CD3D, MALAT1, CD3E
# significant marker in cluster 8: KRLB1, CCL5, GNLY

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
de_genes[de_genes$gene %in% tcell_genes,]
# yes


# Differential expression between groups of cells -------------------------

# select the annotation performed by singleR (to replace the gene selection for the clustering)
seu_int <- Seurat::SetIdent(seu_int, value = "SingleR_annot")

# perform comparison between LT8 vs LT4 cells
deg_cd8_cd4 <- Seurat::FindMarkers(seu_int,
                                   ident.1 = "CD8+ T cells",
                                   ident.2 = "CD4+ T cells",
                                   group.by = seu_int$SingleR_annot,
                                   test.use = "wilcox")

# select the significant DEG
deg_cd8_cd4 <- subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj<0.05)

# Are CD8A, CD8B and CD4 in there? 
deg_cd8_cd4[c("CD8A", "CD8B", "CD4"),]
# yes
# What does the sign (i.e. positive or negative) mean in the log fold change values? 
# if the gene is over-expressed (positive sign) or under-expressed (negative sign) in CD8+ T cells.
# Are they according to the CD8+ and CD4+ annotations?
Seurat::VlnPlot(seu_int, features = c("CD8A", "CD8B", "CD4"),
                idents = c("CD8+ T cells", "CD4+ T cells"))
# idents arguments allows to select the terms to visualize.


# differential expression with limma --------------------------------------

# load the data
proB <- readRDS("course_data/proB.rds")

# display the UMAP
DimPlot(proB, group.by = "orig.ident")

# count the cells
table(proB@meta.data$type)
# ETV6-RUNX1      PBMMC 
#      2000       1021

# display the cell information
head(proB@meta.data)

# prepare the data for the analysis
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

# check the UMAP
Seurat::DimPlot(proB)

# pseudo-bulk analysis preparation 
#taking the proB data 
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

## add the patient id also for paired DGE
proB$patient.id<-gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

## Here we do perform pseudo-bulk:
##first a mandatory column of sample needs to be added to the meta data that is the grouping factor, should be the samples
proB$sample <- factor(proB$orig.ident)

##first an sce object is needed
sce_proB <- as.SingleCellExperiment(proB)

##aggregateAcrossCells here it is only aggregated by sample, one could imagine
##to aggregate by sample and by celltype for instance
summed <- aggregateAcrossCells(sce_proB, 
                               id=colData(sce_proB)[,c("sample")])

##have a look at the counts
counts(summed)[1:3,]

#have a look at the colData of our new object summed, can you see type and 
#patient.id are there
head(colData(summed))

#As in the standard limma analysis generate a DGE object

y <- DGEList(counts(summed), samples=colData(summed)$sample)
summary(keep)

##filter lowly expressed (recommanded for limma)
keep <- filterByExpr(y, group=summed$type)
y <- y[keep,]

##see how many genes were kept 
summary(keep)

## Create the design matrix and include the technology as a covariate:
design <- model.matrix(~0 + summed$type + summed$patient.id)

# Have a look
design

# change column/rownames names to more simple group names: 
colnames(design) <- make.names(c("ETV6-RUNX1", "PBMMC","patient2","patient3"))
rownames(design)<-colData(summed)$sample

# Have a look
design

# specify the contrast to analyze
contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC, levels = design) # tumoral vs ctrl

# perfom the normalization
dge <- edgeR::calcNormFactors(y)  

#Do limma
vm <- limma::voom(dge, design = design, plot = TRUE) # necessary for RNAseq, limma was created for microarray analysis
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)


# Show the top differentially expressed genes:
limma::topTable(fit.contrasts, number = 10, sort.by = "P")
limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))
# 2738

# check result in violin plot
Seurat::VlnPlot(proB, "S100A9", split.by = "type")
Seurat::VlnPlot(proB, "SOCS2", split.by = "type")

# same analysis with Seurat (does not take in account that the samples are paired)
tum_vs_norm <- Seurat::FindMarkers(proB, 
                                   ident.1 = "ETV6-RUNX1", 
                                   ident.2 = "PBMMC", 
                                   group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj<0.05)
# How many genes are significant? 
dim(tum_vs_norm)
# 1893 DEGs
# How does the fold change of these genes compare to the fold change of the top genes found by limma?
limma_de[c("S100A9", "SOCS2"),]
tum_vs_norm[c("S100A9", "SOCS2"),]
# FC are lower in the Seurat analysis
# correct answer:
# If we merge the FindMarkers and the limma results, keep limma‘s most significant genes and plot:
merge_limma_FindMarkers <- merge(tum_vs_norm, limma_de, by="row.names",
                                   all.x=T)
par(mar=c(4,4,4,4))
plot(merge_limma_FindMarkers$avg_log2FC,
     merge_limma_FindMarkers$logFC,
     xlab="log2FC Wilcoxon", ylab="log2FC limma",
     pch=15, cex=0.5)
abline(a=0, b=1, col="red")

