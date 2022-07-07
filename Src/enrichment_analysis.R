# command line for the enrichment analysis,
# next step after differential_gene_expression.R
# source: https://sib-swiss.github.io/single-cell-training/2022.7/day3/enrichment_analysis/
# Wed Jul  6 11:59:28 2022

# loading step ------------------------------------------------------------

# load the libraries
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# download the gene information
BiocManager::install("org.Hs.eg.db", update = FALSE)
AnnotationDbi::keytypes(org.Hs.eg.db)

# selection of down-regulated genes in tumor cells
tum_down <- subset(tum_vs_norm,
                   tum_vs_norm$avg_log2FC < -1 &
                     tum_vs_norm$p_val_adj < 0.05)
tum_down_genes <- rownames(tum_down)

# perform the enrichment GO analysis
tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)

# check the result
View(tum_vs_norm_go@result)

# remove redundant gene set
enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)

# display the enrichment map
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.5)


## perform with another database

# select the hallmark H database
gmt <- msigdbr::msigdbr(species = "human", category = "H")
?clusterProfiler::read.gmt # other possibility to use for the gmt

# perform the corresponding enrichment analysis
tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
                                                universe = rownames(proB),
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])

# check the result
View(tum_vs_norm_enrich@result[which(tum_vs_norm_enrich@result$p.adjust<0.05),])

# clean the environment -------------------------------------------------
rm(list = ls())
gc()
.rs.restartR()
