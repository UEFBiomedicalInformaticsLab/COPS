library(parallel)
library(WGCNA)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 

WGCNA::enableWGCNAThreads(PARALLEL)
gene_filter <- setdiff(rownames(tbrca_norm_deg), rownames(tbrca_norm)[zero_var])
gene_correlation <- WGCNA::cor(t(tbrca_norm[gene_filter,]), method =  "spearman")
adj <- WGCNA::adjacency.fromSimilarity(gene_correlation, power = MODULE_POWER)
TOM <- WGCNA::TOMsimilarity(adj, TOMType = "unsigned")
geneTree <- flashClust::flashClust(as.dist(1 - TOM), method="average")
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 20)
adj_modules <- WGCNA::adjacency.fromSimilarity(gene_correlation[dynamicMods != 0, dynamicMods != 0], power = MODULE_POWER)
MEList <- WGCNA::moduleEigengenes(t(log2(tbrca_norm[gene_filter,][dynamicMods != 0,] + 1)), 
                                  colors = dynamicMods[dynamicMods != 0])
MEs <- MEList$eigengenes

write.csv(MEs, paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"))
write.csv(data.frame(SYMBOL = gene_filter[dynamicMods != 0], module = dynamicMods[dynamicMods != 0]),
          paste0(path_intermediate_results, "/brca/wgcna/module_genes_deg.csv"))