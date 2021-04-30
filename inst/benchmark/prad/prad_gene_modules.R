library(parallel)
library(WGCNA)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 

WGCNA::enableWGCNAThreads(PARALLEL)
gene_filter <- setdiff(rownames(tprad_norm_deg), rownames(tprad_norm)[zero_var])
gene_correlation <- WGCNA::cor(t(tprad_norm[gene_filter,]), method =  "spearman")
adj <- WGCNA::adjacency.fromSimilarity(gene_correlation, power = MODULE_POWER)
TOM <- WGCNA::TOMsimilarity(adj, TOMType = "unsigned")
geneTree <- flashClust::flashClust(as.dist(1 - TOM), method="average")
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 20)
adj_modules <- WGCNA::adjacency.fromSimilarity(gene_correlation[dynamicMods != 0, dynamicMods != 0], power = MODULE_POWER)
MEList <- WGCNA::moduleEigengenes(t(log2(tprad_norm[gene_filter,][dynamicMods != 0,] + 1)), 
                                  colors = dynamicMods[dynamicMods != 0])
MEs <- MEList$eigengenes

write.csv(MEs, paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"))
write.csv(data.frame(SYMBOL = gene_filter[dynamicMods != 0], module = dynamicMods[dynamicMods != 0]),
          paste0(path_intermediate_results, "/prad/wgcna/module_genes_deg.csv"))