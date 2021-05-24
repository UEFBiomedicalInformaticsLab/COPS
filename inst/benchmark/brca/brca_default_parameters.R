# General
PARALLEL <- 10
NRUNS <- 20
NFOLDS <- 5
SUMMARISE <- FALSE

# Data
OTP_CUTOFF <- 0.0

# DR
PCA_DIMS <- c(2,5,10)
UMAP_DIMS <- c(2,5,10)
TSNE_PERPLEXITY <- c(15,30,45)

# Clustering
NCLUSTERS <- 2:6
DR_CLUST_METHODS <- c("hierarchical", "kmeans", "diana", "model")
OG_CLUST_METHODS <- c("hierarchical", "diana")
DR_LINKAGES <- c("complete", "average", "ward")
OG_LINKAGES <- c("complete", "average")
DR_METRIC <- "euclidean"
OG_METRIC <- "correlation"
COR_METHOD <- "pearson"
GMM_MODEL <- "VVV"

# Random walk
RWR_SEEDS <- 25
RWR_RESTART_PROB <- 0.75
RWR_GCN_THRESHOLD <- 0.475

# Gene modules
MODULE_POWER <- 6

# Random seed
set.seed(0)