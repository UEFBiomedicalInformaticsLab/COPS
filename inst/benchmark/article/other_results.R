#library(ggplot2)

## BRCA
# kNN community analysis results
brca_res_knn_comm <- read.csv(paste0(path_intermediate_results, "/brca/knn_communities/scores.csv"), header = TRUE, row.names = 1)
brca_res_knn_comm$Approach <- "Direct clustering"
brca_res_knn_comm$cNMI <- 0.5 + brca_res_knn_comm$NMI.BRCA_Subtype_PAM50 / 2 - brca_res_knn_comm$NMI.tss / 2

# Check the number of communities
#brca_clust_knn_comm <- read.csv(paste0(path_intermediate_results, "/brca/knn_communities/clusters.csv.gz"), header = TRUE, row.names = 1)
#brca_knn_conn_k <- table(tapply(brca_clust_knn_comm$cluster, brca_clust_knn_comm[,c("run", "fold")], function(x) length(unique(x))))

# Spectrum 
brca_res_spectral <- read.csv(paste0(path_intermediate_results, "/brca/spectral/scores.csv"), header = TRUE, row.names = 1)
brca_res_spectral$Approach <- "Direct clustering"
brca_res_spectral$cNMI <- brca_res_spectral$wsum

# SC3
brca_res_sc3 <- read.csv(paste0(path_intermediate_results, "/brca/sc3/scores.csv"), header = TRUE, row.names = 1)
brca_res_sc3$Approach <- "SC3"
brca_res_sc3$drname <- "SC3"
brca_res_sc3$datname <- "SC3"
brca_res_sc3$cNMI <- 0.5 + brca_res_sc3$NMI.BRCA_Subtype_PAM50 / 2 - brca_res_sc3$NMI.tss / 2

## PRAD
# kNN community analysis results
prad_res_knn_comm <- read.csv(paste0(path_intermediate_results, "/prad/knn_communities/scores.csv"), header = TRUE, row.names = 1)
prad_res_knn_comm$Approach <- "Direct clustering"
prad_res_knn_comm$cNMI <- 0.5 + prad_res_knn_comm$NMI.Gleason_category / 2 - prad_res_knn_comm$NMI.tss / 2

# Check the number of communities
#prad_clust_knn_comm <-  read.csv(paste0(path_intermediate_results, "/prad/knn_communities/clusters.csv.gz"), header = TRUE, row.names = 1)
#table(tapply(prad_clust_knn_comm$cluster, prad_clust_knn_comm[,c("run", "fold")], function(x) length(unique(x))))

# Spectrum 
prad_res_spectral <- read.csv(paste0(path_intermediate_results, "/prad/spectral/scores.csv"), header = TRUE, row.names = 1)
prad_res_spectral$Approach <- "Direct clustering"
prad_res_spectral$cNMI <- prad_res_spectral$wsum

# SC3

prad_res_sc3 <- read.csv(paste0(path_intermediate_results, "/prad/sc3/scores.csv"), header = TRUE, row.names = 1)
prad_res_sc3$Approach <- "SC3"
prad_res_sc3$drname <- "SC3"
prad_res_sc3$datname <- "SC3"
prad_res_sc3$cNMI <- 0.5 + prad_res_sc3$NMI.Gleason_category / 2 - prad_res_sc3$NMI.tss / 2

## Tables

brca_res_other <- Reduce("rbind", list(brca_res_knn_comm, brca_res_spectral, brca_res_sc3))
prad_res_other <- Reduce("rbind", list(prad_res_knn_comm, prad_res_spectral, prad_res_sc3))

source("article/plot_renaming.R")
brca_res_other <- rename_methods(brca_res_other)
prad_res_other <- rename_methods(prad_res_other)

summary_by <- c("datname", "drname", "k", "m", "Approach", "Method", "Embedding", "Clustering")

# Calculate SD
summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[paste0(i, "_mean")]] <- mean(x[[i]], na.rm = TRUE)
    out[[paste0(i, "_sd")]] <- sd(x[[i]], na.rm = TRUE)
    out[[paste0(i, "_median")]] <- median(x[[i]], na.rm = TRUE)
    out[[paste0(i, "_Q95")]] <- quantile(x[[i]], probs = 0.95, na.rm = TRUE)
    out[[paste0(i, "_Q05")]] <- quantile(x[[i]], probs = 0.05, na.rm = TRUE)
  }
  return(out)
}

brca_table_other <- plyr::ddply(brca_res_other[brca_res_other$fold != 6,], summary_by, summary_function)
brca_table_other$fold <- NULL
brca_table_other$run <- NULL

prad_table_other <- plyr::ddply(prad_res_other[prad_res_other$fold != 6,], summary_by, summary_function)
prad_table_other$fold <- NULL
prad_table_other$run <- NULL

metrics <- c("Silhouette_mean", "Silhouette_sd", "ClusteringStabilityJaccard_mean", "ClusteringStabilityJaccard_sd", 
             "cNMI_mean", "cNMI_sd", "Module_score_mean", "Module_score_sd", 
             "SurvivalPValue_Q05", "SurvivalPValue_median", "SurvivalPValue_Q95")

# All in one
writexl::write_xlsx(brca_table_other[, c(summary_by, metrics)],
                    path = paste0(path_plots, "/brca_extra_methods.xlsx"), 
                    col_names = TRUE)
writexl::write_xlsx(prad_table_other[, c(summary_by, metrics)],
                    path = paste0(path_plots, "/prad_extra_methods.xlsx"), 
                    col_names = TRUE)

# One table per method
for (clustering_i in c("spectral", "knn_communities", "sc3")) {
  writexl::write_xlsx(brca_table_other[brca_table_other$Clustering == clustering_i, 
                                       c(summary_by, metrics)],
                      path = paste0(path_plots, "/brca_", clustering_i, ".xlsx"), 
                      col_names = TRUE)
}

for (clustering_i in c("spectral", "knn_communities", "sc3")) {
  writexl::write_xlsx(prad_table_other[prad_table_other$Clustering == clustering_i, 
                                       c(summary_by, metrics)],
                      path = paste0(path_plots, "/prad_", clustering_i, ".xlsx"), 
                      col_names = TRUE)
}