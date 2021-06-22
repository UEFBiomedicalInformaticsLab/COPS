source("load_config.R")

# Load results, keep and rename relevant factors for plots
#source("article/plot_renaming.R")

source("extra/load_deg_results.R")
#source("brca/brca_load_results.R")
#brca_scores_sampled <- rbind(brca_scores_sampled, brca_scores_deg_sampled)
#brca_scores <- rbind(brca_scores, brca_scores_deg)
brca_scores_sampled <- brca_scores_deg_sampled
#brca_scores <- brca_scores_deg

brca_scores_sampled <- brca_scores_sampled[brca_scores_sampled$Clustering != "HC_single",]
temp <- brca_scores_sampled[brca_scores_sampled$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                              ! brca_scores_sampled$Approach %in% c("DR*", "DR", "DR DEG") | 
                              brca_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
brca_scores_sampled <- rename_dr_methods(temp)
#brca_scores_sampled$Approach <- factor(brca_scores_sampled$Approach, levels = unique(brca_scores_sampled$Approach)[c(2,3,5,1,4,6)])

#brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]
#temp <- brca_scores[brca_scores$drname %in% c("umap5_20n", "tsne45", "pca2") | 
#                      ! brca_scores$Approach %in% c("DR*", "DR", "DR DEG") | 
#                      brca_scores$drname == "VAE",]
#temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
#brca_scores <- rename_dr_methods(temp)
#brca_scores$Approach <- factor(brca_scores$Approach, levels = unique(brca_scores$Approach)[c(2,3,5,1,4,6)])

#brca_scores$Unsupervised_wsum <- brca_scores$Silhouette + 
#  brca_scores$ClusteringStabilityJaccard + 
#  brca_scores$Module_score + 
#  brca_scores$SurvivalPValue_score + 
#  1 - brca_scores$NMI.tss

brca_scores_sampled$Unsupervised_wsum <- brca_scores_sampled$Silhouette + 
  brca_scores_sampled$ClusteringStabilityJaccard + 
  brca_scores_sampled$Module_score + 
  brca_scores_sampled$SurvivalPValue_score + 
  1 - brca_scores_sampled$NMI.tss

#source("prad/prad_load_results.R")
#prad_scores_sampled <- rbind(prad_scores_sampled, prad_scores_deg_sampled)
#prad_scores <- rbind(prad_scores, prad_scores_deg)
prad_scores_sampled <- prad_scores_deg_sampled
#prad_scores <- prad_scores_deg

prad_scores_sampled <- prad_scores_sampled[prad_scores_sampled$Clustering != "HC_single",]
temp <- prad_scores_sampled[prad_scores_sampled$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                              ! prad_scores_sampled$Approach %in% c("DR*", "DR", "DR DEG") |
                              prad_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
prad_scores_sampled <- rename_dr_methods(temp)
#prad_scores_sampled$Approach <- factor(prad_scores_sampled$Approach, levels = unique(prad_scores_sampled$Approach)[c(2,3,5,1,4,6)])

#prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]
#temp <- prad_scores[prad_scores$drname %in% c("umap2_20n", "tsne45", "pca2") | 
#                      ! prad_scores$Approach %in% c("DR*", "DR", "DR DEG") |
#                      prad_scores$drname == "VAE",]
#temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
#prad_scores <- rename_dr_methods(temp)
#prad_scores$Approach <- factor(prad_scores$Approach, levels = unique(prad_scores$Approach)[c(2,3,5,1,4,6)])

#prad_scores$Unsupervised_wsum <- prad_scores$Silhouette + 
#  prad_scores$ClusteringStabilityJaccard + 
#  prad_scores$Module_score + 
#  prad_scores$SurvivalPValue_score + 
#  1 - prad_scores$NMI.tss

prad_scores_sampled$Unsupervised_wsum <- prad_scores_sampled$Silhouette + 
  prad_scores_sampled$ClusteringStabilityJaccard + 
  prad_scores_sampled$Module_score + 
  prad_scores_sampled$SurvivalPValue_score + 
  1 - prad_scores_sampled$NMI.tss

# Summarise
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

brca_table <- plyr::ddply(brca_scores_sampled[brca_scores_sampled$fold != 6,], summary_by, summary_function)
brca_table$fold <- NULL
brca_table$run <- NULL

prad_table <- plyr::ddply(prad_scores_sampled[prad_scores_sampled$fold != 6,], summary_by, summary_function)
prad_table$fold <- NULL
prad_table$run <- NULL

metrics <- c("Silhouette_mean", "Silhouette_sd", "ClusteringStabilityJaccard_mean", "ClusteringStabilityJaccard_sd", 
             "cNMI_mean", "cNMI_sd", "Module_score_mean", "Module_score_sd", 
             "SurvivalPValue_Q05", "SurvivalPValue_median", "SurvivalPValue_Q95")

# All in one
writexl::write_xlsx(brca_table[, c(summary_by, metrics)],
                    path = paste0(path_plots, "/brca_deg_scores.xlsx"), 
                    col_names = TRUE)
writexl::write_xlsx(prad_table[, c(summary_by, metrics)],
                    path = paste0(path_plots, "/prad_deg_scores.xlsx"), 
                    col_names = TRUE)

# DR-CL
writexl::write_xlsx(brca_table[brca_table$Approach %in% c("DR DEG"), c(summary_by, metrics)],
                    path = paste0(path_plots, "/brca_drcl_deg_scores.xlsx"), 
                    col_names = TRUE)
writexl::write_xlsx(prad_table[prad_table$Approach %in% c("DR DEG"), c(summary_by, metrics)],
                    path = paste0(path_plots, "/prad_drcl_deg_scores.xlsx"), 
                    col_names = TRUE)



