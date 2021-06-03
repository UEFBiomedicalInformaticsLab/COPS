# Load results, keep and rename relevant factors for plots
#source("article/plot_renaming.R")

source("extra/load_deg_results.R")
source("brca/brca_load_results.R")
brca_scores_sampled <- rbind(brca_scores_sampled, brca_scores_deg_sampled)
brca_scores <- rbind(brca_scores, brca_scores_deg)

brca_scores_sampled <- brca_scores_sampled[brca_scores_sampled$Clustering != "HC_single",]
temp <- brca_scores_sampled[brca_scores_sampled$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                              ! brca_scores_sampled$Approach %in% c("DR*", "DR", "DR, DEG") | 
                              brca_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
brca_scores_sampled <- rename_dr_methods(temp)
#brca_scores_sampled$Approach <- factor(brca_scores_sampled$Approach, levels = unique(brca_scores_sampled$Approach)[c(2,3,5,1,4,6)])

brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]
temp <- brca_scores[brca_scores$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                      ! brca_scores$Approach %in% c("DR*", "DR", "DR DEG") | 
                      brca_scores$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
brca_scores <- rename_dr_methods(temp)
#brca_scores$Approach <- factor(brca_scores$Approach, levels = unique(brca_scores$Approach)[c(2,3,5,1,4,6)])

brca_scores$Unsupervised_wsum <- brca_scores$Silhouette + 
  brca_scores$ClusteringStabilityJaccard + 
  brca_scores$Module_score + 
  brca_scores$SurvivalPValue_score + 
  1 - brca_scores$NMI.tss

brca_scores_sampled$Unsupervised_wsum <- brca_scores_sampled$Silhouette + 
  brca_scores_sampled$ClusteringStabilityJaccard + 
  brca_scores_sampled$Module_score + 
  brca_scores_sampled$SurvivalPValue_score + 
  1 - brca_scores_sampled$NMI.tss

source("prad/prad_load_results.R")
prad_scores_sampled <- rbind(prad_scores_sampled, prad_scores_deg_sampled)
prad_scores <- rbind(prad_scores, prad_scores_deg)

prad_scores_sampled <- prad_scores_sampled[prad_scores_sampled$Clustering != "HC_single",]
temp <- prad_scores_sampled[prad_scores_sampled$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                              ! prad_scores_sampled$Approach %in% c("DR*", "DR", "DR DEG") |
                              prad_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
prad_scores_sampled <- rename_dr_methods(temp)
#prad_scores_sampled$Approach <- factor(prad_scores_sampled$Approach, levels = unique(prad_scores_sampled$Approach)[c(2,3,5,1,4,6)])

prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]
temp <- prad_scores[prad_scores$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                      ! prad_scores$Approach %in% c("DR*", "DR", "DR DEG") |
                      prad_scores$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
prad_scores <- rename_dr_methods(temp)
#prad_scores$Approach <- factor(prad_scores$Approach, levels = unique(prad_scores$Approach)[c(2,3,5,1,4,6)])

prad_scores$Unsupervised_wsum <- prad_scores$Silhouette + 
  prad_scores$ClusteringStabilityJaccard + 
  prad_scores$Module_score + 
  prad_scores$SurvivalPValue_score + 
  1 - prad_scores$NMI.tss

prad_scores_sampled$Unsupervised_wsum <- prad_scores_sampled$Silhouette + 
  prad_scores_sampled$ClusteringStabilityJaccard + 
  prad_scores_sampled$Module_score + 
  prad_scores_sampled$SurvivalPValue_score + 
  1 - prad_scores_sampled$NMI.tss

## BRCA
# Calculate SD
sd_summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[paste0(i, "_sd")]] <- sd(x[[i]], na.rm = TRUE)
  }
  return(out)
}


brca_scores_sd <- plyr::ddply(brca_scores_sampled, summary_by, sd_summary_function)
brca_scores_sd$fold <- NULL
brca_scores_sd$run <- NULL

brca_table <- plyr::join(brca_scores, brca_scores_sd, by = summary_by)

prad_scores_sd <- plyr::ddply(prad_scores_sampled, summary_by, sd_summary_function)
prad_scores_sd$fold <- NULL
prad_scores_sd$run <- NULL

prad_table <- plyr::join(prad_scores, prad_scores_sd, by = summary_by)

# Selections
brca_means_nk <- plyr::ddply(brca_table, 
                             summary_by[summary_by != "k"],
                             function(x) data.frame(nk_mean = mean(x$Unsupervised_wsum, na.rm = TRUE)))
brca_nk_best <- plyr::ddply(brca_means_nk, 
                            c("Approach", "Embedding"),
                            function(x) data.frame(Clustering = x$Clustering[which.max(x$nk_mean)]))
brca_table <- plyr::join(brca_nk_best, brca_table)

# Create tables
brca_wsum_table <- tidyr::pivot_wider(brca_table, 
                                         id_cols = c("Approach", "Embedding", "Clustering"),
                                         names_from = c("k"),
                                         names_sep = "_",
                                         values_from = c("Unsupervised_wsum", "Unsupervised_wsum_sd"))
#ind <- (2:(ncol(brca_wsum_table) %/% 2)) * 2

#brca_wsum_table$DR_method <- sapply(strsplit(brca_wsum_table$Embedding, " "), function(x) x[[1]])
#brca_wsum_table <- dplyr::arrange(brca_wsum_table, DR_method, Approach)

#write.table(brca_wsum_table[,c(1,2, 2 + order(c(ind - 1, ind)))], 
#            file = paste0(path_plots, "/brca_silh.csv"), 
#            sep = ";", dec = ",", row.names = FALSE, quote = FALSE)

temp <- as.data.frame(t(1:ncol(brca_wsum_table)))
colnames(temp) <- colnames(brca_wsum_table)

brca_wsum_table <- rbind(temp, brca_wsum_table)
brca_wsum_table <- rbind(temp, brca_wsum_table)

brca_wsum_table[1,] <- c("", "", "", sapply(strsplit(colnames(brca_wsum_table), "_")[-c(1:3)], 
                                               function(x) paste0("k=", x[[length(x)]])))
brca_wsum_table[2,] <- c("", "", "", sapply(strsplit(colnames(brca_wsum_table), "_")[-c(1:3)], 
                                               function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")))

#ind <- order(brca_wsum_table[1,-(1:2)], brca_wsum_table[2,-(1:2)], brca_wsum_table[3,-(1:2)])

ind <- c(1:3, 3 + order(brca_wsum_table[1, -c(1:3)], 
                        brca_wsum_table[2, -c(1:3)]))

write.table(brca_wsum_table[,ind], 
            file = paste0(path_plots, "/brca_wsum_best.csv"), 
            sep = ";", dec = ",", row.names = FALSE, col.names = FALSE,
            quote = FALSE)







