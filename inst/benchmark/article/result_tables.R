# Load results, keep and rename relevant factors for plots
#source("article/plot_renaming.R")

source("brca/brca_load_results.R")
brca_scores_sampled <- brca_scores_sampled[brca_scores_sampled$Clustering != "HC_single",]
temp <- brca_scores_sampled[brca_scores_sampled$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                      ! brca_scores_sampled$Approach %in% c("DR*", "DR") | 
                      brca_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
brca_scores_sampled <- rename_dr_methods(temp)
brca_scores_sampled$Approach <- factor(brca_scores_sampled$Approach, levels = unique(brca_scores_sampled$Approach)[c(2,3,5,1,4,6)])

brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]
temp <- brca_scores[brca_scores$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                      ! brca_scores$Approach %in% c("DR*", "DR") | 
                      brca_scores$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
brca_scores <- rename_dr_methods(temp)
brca_scores$Approach <- factor(brca_scores$Approach, levels = unique(brca_scores$Approach)[c(2,3,5,1,4,6)])

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
prad_scores_sampled <- prad_scores_sampled[prad_scores_sampled$Clustering != "HC_single",]
temp <- prad_scores_sampled[prad_scores_sampled$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                      ! prad_scores_sampled$Approach %in% c("DR*", "DR") |
                      prad_scores_sampled$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
prad_scores_sampled <- rename_dr_methods(temp)
prad_scores_sampled$Approach <- factor(prad_scores_sampled$Approach, levels = unique(prad_scores_sampled$Approach)[c(2,3,5,1,4,6)])

prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]
temp <- prad_scores[prad_scores$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                      ! prad_scores$Approach %in% c("DR*", "DR") |
                      prad_scores$drname == "VAE",]
temp <- temp[temp$k %in% 2:6,]
#temp <- rename_methods(temp)
prad_scores <- rename_dr_methods(temp)
prad_scores$Approach <- factor(prad_scores$Approach, levels = unique(prad_scores$Approach)[c(2,3,5,1,4,6)])

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

for (metric_i in c("Silhouette", "ClusteringStabilityJaccard", "cNMI", "Module_score")) {
  # Create tables
  brca_table_i <- tidyr::pivot_wider(brca_table, 
                                     id_cols = c("Approach", "Embedding", "Clustering"),
                                     names_from = c("k"),
                                     #names_prefix = "k",
                                     names_sep = "_",
                                     values_from = paste0(metric_i, c("", "_sd")))
  brca_table_i$Approach <- factor(brca_table_i$Approach, unique(brca_table_i$Approach))
  brca_table_i <- dplyr::arrange(brca_table_i, Approach, Embedding, Clustering)
  
  temp <- as.data.frame(t(1:ncol(brca_table_i)))
  colnames(temp) <- colnames(brca_table_i)
  
  brca_table_i <- rbind(temp, brca_table_i)
  brca_table_i <- rbind(temp, brca_table_i)
  
  brca_table_i[1,] <- c("", "", "", sapply(strsplit(colnames(brca_table_i), "_")[-c(1:3)], 
                                              function(x) paste0("k=", x[[length(x)]])))
  brca_table_i[2,] <- c("", "", "", sapply(strsplit(colnames(brca_table_i), "_")[-c(1:3)], 
                                              function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")))
  
  #ind <- order(brca_table_i[1,-(1:2)], brca_table_i[2,-(1:2)], brca_table_i[3,-(1:2)])
  
  ind <- c(1:3, 3 + order(brca_table_i[1, -c(1:3)], 
                          brca_table_i[2, -c(1:3)]))
  
  writexl::write_xlsx(brca_table_i[,ind], path = paste0(path_plots, "/brca_", metric_i, ".xlsx"), 
              #row.names = FALSE, 
              col_names = FALSE)#, sheetName = "BRCA")
}

for (metric_i in c("Silhouette", "ClusteringStabilityJaccard", "cNMI", "Module_score")) {
  # Create tables
  prad_table_i <- tidyr::pivot_wider(prad_table, 
                                     id_cols = c("Approach", "Embedding", "Clustering"),
                                     names_from = c("k"),
                                     #names_prefix = "k",
                                     names_sep = "_",
                                     values_from = paste0(metric_i, c("", "_sd")))
  prad_table_i$Approach <- factor(prad_table_i$Approach, unique(prad_table_i$Approach))
  prad_table_i <- dplyr::arrange(prad_table_i, Approach, Embedding, Clustering)
  
  temp <- as.data.frame(t(1:ncol(prad_table_i)))
  colnames(temp) <- colnames(prad_table_i)
  
  prad_table_i <- rbind(temp, prad_table_i)
  prad_table_i <- rbind(temp, prad_table_i)
  
  prad_table_i[1,] <- c("", "", "", sapply(strsplit(colnames(prad_table_i), "_")[-c(1:3)], 
                                           function(x) paste0("k=", x[[length(x)]])))
  prad_table_i[2,] <- c("", "", "", sapply(strsplit(colnames(prad_table_i), "_")[-c(1:3)], 
                                           function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")))
  
  #ind <- order(prad_table_i[1,-(1:2)], prad_table_i[2,-(1:2)], prad_table_i[3,-(1:2)])
  
  ind <- c(1:3, 3 + order(prad_table_i[1, -c(1:3)], 
                          prad_table_i[2, -c(1:3)]))
  
  writexl::write_xlsx(prad_table_i[,ind], path = paste0(path_plots, "/prad_", metric_i, ".xlsx"), 
                      #row.names = FALSE, 
                      col_names = FALSE)#, sheetName = "PRAD")
}

# Best clustering result only 
brca_means_nk <- plyr::ddply(brca_table, 
                             summary_by[summary_by != "k"],
                             function(x) data.frame(nk_mean = mean(x$Unsupervised_wsum, na.rm = TRUE)))
brca_nk_best <- plyr::ddply(brca_means_nk, 
                            c("Approach", "Embedding"),
                            function(x) data.frame(Clustering = x$Clustering[which.max(x$nk_mean)]))
brca_table_best <- plyr::join(brca_nk_best, brca_table)


prad_means_nk <- plyr::ddply(prad_table, 
                             summary_by[summary_by != "k"],
                             function(x) data.frame(nk_mean = mean(x$Unsupervised_wsum, na.rm = TRUE)))
prad_nk_best <- plyr::ddply(prad_means_nk, 
                            c("Approach", "Embedding"),
                            function(x) data.frame(Clustering = x$Clustering[which.max(x$nk_mean)]))
prad_table_best <- plyr::join(prad_nk_best, prad_table)

for (metric_i in c("Silhouette", "ClusteringStabilityJaccard", "cNMI", "Module_score")) {
  # Create tables
  brca_table_best_i <- tidyr::pivot_wider(brca_table_best, 
                                     id_cols = c("Approach", "Embedding", "Clustering"),
                                     names_from = c("k"),
                                     #names_prefix = "k",
                                     names_sep = "_",
                                     values_from = paste0(metric_i, c("", "_sd")))
  brca_table_best_i$Approach <- factor(brca_table_best_i$Approach, unique(brca_table_best_i$Approach))
  brca_table_best_i <- dplyr::arrange(brca_table_best_i, Approach, Embedding, Clustering)
  
  temp <- as.data.frame(t(1:ncol(brca_table_best_i)))
  colnames(temp) <- colnames(brca_table_best_i)
  
  brca_table_best_i <- rbind(temp, brca_table_best_i)
  brca_table_best_i <- rbind(temp, brca_table_best_i)
  
  brca_table_best_i[1,] <- c("", "", "", sapply(strsplit(colnames(brca_table_best_i), "_")[-c(1:3)], 
                                           function(x) paste0("k=", x[[length(x)]])))
  brca_table_best_i[2,] <- c("", "", "", sapply(strsplit(colnames(brca_table_best_i), "_")[-c(1:3)], 
                                           function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")))
  
  #ind <- order(brca_table_best_i[1,-(1:2)], brca_table_best_i[2,-(1:2)], brca_table_best_i[3,-(1:2)])
  
  ind <- c(1:3, 3 + order(brca_table_best_i[1, -c(1:3)], 
                          brca_table_best_i[2, -c(1:3)]))
  
  writexl::write_xlsx(brca_table_best_i[,ind], path = paste0(path_plots, "/brca_", metric_i, "_best.xlsx"), 
                      #row.names = FALSE, 
                      col_names = FALSE)#, sheetName = "BRCA")
}

for (metric_i in c("Silhouette", "ClusteringStabilityJaccard", "cNMI", "Module_score")) {
  # Create tables
  prad_table_best_i <- tidyr::pivot_wider(prad_table_best, 
                                     id_cols = c("Approach", "Embedding", "Clustering"),
                                     names_from = c("k"),
                                     #names_prefix = "k",
                                     names_sep = "_",
                                     values_from = paste0(metric_i, c("", "_sd")))
  prad_table_best_i$Approach <- factor(prad_table_best_i$Approach, unique(prad_table_best_i$Approach))
  prad_table_best_i <- dplyr::arrange(prad_table_best_i, Approach, Embedding, Clustering)
  
  temp <- as.data.frame(t(1:ncol(prad_table_best_i)))
  colnames(temp) <- colnames(prad_table_best_i)
  
  prad_table_best_i <- rbind(temp, prad_table_best_i)
  prad_table_best_i <- rbind(temp, prad_table_best_i)
  
  prad_table_best_i[1,] <- c("", "", "", sapply(strsplit(colnames(prad_table_best_i), "_")[-c(1:3)], 
                                           function(x) paste0("k=", x[[length(x)]])))
  prad_table_best_i[2,] <- c("", "", "", sapply(strsplit(colnames(prad_table_best_i), "_")[-c(1:3)], 
                                           function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")))
  
  #ind <- order(prad_table_best_i[1,-(1:2)], prad_table_best_i[2,-(1:2)], prad_table_best_i[3,-(1:2)])
  
  ind <- c(1:3, 3 + order(prad_table_best_i[1, -c(1:3)], 
                          prad_table_best_i[2, -c(1:3)]))
  
  writexl::write_xlsx(prad_table_best_i[,ind], path = paste0(path_plots, "/prad_", metric_i, "_best.xlsx"), 
                      #row.names = FALSE, 
                      col_names = FALSE)#, sheetName = "PRAD")
}

## xtables
# BRCA silhouette
brca_table_best_silh <- tidyr::pivot_wider(brca_table_best, 
                                        id_cols = c("Approach", "Embedding", "Clustering"),
                                        names_from = c("k"),
                                        names_sep = "_",
                                        values_from = "Silhouette")

brca_table_best_silh$Approach <- factor(brca_table_best_silh$Approach, unique(brca_table_best_silh$Approach))
brca_table_best_silh <- dplyr::arrange(brca_table_best_silh, Approach, Embedding, Clustering)

temp <- as.data.frame(t(1:ncol(brca_table_best_silh)))
colnames(temp) <- colnames(brca_table_best_silh)


ind <- c(1:3, 3 + order(colnames(brca_table_best_silh)[-c(1:3)]))

print(xtable::xtable(brca_table_best_silh[,ind[-4]], digits = 2), include.rownames = FALSE)

# BRCA stability
brca_table_best_stab <- tidyr::pivot_wider(brca_table_best, 
                                           id_cols = c("Approach", "Embedding", "Clustering"),
                                           names_from = c("k"),
                                           names_sep = "_",
                                           values_from = "ClusteringStabilityJaccard")

brca_table_best_stab$Approach <- factor(brca_table_best_stab$Approach, unique(brca_table_best_stab$Approach))
brca_table_best_stab <- dplyr::arrange(brca_table_best_stab, Approach, Embedding, Clustering)

temp <- as.data.frame(t(1:ncol(brca_table_best_stab)))
colnames(temp) <- colnames(brca_table_best_stab)

#ind <- c(1:3, 3 + order(brca_table_best_stab[1, -c(1:3)], 
#                        brca_table_best_stab[2, -c(1:3)]))

ind <- c(1:3, 3 + order(colnames(brca_table_best_stab)[-c(1:3)]))

print(xtable::xtable(brca_table_best_stab[,ind[-4]], digits = 2), include.rownames = FALSE)

# PRAD silhouette
prad_table_best_silh <- tidyr::pivot_wider(prad_table_best, 
                                           id_cols = c("Approach", "Embedding", "Clustering"),
                                           names_from = c("k"),
                                           names_sep = "_",
                                           values_from = "Silhouette")

prad_table_best_silh$Approach <- factor(prad_table_best_silh$Approach, unique(prad_table_best_silh$Approach))
prad_table_best_silh <- dplyr::arrange(prad_table_best_silh, Approach, Embedding, Clustering)

temp <- as.data.frame(t(1:ncol(prad_table_best_silh)))
colnames(temp) <- colnames(prad_table_best_silh)


ind <- c(1:3, 3 + order(colnames(prad_table_best_silh)[-c(1:3)]))

print(xtable::xtable(prad_table_best_silh[,ind[-4]], digits = 2), include.rownames = FALSE)

# PRAD stability
prad_table_best_stab <- tidyr::pivot_wider(prad_table_best, 
                                           id_cols = c("Approach", "Embedding", "Clustering"),
                                           names_from = c("k"),
                                           names_sep = "_",
                                           values_from = "ClusteringStabilityJaccard")

prad_table_best_stab$Approach <- factor(prad_table_best_stab$Approach, unique(prad_table_best_stab$Approach))
prad_table_best_stab <- dplyr::arrange(prad_table_best_stab, Approach, Embedding, Clustering)

temp <- as.data.frame(t(1:ncol(prad_table_best_stab)))
colnames(temp) <- colnames(prad_table_best_stab)

#ind <- c(1:3, 3 + order(prad_table_best_stab[1, -c(1:3)], 
#                        prad_table_best_stab[2, -c(1:3)]))

ind <- c(1:3, 3 + order(colnames(prad_table_best_stab)[-c(1:3)]))

print(xtable::xtable(prad_table_best_stab[,ind[-4]], digits = 2), include.rownames = FALSE)


