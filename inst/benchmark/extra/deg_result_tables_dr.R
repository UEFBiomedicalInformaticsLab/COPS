# Load results, keep and rename relevant factors for plots
#source("article/plot_renaming.R")

source("extra/load_deg_results.R")
source("brca/brca_load_results.R")
brca_scores_sampled <- rbind(brca_scores_sampled, brca_scores_deg_sampled)
brca_scores <- rbind(brca_scores, brca_scores_deg)

brca_scores_sampled <- brca_scores_sampled[brca_scores_sampled$Clustering != "HC_single",]
temp <- brca_scores_sampled <- brca_scores_sampled[brca_scores_sampled$Approach %in% c("DR*", "DR", "DR DEG"),]
#temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
brca_scores_sampled <- temp
brca_scores_sampled$Embedding <- brca_scores_sampled$drname

brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]
temp <- brca_scores <- brca_scores[brca_scores$Approach %in% c("DR*", "DR", "DR DEG"),]
#temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
brca_scores <- temp
brca_scores$Embedding <- brca_scores$drname

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
temp <- prad_scores_sampled <- prad_scores_sampled[prad_scores_sampled$Approach %in% c("DR*", "DR", "DR DEG"),]
#temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
prad_scores_sampled <- temp
prad_scores_sampled$Embedding <- prad_scores_sampled$drname

prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]
temp <- prad_scores <- prad_scores[prad_scores$Approach %in% c("DR*", "DR", "DR DEG"),]
#temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
prad_scores <- temp
prad_scores$Embedding <- prad_scores$drname

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

# Selections
brca_means_nk <- plyr::ddply(brca_table, 
                             summary_by[summary_by != "k"],
                             function(x) data.frame(nk_mean = mean(x$Unsupervised_wsum, na.rm = TRUE)))
brca_nk_best <- plyr::ddply(brca_means_nk, 
                            c("Approach", "Embedding"),
                            function(x) data.frame(Clustering = x$Clustering[which.max(x$nk_mean)]))
#brca_table <- plyr::join(brca_nk_best, brca_table)

test <- plyr::join(brca_nk_best, brca_means_nk)
test$DR_method <- sapply(strsplit(test$Embedding, " "), function(x) x[[1]])
test <- plyr::ddply(test, 
                    c("Approach", "DR_method"),
                    function(x) data.frame(Embedding = x$Embedding[which.max(x$nk_mean)]))

brca_table <- plyr::join(test, brca_table)

# Create tables
brca_dr_wsum_table <- tidyr::pivot_wider(brca_table, 
                                      id_cols = c("Embedding", "Clustering"),
                                      names_from = c("Approach", "k"),
                                      #names_prefix = "k",
                                      names_sep = "_",
                                      values_from = c("Unsupervised_wsum", "Unsupervised_wsum_sd"))
#ind <- (2:(ncol(brca_dr_wsum_table) %/% 2)) * 2

brca_dr_wsum_table <- dplyr::arrange(brca_dr_wsum_table, Embedding, Clustering)

#write.table(brca_dr_wsum_table[,c(1,2, 2 + order(c(ind - 1, ind)))], 
#            file = paste0(path_plots, "/brca_silh.csv"), 
#            sep = ";", dec = ",", row.names = FALSE, quote = FALSE)

temp <- as.data.frame(t(1:ncol(brca_dr_wsum_table)))
colnames(temp) <- colnames(brca_dr_wsum_table)

brca_dr_wsum_table <- rbind(temp, brca_dr_wsum_table)
brca_dr_wsum_table <- rbind(temp, brca_dr_wsum_table)
brca_dr_wsum_table <- rbind(temp, brca_dr_wsum_table)

brca_dr_wsum_table[1,] <- c("", "", sapply(strsplit(colnames(brca_dr_wsum_table), "_")[-(1:2)], 
                                              function(x) x[[length(x) - 1]]))
brca_dr_wsum_table[2,] <- c("", "", sapply(strsplit(colnames(brca_dr_wsum_table), "_")[-(1:2)], 
                                              function(x) paste0("k=", x[[length(x)]])))
brca_dr_wsum_table[3,] <- c("", "", sapply(strsplit(colnames(brca_dr_wsum_table), "_")[-(1:2)], 
                                              function(x) ifelse(x[[length(x) - 2]] == "sd", "sd", "mean")))

#ind <- order(brca_dr_wsum_table[1,-(1:2)], brca_dr_wsum_table[2,-(1:2)], brca_dr_wsum_table[3,-(1:2)])

ind<- c(1, 2, 2 + order(brca_dr_wsum_table[2, -c(1:2)], 
                        brca_dr_wsum_table[1, -c(1:2)], 
                        brca_dr_wsum_table[3, -c(1:2)]))

write.table(brca_dr_wsum_table[,ind], 
            file = paste0(path_plots, "/brca_dr_wsum.csv"), 
            sep = ";", dec = ",", row.names = FALSE, col.names = FALSE,
            quote = FALSE)



# Create tables
brca_dr_wsum_table <- tidyr::pivot_wider(brca_table, 
                                         id_cols = c("Approach", "Embedding", "Clustering"),
                                         names_from = c("k"),
                                         names_sep = "_",
                                         values_from = c("Unsupervised_wsum", "Unsupervised_wsum_sd"))
#ind <- (2:(ncol(brca_dr_wsum_table) %/% 2)) * 2

brca_dr_wsum_table$DR_method <- sapply(strsplit(brca_dr_wsum_table$Embedding, " "), function(x) x[[1]])
brca_dr_wsum_table <- dplyr::arrange(brca_dr_wsum_table, DR_method, Approach)

#write.table(brca_dr_wsum_table[,c(1,2, 2 + order(c(ind - 1, ind)))], 
#            file = paste0(path_plots, "/brca_silh.csv"), 
#            sep = ";", dec = ",", row.names = FALSE, quote = FALSE)

temp <- as.data.frame(t(1:ncol(brca_dr_wsum_table)))
colnames(temp) <- colnames(brca_dr_wsum_table)

brca_dr_wsum_table <- rbind(temp, brca_dr_wsum_table)
brca_dr_wsum_table <- rbind(temp, brca_dr_wsum_table)

brca_dr_wsum_table[1,] <- c("", "", "", sapply(strsplit(colnames(brca_dr_wsum_table), "_")[-c(1:3, ncol(brca_dr_wsum_table))], 
                                           function(x) paste0("k=", x[[length(x)]])),
                            "")
brca_dr_wsum_table[2,] <- c("", "", "", sapply(strsplit(colnames(brca_dr_wsum_table), "_")[-c(1:3, ncol(brca_dr_wsum_table))], 
                                           function(x) ifelse(x[[length(x) - 1]] == "sd", "sd", "mean")),
                            "")

#ind <- order(brca_dr_wsum_table[1,-(1:2)], brca_dr_wsum_table[2,-(1:2)], brca_dr_wsum_table[3,-(1:2)])

ind <- c(1:3, 3 + order(brca_dr_wsum_table[1, -c(1:3, ncol(brca_dr_wsum_table))], 
                        brca_dr_wsum_table[2, -c(1:3, ncol(brca_dr_wsum_table))]))

write.table(brca_dr_wsum_table[,ind], 
            file = paste0(path_plots, "/brca_dr_wsum_best.csv"), 
            sep = ";", dec = ",", row.names = FALSE, col.names = FALSE,
            quote = FALSE)






