## BRCA
# kNN community analysis results
brca_res_knn_comm <- read.csv(paste0(path_intermediate_results, "/brca/knn_communities/scores.csv"), header = TRUE, row.names = 1)
brca_res_knn_comm$Approach <- "Direct clustering"

brca_res_knn_comm$cNMI <- 0.5 + brca_res_knn_comm$NMI.BRCA_Subtype_PAM50 / 2 - brca_res_knn_comm$NMI.tss / 2
sapply(brca_res_knn_comm[brca_res_knn_comm$fold != 6,-which(colnames(brca_res_knn_comm) %in% summary_by)], mean)

# Survival
ggplot(brca_res_knn_comm[brca_res_knn_comm$fold != 6,], aes(x = k, y = SurvivalPValue, fill = m)) + geom_boxplot() + theme_bw() +
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8)) + 
  ggtitle("Survival p-value of direct clustering methods in BRCA")

# Check the number of communities
clust_knn_comm <-  read.csv(paste0(path_intermediate_results, "/brca/knn_communities/clusters.csv.gz"), header = TRUE, row.names = 1)
table(tapply(clust_knn_comm$cluster, clust_knn_comm[,c("run", "fold")], function(x) length(unique(x))))

# Spectrum 
brca_res_spectral <- read.csv(paste0(path_intermediate_results, "/brca/spectral/scores.csv"), header = TRUE, row.names = 1)
brca_res_spectral$Approach <- "Direct clustering"
brca_res_spectral$cNMI <- brca_res_spectral$wsum

ggplot(brca_res_spectral[brca_res_spectral$fold != 6,], aes(x = k, y = Silhouette, group = k, fill = m)) + geom_boxplot()
ggplot(brca_res_spectral[brca_res_spectral$fold != 6,], aes(x = k, y = TrainStabilityJaccard, group = k, fill = m)) + geom_boxplot()
ggplot(brca_res_spectral[brca_res_spectral$fold != 6,], aes(x = k, y = cNMI, group = k, fill = m)) + geom_boxplot()
ggplot(brca_res_spectral[brca_res_spectral$fold != 6,], aes(x = k, y = Module_score, group = k, fill = m)) + geom_boxplot()
ggplot(brca_res_spectral[brca_res_spectral$fold != 6,], aes(x = k, y = SurvivalPValue, group = k, fill = m)) + geom_boxplot() + theme_bw() +
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8)) + 
  ggtitle("Survival p-value of direct clustering methods in BRCA")

## PRAD
# kNN community analysis results
prad_res_knn_comm <- read.csv(paste0(path_intermediate_results, "/prad/knn_communities/scores.csv"), header = TRUE, row.names = 1)
prad_res_knn_comm$Approach <- "Direct clustering"

prad_res_knn_comm$cNMI <- 0.5 + prad_res_knn_comm$NMI.Gleason_category / 2 - prad_res_knn_comm$NMI.tss / 2
sapply(prad_res_knn_comm[prad_res_knn_comm$fold != 6,-which(colnames(prad_res_knn_comm) %in% summary_by)], mean)

# Survival
ggplot(prad_res_knn_comm[prad_res_knn_comm$fold != 6,], aes(x = k, y = SurvivalPValue, fill = m)) + geom_boxplot() + theme_bw() +
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8)) + 
  ggtitle("Survival p-value of direct clustering methods in PRAD")

# Check the number of communities
clust_knn_comm <-  read.csv(paste0(path_intermediate_results, "/prad/knn_communities/clusters.csv.gz"), header = TRUE, row.names = 1)
table(tapply(clust_knn_comm$cluster, clust_knn_comm[,c("run", "fold")], function(x) length(unique(x))))

