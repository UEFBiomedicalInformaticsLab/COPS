# Load results, keep and rename relevant factors for plots
#source("article/plot_renaming.R")

source("brca/brca_load_results.R")
brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]
temp <- brca_scores[brca_scores$drname %in% c("umap5_20n", "tsne45", "pca2") | ! brca_scores$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6,]
#temp <- rename_methods(temp)
brca_scores <- rename_dr_methods(temp)
brca_scores$Approach <- factor(brca_scores$Approach, levels = unique(brca_scores$Approach)[c(2,3,5,1,4,6)])

brca_scores_survival_sampled <- brca_scores_survival_sampled[brca_scores_survival_sampled$Clustering != "HC_single",]
temp <- brca_scores_survival_sampled[brca_scores_survival_sampled$drname %in% c("umap5_20n", "tsne45", "pca2") | 
                                       ! brca_scores_survival_sampled$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6,]
#temp <- rename_methods(temp)
brca_scores_survival_sampled <- rename_dr_methods(temp)
brca_scores_survival_sampled$Approach <- factor(brca_scores_survival_sampled$Approach, 
                                                levels = unique(brca_scores_survival_sampled$Approach)[c(2,3,5,1,4,6)])

source("prad/prad_load_results.R")
prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]
temp <- prad_scores[prad_scores$drname %in% c("umap2_20n", "tsne45", "pca2") | ! prad_scores$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6,]
#temp <- rename_methods(temp)
prad_scores <- rename_dr_methods(temp)
prad_scores$Approach <- factor(prad_scores$Approach, levels = unique(prad_scores$Approach)[c(2,3,5,1,4,6)])

prad_scores_survival_sampled <- prad_scores_survival_sampled[prad_scores_survival_sampled$Clustering != "HC_single",]
temp <- prad_scores_survival_sampled[prad_scores_survival_sampled$drname %in% c("umap2_20n", "tsne45", "pca2") | 
                                       ! prad_scores_survival_sampled$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6,]
#temp <- rename_methods(temp)
prad_scores_survival_sampled <- rename_dr_methods(temp)
prad_scores_survival_sampled$Approach <- factor(prad_scores_survival_sampled$Approach, 
                                                levels = unique(prad_scores_survival_sampled$Approach)[c(2,3,5,1,4,6)])

# Pareto
brca_pareto_score_subset <- brca_scores

brca_pareto_score_subset <- brca_pareto_score_subset[brca_pareto_score_subset$ClusteringStabilityJaccard > 0.8,]

brca_pareto_score_subset$Approach_unmodified <- brca_pareto_score_subset$Approach
brca_pareto_score_subset$Approach <- paste0(brca_pareto_score_subset$Approach, " (", brca_pareto_score_subset$Embedding, ")")

pref <- rPref::high(Silhouette) * rPref::high(cNMI) * rPref::high(SurvivalPValue_score) * 
  rPref::high(Module_score) * rPref::high(ClusteringStabilityJaccard)
brca_pres <- rPref::psel(brca_pareto_score_subset, pref, top_level = 1)

prad_pareto_score_subset <- prad_scores

prad_pareto_score_subset <- prad_pareto_score_subset[prad_pareto_score_subset$ClusteringStabilityJaccard > 0.8,]

prad_pareto_score_subset$Approach_unmodified <- prad_pareto_score_subset$Approach
prad_pareto_score_subset$Approach <- paste0(prad_pareto_score_subset$Approach, " (", prad_pareto_score_subset$Embedding, ")")

pref <- rPref::high(Silhouette) * rPref::high(cNMI) * 
  rPref::high(SurvivalPValue_score) * 
  rPref::high(Module_score) * rPref::high(ClusteringStabilityJaccard)
prad_pres <- rPref::psel(prad_pareto_score_subset, pref, top_level = 1)


