#source("load_config.R")

# Read BRCA DEG scores from files

brca_res_dimred_deg <- read.csv(paste0(path_intermediate_results, "/brca/dimred_deg/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_deg$Approach <- "DR DEG"
brca_res_dimred_deg_umap10n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_deg_umap10n$Approach <- "DR DEG"
brca_res_dimred_deg_umap10n$drname <- paste0(brca_res_dimred_deg_umap10n$drname, "_10n")
brca_res_dimred_deg_umap20n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_deg_umap20n$Approach <- "DR DEG"
brca_res_dimred_deg_umap20n$drname <- paste0(brca_res_dimred_deg_umap20n$drname, "_20n")
brca_res_dimred_deg_umap30n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_deg_umap30n$Approach <- "DR DEG"
brca_res_dimred_deg_umap30n$drname <- paste0(brca_res_dimred_deg_umap30n$drname, "_30n")

brca_scores_deg <- plyr::rbind.fill(brca_res_dimred_deg,
                                    brca_res_dimred_deg_umap10n,
                                    brca_res_dimred_deg_umap20n,
                                    brca_res_dimred_deg_umap30n)
brca_scores_deg <- brca_scores_deg[brca_scores_deg$fold != 6, ]

source("article/plot_renaming.R")

brca_scores_deg <- rename_methods(brca_scores_deg)
#brca_scores_deg$Approach <- factor(brca_scores_deg$Approach, unique(brca_scores_deg$Approach)[c(1,3,5,2,4,6)])
brca_scores_deg$cNMI <- (brca_scores_deg$NMI.BRCA_Subtype_PAM50 + 1 - brca_scores_deg$NMI.tss) / 2

summary_by <- c("datname", "drname", "k", "m", "Approach", "Method", "Embedding", "Clustering")

# Separate survival data for box-plots
brca_scores_deg_survival_sampled <- brca_scores_deg[c(summary_by, "SurvivalPValue")]

# Prepare averages
summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[i]] <- mean(x[[i]], na.rm = TRUE)
  }
  return(out)
}

brca_scores_deg_sampled <- brca_scores_deg

brca_scores_deg <- plyr::ddply(brca_scores_deg, summary_by, summary_function)
brca_scores_deg$fold <- NULL
brca_scores_deg$run <- NULL

# Read PRAD DEG scores from files

prad_res_dimred_deg <- read.csv(paste0(path_intermediate_results, "/prad/dimred_deg/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_deg$Approach <- "DR DEG"
prad_res_dimred_deg_umap10n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_deg_umap10n$Approach <- "DR DEG"
prad_res_dimred_deg_umap10n$drname <- paste0(prad_res_dimred_deg_umap10n$drname, "_10n")
prad_res_dimred_deg_umap20n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_deg_umap20n$Approach <- "DR DEG"
prad_res_dimred_deg_umap20n$drname <- paste0(prad_res_dimred_deg_umap20n$drname, "_20n")
prad_res_dimred_deg_umap30n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_deg_umap30n$Approach <- "DR DEG"
prad_res_dimred_deg_umap30n$drname <- paste0(prad_res_dimred_deg_umap30n$drname, "_30n")

prad_scores_deg <- plyr::rbind.fill(prad_res_dimred_deg,
                                    prad_res_dimred_deg_umap10n,
                                    prad_res_dimred_deg_umap20n,
                                    prad_res_dimred_deg_umap30n)
prad_scores_deg <- prad_scores_deg[prad_scores_deg$fold != 6, ]

prad_scores_deg <- rename_methods(prad_scores_deg)
#prad_scores_deg$Approach <- factor(prad_scores_deg$Approach, unique(prad_scores_deg$Approach)[c(1,3,5,2,4,6)])
prad_scores_deg$cNMI <- (prad_scores_deg$NMI.Gleason_category + 1 - prad_scores_deg$NMI.tss) / 2

summary_by <- c("datname", "drname", "k", "m", "Approach", "Method", "Embedding", "Clustering")

# Separate survival data for box-plots
prad_scores_deg_survival_sampled <- prad_scores_deg[c(summary_by, "SurvivalPValue")]

# Prepare averages
summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[i]] <- mean(x[[i]], na.rm = TRUE)
  }
  return(out)
}

prad_scores_deg_sampled <- prad_scores_deg

prad_scores_deg <- plyr::ddply(prad_scores_deg, summary_by, summary_function)
prad_scores_deg$fold <- NULL
prad_scores_deg$run <- NULL
