# Renames default method names to something more understandable. 
rename_methods <- function(x) {
  # Factor for grouping lines in plots
  x$Method <- paste0(x$datname, "+", x$drname, "+", x$m)
  
  # Embedding
  x$Embedding <- gsub("\\+original$", "", paste0(x$datname, "+", x$drname))
  x$Embedding <- gsub("^expr\\+", "", x$Embedding)
  x$Embedding <- gsub("^GENE_RWR\\+", "", x$Embedding)
  x$Embedding <- gsub("^KEGG_.*", "KEGG", x$Embedding)
  x$Embedding <- gsub("^GO_.*", "GO", x$Embedding)
  x$Embedding <- gsub("^REACTOME_.*", "REACTOME", x$Embedding)
  x$Embedding <- gsub("^HALLMARK_.*", "Hallmark", x$Embedding)
  pw_ind <- match(c("GO", "KEGG", "REACTOME", "Hallmark"), unique(x$Embedding))
  pw_ind <- pw_ind[!is.na(pw_ind)]
  if (length(pw_ind) > 0) {
    x$Embedding <- factor(x$Embedding, c(unique(x$Embedding)[c(pw_ind[-4], (1:length(unique(x$Embedding)))[-pw_ind])], "Hallmark"))
  } else {
    x$Embedding <- factor(x$Embedding, unique(x$Embedding))
  }
  
  # Clustering method
  x$Clustering <- x$m
  x$Clustering <- gsub("model", "GMM", x$Clustering)
  x$Clustering <- gsub("kmeans", "k-means", x$Clustering)
  x$Clustering <- gsub("hierarchical", "HC", x$Clustering)
  x$Clustering <- gsub("_average$", " (average)", x$Clustering)
  x$Clustering <- gsub("_ward$", " (Ward)", x$Clustering)
  x$Clustering <- gsub("_complete$", " (complete)", x$Clustering)
  x$Clustering <- gsub("^diana$", "DIANA", x$Clustering)
  
  # Survival
  colnames(x)[colnames(x) == "cluster_significance"] <- "SurvivalPValue"
  if ("SurvivalPValue" %in% colnames(x)) {
    x$SurvivalPValue_score <- (0.25 * (x$SurvivalPValue < 0.05) + 0.25 * (x$SurvivalPValue < 0.01) + 0.5 * (x$SurvivalPValue < 0.001))
  }
  
  # Other
  x$k <- factor(x$k)
  colnames(x) <- gsub("TrainStability", "ClusteringStability", colnames(x))
  colnames(x) <- gsub("TestStability", "ProjectionClusteringStability", colnames(x))
  
  return(x)
}

rename_dr_methods <- function(x) {
  x$Embedding <- gsub("^pca.*", "PCA", x$Embedding)
  x$Embedding <- gsub("^tsne.*", "t-SNE", x$Embedding)
  x$Embedding <- gsub("^umap.*", "UMAP", x$Embedding)
  pw_ind <- match(c("GO", "KEGG", "REACTOME", "Hallmark"), unique(x$Embedding))
  pw_ind <- pw_ind[!is.na(pw_ind)]
  if (length(pw_ind) > 0) {
    x$Embedding <- factor(x$Embedding, c(unique(x$Embedding)[c(pw_ind[-4], (1:length(unique(x$Embedding)))[-pw_ind])], "Hallmark"))
  } else {
    x$Embedding <- factor(x$Embedding, unique(x$Embedding))
  }
  return(x)
}