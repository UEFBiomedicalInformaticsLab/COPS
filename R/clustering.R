clustering_only <- function(dat,
                            n_clusters = 2:5,
                            cluster_methods = c("hierarchical","diana","kmeans"),
                            clustering_dissimilarity = NULL, 
                            distance_metric = "euclidean",
                            correlation_method = "spearman",
                            hierarchical_linkage = "complete",
                            kmeans_num_init = 100,
                            kmeans_max_iters = 100,
                            kmeans_tol = 0.0001,
                            gmm_modelNames = NULL,  
                            gmm_shrinkage = 0.01, 
                            knn_neighbours = 30, 
                            knn_jaccard = TRUE, 
                            ...) {
  temp <- dat[grepl("^dim[0-9]+$", colnames(dat))]
  temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
  rownames(temp) <- dat$id
  
  # Create dissimilarity matrix for Silhouette computation and HC
  if (!is.null(clustering_dissimilarity)) {
    diss <- clustering_dissimilarity
  } else {
    clustering_dissimilarity(temp, distance_metric, correlation_method)
  }
  
  # Prepare case with multiple linkage methods
  cluster_methods_expanded <- cluster_methods
  multiple_linkages <- "hierarchical" %in% cluster_methods & length(hierarchical_linkage) > 1
  if (multiple_linkages) {
    cluster_methods_expanded <- c(cluster_methods[cluster_methods != "hierarchical"], 
                                  paste0("hierarchical_", hierarchical_linkage))
  }
  
  clusters <- data.frame()
  
  for (i in 1:length(cluster_methods_expanded)) {
    if (grepl("^hierarchical", cluster_methods_expanded[i])) {
      if (!multiple_linkages) {
        if (hierarchical_linkage == "ward") {
          clust_i <- flashClust::flashClust(diss**2, method = hierarchical_linkage)
        } else {
          clust_i <- flashClust::flashClust(diss, method = hierarchical_linkage)
        }
      } else {
        linkage_i <- gsub("^hierarchical_", "", cluster_methods_expanded[i])
        if (linkage_i == "ward") {
          if (distance_metric != "euclidean") stop("Ward linkage should be used with euclidean distance")
          clust_i <- flashClust::flashClust(diss**2, method = linkage_i)
        } else {
          clust_i <- flashClust::flashClust(diss, method = linkage_i)
        }
      }
      for (j in 1:length(n_clusters)) {
        temp_k <- cutree(clust_i, n_clusters[j])
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = temp_k))
      }
    } else if (cluster_methods_expanded[i] == "diana") {
      clust_i <- cluster::diana(x = diss, diss = TRUE)
      for (j in 1:length(n_clusters)) {
        temp_k <- cutree(clust_i, n_clusters[j])
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = temp_k))
      }
      #} else if (cluster_methods_expanded[i] == "agnes") {
      #  
      #} else if (cluster_methods_expanded[i] == "pam") {
      #  
    } else if (cluster_methods_expanded[i] == "kmeans") {
      if (distance_metric != "euclidean") stop("Only euclidean distance is supported by k-means")
      for (j in 1:length(n_clusters)) {
        clust_k <- ClusterR::KMeans_rcpp(data = temp, clusters = n_clusters[j], 
                                         num_init = kmeans_num_init, max_iters = kmeans_max_iters,
                                         tol = kmeans_tol)
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = clust_k$clusters))
      }
      #} else if (cluster_methods_expanded[i] == "sota") {
      #  
    } else if (cluster_methods_expanded[i] == "model") {
      if (distance_metric != "euclidean") stop("Only euclidean distance is supported by GMM")
      for (j in 1:length(n_clusters)) {
        clust_k <- mclust::Mclust(data = temp, G = n_clusters[j], modelNames = gmm_modelNames, 
                                  prior = priorControl(functionName="defaultPrior", shrinkage = gmm_shrinkage), 
                                  verbose = FALSE)
        if (is.null(clust_k)) {
          # Fail
          warning(paste0("GMM fitting failed (model: ", gmm_modelNames, ", samples: ", 
                         dim(temp)[1], ", dimensions = ", dim(temp)[2], ", clusters: ",
                         n_clusters[j], ")"))
          clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                                 k = n_clusters[j], cluster = NA))
        } else {
          # Success
          clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                                 k = n_clusters[j], cluster = clust_k$classification))
        }
      }
    } else if (cluster_methods_expanded[i] == "knn_communities") {
      clust_k <- knn_communities(t(temp), k = knn_neighbours, jaccard_kernel = knn_jaccard)
      clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                             k = "variable", cluster = clust_k))
    } else if (cluster_methods_expanded[i] == "spectral") {
      for (j in 1:length(n_clusters)) {
        clust_k <- Spectrum::Spectrum(t(temp), method = 3, fixk = n_clusters[j])
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = clust_k$assignments))
      }
    } else if (cluster_methods_expanded[i] == "SC3") {
      # Multi-thread errors that seemingly cannot be avoided
      for (j in 1:length(n_clusters)) {
        # SC3 only accepts input in the form of SingleCellExperiment 
        hack <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(temp)))
        SummarizedExperiment::rowData(hack)$feature_symbol <- colnames(temp)
        hack <- SC3::sc3(hack, ks = n_clusters[j], gene_filter = FALSE, n_cores = NULL)
        clust_k <- cutree(hack@metadata$sc3$consensus[[1]]$hc, n_clusters[j])
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = clust_k))
      }
    }else {
      stop(paste("Unsupported method:", cluster_methods_expanded[i]))
    }
  }
  
  # Combine outputs with metadata
  out <- plyr::join(dat[!grepl("^dim[0-9]+$", colnames(dat))], 
                                  clusters, 
                                  by = "id")
  out <- out[!is.na(out$cluster),] # potential issue with reference fold missing later
  return(out)
}

temp_removed <- function(){
  out_list$metrics <- metrics
  
  # Insert meta data (if present)
  out_list$metrics$run <- dat$run[1]
  out_list$metrics$fold <- dat$fold[1]
  out_list$metrics$datname <- dat$datname[1]
  out_list$metrics$drname <- dat$drname[1]
  
  # Pearson's chi-squared test
  f1 <- function(x, c1, c2) {
    temp <- tryCatch(suppressWarnings(chisq.test(x[[c1]], x[[c2]])), error = function(e) NULL)
    if (is.null(temp)) {
      temp <- data.frame(p = NA)
    } else {
      temp <- data.frame(p = temp$p.value)
    }
    temp$run <- x$run[1]
    temp$fold <- x$fold[1]
    temp$datname <- x$datname[1]
    temp$drname <- x$drname[1]
    temp$k <- x$k[1]
    temp$m <- x$m[1]
    temp$label <- c2
    return(temp)
  }
  # Other contingency table based metrics
  f2 <- function(x, c1, c2) {
    temp <- data.frame(nmi = aricode::NMI(x[[c1]], x[[c2]]),
                       ari = aricode::ARI(x[[c1]], x[[c2]]))
    temp$run <- x$run[1]
    temp$fold <- x$fold[1]
    temp$datname <- x$datname[1]
    temp$drname <- x$drname[1]
    temp$k <- x$k[1]
    temp$m <- x$m[1]
    temp$label <- c2
    return(temp)
  }
  
  if (!is.null(batch_label_names)) {
    batch_label_chisq <- list()
    batch_label_assoc <- list()
    for (i in 1:length(batch_label_names)) {
      batch_label_chisq[[i]] <- Reduce("rbind", 
                                       lapply(split(out_list$clusters, 
                                                    out_list$clusters[c("k", "m")]), 
                                              f1, c1 = "cluster", c2 = batch_label_names[i]))
      batch_label_assoc[[i]] <- Reduce("rbind", 
                                       lapply(split(out_list$clusters, 
                                                    out_list$clusters[c("k", "m")]), 
                                              f2, c1 = "cluster", c2 = batch_label_names[i]))
    }
    out_list$chisq_pval <- Reduce("rbind", batch_label_chisq)
    out_list$batch_association <- Reduce("rbind", batch_label_assoc)
  }
  
  if (!is.null(subtype_label_names)) {
    subtype_label_assoc <- list()
    for (i in 1:length(subtype_label_names)) {
      temp <- out_list$clusters[!is.na(out_list$clusters[[subtype_label_names[i]]]),]
      subtype_label_assoc[[i]] <- Reduce("rbind", 
                                         lapply(split(temp, 
                                                      temp[c("k", "m")]), 
                                                f2, c1 = "cluster", c2 = subtype_label_names[i]))
    }
    out_list$subtype_association <- Reduce("rbind", subtype_label_assoc)
  }
  
  f3 <- function(x, c1) {
    temp <- as.data.frame(t(as.matrix(table(x[[c1]]))))
    temp$run <- x$run[1]
    temp$fold <- x$fold[1]
    temp$datname <- x$datname[1]
    temp$drname <- x$drname[1]
    temp$k <- x$k[1]
    temp$m <- x$m[1]
    return(temp)
  }
  
  if (cluster_size_table) {
    out_list$cluster_sizes <-  Reduce(rbind_fill, 
                                      lapply(split(out_list$clusters, 
                                                   out_list$clusters[c("k", "m")]), 
                                             f3, c1 = "cluster"))
  }
  
  return(out_list)
}

clustering_dissimilarity <- function(x, 
                                     distance_metric = "euclidean", 
                                     correlation_method = "spearman", 
                                     ...) {
  temp <- x[grepl("^dim[0-9]+$", colnames(x))]
  temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
  rownames(temp) <- x$id
  if (distance_metric == "euclidean") {
    diss <- dist(x)
  } else if(distance_metric == "correlation") {
    diss <- as.dist(0.5 - cor(t(x), method = correlation_method)/2)
  } else {
    stop(paste("Unsupported distance metric:", distance_metric))
  }
  return(diss)
}

clustering_metrics <- function(x, 
                               dat = NULL, 
                               by = c("k", "m"),
                               clustering_dissimilarity = NULL, 
                               distance_metric = "euclidean",
                               correlation_method = "spearman",
                               cluster_size_table = TRUE, 
                               silhouette_min_cluster_size = 0.0,
                               ...) {
  # Create dissimilarity matrix for Silhouette computation and HC
  if (!is.null(clustering_dissimilarity)) {
    diss <- clustering_dissimilarity
  } else {
    if (is.null(dat)) {
      stop("Cannot evaluate metrics, both dissimilarity and data are missing.")
    }
    temp <- dat[grepl("^dim[0-9]+$", colnames(dat))]
    temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
    rownames(temp) <- dat$id
    clustering_dissimilarity(temp, distance_metric, correlation_method)
  }
  
  metrics <- data.frame()
  csize <- list()
  clusters <- split(x, f = x[by])
  for (i in 1:length(clusters)) {
    # Silhouette
    silh_i <- silhouette_adjusted(clusters[[i]]$cluster, 
                                  diss = diss, 
                                  min_size = silhouette_min_cluster_size)
    metrics <- rbind(metrics, data.frame(m = clusters[[i]]$m[1], 
                                         k = clusters[[i]]$k[1],
                                         metric = "Silhouette",
                                         value = mean(silh_i[,"sil_width"])))
    if (cluster_size_table) {
      csize[[i]] <- as.data.frame(t(as.matrix(table(clusters[[i]]$cluster))),
                                  m = clusters[[i]]$m[1], 
                                  k = clusters[[i]]$k[1])
    }
  }
  out_list <- list(metrics = metrics)
  if (length(csize) > 0) out_list$cluster_sizes <- Reduce(rbind_fill, csize)
  return(out_list)
}


