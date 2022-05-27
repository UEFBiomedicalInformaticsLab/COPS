#' Clustering, internal evaluation and batch effect estimation
#'
#' Single embedding or dataset evaluation
#'
#' Supported clustering methods are:
#' \itemize{
#' \item "hierarchical" - agglomerative hierarchical clustering
#' \item "diana" - divisive hierarchical clustering analysis
#' \item "kmeans" - k-means++
#' \item "model" - Gaussian Mixture Models
#' \item "knn_communities" - Louvain community detection on shared k nearest neighbour graphs
#' \item "spectral" - spectral clustering
#' }
#'
#' @param dat A data matrix with samples on columns.
#' @param n_clusters A vector of integers, numbers of clusters to be generated.
#' @param cluster_methods A vector of clustering method names, see details for options.
#' @param clustering_dissimilarity A dissimilarity matrix used in some methods such as hierarchical clustering. Computed with \code{\link{clustering_dissimilarity_from_data}} if missing. 
#' @param distance_metric Either "euclidean" or "correlation".
#' @param correlation_method Method for \code{\link[stats]{cor}}.
#' @param hierarchical_linkage See \code{\link[flashClust]{flashClust}}.
#' @param kmeans_num_init See \code{\link[ClusterR]{KMeans_rcpp}}.
#' @param kmeans_max_iters See \code{\link[ClusterR]{KMeans_rcpp}}.
#' @param kmeans_tol See \code{\link[ClusterR]{KMeans_rcpp}}.
#' @param gmm_modelNames Sepcifies model type for \code{\link[mclust]{Mclust}}
#' @param gmm_shrinkage Shrinkage parameter for \code{\link[mclust]{priorControl}}. 
#' @param knn_neighbours number of nearest neighbours for community detection.
#' @param knn_jaccard computes shared neighbour weights with Jaccard ubdex if \code{TRUE}. 
#' @param ... extra arguments are ignored currently
#'
#' @return Returns a \code{list} containing clusters, metrics, and
#'         \code{\link[stats]{chisq.test}} p-values
#'         if batch_label was supplied
#' @export
#' @importFrom stats chisq.test cutree cor dist as.dist
#' @importFrom aricode NMI ARI
#' @importFrom flashClust flashClust
#' @importFrom ClusterR KMeans_rcpp
#' @import mclust
#' @importFrom Spectrum Spectrum
#' @importFrom cluster silhouette
clustering_analysis <- function(dat,
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
  temp <- dat[,grepl("^dim[0-9]+$", colnames(dat))]
  temp <- temp[,sapply(temp, function(x) all(!is.na(x)))]
  rownames(temp) <- dat$id
  
  # Create dissimilarity matrix for Silhouette computation and HC
  if (!is.null(clustering_dissimilarity)) {
    diss <- clustering_dissimilarity
  } else {
    diss <- clustering_dissimilarity_from_data(temp, distance_metric, correlation_method)
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

#' Dissimilarity matrix
#'
#' @param x 
#' @param distance_metric 
#' @param correlation_method 
#' @param ... 
#'
#' @return
#' @export
clustering_dissimilarity_from_data <- function(x, 
                                               distance_metric = "euclidean", 
                                               correlation_method = "spearman", 
                                               preprocess = FALSE, 
                                               ...) {
  if (preprocess) {
    temp <- x[,grepl("^dim[0-9]+$", colnames(x))]
    temp <- temp[,sapply(temp, function(x) all(!is.na(x)))]
    rownames(temp) <- x$id
    x <- temp
  }
  if (distance_metric == "euclidean") {
    diss <- dist(x)
  } else if(distance_metric == "correlation") {
    diss <- as.dist(0.5 - cor(t(x), method = correlation_method)/2)
  } else {
    stop(paste("Unsupported distance metric:", distance_metric))
  }
  return(diss)
}

#' Clustering, internal evaluation and batch effect estimation
#'
#' Single embedding or dataset evaluation
#'
#' Supported clustering methods are:
#' \itemize{
#' \item "hierarchical" - agglomerative hierarchical clustering
#' \item "diana" - divisive hierarchical clustering analysis
#' \item "kmeans" - k-means++
#' \item "model" - Gaussian Mixture Models
#' \item "knn_communities" - Louvain community detection on shared k nearest neighbour graphs
#' \item "spectral" - spectral clustering
#' }
#'
#' @param x A data.frame of clustering results. 
#' @param dat A data.frame with data columns identified with "dim". Not required if clustering_dissimilarity is defined. 
#' @param by vector of variable names to split by
#' @param clustering_dissimilarity a \code{dist} object to use in silhouette calculation, defaults to euclidean distance matrix if left \code{NULL}
#' @param cluster_size_table return cluster sizes if \code{TRUE}.
#' @param silhouette_min_cluster_size proportional cluster size threshold for merging into nearest neighbours' cluster for silhouette computation.
#' @param ... extra arguments are passed to \code{\link{clustering_dissimilarity_from_data}}
#'
#' @return Returns a \code{list} containing metrics and cluster sizes
#' @export
#' @importFrom cluster silhouette
clustering_metrics <- function(x, 
                               dat = NULL, 
                               by = c("k", "m"), # must have k and m # TODO: fix?
                               clustering_dissimilarity = NULL, 
                               cluster_size_table = TRUE, 
                               silhouette_min_cluster_size = 0.0,
                               ...) {
  x <- as.data.frame(x)
  by <- by[by %in% colnames(x)]
  # Create dissimilarity matrix for Silhouette computation and HC
  if (!is.null(clustering_dissimilarity)) {
    diss <- clustering_dissimilarity
  } else {
    if (is.null(dat)) {
      stop("Cannot evaluate metrics, both dissimilarity and data are missing.")
    }
    temp <- dat[,grepl("^dim[0-9]+$", colnames(dat))]#, with = FALSE]
    temp <- temp[,apply(!is.na(temp), 2, all)]#, with = FALSE]
    rownames(temp) <- dat$id
    diss <- clustering_dissimilarity_from_data(temp, ...)
  }
  
  metrics <- data.frame()
  csize <- list()
  clusters <- split_by_safe(x, by)
  for (i in 1:length(clusters)) {
    # Silhouette
    matched_ind <- match(clusters[[i]]$id, labels(diss))
    silh_i <- silhouette_adjusted(clusters[[i]]$cluster, 
                                  diss = as.matrix(diss)[matched_ind, matched_ind], 
                                  min_size = silhouette_min_cluster_size)
    metrics <- rbind(metrics, data.frame(clusters[[i]][1, by], 
                                         metric = "Silhouette",
                                         value = mean(silh_i[,"sil_width"])))
    if (cluster_size_table) {
      csize[[i]] <- data.frame(clusters[[i]][1, by], 
                               as.data.frame(t(as.matrix(table(clusters[[i]]$cluster)))))
    }
  }
  out_list <- list(metrics = metrics)
  if (length(csize) > 0) out_list$cluster_sizes <- Reduce(rbind_fill, csize)
  return(out_list)
}

#' Clustering analysis on cross-validated data sets
#'
#' Performs clustering analysis on each fold of an external cross validation.
#'
#' Produces clusterings using multiple methods and settings while computing internal validation metrics 
#' such as Connectivity, Dunn and Silhouette scores. Also computes chi-squared tests with respect to 
#' a batch label if one is provided. 
#'
#' @param dat_embedded list of \code{data.frame}s
#' @param parallel number of threads
#' @param by variables to split input data by
#' @param ... extra arguments are passed through to clustering_evaluation
#'
#' @return Returns a \code{list} of \code{data.frames} containing \code{\link{clustering_analysis}} and \code{\link{clustering_metrics}} outputs for every
#'         combination of CV run, CV fold, clustering method, number of clusters as well as all combinations of
#'         data sets and dimensionality reduction techniques found in the input \code{data.frame}.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table rbindlist
cv_clusteval <- function(dat_embedded, 
                         parallel = 1, 
                         by = c("run", "fold"),
                         ...) {
  temp_list <- list()
  for (i in 1:length(dat_embedded)) {
    temp <- dat_embedded[[i]]
    drname <- names(dat_embedded)[i]
    if (is.null(drname)) drname <- i
    temp$drname <- drname
    temp <- split_by_safe(temp, by)
    temp_list <- c(temp_list, temp)
  }
  # Binding function that concatenates relevant list components
  cfun <- function(...){
    bound_list <- list()
    bound_list$clusters <- rbindlist(lapply(list(...), function(x) x$clusters))
    bound_list$metrics <- rbindlist(lapply(list(...), function(x) x$metrics))
    #bound_list$chisq_pval <- rbindlist(lapply(list(...), function(x) x$chisq_pval))
    #bound_list$batch_association <- rbindlist(lapply(list(...), function(x) x$batch_association))
    #bound_list$subtype_association <- rbindlist(lapply(list(...), function(x) x$subtype_association))
    bound_list$cluster_sizes <- rbindlist(lapply(list(...), function(x) x$cluster_sizes), fill = TRUE)
    return(bound_list)
  }
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(foreach(temp = temp_list,
                          .combine = cfun,
                          .export = c("clustering_analysis", "clustering_metrics"),
                          .packages = c("reshape2", "mclust", "cluster", "flashClust", "ClusterR"),
                          .multicombine = TRUE,
                          .maxcombine = max(length(temp_list), 2)) %dopar% {
                            temp_diss <- clustering_dissimilarity_from_data(temp, ..., preprocess = TRUE)
                            temp_clust <- clustering_analysis(temp, clustering_dissimilarity = temp_diss, ...)
                            temp_metrics <- clustering_metrics(temp_clust, clustering_dissimilarity = temp_diss, ...)
                            res <- list(clusters = temp_clust, metrics = temp_metrics$metrics, cluster_sizes = temp_metrics$cluster_sizes)
                            res
                          }, finally = close_parallel_cluster(parallel_clust))
  return(out)
}
