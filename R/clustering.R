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
#' \item "SC3" - consensus clustering http://dx.doi.org/10.1038/nmeth.4236, note that this requires SC3 installation which is not required by default
#' \item "kkmeans" - kernelized k-means initialized by a spectral approximation
#' \item "kkmeanspp" - kernelized k-means++ with random initializations
#' }
#'
#' @param dat A data.frame with features on columns labeled as "dim[0-9]+", must also contain "id" column.
#' @param n_clusters A vector of integers defining the number of clusters.
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
#' @param kernel kernel for kernel k-means, options: "linear", "gaussian", "rbf", "jaccard", "tanimoto"
#' @param kernel_gamma gamma for the Gaussian/RBF kernel, higher values correspond to more complicated boundaries
#' @param kernel_center center kernels if TRUE
#' @param kernel_normalize normalize kernels to L2 norm 1 if TRUE
#' @param kkmeans_algorithm See \code{\link{kernel_kmeans}} options. 
#' @param kkmeans_refine See \code{\link{kernel_kmeans}}. 
#' @param kkmeans_maxiter maximum number of iterations for kernel k-means
#' @param kkmeans_n_init number of random initializations for kernel k-means++
#' @param kkmeans_tol delta error convergence threshold for spectral clustering
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{list} containing clusters, metrics, and
#'         \code{\link[stats]{chisq.test}} p-values
#'         if batch_label was supplied
#' @export
#' @importFrom stats chisq.test cutree cor dist as.dist
#' @importFrom flashClust flashClust
#' @importFrom ClusterR KMeans_rcpp
#' @import mclust
#' @importFrom cluster silhouette
clustering_analysis <- function(
    dat,
    n_clusters = 2:5,
    cluster_methods = c("hierarchical","diana","kmeans"),
    clustering_dissimilarity = NULL, 
    distance_metric = "euclidean",
    correlation_method = "spearman",
    hierarchical_linkage = "complete",
    kmeans_num_init = 100,
    kmeans_max_iters = 100,
    kmeans_tol = 1e-8,
    gmm_modelNames = NULL,  
    gmm_shrinkage = 0.01, 
    knn_neighbours = 30, 
    knn_jaccard = TRUE, 
    kernel = "linear",
    kernel_gamma = 1,
    kernel_center = TRUE,
    kernel_normalize = TRUE,
    kkmeans_algorithm = "spectral_qr", 
    kkmeans_refine = TRUE, 
    kkmeans_maxiter = 100,
    kkmeans_n_init = 100,
    kkmeans_tol = 1e-8,
    ...
) {
  temp <- dat[,grepl("^dim[0-9]+$", colnames(dat))]
  temp <- temp[,sapply(temp, function(x) all(!is.na(x)))]
  rownames(temp) <- dat$id
  
  # Create dissimilarity matrix for Silhouette computation and HC
  if (!is.null(clustering_dissimilarity)) {
    diss <- clustering_dissimilarity
  } else {
    diss <- clustering_dissimilarity_from_data(
      temp, 
      distance_metric, 
      correlation_method)
  }
  
  # Prepare case with multiple linkage methods
  cluster_methods_expanded <- cluster_methods
  multiple_linkages <- "hierarchical" %in% cluster_methods & 
    length(hierarchical_linkage) > 1
  if (multiple_linkages) {
    cluster_methods_expanded <- c(
      cluster_methods[cluster_methods != "hierarchical"], 
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
          if (distance_metric != "euclidean") {
            stop("Ward linkage should be used with euclidean distance")
          }
          clust_i <- flashClust::flashClust(diss**2, method = linkage_i)
        } else {
          clust_i <- flashClust::flashClust(diss, method = linkage_i)
        }
      }
      for (j in 1:length(n_clusters)) {
        temp_k <- cutree(clust_i, n_clusters[j])
        clusters <- rbind(
          clusters, 
          data.frame(
            id = rownames(temp), 
            m = cluster_methods_expanded[i], 
            k = n_clusters[j], 
            cluster = temp_k)
        )
      }
    } else if (cluster_methods_expanded[i] == "diana") {
      clust_i <- cluster::diana(x = diss, diss = TRUE)
      for (j in 1:length(n_clusters)) {
        temp_k <- cutree(clust_i, n_clusters[j])
        clusters <- rbind(
          clusters, 
          data.frame(
            id = rownames(temp), 
            m = cluster_methods_expanded[i], 
            k = n_clusters[j], 
            cluster = temp_k)
          )
      }
      #} else if (cluster_methods_expanded[i] == "agnes") {
      #  
      #} else if (cluster_methods_expanded[i] == "pam") {
      #  
    } else if (cluster_methods_expanded[i] == "kmeans") {
      if (distance_metric != "euclidean") {
        stop("Only euclidean distance is supported by k-means")
      }
      for (j in 1:length(n_clusters)) {
        clust_k <- ClusterR::KMeans_rcpp(
          data = temp, 
          clusters = n_clusters[j], 
          num_init = kmeans_num_init, 
          max_iters = kmeans_max_iters,
          tol = kmeans_tol)
        clusters <- rbind(
          clusters, 
          data.frame(
            id = rownames(temp), 
            m = cluster_methods_expanded[i], 
            k = n_clusters[j], 
            cluster = clust_k$clusters))
      }
      #} else if (cluster_methods_expanded[i] == "sota") {
      #  
    } else if (cluster_methods_expanded[i] == "model") {
      if (distance_metric != "euclidean") stop("Only euclidean distance is supported by GMM")
      for (j in 1:length(n_clusters)) {
        clust_k <- mclust::Mclust(
          data = temp, 
          G = n_clusters[j], 
          modelNames = gmm_modelNames, 
          prior = priorControl(
            functionName="defaultPrior", 
            shrinkage = gmm_shrinkage), 
          verbose = FALSE)
        if (is.null(clust_k)) {
          # Fail
          warning(paste0(
            "GMM fitting failed (model: ", gmm_modelNames, ", samples: ", 
            dim(temp)[1], ", dimensions = ", dim(temp)[2], ", clusters: ",
            n_clusters[j], ")"))
          clusters <- rbind(
            clusters, 
            data.frame(
              id = rownames(temp), 
              m = cluster_methods_expanded[i], 
              k = n_clusters[j], 
              cluster = NA)
            )
        } else {
          # Success
          clusters <- rbind(
            clusters, 
            data.frame(
              id = rownames(temp), 
              m = cluster_methods_expanded[i], 
              k = n_clusters[j], 
              cluster = clust_k$classification)
            )
        }
      }
    } else if (cluster_methods_expanded[i] == "knn_communities") {
      clust_k <- knn_communities(
        t(temp), 
        k = knn_neighbours, 
        jaccard_kernel = knn_jaccard)
      clusters <- rbind(
        clusters, 
        data.frame(
          id = rownames(temp), 
          m = cluster_methods_expanded[i], 
          k = "variable", 
          cluster = clust_k)
        )
    } else if (cluster_methods_expanded[i] == "spectral") {
      if (!requireNamespace("Spectrum", quietly = TRUE)) {
        stop("Please install Spectrum to enable spectral clustering.")
      }
      for (j in 1:length(n_clusters)) {
        clust_k <- Spectrum::Spectrum(
          t(temp), 
          method = 3, 
          fixk = n_clusters[j])
        clusters <- rbind(
          clusters, 
          data.frame(
            id = rownames(temp), 
            m = cluster_methods_expanded[i], 
            k = n_clusters[j], 
            cluster = clust_k$assignments)
          )
      }
    } else if (cluster_methods_expanded[i] == "SC3") {
      sc3 <- requireNamespace("SC3", quietly = TRUE)
      sce <- requireNamespace("SingleCellExperiment", quietly = TRUE)
      sc <- requireNamespace("SummarizedExperiment", quietly = TRUE)
      if (!(sc3 & sce & sc)) {
        stop("To run the SC3 clustering method, please install SC3 and required packages.")
      }
      # Multi-thread errors that seemingly cannot be avoided
      for (j in 1:length(n_clusters)) {
        # SC3 only accepts input in the form of SingleCellExperiment 
        sc3in <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(temp)))
        SummarizedExperiment::rowData(sc3in)$feature_symbol <- colnames(temp)
        sc3in <- SC3::sc3(sc3in, ks = n_clusters[j], gene_filter = FALSE, n_cores = NULL)
        clust_k <- cutree(sc3in@metadata$sc3$consensus[[1]]$hc, n_clusters[j])
        clusters <- rbind(
          clusters, 
          data.frame(
            id = rownames(temp), 
            m = cluster_methods_expanded[i], 
            k = n_clusters[j], 
            cluster = clust_k)
          )
      }
    } else if (cluster_methods_expanded[i] == "kkmeans") {
      for (kernel_i in kernel) {
        if (kernel_i == "linear") {
          temp_kernel <- temp %*% t(temp)
          if (kernel_center) temp_kernel <- center_kernel(temp_kernel)
          if (kernel_normalize) temp_kernel <- normalize_kernel(temp_kernel)
          temp_kernel <- list(temp_kernel)
        } else if (kernel_i %in% c("gaussian", "rbf")) {
          temp_kernel <- lapply(
            kernel_gamma, 
            function(kg) exp(- kg * as.matrix(dist(temp))**2))
        } else if (kernel_i %in% c("jaccard", "tanimoto")) {
          temp_kernel <- jaccard_matrix(t(temp))
          temp_kernel[is.nan(temp_kernel)] <- 0
          diag(temp_kernel) <- 1
          if (kernels_center[i]) temp_kernel <- center_kernel(temp_kernel)
          if (kernels_normalize[i]) temp_kernel <- normalize_kernel(temp_kernel)
          temp_kernel <- list(temp_kernel)
        } else {
          stop(paste0("Kernel \"", kernel_i, "\" is not supported."))
        }
        for (k_i in 1:length(temp_kernel)) {
          if (kkmeans_algorithm %in% c("spectral", "spectral_qr")) {
            eigs <- eigen(temp_kernel[[k_i]], symmetric = TRUE)$vectors
          } else {
            eigs <- NULL
          }
          for (k in n_clusters) {
            temp_res <- kernel_kmeans(
              K = temp_kernel[[k_i]], 
              n_k = k, 
              algorithm = kkmeans_algorithm, 
              spectral_qr_refine = kkmeans_refine, 
              kernel_eigen_vectors = eigs, 
              max_iter = kkmeans_maxiter, 
              num_init = kkmeans_n_init, 
              tol = kkmeans_tol, 
              parallel = 1
            )
            method_string <- paste0(
              kernel_i, 
              ifelse(kernel_i %in% c("gaussian", "rbf"), 
                     paste0("_", kernel_gamma[k_i]), ""), 
              "_", cluster_methods_expanded[i])
            clusters <- rbind(
              clusters, 
              data.frame(
                id = rownames(temp), 
                m = method_string, 
                k = k, 
                cluster = temp_res$clusters)
            )
          }
        }
      } 
    } else {
      stop(paste("Unsupported method:", cluster_methods_expanded[i]))
    }
  }
  
  # Combine outputs with metadata
  out <- plyr::join(
    dat[!grepl("^dim[0-9]+$", colnames(dat))], 
    clusters, 
    by = "id"
  )
  out <- out[!is.na(out$cluster),] # potential issue with reference fold missing later
  return(out)
}

#' Dissimilarity matrix
#'
#' @param x input data \code{matrix}
#' @param distance_metric Either "euclidean" or "correlation"
#' @param correlation_method method for \code{\link[stats]{cor}}
#' @param preprocess set to \code{TRUE}if \code{x} is a \code{data.frame} that 
#'   is processed within the COPS pipeline with data columns named "dim1", "dim2" and so on
#' @param ... ignored
#'
#' @return \code{dist}-object
#' @export
clustering_dissimilarity_from_data <- function(
    x, 
    distance_metric = "euclidean", 
    correlation_method = "spearman", 
    preprocess = FALSE, 
    ...
) {
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

#' Internal evaluation of clustering results
#'
#' Single embedding or dataset evaluation
#'
#' @param x A data.frame of clustering results. 
#' @param dat A data.frame with data columns identified with "dim". Not required if dissimilarity is defined. 
#' @param by vector of variable names to split by
#' @param dissimilarity a \code{dist} object to use in silhouette calculation, defaults to euclidean distance matrix if left \code{NULL}
#' @param cluster_size_table return cluster sizes if \code{TRUE}.
#' @param silhouette_min_cluster_size proportional cluster size threshold for merging into nearest neighbours' cluster for silhouette computation.
#' @param internal_metrics Internal metric names passed to \code{\link[clusterCrit]{intCriteria}}. This will slow the pipeline down considerably. 
#' @param ... extra arguments are passed to \code{\link{clustering_dissimilarity_from_data}}
#'
#' @return Returns a \code{list} containing metrics and cluster sizes
#' @export
#' @importFrom cluster silhouette
clustering_metrics <- function(
    x, 
    dat = NULL, 
    by = c("k", "m"), # must have k and m # TODO: fix?
    dissimilarity = NULL, 
    cluster_size_table = TRUE, 
    silhouette_min_cluster_size = 0.0,
    internal_metrics = NULL, 
    ...
) {
  x <- as.data.frame(x)
  by <- by[by %in% colnames(x)]
  # If internal_metrics is not set, we use the given dissimilarity matrix
  if (is.null(internal_metrics)) {
    # Create dissimilarity matrix for Silhouette computation and HC
    if (!is.null(dissimilarity)) {
      diss <- dissimilarity
    } else {
      if (is.null(dat)) {
        stop("Cannot evaluate metrics, both dissimilarity and data are missing.")
      }
      if (length(class(dat)) == 1 && class(dat) == "list") {
        if (length(dat) > 1) warning(paste(
          "More than one input dataset was found for internal evaluation.", 
          "Using only the first one."))
        dat <- dat[[1]]
      }
      temp <- dat[,grepl("^dim[0-9]+$", colnames(dat))]#, with = FALSE]
      temp <- temp[,apply(!is.na(temp), 2, all)]#, with = FALSE]
      rownames(temp) <- dat$id
      diss <- clustering_dissimilarity_from_data(temp, ...)
    }
  }
  
  metrics <- data.frame()
  csize <- list()
  clusters <- split_by_safe(x, by)
  for (i in 1:length(clusters)) { # TODO: add foreach and prevent nested parallelization
    if (!is.null(internal_metrics)) {
      if (!requireNamespace("clusterCrit", quietly = TRUE)) {
        stop("Please install clusterCrit to enable advanced internal metrics.")
      }
      if ("data.frame" %in% class(dat)) dat <- list(dat) # make list of dataframes
      if (length(dat) > 1 & is.null(names(dat))) names(dat) <- paste0("input", 1:length(dat))
      if (length(dat) == 0) stop("No data for internal metric calculations.")
      for (j in 1:length(dat)) {
        temp <- as.data.frame(dat[[j]])
        rownames(temp) <- temp[["id"]]
        temp <- temp[,grepl("^dim[0-9]+$", colnames(temp))]
        #temp <- temp[,apply(!is.na(temp), 2, all)]
        temp <- as.matrix(temp)
        
        if(any(!clusters[[i]]$id %in% rownames(temp))) stop("ID mismatch in internal indices data.")
        if(any(is.na(temp))) stop("NAs in internal indices data.")
        
        ccrit <- clusterCrit::intCriteria(
          temp[clusters[[i]]$id,], 
          as.integer(clusters[[i]]$cluster), 
          internal_metrics)
        
        ccrit <- data.frame(
          clusters[[i]][1, by], 
          metric = paste0(names(dat)[j], names(ccrit)), 
          value = unlist(ccrit),
          row.names = 1:length(ccrit))
        metrics <- rbind(metrics, ccrit)
      }
    } else {
      # Silhouette using given distance matrix
      matched_ind <- match(clusters[[i]]$id, labels(diss))
      silh_i <- silhouette_adjusted(
        clusters[[i]]$cluster, 
        diss = as.matrix(diss)[matched_ind, matched_ind], 
        min_size = silhouette_min_cluster_size)
      metrics <- rbind(
        metrics, 
        data.frame(
          clusters[[i]][1, by], 
          metric = "Silhouette",
          value = mean(silh_i[,"sil_width"]))
      )
    }
    
    if (cluster_size_table) {
      csize[[i]] <- data.frame(
        clusters[[i]][1, by], 
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
#' @param silhouette_dissimilarity dissimilarity matrix used for silhouette evaluation
#' @param dat_list list of input data matrices used for calculating clustering indices
#' @param ... extra arguments are passed through to \code{\link{clustering_dissimilarity_from_data}}, \code{\link{clustering_analysis}} and \code{\link{clustering_metrics}}
#'
#' @return Returns a \code{list} of \code{data.frames} containing \code{\link{clustering_analysis}} and \code{\link{clustering_metrics}} outputs for every
#'         combination of CV run, CV fold, clustering method, number of clusters as well as all combinations of
#'         data sets and dimensionality reduction techniques found in the input \code{data.frame}.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table rbindlist
subsample_clustering_evaluation <- function(
    dat_embedded, 
    parallel = 1, 
    by = c("datname", "drname", "run", "fold"), 
    silhouette_dissimilarity = NULL, 
    dat_list = NULL, 
    ...
) {
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
    bound_list$clusters <- data.table::rbindlist(
      lapply(list(...), function(x) x$clusters))
    bound_list$metrics <- data.table::rbindlist(
      lapply(list(...), function(x) x$metrics))
    bound_list$cluster_sizes <- data.table::rbindlist(
      lapply(list(...), function(x) x$cluster_sizes), fill = TRUE)
    return(bound_list)
  }
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(
    foreach(
    temp = temp_list,
    .combine = cfun,
    .export = c("clustering_analysis", "clustering_metrics"),
    .packages = c("reshape2", "mclust", "cluster", "flashClust", "ClusterR"),
    .multicombine = TRUE,
    .maxcombine = max(length(temp_list), 2)) %dopar% 
    {
      temp_diss <- clustering_dissimilarity_from_data(
        temp, 
        ..., 
        preprocess = TRUE
      )
      temp_clust <- clustering_analysis(
        temp, 
        clustering_dissimilarity = temp_diss, 
        ...
      )
      if (is.null(silhouette_dissimilarity)) silhouette_dissimilarity <- temp_diss
      temp_metrics <- clustering_metrics(
        temp_clust, 
        dat = dat_list, 
        dissimilarity = silhouette_dissimilarity, 
        by = c(by, "k", "m"), 
        ...
      )
      res <- list(
        clusters = temp_clust, 
        metrics = temp_metrics$metrics, 
        cluster_sizes = temp_metrics$cluster_sizes
      )
      res
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(out)
}

