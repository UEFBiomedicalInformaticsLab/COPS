#' Dispersion Separability Criterion
#'
#' Used by TCGA Batch Effects Viewer \url{https://bioinformatics.mdanderson.org/public-software/tcga-batch-effects/}.
#' Based on \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.1128}.
#' Can be used to measure batch effect. 
#'
#' @param data_matrix numeric matrix, samples on columns
#' @param batch_label categorical variable, must be vector
#'
#' @return Returns the DSC of \code{data_matrix} with respect to \code{batch_label} as a scalar value
#' @export
DSC <- function(data_matrix, batch_label) {
  if (is.null(names(batch_label))) names(batch_label) <- colnames(data_matrix)
  n <- table(batch_label)
  p <- n / sum(n)
  
  # traces of covariance matrices
  Sw <- sapply(split(names(batch_label), batch_label), function(x) sum(apply(data_matrix[, x, drop = FALSE], 1, var)))
  
  Dw <- sqrt(sum(Sw * p[names(Sw)], na.rm = TRUE))
  
  mu <-  lapply(split(names(batch_label), batch_label), 
                function(x) apply(data_matrix[, x, drop = FALSE], 1, mean))
  
  # Global mean
  #M <- apply(data_matrix, 1, mean)
  # This is actually more efficient (not much though)
  M <- rep(0, nrow(data_matrix))
  for (i in names(mu)) {
    M <- M + p[[i]] * mu[[i]]
  }
  
  Sb <- c()
  for (i in names(mu)) {
    Sb[i] <- sum((mu[[i]] - M)**2)
  }
  Db <- sqrt(sum(Sb * p[names(Sb)]))
  
  return(Db/Dw)
}

#' Jaccard similarity coefficient between two categorical vectors
#' 
#' R implementation of clusteval::cluster_similarity to avoid crashing.
#' 
#' @param a categorical vector
#' @param b categorical vector
#' 
#' @return Scalar, Jaccard coefficient 
#' @export
JaccardSimCoef <- function(a,b) {
  na <- is.na(a) | is.na(b)
  a <- a[!na]
  b <- b[!na]
  # flattened matrix of observation pairs (only need upper triangular, but its not efficient in R?)
  A <- rep(a, length(a)) == rep(a, each = length(a))
  B <- rep(b, length(b)) == rep(b, each = length(b))
  # compare matrices (remove diagonal)
  s1 <- sum(A & B) - length(a)
  s2 <- sum(A | B) - length(a)
  # return (first remove diagonal from result)
  return(s1 / s2)
}



#' Average Jaccard dissimilarity coefficient between multiple clustering results
#'
#' Measures the stability of a given clustering method when combined with
#' cross-validation or bootstrap by estimating the repeatability of clustering results.
#'
#' @param clustering_list \code{list} of clustering vectors
#'
#' @return Returns a scalar value giving the mean jaccard distance between clustering vectors
#' @importFrom clusteval cluster_similarity
#' @export
jdist <- function(clustering_list) {
  n_runs <- length(clustering_list)
  jdist <- array(NA, dim = c(n_runs, n_runs))
  for (i in 1:(n_runs-1)) {
    for (j in (i+1):n_runs) {
      jdist[i,j] <- 1 - clusteval::cluster_similarity(clustering_list[[i]],
                                                      clustering_list[[j]],
                                                      similarity = "jaccard")
    }
  }
  mean_dist <- mean(jdist, na.rm = TRUE)
  return(mean_dist)
}

#' Average Jaccard dissimilarity coefficient between clusterings and references
#'
#' Measures the stability of a given clustering method when combined with
#' cross-validation or bootstrap by estimating the repeatability of clustering results.
#'
#' @param clustering_list \code{list} of clustering vectors
#' @param clustering_reference_list \code{list} of reference clustering vectors
#'
#' @return Returns a scalar value giving the mean jaccard distance between clustering vectors
#' @importFrom clusteval cluster_similarity
#' @export
jdist_ref <- function(clustering_list, clustering_reference_list) {
  if (length(clustering_list) != length(clustering_reference_list)) {
    stop("Number of given references does not match number of given inputs.")
  }
  if (length(clustering_list) > 0) {
    jdist <- c()
    for (i in 1:length(clustering_list)) {
      jdist[i] <- 1 - clusteval::cluster_similarity(clustering_list[[i]],
                                                    clustering_reference_list[[i]],
                                                    similarity = "jaccard")
    }
    mean_dist <- mean(jdist, na.rm = TRUE)
  } else {
    mean_dist <- NA
  }
  return(mean_dist)
}

#' Clustering stability evaluation
#'
#' Performs stability analysis on cross-validated clusterings using \code{\link{jdist}}.
#'
#' Default settings work with \code{\link{cv_clusteval}} output 'clusters'.
#'
#' @param clust clustering \code{data.frame} such as returned by \code{\link{cv_clusteval}}
#' @param by vector of column names to keep
# @param by2 vector of column names to split by, must contain "fold"
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{data.frame} where each row corresponds to clustering stability
#'         with respect to kept column variables
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table data.table is.data.table rbindlist setDT setDTthreads
#'
stability_eval <- function(clust,
                           by = c("datname", "drname", "run", "k", "m"),
                           #by2 = c("fold"),
                           ...)
{
  by2 = c("fold")
  # Function to be applied for each
  f1 <- function(clust, clustref) {
    if (length(clust) != length(clustref)) {
      stop("Number of given references does not match number of given inputs.")
    }
    if (length(clust) > 0) {
      jsc <- c()
      nmi <- c()
      ari <- c()
      for (i in 1:length(clust)) {
        jsc[i] <- clusteval::cluster_similarity(clust[[i]],
                                                clustref[[i]],
                                                similarity = "jaccard")
        nmi[i] <- aricode::NMI(clust[[i]], clustref[[i]])
        ari[i] <- aricode::ARI(clust[[i]], clustref[[i]])
      }
      #mean_jsc <- mean(jsc, na.rm = TRUE)
      #mean_nmi <- mean(nmi, na.rm = TRUE)
      #mean_ari <- mean(ari, na.rm = TRUE)
    } else {
      #mean_jsc <- NA
      #mean_nmi <- NA
      #mean_ari <- NA
      jsc <- NA
      nmi <- NA
      ari <- NA
    }
    #return(list(mean_jsc = mean_jsc, mean_nmi = mean_nmi, mean_ari = mean_ari))
    return(list(jsc = jsc, nmi = nmi, ari = ari))
  }
  # Function to be applied for method combinations in clustering table
  f2 <- function(x) {
    data.table::setDTthreads(1)
    ref_i <- unique(x$fold)
    ref_i <- ref_i[!ref_i %in% unique(x$cv_index)]
    
    ref <- x[fold == ref_i,]
    colnames(ref)[colnames(ref) == "cluster"] <- "reference_cluster"
    ref$fold <- NULL
    
    nonref <- x[fold != ref_i,]
    nonref$test_ind <- nonref$cv_index == nonref$fold
    
    ref_cols <- c("id", by2[by2 != "fold"], "reference_cluster")
    nonref <- merge(nonref, ref[, ..ref_cols], by = c("id", by2[by2 != "fold"]))
    
    train_nonref <- split(nonref$cluster[!nonref$test_ind], nonref[!nonref$test_ind, ..by2])
    train_ref <- split(nonref$reference_cluster[!nonref$test_ind], nonref[!nonref$test_ind, ..by2])
    train_res <- f1(train_nonref, train_ref)
    
    test_ref <- split(nonref$cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
    test_nonref <- split(nonref$reference_cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
    test_res <- f1(test_nonref, test_ref)
    
    out_f2 <- data.table::data.table(fold = names(train_ref), 
                                     train_jsc = train_res$jsc,
                                     train_nmi = train_res$nmi,
                                     train_ari = train_res$ari,
                                     test_jsc = test_res$jsc,
                                     test_nmi = test_res$nmi,
                                     test_ari = test_res$ari)
    
    #out_f2 <- data.table::data.table(train_jsc = train_res$mean_jsc, 
    #                                 train_nmi = train_res$mean_nmi, 
    #                                 train_ari = train_res$mean_ari, 
    #                                 test_jsc = test_res$mean_jsc, 
    #                                 test_nmi = test_res$mean_nmi, 
    #                                 test_ari = test_res$mean_ari)
    return(out_f2)
  }
  by <- by[by %in% colnames(clust)]
  
  if (!is.data.table(clust)) data.table::setDT(clust)
  
  temp_list <- split(clust, by = by)
  #temp_list <- split(clust, clust[by])
  stability <- foreach(temp = temp_list,
                      .combine = function(...) data.table::rbindlist(list(...)),
                      .export = c(),
                      .packages = c("clusteval", "data.table", "aricode"),
                      .multicombine = TRUE,
                      .maxcombine = length(temp_list)) %dopar% {
    out <- f2(temp)
    for (j in by) {
      out[[j]] <- temp[[j]][1]
    }
    out
  }
  return(as.data.frame(stability))
}

#' Continuous variable to categorical variable association analysis
#'
#' Three non-clustering batch effect estimators put together. \cr
#' Can be used for other purposes as well. 
#'
#' @param dat data matrix, samples on columns
#' @param class vector or matrix where columns are categorical variables
#' @param n_pc_max maximum number of principal components to analyze
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{list} containing \code{\link[kBET]{batch_sil}},
#'         \code{\link[kBET]{pcRegression}} and \code{\link{DSC}}
#'         computed for each column of \code{class}
#' @export
#' @importFrom FactoMineR PCA
#' @importFrom kBET batch_sil
#' @importFrom kBET pcRegression
class_associations <-  function(dat, class, n_pc_max = 10, ...){
  out <- list()
  if (is.null(dim(class))) class <- cbind(as.character(class), c())

  pca_silh <- list()
  pca_reg <- list()
  #DSC_res <- c()
  for (i in 1:ncol(class)) {
    # PCA based
    dat_pca <- FactoMineR::PCA(t(dat),
                               scale.unit = FALSE,
                               ncp = min(n_pc_max, nrow(dat)),
                               graph = FALSE)
    pca_silh[[i]] <- kBET::batch_sil(list(x = dat_pca$ind$coord),
                                     class[,i],
                                     nPCs = min(n_pc_max, nrow(dat)))
    pca_reg[[i]] <- suppressWarnings(kBET::pcRegression(list(x = dat_pca$ind$coord),
                                       class[,i],
                                       n_top = min(n_pc_max, nrow(dat))))

    # Other
    #DSC_res[i] <- DSC(dat, class[,i])
  }
  if (!is.null(dimnames(class))) {
    names(pca_silh) <- dimnames(class)[[2]]
    names(pca_reg) <- dimnames(class)[[2]]
    names(DSC_res) <- dimnames(class)[[2]]
  }
  out <- list(PCA_silhouette = pca_silh,
              PCA_regression = pca_reg)#,
              #DSC = DSC_res)
  return(out)
}


#' Clustering, internal evaluation and batch effect estimation
#'
#' Single embedding or dataset evaluation
#'
#' Supported clustering methods are:
#' \itemize{
#' \item "hierarchical" -
#' \item "diana" -
#' \item "kmeans" -
#' \item "model" -
#' }
#'
#' @param dat A data matrix with samples on columns.
#' @param batch_label_names A character vector containing column names corresponding to batch labels.
#' @param n_clusters A vector of integers, numbers of clusters to be generated.
#' @param cluster_methods A vector of clustering method names, see details for options.
#' @param subtype_label_names A character vector containing column names corresponding to batch labels.
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
#' @importFrom mclust Mclust
#' @importFrom Spectrum Spectrum
#' @importFrom cluster silhouette
clustering_evaluation <- function(dat,
                                  batch_label_names = NULL,
                                  subtype_label_names = NULL,
                                  n_clusters = 2:5,
                                  cluster_methods = c("hierarchical","diana","kmeans"),
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
  if (distance_metric == "euclidean") {
    diss <- dist(temp)
  } else if(distance_metric == "correlation") {
    diss <- as.dist(0.5 - cor(t(temp), method = correlation_method)/2)
  } else {
    stop(paste("Unsupported distance metric:", distance_metric))
  }
  
  # Prepare case with multiple linkage methods
  cluster_methods_expanded <- cluster_methods
  multiple_linkages <- "hierarchical" %in% cluster_methods & length(hierarchical_linkage) > 1
  if (multiple_linkages) {
    cluster_methods_expanded <- c(cluster_methods[cluster_methods != "hierarchical"], 
                                  paste0("hierarchical_", hierarchical_linkage))
  }
  
  # Allocate array for clustering results
  #clusters <- array(dim = c(dim(temp)[1], length(n_clusters), length(cluster_methods_expanded)))
  #dimnames(clusters) <- list(id = rownames(temp), k = n_clusters, m = cluster_methods_expanded)
  
  # Alternative
  clusters <- data.frame()
  
  # Metrics collected to data.frame
  metrics <- data.frame()
  
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
        silh_k <- cluster::silhouette(x = temp_k, dist = diss)
        metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                             metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = temp_k))
      }
    } else if (cluster_methods_expanded[i] == "diana") {
      clust_i <- cluster::diana(x = diss, diss = TRUE)
      for (j in 1:length(n_clusters)) {
        temp_k <- cutree(clust_i, n_clusters[j])
        silh_k <- cluster::silhouette(x = temp_k, dist = diss)
        metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                             metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
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
        silh_k <- cluster::silhouette(x = clust_k$clusters, dist = diss)
        metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                             metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
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
          silh_k <- cluster::silhouette(x = clust_k$classification, dist = diss)
          metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                               metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
          clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                                 k = n_clusters[j], cluster = clust_k$classification))
        }
      }
    } else if (cluster_methods_expanded[i] == "knn_communities") {
      clust_k <- knn_communities(t(temp), k = knn_neighbours, jaccard_kernel = knn_jaccard)
      silh_k <- cluster::silhouette(x = clust_k, dist = diss)
      metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = "variable", 
                                           metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
      clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                             k = "variable", cluster = clust_k))
    } else if (cluster_methods_expanded[i] == "spectral") {
      for (j in 1:length(n_clusters)) {
        clust_k <- Spectrum::Spectrum(t(temp), method = 3, fixk = n_clusters[j])
        silh_k <- cluster::silhouette(x = clust_k$assignments, dist = diss)
        metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                             metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = clust_k$assignments))
      }
    } else if (cluster_methods_expanded[i] == "SC3") {
      for (j in 1:length(n_clusters)) {
        # SC3 only accepts input in the form of SingleCellExperiment 
        hack <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(temp)))
        SummarizedExperiment::rowData(hack)$feature_symbol <- colnames(temp)
        hack <- SC3::sc3(hack, ks = n_clusters[j], gene_filter = FALSE, n_cores = NULL)
        clust_k <- cutree(hack@metadata$sc3$consensus[[1]]$hc, n_clusters[j])
        silh_k <- cluster::silhouette(x = clust_k, dist = diss)
        metrics <- rbind(metrics, data.frame(m = cluster_methods_expanded[i], k = n_clusters[j], 
                                             metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
        clusters <- rbind(clusters, data.frame(id = rownames(temp), m = cluster_methods_expanded[i], 
                                               k = n_clusters[j], cluster = clust_k))
      }
    }else {
      stop(paste("Unsupported method:", cluster_methods_expanded[i]))
    }
  }
  
  out_list <- list()
  
  # Combine outputs with metadata
  out_list$clusters <- plyr::join(dat[!grepl("^dim[0-9]+$", colnames(dat))], 
                                  clusters, 
                                  by = "id")
  out_list$clusters <- out_list$clusters[!is.na(out_list$clusters$cluster),] # potential issue with reference fold missing later
  
  out_list$metrics <- metrics
  
  # Insert meta data (if present)
  out_list$metrics$run <- dat$run[1]
  out_list$metrics$fold <- dat$fold[1]
  out_list$metrics$datname <- dat$datname[1]
  out_list$metrics$drname <- dat$drname[1]
  
  # Pearson's chi-squared test
  f1 <- function(x, c1, c2) {
    temp <- tryCatch(suppressWarnings(chisq.test(x[[c1]], x[[c2]])), error = function(e) NULL)
    if (!is.null(temp)) {
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
  return(out_list)
}

#' Gene module score
#' 
#' Metric based on gene module eigen gene correlation
#' 
#' @param clust A \code{data.frame} with columns "id" and "cluster".
#' @param module_eigs Gene module eigen-genes for each sample (samples x modules).
#' @param module_cor_threshold Threshold for counting correlations.
#' @param module_nan.substitute Substituted value when dividing by zero when there are no correlated clusters for a module.
#' 
#' @export
#' @examples library(COPS)
#' library(WGCNA)
#' library(parallel)
#' 
#' # Generate module eigen genes with WGCNA
#' gene_correlation <- cor(t(ad_ge_micro_zscore), method =  "spearman")
#' adj <- WGCNA::adjacency.fromSimilarity(gene_correlation, power = 2)
#' TOM <- WGCNA::TOMsimilarity(adj, TOMType = "unsigned")
#' geneTree <- flashClust::flashClust(as.dist(1 - TOM), method="average")
#' dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 20)
#' adj_modules <- WGCNA::adjacency.fromSimilarity(gene_correlation[dynamicMods != 0, dynamicMods != 0], power = 2)
#' MEList <- WGCNA::moduleEigengenes(t(ad_ge_micro_zscore[dynamicMods != 0,]), colors = dynamicMods[dynamicMods != 0])
#' MEs <- MEList$eigengenes
#' 
#' # Compute the module score
#' clust <- cutree(hclust(as.dist(1-cor(ad_ge_micro_zscore, method = "spearman")), method = "average"), k = 3)
#' clust <- data.frame(id = names(clust), cluster = clust)
#' 
#' score <- gene_module_score(clust, MEs)
#' 
#' # Within full pipeline
#' res <- dimred_clusteval_pipeline(ad_ge_micro_zscore, 
#' batch_label = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap", "tsne"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' metric = "euclidean", 
#' n_clusters = 2:4, 
#' module_eigs = MEs)
#' 
#' scores <- clusteval_scoring(res, wsum = Silhouette - Module_score, summarise = TRUE)
gene_module_score <- function(clust, 
                              module_eigs, 
                              module_cor_threshold = 0.3, 
                              module_nan.substitute = 0) {
  clust_cor <- lapply(as.data.frame(module_eigs[clust$id,]), 
                                    function(x) sapply(unique(clust$cluster), 
                                                       function(y) cor(x, clust$cluster == y)))
  clust_cor_mat <- Reduce("rbind", clust_cor)
  clust_cor_mat_pos <- clust_cor_mat > module_cor_threshold
  clust_cor_mat_neg <- clust_cor_mat < -module_cor_threshold
  
  score <- apply(clust_cor_mat_pos, 1, function(x) min(1, sum(x)))
  score <- score + apply(clust_cor_mat_neg, 1, function(x) min(1, sum(x)))
  score <- score / apply(clust_cor_mat_pos + clust_cor_mat_neg, 1, sum)
  score[is.nan(score)] <- module_nan.substitute
  score <- mean(score, na.rm = TRUE)
  return(score)
}

#' Gene module correlation score for multiple clustering results
#' 
#' Runs \code{\link{gene_module_score}} in parallel. 
#' 
#' Assumes that a parallel backend has been registered for \code{foreach}.
#'
#' @param clusters A data.table or data.frame with clustering information. 
#' @param module_eigs See \code{\link{gene_module_score}}.
#' @param module_cor_threshold See \code{\link{gene_module_score}}.
#' @param module_nan.substitute See \code{\link{gene_module_score}}.
#' @param by Column names that identify a single clustering result.
#' @param ... Extra arguments are ignored.
#'
#' @return
#' @export
module_evaluation <- function(clusters, 
                              module_eigs, 
                              module_cor_threshold = 0.3, 
                              module_nan.substitute = 0, 
                              by = c("run", "fold", "datname", "drname", "k", "m"), 
                              ...) {
  if (data.table::is.data.table(clusters)) {
    clust_list <- split(clusters, by = by)
  } else {
    clust_list <- split(clusters, clusters[, by])
  }
  
  out <- foreach(clust = clust_list,
                 .combine = function(...) data.table::rbindlist(list(...)),
                 .packages = c(),
                 .multicombine = TRUE,
                 .maxcombine = length(clust_list)) %dopar% {
                   gm_score <- gene_module_score(clust, module_eigs, module_cor_threshold, module_nan.substitute)
                   gm_score <- data.frame(as.data.frame(clust)[1,by], Module_score = gm_score)
                   gm_score
                 }
  return(out)
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
#' @param ... extra arguments are passed through to clustering_evaluation
#'
#' @return Returns a \code{list} of \code{data.frames} containing \code{\link{clustering_evaluation}} outputs for every
#'         combination of CV run, CV fold, clustering method, number of clusters as well as all combinations of
#'         data sets and dimensionality reduction techniques found in the input \code{data.frame}.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table rbindlist
cv_clusteval <- function(dat_embedded, ...) {
  temp_list <- list()
  for (i in 1:length(dat_embedded)) {
    temp <- dat_embedded[[i]]
    drname <- names(dat_embedded)[i]
    if (is.null(drname)) drname <- i
    temp$drname <- drname
    temp <- split(temp, temp[,c("run", "fold")])
    temp_list <- c(temp_list, temp)
  }
  # Binding function that concatenates relevant list components
  cfun <- function(...){
    bound_list <- list()
    bound_list$clusters <- rbindlist(lapply(list(...), function(x) x$clusters))
    bound_list$metrics <- rbindlist(lapply(list(...), function(x) x$metrics))
    bound_list$chisq_pval <- rbindlist(lapply(list(...), function(x) x$chisq_pval))
    bound_list$batch_association <- rbindlist(lapply(list(...), function(x) x$batch_association))
    bound_list$subtype_association <- rbindlist(lapply(list(...), function(x) x$subtype_association))
    return(bound_list)
  }
  
  out <- foreach(temp = temp_list,
                 .combine = cfun,
                 .export = c("clustering_evaluation"),
                 .packages = c("reshape2", "mclust", "cluster", "flashClust", "ClusterR"),
                 .multicombine = TRUE,
                 .maxcombine = length(temp_list)) %dopar% {
    temp <- clustering_evaluation(temp, ...)
    temp
  }
  return(out)
}





