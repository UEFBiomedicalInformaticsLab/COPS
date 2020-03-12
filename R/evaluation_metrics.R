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
  # Dispersion within batches
  Sw <- scale(t(data_matrix), center = TRUE, scale = FALSE)**2 / dim(data_matrix)[2] # vectorized
  Dw <- sqrt(sum(Sw))
  #Dw <- sqrt(sum(sapply(split(data.frame(t(data_matrix)), batch_label), function(x) sum(diag(cov(x))))))
  # Dispersion between batches
  M = apply(data_matrix, 1, mean)
  mean_var <- function(x) {o <- apply(x, 2, mean) - M; return(sum(o**2))} # vectorized
  #mean_var <- function(x) {o <- apply(x, 2, mean) - M; return(sum(diag(o%*%t(o))))}
  Sb <- sapply(split(data.frame(t(data_matrix)), batch_label), mean_var)
  Db = sqrt(sum(Sb))
  return(Db/Dw)
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
#' @param by2 vector of column names to split by, "run" and "fold" by default
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{data.frame} where each row corresponds to clustering stability
#'         with respect to kept column variables
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table data.table setDT rbindlist
#'
stability_eval <- function(clust,
                           by = c("datname", "drname", "k", "m"),
                           by2 = c("run", "fold"),
                           ...)
{
  f <- function(x) {
    data.table::setDTthreads(1)
    ref_i <- unique(x$fold)
    ref_i <- ref_i[!ref_i %in% unique(x$cv_index)]
    
    ref <- x[fold == ref_i,]
    colnames(ref)[colnames(ref) == "cluster"] <- "reference_cluster"
    ref$fold <- NULL
    
    nonref <- x[fold != ref_i,]
    nonref$test_ind <- nonref$cv_index == nonref$fold
    
    #nonref <- plyr::join(nonref, ref[c("id", by2, "reference_cluster")], 
    #                     by = c("id", by2[by2 != "fold"]), type = "inner")
    ref_cols <- c("id", by2[by2 != "fold"], "reference_cluster")
    nonref <- merge(nonref, ref[, ..ref_cols], by = c("id", by2[by2 != "fold"]))
    
    #train_nonref <- nonref[test_ind == FALSE, list(cluster), by = by2]
    train_nonref <- split(nonref$cluster[!nonref$test_ind], nonref[!nonref$test_ind, ..by2])
    train_ref <- split(nonref$reference_cluster[!nonref$test_ind], nonref[!nonref$test_ind, ..by2])
    train_jdist <- jdist_ref(train_nonref, train_ref)
    
    test_ref <- split(nonref$cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
    test_nonref <- split(nonref$reference_cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
    test_jdist <- jdist_ref(test_nonref, test_ref)
    
    return(data.table::data.table(jdist_train = train_jdist, jdist_test = test_jdist))
  }
  by <- by[by %in% colnames(clust)]
  
  data.table::setDT(clust)
  
  temp_list <- split(clust, by = by)
  #temp_list <- split(clust, clust[by])
  stability <- foreach(temp = temp_list,
                      .combine = function(...) data.table::rbindlist(list(...)),
                      .export = c("jdist_ref"),
                      .packages = c("clusteval", "data.table"),
                      .multicombine = TRUE,
                      .maxcombine = length(temp_list)) %dopar% {
    out <- f(temp)
    for (j in by) {
      out[[j]] <- temp[[j]][1]
    }
    out
  }
  return(as.data.frame(stability))
}

#' Categorical variable association estimates
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
#'         \code{\link[KBET]{pcRegression}} and \code{\link{DSC}}
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
  DSC_res <- c()
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
    DSC_res[i] <- DSC(dat, class[,i])
  }
  if (!is.null(dimnames(class))) {
    names(pca_silh) <- dimnames(class)[[2]]
    names(pca_reg) <- dimnames(class)[[2]]
    names(DSC_res) <- dimnames(class)[[2]]
  }
  out <- list(PCA_silhouette = pca_silh,
              PCA_regression = pca_reg,
              DSC = DSC_res)
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
#' \item "agnes" -
#' \item "kmeans" -
#' \item "pam" -
#' \item "sota" -
#' }
#' The clustering methods are interfaced through \code{\link[clValid]{clValid}} so please refer to
#' their documentation for options on metrics.
#'
#' @param dat data matrix with samples on columns
#' @param batch_label_names character vector containing column names corresponding to batch labels
#' @param n_clusters vector of integers, numbers of clusters to be generated
#' @param cluster_methods vector of clustering method names, see details for options
#' @param metric distance metric used for clustering, see details for options
#' @param ... extra arguments are ignored currently
#'
#' @return Returns a \code{list} containing clusters, metrics, and
#'         \code{\link[stats]{chisq.test}} p-values
#'         if batch_label was supplied
#' @export
#' @importFrom stats chisq.test cutree
clustering_evaluation <- function(dat,
                                  batch_label_names = NULL,
                                  n_clusters = 2:5,
                                  cluster_methods = c("hierarchical","pam","diana","kmeans"),
                                  metric = "euclidean",
                                  ...) {
  temp <- dat[grepl("^dim[0-9]+$", colnames(dat))]
  temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
  rownames(temp) <- dat$id
  out <- clValid(as.matrix(temp),
                 n_clusters,
                 clMethods = cluster_methods,
                 metric = metric,
                 validation="internal")
  names(dimnames(out@measures)) <- c("metric", "k", "m")
  # Extract clusters into array
  clusters <- array(dim = c(dim(temp)[1], length(n_clusters), length(cluster_methods)))
  dimnames(clusters) <- list(id = rownames(temp), k = n_clusters, m = cluster_methods)
  for (j in 1:length(cluster_methods)) {
    temp <- clusters(out, cluster_methods[j])
    for (i in 1:length(n_clusters)) {
      if (cluster_methods[j] %in% c("hierarchical", "diana", "agnes")) {
        temp_k <- cutree(temp, n_clusters[i])
      } else if (cluster_methods[j] == "pam") {
        temp_k <- temp[[i]]$clustering
      } else if (cluster_methods[j] == "kmeans") {
        temp_k <- temp[[i]]$cluster
      } else if (cluster_methods[j] == "sota") {
        temp_k <- temp[[i]]$clust
      } else {
        stop(paste("Unsupported method:", cluster_methods[j]))
      }
      clusters[,i,j] <- temp_k
    }
  }
  out_list <- list()
  # Combine outputs with metadata
  out_list$clusters <- plyr::join(dat[!grepl("^dim[0-9]+$", colnames(dat))], 
                                  reshape2::melt(clusters, value.name = "cluster"), 
                                  by = "id")
  out_list$metrics <- reshape2::melt(out@measures, value.name = "value")
  # Insert meta data (if present)
  out_list$metrics$run <- dat$run[1]
  out_list$metrics$fold <- dat$fold[1]
  out_list$metrics$datname <- dat$datname[1]
  out_list$metrics$drname <- dat$drname[1]
  
  if (!is.null(batch_label_names)) {
    f <- function(x) {
      out <- data.frame()
      for (i in 1:length(batch_label_names)) {
        temp <- suppressWarnings(chisq.test(table(x[c(batch_label_names[i], "cluster")])))
        temp <- data.frame(p = temp$p.value)
        temp$run <- x$run[1]
        temp$fold <- x$fold[1]
        temp$datname <- x$datname[1]
        temp$drname <- x$drname[1]
        temp$k <- x$k[1]
        temp$m <- x$m[1]
        temp$batch_label <- batch_label_names[i]
        out <- rbind(out, temp)
      }
      return(out)
    }
    
    #out_list$chisq_pval <- plyr::ddply(out_list$clusters, c("k", "m"), f)
    out_list$chisq_pval <- Reduce(rbind, lapply(split(out_list$clusters, out_list$clusters[c("k", "m")]), f))
  }
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
#' @param dat_folded list of \code{data.frame}s
#' @param ... extra arguments are passed through to clustering_evaluation
#'
#' @return Returns a \code{list} of \code{data.frames} containing \code{\link{clustering_evaluation}} outputs for every
#'         combination of CV run, CV fold, clustering method, number of clusters as well as all combinations of
#'         data sets and dimensionality reduction techniques found in the input \code{data.frame}.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table rbindlist
cv_clusteval <- function(dat_folded, ...) {
  # If runs and folds are already separated, this produces a list of length 1
  temp_list <- list()
  for (i in 1:length(dat_folded)) {
    temp <- dat_folded[[i]]
    temp$drname <- names(dat_folded)[i]
    if (is.null(temp$drname)) temp$drname <- i
    #temp <- plyr::dlply(temp, c("run", "fold"), function(x) x)
    temp <- split(temp, temp[c("run", "fold")])
    temp_list <- c(temp_list, temp)
  }
  # Binding function that concatenates relevant list components
  cfun <- function(...){
    bound_list <- list()
    bound_list$clusters <- as.data.frame(rbindlist(lapply(list(...), function(x) x$clusters)))
    bound_list$metrics <- as.data.frame(rbindlist(lapply(list(...), function(x) x$metrics)))
    bound_list$chisq_pval <- as.data.frame(rbindlist(lapply(list(...), function(x) x$chisq_pval)))
    return(bound_list)
  }
  
  out <- foreach(temp = temp_list,
                 .combine = cfun,
                 .export = c("clustering_evaluation"),
                 .packages = c("clValid", "reshape2"),
                 .multicombine = TRUE,
                 .maxcombine = length(temp_list)) %dopar% {
    temp <- clustering_evaluation(temp, ...)
    temp
  }
  return(out)
}





