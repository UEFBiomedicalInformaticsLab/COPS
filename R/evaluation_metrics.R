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

silhouette_adjusted <- function(x, diss, min_size = 0.05) {
  csize <- table(x)
  under_size <- which(csize/length(x) < min_size)
  diss_main <- as.matrix(diss)[,!x %in% names(csize)[under_size]]
  x_main <- x[!x %in% names(csize)[under_size]]
  out <- x
  for (i in under_size) {
    out[x == names(csize)[i]] <- x_main[apply(diss_main[x == names(csize)[i], ], 1, which.min)]
  }
  silh <- cluster::silhouette(out, diss)
  if (is.na(silh)) {
    silh <- data.frame(sil_width = NA)
  }
  return(silh)
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
                           parallel = 1, 
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
  
  temp_list <- split_by_safe(clust, by)
  
  parallel_clust <- setup_parallelization(parallel)
  
  stability <- tryCatch(foreach(temp = temp_list,
                      .combine = function(...) data.table::rbindlist(list(...)),
                      .export = c(),
                      .packages = c("clusteval", "data.table", "aricode"),
                      .multicombine = TRUE,
                      .maxcombine = max(length(temp_list), 2)) %dopar% {
    out <- f2(temp)
    for (j in by) {
      out[[j]] <- temp[[j]][1]
    }
    out
  }, finally = if(parallel > 1) parallel::stopCluster(parallel_clust))
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
embedding_associations <-  function(dat, 
                                    class, 
                                    n_pc_max = 10, 
                                    ...){
  out <- list()
  if (is.null(dim(class))) class <- cbind(as.character(class), c())

  pca_silh <- list()
  pca_reg <- list()
  #DSC_res <- c()
  for (i in 1:ncol(class)) {
    # PCA based
    dat_pca <- FactoMineR::PCA(t(dat),
                               scale.unit = FALSE,
                               ncp = min(c(n_pc_max, dim(dat))),
                               graph = FALSE)
    pca_silh[[i]] <- kBET::batch_sil(list(x = dat_pca$ind$coord),
                                     class[,i],
                                     nPCs = min(c(n_pc_max, dim(dat))))
    pca_reg[[i]] <- suppressWarnings(kBET::pcRegression(list(x = dat_pca$ind$coord),
                                       class[,i],
                                       n_top = min(c(n_pc_max, dim(dat)))))

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

#' Association analysis for clustering results
#'
#' @param clust data.frame with columns id and cluster.
#' @param association_data data.frame with association variables, rownames should match with id in \code{clust}.
#'
#' @return A data.frame of NMI, ARI and Chi-squared test p-values for categorical variables 
#'         as well as analysis of variance and kruskal-wallis test p-values for 
#'         continuous variables. 
#' @export
cluster_associations <- function(clust, association_data) {
  association_data_matched <- association_data[match(clust$id, 
                                               rownames(association_data)), 
                                         , 
                                         drop = FALSE]
  out <- list()
  for (i in 1:ncol(association_data)) {
    association_var <- association_data_matched[,i]
    nna_ind <- which(!is.na(association_var))
    if (length(nna_ind) > 1) {
      if (length(unique(association_var[nna_ind])) > 1 &
          length(unique(clust$cluster[nna_ind])) > 1) {
        valid <- TRUE
      } else {
        valid <- FALSE
      }
    } else {
      valid <- FALSE
    }
    if(valid) {
      if (class(association_var) %in% c("character", "factor")) {
        out[[i]] <- data.frame(nmi = aricode::NMI(clust$cluster[nna_ind], 
                                                  association_var[nna_ind]), 
                               ari = aricode::ARI(clust$cluster[nna_ind], 
                                                  association_var[nna_ind]))
        temp <- tryCatch(suppressWarnings(chisq.test(clust$cluster[nna_ind], 
                                                     association_var[nna_ind])), 
                         error = function(e) NULL)
        if (is.null(temp)) {
          out[[i]]$chisq.p <- NA
        } else {
          out[[i]]$chisq.p <- temp$p.value
        }
        colnames(out[[i]]) <- paste0(colnames(association_data)[i], ".", colnames(out[[i]]))
      } else if (class(association_var) %in% c("numeric", "integer")) {
        anova_res <- stats::aov(association_var[nna_ind] ~ clust$cluster[nna_ind])
        kruskal_res <- stats::kruskal.test(association_var[nna_ind], 
                                           clust$cluster[nna_ind])
        out[[i]] <- data.frame(aov.p = summary(anova_res)[[1]][1,"Pr(>F)"],
                               kruskal.p = kruskal_res$p.value)
        colnames(out[[i]]) <- paste0(colnames(association_data)[i], ".", colnames(out[[i]]))
      } else {
        warning(paste0("Unknown data type in association data column ", 
                       colnames(association_data)[i], "."))
      }
    }
  }
  return(Reduce("cbind", out[!sapply(out, is.null)]))
}


#' 
#' 
#' Runs \code{\link{cluster_associations}} in parallel. 
#' 
#' Assumes that a parallel backend has been registered for \code{foreach}.
#'
#' @param clusters A data.table or data.frame with clustering information. 
#' @param association_data A data.frame with association variables (categoric or numeric).
#' @param by Column names that identify a single clustering result.
#' @param ... Extra arguments are ignored.
#'
#' @return \code{data.table} containing association variables
#' @export
association_analysis_cv <- function(clusters, 
                                association_data, 
                                by = c("run", "fold", "datname", "drname", "k", "m"), 
                                parallel = 1, 
                                ...) {
  clust_list <- split_by_safe(clusters, by)
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(foreach(clust = clust_list,
                 .combine = function(...) data.table::rbindlist(list(...), fill = TRUE),
                 .export = c("cluster_associations"),
                 .packages = c(),
                 .multicombine = TRUE,
                 .maxcombine = max(length(clust_list), 2)) %dopar% {
                   cluster_assoc <- cluster_associations(clust, association_data)
                   if (!is.null(cluster_assoc)) {
                     cluster_assoc <- data.frame(as.data.frame(clust)[1,by], cluster_assoc)
                   }
                   cluster_assoc
                 }, finally = if(parallel > 1) parallel::stopCluster(parallel_clust))
  return(out)
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
                              parallel = 1, 
                              ...) {
  clust_list <- split_by_safe(clusters, by)
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(foreach(clust = clust_list,
                 .combine = function(...) data.table::rbindlist(list(...)),
                 .export = c("gene_module_score"),
                 .packages = c(),
                 .multicombine = TRUE,
                 .maxcombine = max(length(clust_list), 2)) %dopar% {
                   gm_score <- gene_module_score(clust, module_eigs, module_cor_threshold, module_nan.substitute)
                   gm_score <- data.frame(as.data.frame(clust)[1,by], Module_score = gm_score)
                   gm_score
                 }, finally = if(parallel > 1) parallel::stopCluster(parallel_clust))
  return(out)
}
