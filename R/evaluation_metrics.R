#' Dispersion Separability Criterion
#'
#' Used by TCGA Batch Effects Viewer \url{https://bioinformatics.mdanderson.org/public-software/tcga-batch-effects/}.
#' Based on \url{https://www.jmlr.org/papers/v5/dy04a.html}.
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
  if (length(silh) == 1 && is.na(silh)) {
    silh <- data.frame(sil_width = NA)
  }
  return(silh)
}

#' Clustering stability evaluation
#'
#' Performs stability analysis on cross-validated clusterings.
#'
#' Default settings work with \code{\link{cv_clusteval}} output 'clusters'.
#'
#' @param clusters clustering \code{data.frame} such as returned by \code{\link{cv_clusteval}}
#' @param by vector of column names to keep
#' @param parallel number of threads
#' @param reference_fold fold number that corresponds to reference which other folds are compared against, inferred from input by default
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{data.frame} where each row corresponds to clustering stability
#'         with respect to kept column variables
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table data.table is.data.table rbindlist setDT setDTthreads
#'
stability_eval <- function(clusters,
                           by = c("datname", "drname", "run", "k", "m"),
                           #by2 = c("fold"),
                           parallel = 1, 
                           reference_fold = NULL,
                           ...)
{
  if (is.null(reference_fold)) {
    if (!"cv_index" %in% colnames(clusters)) {
      stop("Please define the reference fold number.")
    } else {
      reference_fold <- unique(clusters$fold)
      reference_fold <- reference_fold[!reference_fold %in% unique(clusters$cv_index)]
    }
  }
  by2 = c("fold")
  # Function to be applied for each result
  f1 <- function(clust, clustref) {
    if (length(clust) != length(clustref)) {
      stop("Number of given references does not match number of given inputs.")
    }
    if (length(clust) > 0) {
      jsc <- c()
      nmi <- c()
      ari <- c()
      for (i in 1:length(clust)) {
        jsc[i] <- jaccard_similarity(clust[[i]], clustref[[i]])
        nmi[i] <- igraph::compare(clust[[i]], clustref[[i]], method = "nmi")
        ari[i] <- igraph::compare(clust[[i]], clustref[[i]], method = "adjusted.rand")
      }
    } else {
      jsc <- NA
      nmi <- NA
      ari <- NA
    }
    return(list(jsc = jsc, nmi = nmi, ari = ari))
  }
  # Function to be applied for method combinations in clustering table
  f2 <- function(x, ref_i) {
    data.table::setDTthreads(1)
    
    ref <- x[fold == ref_i,]
    colnames(ref)[colnames(ref) == "cluster"] <- "reference_cluster"
    ref$fold <- NULL
    
    nonref <- x[fold != ref_i,]
    if ("cv_index" %in% colnames(x)) {
      nonref$test_ind <- nonref$cv_index == nonref$fold
    } else {
      # Assume all samples are from the training set
      nonref$test_ind <- FALSE
    }
    
    ref_cols <- c("id", by2[by2 != "fold"], "reference_cluster")
    if ("data.table" %in% class(ref)) {
      nonref <- merge(nonref, ref[, ..ref_cols], by = c("id", by2[by2 != "fold"]))
    } else {
      nonref <- merge(nonref, ref[, ref_cols], by = c("id", by2[by2 != "fold"]))
    }
    
    if (any(!nonref$test_ind)) {
      if("data.table" %in% class(nonref)) {
        train_nonref <- split(nonref$cluster[!nonref$test_ind], nonref[!nonref$test_ind, ..by2])
        train_ref <- split(nonref$reference_cluster[!nonref$test_ind], nonref[!nonref$test_ind,  ..by2])
        train_res <- f1(train_nonref, train_ref)
      } else {
        train_nonref <- split(nonref$cluster[!nonref$test_ind], nonref[!nonref$test_ind, by2])
        train_ref <- split(nonref$reference_cluster[!nonref$test_ind], nonref[!nonref$test_ind, by2])
        train_res <- f1(train_nonref, train_ref)
      }
    } else {
      train_res <- list()
    }
    
    
    if (any(nonref$test_ind)) {
      if("data.table" %in% class(nonref)) {
        test_ref <- split(nonref$cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
        test_nonref <- split(nonref$reference_cluster[nonref$test_ind], nonref[nonref$test_ind, ..by2])
        test_res <- f1(test_nonref, test_ref)
      } else {
        test_ref <- split(nonref$cluster[nonref$test_ind], nonref[nonref$test_ind, by2])
        test_nonref <- split(nonref$reference_cluster[nonref$test_ind], nonref[nonref$test_ind, by2])
        test_res <- f1(test_nonref, test_ref)
      }
    } else {
      test_res <- list()
    }
    
    
    out_f2 <- data.table::data.table(fold = names(train_ref), 
                                     train_jsc = train_res$jsc,
                                     train_nmi = train_res$nmi,
                                     train_ari = train_res$ari,
                                     test_jsc = test_res$jsc,
                                     test_nmi = test_res$nmi,
                                     test_ari = test_res$ari)
    return(out_f2)
  }
  by <- by[by %in% colnames(clusters)]
  
  temp_list <- split_by_safe(clusters, by)
  
  parallel_clust <- setup_parallelization(parallel)
  
  stability <- tryCatch(foreach(temp = temp_list,
                        .combine = function(...) data.table::rbindlist(list(...)),
                        .export = c(),
                        .packages = c("data.table", "igraph"),
                        .multicombine = TRUE,
                        .maxcombine = max(length(temp_list), 2)) %dopar% {
    out <- tryCatch({
      out <- f2(temp, reference_fold)
      for (j in by) {
        out[[j]] <- temp[[j]][1]
      }
      out}, error = function(e) return(NULL))
    out
  }, finally = close_parallel_cluster(parallel_clust))
  return(as.data.frame(stability))
}

#' Feature batch-effect analysis
#'
#' Three batch effect estimators for numeric features put together. \cr
#' Can be used for other purposes, e.g., as a multivariate association index. 
#'
#' @param x data matrix, samples on columns
#' @param class vector, data.frame or list where columns are categorical variables
#' @param n_pc_max maximum number of principal components to analyze
#' @param ... extra arguments are ignored
#'
#' @return Returns a \code{list} containing the following objects: 
#' \itemize{
#'   \item{associations}{\code{data.frame} summary of \code{\link[kBET]{batch_sil}},
#'                       \code{\link[kBET]{pcRegression}} (\code{$maxR2}) and 
#'                       \code{\link{DSC}} computed for each column of \code{class}}
#'   \item{eigenvalues}{PC eigenvalues and percentage of explained variance from 
#'                      \code{\link[FactoMineR]{PCA}}}
#' }
#' @export
#' @importFrom FactoMineR PCA
#' @importFrom kBET batch_sil
#' @importFrom kBET pcRegression
feature_associations <- function(x, 
                                 class, 
                                 n_pc_max = 10, 
                                 ...){
  if (class(class) != "list" & is.null(dim(class))) class <- data.frame(class = as.character(class))
  class <- lapply(class, as.factor)
  class <- class[sapply(class, function(x) length(levels(x)))>1]
  
  pca_silh <- c()
  pca_reg <- list()
  DSC_res <- c()
  
  dat_pca <- FactoMineR::PCA(t(x),
                             scale.unit = FALSE,
                             ncp = min(c(n_pc_max, dim(x))),
                             graph = FALSE)
  
  for (i in 1:length(class)) {
    # PCA based
    pca_silh[i] <- kBET::batch_sil(list(x = dat_pca$ind$coord[!is.na(class[[i]]),]),
                                     class[[i]][!is.na(class[[i]])],
                                     nPCs = min(c(n_pc_max, dim(x))))
    pca_reg[[i]] <- kBET::pcRegression(list(x = dat_pca$ind$coord[!is.na(class[[i]]),], 
                                          sdev = sqrt(dat_pca$eig[,"eigenvalue"])),
                                     class[[i]][!is.na(class[[i]])],
                                     n_top = min(c(n_pc_max, dim(x))))

    # Other
    DSC_res[i] <- DSC(x[,!is.na(class[[i]])], class[[i]][!is.na(class[[i]])])
  }
  if (is.null(names(class))) names(class) <- 1:length(class)
  out <- data.frame(class = names(class),
                    PC_silhouette = pca_silh,
                    PC_R2_max = sapply(pca_reg, function(a) a[["maxR2"]]), 
                    DSC = DSC_res)
  
  npc <- min(c(n_pc_max, dim(x)))
  pc_wise_silhouettes <- list()
  for (i in names(class)) {
    pc_wise_silhouettes[[i]] <- sapply(1:npc, function(j) {
      mean(cluster::silhouette(
        as.integer(class[[i]][!is.na(class[[i]])]), 
        dist(dat_pca$ind$coord[!is.na(class[[i]]), j, drop = FALSE]))[,"sil_width"])
      })
  }
  
  return(list(associations = out, eigenvalues = dat_pca$eig, pc_wise_silhouettes = pc_wise_silhouettes))
}

#' Association analysis for clustering results
#'
#' @param clusters data.frame with columns id and cluster.
#' @param association_data data.frame with association variables, rownames should match with id in \code{clusters}.
#'
#' @return A data.frame of NMI, ARI and Chi-squared test p-values for categorical variables 
#'         as well as analysis of variance and kruskal-wallis test p-values for 
#'         continuous variables. 
#' @export
cluster_associations <- function(clusters, 
                                 association_data) {
  association_data_matched <- association_data[match(clusters$id, 
                                               rownames(association_data)), 
                                         , 
                                         drop = FALSE]
  out <- list()
  for (i in 1:ncol(association_data)) {
    association_var <- association_data_matched[,i]
    nna_ind <- which(!is.na(association_var))
    n_samples <- length(nna_ind)
    valid <- n_samples > 1
    if (valid) {
      n_categories <- length(unique(association_var[nna_ind]))
      n_clusters <- length(unique(clusters$cluster[nna_ind]))
      valid <- n_categories > 1 &
        n_clusters > 1
    }
    if(valid) {
      if (class(association_var) %in% c("character", "factor")) {
        valid <- n_categories < n_samples & 
          n_clusters < n_samples
        if (valid) {
          if (class(association_var) == "factor") {
            clean_var <- droplevels(association_var[nna_ind])
          } else {
            clean_var <- association_var[nna_ind]
          }
          out[[i]] <- data.frame(
            nmi = igraph::compare(
              clusters$cluster[nna_ind], 
              clean_var,
              method = "nmi"), 
            ari = igraph::compare(
              clusters$cluster[nna_ind], 
              clean_var, 
              method = "adjusted.rand"))
          temp <- tryCatch(
            suppressWarnings(
              chisq.test(clusters$cluster[nna_ind], clean_var)
            ), error = function(e) NULL)
          if (is.null(temp)) {
            out[[i]]$chisq.p <- NA
          } else {
            out[[i]]$chisq.p <- temp$p.value
          }
          colnames(out[[i]]) <- paste0(colnames(association_data)[i], ".", colnames(out[[i]]))
        }
      } else if (class(association_var) %in% c("numeric", "integer")) {
        anova_res <- stats::aov(association_var[nna_ind] ~ clusters$cluster[nna_ind])
        kruskal_res <- stats::kruskal.test(association_var[nna_ind], 
                                           clusters$cluster[nna_ind])
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


#' Cross-validated cluster association analysis
#' 
#' Runs \code{\link{cluster_associations}} in parallel. 
#' 
#' Assumes that a parallel backend has been registered for \code{foreach}.
#'
#' @param clusters A data.table or data.frame with clustering information. 
#' @param association_data A data.frame with association variables (categoric or numeric).
#' @param by Column names that identify a single clustering result.
#' @param parallel Number of parallel threads.
#' @param ... Extra arguments are ignored.
#'
#' @return \code{data.table} containing association variables
#' @export
cv_association_analysis <- function(clusters, 
                                association_data, 
                                by = c("run", "fold", "datname", "drname", "k", "m"), 
                                parallel = 1, 
                                ...) {
  by <- by[by %in% colnames(clusters)]
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
                 }, finally = close_parallel_cluster(parallel_clust))
  return(out)
}

#' Gene module score
#' 
#' Metric based on gene module eigen gene correlation
#' 
#' @param x A \code{data.frame} with columns "id" and "cluster".
#' @param module_eigs Gene module eigen-genes for each sample (samples x modules).
#' @param module_cor_threshold Threshold for counting correlations.
#' @param module_nan.substitute Substituted value when dividing by zero when there are no correlated clusters for a module.
#' @param ... Extra arguments are ignored.
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
#' # Compute the module score for a given clustering result
#' clust <- cutree(hclust(as.dist(1-cor(ad_ge_micro_zscore, method = "spearman")), method = "average"), k = 3)
#' clust <- data.frame(id = names(clust), cluster = clust)
#' 
#' score <- gene_module_score(clust, MEs)
#' 
#' # Within full pipeline
#' res <- COPS(ad_ge_micro_zscore, 
#' association_data = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' distance_metric = "euclidean", 
#' n_clusters = 2:4, 
#' module_eigs = MEs)
#' 
#' scores <- scoring(res, wsum = Silhouette - Module_score, summarise = TRUE)
gene_module_score <- function(x, 
                              module_eigs, 
                              module_cor_threshold = 0.3, 
                              module_nan.substitute = 0, 
                              ...) {
  clust_cor <- lapply(as.data.frame(module_eigs[x$id,]), 
                                    function(a) sapply(unique(x$cluster), 
                                                       function(b) cor(a, x$cluster == b)))
  clust_cor_mat <- Reduce("rbind", clust_cor)
  clust_cor_mat_pos <- clust_cor_mat > module_cor_threshold
  clust_cor_mat_neg <- clust_cor_mat < -module_cor_threshold
  
  score <- apply(clust_cor_mat_pos, 1, function(a) min(1, sum(a)))
  score <- score + apply(clust_cor_mat_neg, 1, function(a) min(1, sum(a)))
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
#' @param by Column names that identify a single clustering result.
#' @param parallel Number of parallel threads.
#' @param ... Extra arguments are passed to \code{\link{gene_module_score}}.
#' 
#' @return \code{data.table}
#' @export
cv_module_evaluation <- function(clusters, 
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
                   gm_score <- gene_module_score(x = clust, ...)
                   gm_score <- data.frame(as.data.frame(clust)[1,by], Module_score = gm_score)
                   gm_score
                 }, finally = close_parallel_cluster(parallel_clust))
  return(out)
}
