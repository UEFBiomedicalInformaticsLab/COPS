#' Pipeline combining dimensionality reduction, clustering, cross-validation and evaluation
#'
#' Old deprecated version of the pipeline. See \code{\link{COPS}} for the newer function. 
#'
#' @param dat A single matrix or list of data matrices corresponding to the same data but different pre-processing.  
#' @param nfolds Number of cross-validation folds for stability evaluation and metric estimates.
#' @param nruns Number of cross-validation replicates for stability evaluation and metric estimates.
#' @param batch_label Vector or matrix with categorical variables on columns.
#' @param subtype_label Vector or matrix with categorical variables on columns.
#' @param survival_data Data for survival analysis, see \code{\link{survival_preprocess}} for details.
#' @param module_eigs Module eigen genes for gene module correlation evaluation (see \code{\link{gene_module_score}} for details).
#' @param verbose If \code{TRUE}, prints progress notifications.
#' @param parallel Sets up and registers \code{parallel} number of threads for supported operations.
#' @param pathway_enrichment_method \code{enrichment_method} for \code{\link{genes_to_pathways}}
#' @param ... extra arguments are passed to pipeline components where appropriate
#'
#' @return Returns a \code{list} of pipeline component outputs for given settings and input data sets
#' @export
#' 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom utils flush.console
#' @importFrom data.table as.data.table setDT setkey
#' @importFrom plyr join
dimred_clusteval_pipeline <- function(dat, 
                                      nfolds, 
                                      nruns, 
                                      batch_label = NULL,
                                      subtype_label = NULL,
                                      association_data = NULL,
                                      survival_data = NULL,
                                      module_eigs = NULL,
                                      verbose = TRUE,
                                      parallel = 1,
                                      pathway_enrichment_method = "none",
                                      ...) {
  if ("list" %in% class(dat)) {
    dat_list <- dat
  } else {
    dat_list <- list(dat)
  }
  
  id <- colnames(dat_list[[1]])
  
  # Add batch label(s) as separate column(s)
  if (!is.null(batch_label)) {
    if (is.null(dim(batch_label))) {
      # 1d array given as batch label
      if (!is.null(names(batch_label))) {
        batch_id <- names(batch_label)
      } else {
        # Assume data and batch labels are in the same order
        batch_id <- id
        if (length(dat_list) > 1) {
          if (!Reduce("&", lapply(dat_list, function(x) colnames == id))) {
            stop("Ambiguous batch label, no names given.")
          }
        }
        if (length(id) != length(batch_label)) {
          stop("Batch label dimensions do not match with data and no names were given.")
        }
      }
      batch_table <- cbind(batch_label = as.character(batch_label), c())
    } else { 
      # 2d batch_label
      if (!is.null(rownames(batch_label))) {
        batch_id <- rownames(batch_label)
        batch_table <- batch_label
      } else {
        batch_id <- id
        if (length(dat_list) > 1) {
          if (!Reduce("&", lapply(dat_list, function(x) colnames == id))) {
            stop("Ambiguous batch label, rownames missing.")
          }
        }
        if (length(id) != nrow(batch_label)) {
          stop("Batch label dimensions do not match with data and no rownames were given.")
        }
      }
    }
    batch_table <- data.table::as.data.table(batch_table)
    batch_table$id <- batch_id
  }
  
  if (!is.null(subtype_label)) {
    if (is.null(dim(subtype_label))) {
      if (!is.null(names(subtype_label))) {
        subtype_id <- names(subtype_label)
      } else {
        # Assume data and labels are in the same order
        subtype_id <- id
        if (length(dat_list) > 1) {
          if (!Reduce("&", lapply(dat_list, function(x) colnames == id))) {
            stop("Ambiguous subtype label, no names given.")
          }
        }
        if (length(id) != length(subtype_label)) {
          stop("Subtype label dimensions do not match with data and no names were given.")
        }
      }
      subtype_table <- cbind(subtype_label = as.character(subtype_label), c())
    } else { 
      if (!is.null(rownames(subtype_label))) {
        subtype_id <- rownames(subtype_label)
        subtype_table <- subtype_label
      } else {
        subtype_id <- id
        if (length(dat_list) > 1) {
          if (!Reduce("&", lapply(dat_list, function(x) colnames == id))) {
            stop("Ambiguous subtype label, rownames missing.")
          }
        }
        if (length(id) != nrow(subtype_label)) {
          stop("Subtype label dimensions do not match with data and no rownames were given.")
        }
      }
      
    }
    subtype_table <- data.table::as.data.table(subtype_table)
    subtype_table$id <- subtype_id
  }
  
  
  if (is.null(association_data)) {
    if(!is.null(batch_label)) {
      association_data <- batch_table
    }
    if(!is.null(subtype_label)) {
      if (is.null(association_data)) {
        association_data <- subtype_table
      } else {
        association_data <- merge(association_data, subtype_table, by = "id", all = TRUE)
      }
    }
    if (!is.null(association_data)) {
      association_data <- as.data.frame(association_data)
      rownames(association_data) <- association_data$id
      association_data$id <- NULL
    }
  }
  
  out <- COPS(dat, 
              nfolds, 
              nruns, 
              association_data = association_data,
              survival_data = survival_data,
              module_eigs = module_eigs,
              verbose = verbose,
              parallel = parallel,
              pathway_enrichment_method = pathway_enrichment_method,
              ...)
  
  return(out)
}

#' @describeIn clusteval_scoring Retrieves best clustering from CV results based on scores. 
#' In practice retrieves reference fold result from first run matching the best results. 
#' 
#' @param res result from \code{\link{COPS}}
#' @param scores scores from \code{\link{clusteval_scoring}}
#' 
#' @export
get_best_result <- function(res,
                            scores) {
  if (all(c("run", "fold") %in% colnames(res$best))) {
    stop("Please summarise CV results in order to avoid ambiguity.")
  }
  reference_fold <- unique(res$clusters$fold)[which(!unique(res$clusters$fold) %in% unique(res$clusters$cv_index))]
  clust_ref <- as.data.frame(res$clusters[res$clusters$fold == reference_fold & res$clusters$run == 1,])
  clust_ref <- clust_ref[,-which(colnames(res$clusters) %in% c("fold", "run", "cv_index"))]
  
  out <- list()
  out$clusters <- plyr::join(scores$best[c("datname", "drname", "k", "m")], 
                             clust_ref, 
                             by = c("datname", "drname", "k", "m"))
  
  for (i in which(names(res$embedding) == scores$best$drname)) {
    if(res$embedding[[i]]$run[1] == 1 & 
       res$embedding[[i]]$fold[1] == reference_fold &
       res$embedding[[i]]$datname[1] == scores$best$datname) {
      out$embedding <- res$embedding[[i]][,-which(colnames(res$embedding[[i]]) %in% c("fold", "run", "cv_index"))]
    }
  }
  return(out)
}

#' Clustering algorithms for Omics based Patient Stratification
#'
#' Combines \code{\link{cv_fold}}, \code{\link{cv_pathway_enrichment}},
#' \code{\link{cv_dimred}}, \code{\link{cv_clusteval}}, 
#' \code{\link{stability_eval}}, \code{\link{survival_evaluation}}, 
#' \code{\link{module_evaluation}} and \code{\link{association_analysis_cv}} 
#' to conveniently and comprehensively test clustering algorithms on a given set of input data. 
#'
#' @param dat A single matrix or list of matrices, patients on columns and features on rows. 
#' @param nfolds Number of cross-validation folds for stability evaluation and metric estimates.
#' @param nruns Number of cross-validation replicates for stability evaluation and metric estimates.
#' @param association_data Data for association tests, see \code{\link{cluster_associations}} for details.
#' @param survival_data Data for survival analysis, see \code{\link{survival_preprocess}} for details.
#' @param module_eigs Data for gene module correlation analysis, see \code{\link{gene_module_score}} for details.
#' @param verbose Prints progress messages and time taken. 
#' @param parallel Number of parallel threads for supported operations.
#' @param pathway_enrichment_method \code{enrichment_method} for \code{\link{genes_to_pathways}}.
#' @param multi_omic_methods Character vector of multi-view clustering method names for \code{\link{multi_omic_clustering}}.
#' @param vertical_parallelization (Experimental) if set, all pipeline steps are evaluated in succession within each fold (instead of evaluating each step for all folds before moving on). Always true for multi-view methods. 
#' @param ... Extra arguments are passed to pipeline components where appropriate.
#' 
#' @details
#' If multi_omic_methods is set, then the input matrices are treated as 
#' different views of the same patients. Available methods are listed in the
#' documentation for \code{\link{multi_omic_clustering}}. 
#'
#' @return 
#' Returns a \code{list} of pipeline component outputs for each run, fold and 
#' method given different settings and input data sets. 
#' \itemize{
#'   \item{clusters} \code{data.frame} defining clusters
#'   \item{internal_metrics} \code{data.frame} of internal metrics
#'   \item{stability} \code{data.frame} of stability scores
#'   \item{survival} \code{data.frame} of survival analysis results
#'   \item{modules} \code{data.frame} of gene module association scores
#'   \item{association} \code{data.frame} of association results to variables of interest
#'   \item{cluster_sizes} \code{data.frame} giving the sizes of clusters
#' }
#' 
#' @export
#' @examples library(parallel)
#' library(COPS)
#' 
#' # Dimensionality reduction and clustering (DR-CL)
#' res <- COPS(ad_ge_micro_zscore, 
#' association_data = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap", "tsne"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' metric = "euclidean", 
#' n_clusters = 2:4)
#' 
#' # Clustering (CL)
#' res <- COPS(ad_ge_micro_zscore, 
#' association_data = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("none"), 
#' cluster_methods = c("hierarchical"), 
#' metric = "correlation", 
#' n_clusters = 2:4)
#' 
#' # Biological knowledge integration and clustering (BK-CL)
#' res <- COPS(ad_ge_micro_zscore, 
#' association_data = ad_studies, 
#' pathway_enrichment_method = "DiffRank", 
#' gene_key_expr = "ENSEMBL", 
#' gs_subcats = "CP:KEGG", 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("none"), 
#' cluster_methods = c("hierarchical"), 
#' metric = "correlation", 
#' n_clusters = 2:4)
#' 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom utils flush.console
#' @importFrom data.table as.data.table setDT setkey
#' @importFrom plyr join
COPS <- function(dat, 
                 nfolds, 
                 nruns, 
                 association_data = NULL,
                 survival_data = NULL,
                 module_eigs = NULL,
                 verbose = TRUE,
                 parallel = 1,
                 pathway_enrichment_method = "none",
                 multi_omic_methods = NULL,
                 vertical_parallelization = FALSE, 
                 ...) {
  pipeline_start <- Sys.time()
  
  dat <- data_preprocess(dat, verbose)
  
  if (length(multi_omic_methods) > 0 | vertical_parallelization) {
    out <- vertical_pipeline(dat$dat_list, 
                             nruns = nruns,
                             nfolds = nfolds,
                             survival_data = survival_data, 
                             association_data = association_data,
                             multi_omic_methods = multi_omic_methods, 
                             parallel = parallel,
                             gene_id_list = dat$gene_id_list,
                             ...)
    return(out)
  }
  
  # Create cross validation folds
  cv_index <- cv_fold(dat_list = dat$dat_list, 
                      nfolds = nfolds, 
                      nruns = nruns, 
                      ...)
  
  # Pathway enrichment
  if (pathway_enrichment_method != "none") {
    if (length(dat$dat_list) > 1) {
      stop("Multiple data sets with pathway enrichment enabled is not implemented.")
      # TODO: make it compatible. Currently datname is being used as pwname. 
    }
    pw_start <- Sys.time()
    if(verbose) print("Starting pathway enrichment ..."); flush.console()
    dat_pw <- cv_pathway_enrichment(dat_list = dat$dat_list, 
                                    cv_index = cv_index, 
                                    gene_id_list = dat$gene_id_list, 
                                    enrichment_method = pathway_enrichment_method, 
                                    parallel = parallel, 
                                    ...)
    if(verbose) print(paste("Finished pathway enrichment in",
                            time_taken_string(pw_start))); flush.console()
    
    # Dimensionality reduction for pathway enriched features
    dimred_start <- Sys.time()
    if(verbose) print("Starting dimensionality reduction ..."); flush.console()
    dat_embedded <- cv_dimred(dat_pw, 
                              cv_index, 
                              cv_split_data = FALSE, 
                              parallel = parallel, 
                              ...)
    if(verbose) print(paste("Finished dimensionality reduction in",
                            time_taken_string(dimred_start))); flush.console()
  } else {
    # Dimensionality reduction
    dimred_start <- Sys.time()
    if(verbose) print("Starting dimensionality reduction ..."); flush.console()
    dat_embedded <- cv_dimred(dat$dat_list, 
                              cv_index, 
                              parallel = parallel, 
                              ...)
    if(verbose) print(paste("Finished dimensionality reduction in",
                            time_taken_string(dimred_start))); flush.console()
  }
  
  # Clustering evaluation
  clusteval_start <- Sys.time()
  if(verbose) print("Starting clustering analysis ..."); flush.console()
  dat_clustered <- cv_clusteval(dat_embedded, 
                                parallel = parallel, 
                                ...)
  if(verbose) print(paste("Finished clustering analysis in",
                          time_taken_string(clusteval_start))); flush.console()
  
  # Clustering stability evaluation
  stability_test_start <- Sys.time()
  if(verbose) print("Starting clustering stability analysis ..."); flush.console()
  dat_stability <- stability_eval(dat_clustered$clusters, 
                                  parallel = parallel, 
                                  ...)
  if(verbose) print(paste("Finished clustering stability analysis in",
                          time_taken_string(stability_test_start))); flush.console()
  
  # Survival evaluation
  if (!is.null(survival_data)) {
    survival_analysis_start <- Sys.time()
    if(verbose) print("Starting survival analysis ..."); flush.console()
    dat_survival <- survival_evaluation(survival_data, 
                                        dat_clustered$clusters, 
                                        parallel = parallel, 
                                        ...)
    if(verbose) print(paste("Finished survival analysis in",
                            time_taken_string(survival_analysis_start))); flush.console()
  }
  
  # Gene module correlation evaluation
  if (!is.null(module_eigs)) {
    module_analysis_start <- Sys.time()
    if(verbose) print("Starting gene module correlation analysis ..."); flush.console()
    dat_gm_score <- module_evaluation(dat_clustered$clusters, 
                                      module_eigs, 
                                      parallel = parallel, 
                                      ...)
    if(verbose) print(paste("Finished gene module correlation analysis in",
                            time_taken_string(module_analysis_start))); flush.console()
  }
  
  if (!is.null(association_data)) {
    association_analysis_start <- Sys.time()
    if(verbose) print("Starting variable association analysis ..."); flush.console()
    dat_association_score <- association_analysis_cv(dat_clustered$clusters, 
                                              association_data, 
                                              parallel = parallel, 
                                              ...)
    if(verbose) print(paste("Finished variable association analysis in",
                            time_taken_string(association_analysis_start))); flush.console()
  }
  
  # Return
  out <- list(embedding = dat_embedded, 
              clusters = dat_clustered$clusters, 
              internal_metrics = dat_clustered$metrics,
              stability = dat_stability)
  if (!is.null(survival_data)) out$survival <- dat_survival
  if (!is.null(module_eigs)) out$modules <- dat_gm_score
  if (!is.null(association_data)) out$association <- dat_association_score
  out$cluster_sizes <- dat_clustered$cluster_sizes
  
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}

#' COPS pipeline vertical parallelization
#'
#' @param dat_list 
#' @param nfolds 
#' @param nruns 
#' @param survival_data 
#' @param association_data 
#' @param multi_omic_methods 
#' @param parallel 
#' @param data_is_kernels 
#' @param silhouette_dissimilarities 
#' @param ... 
#'
#' @return
#' @importFrom foreach foreach %dopar% %:%
vertical_pipeline <- function(dat_list, 
                              nfolds = 5, 
                              nruns = 1,
                              survival_data = NULL,
                              association_data = NULL, 
                              multi_omic_methods = NULL, 
                              parallel = 1, 
                              data_is_kernels = FALSE, 
                              silhouette_dissimilarities = NULL,
                              ...) {
  if (length(multi_omic_methods) > 0) {
    cv_index <- multi_view_cv_fold(dat_list = dat_list, nfolds = nfolds, nruns = nruns, ...)
    cv_index_split <- split(cv_index[[1]], by = c("run", "fold"))
    
    cfun <- function(x, y) {
      for (i in 1:length(x)) {
        x[[i]] <- plyr::rbind.fill(x[[i]], y[[i]])
      }
      return(x)
    }
    
    f_args <- list(...)
    
    parallel_clust <- setup_parallelization(parallel)
    out <- tryCatch(foreach(i = 1:length(cv_index_split), 
                    .combine = cfun, 
                    .inorder = FALSE) %:%
    foreach(mvc = multi_omic_methods, 
            .combine = cfun, 
            .export = c("dat_list", "cv_index_split", "f_args", "silhouette_dissimilarities"), 
            .packages = c("COPS"), #"iClusterPlus", "IntNMF", "MOFA2"), 
            .inorder = FALSE) %dopar% {
      dat_i <- list()
      non_data_cols <- list()
      by <- c("run", "fold", "m", "k", "mkkm_mr_lambda")
      
      if(data_is_kernels) {
        for (j in 1:length(dat_list)) {
          if (sum(grepl("^dim[0-9]+$", colnames(dat_list[[j]]))) > nrow(dat_list[[j]])) {
            stop("Input kernels are not square!")
          }
          ij_ind <- match(cv_index_split[[i]]$id, dat_list[[j]]$id)
          dat_i[[j]] <- as.matrix(as.data.frame(dat_list[[j]])[ij_ind, paste0("dim", ij_ind)])
          
          temp <- merge(cv_index_split[[i]], dat_list[[j]], by = "id")
          sel <- grep("^dim[0-9]+$", colnames(temp))
          if ("data.table" %in% class(temp)) {
            non_data_cols[[j]] <- temp[,-..sel]
          } else {
            non_data_cols[[j]] <- temp[,-sel]
          }
        }
      } else {
        for (j in 1:length(dat_list)) {
          dat_i[[j]] <- merge(cv_index_split[[i]], dat_list[[j]], by = "id")
          sel <- grep("^dim[0-9]+$", colnames(dat_i[[j]]))
          if ("data.table" %in% class(dat_i[[j]])) {
            non_data_cols[[j]] <- dat_i[[j]][,-..sel]
            dat_i[[j]] <- as.matrix(dat_i[[j]][,..sel])
          } else {
            non_data_cols[[j]] <- dat_i[[j]][,-sel]
            dat_i[[j]] <- as.matrix(dat_i[[j]][,sel])
          }
        }
      }
      
      # multi-omic clustering
      # 1) multi-view clustering
      # 2) multi-view integration/embedding + clustering
      # Since we are running different methods in parallel we take the first result only
      temp_args <- c(list(dat_list_clust = dat_i, 
                          non_data_cols = non_data_cols, 
                          multi_view_methods = mvc,
                          data_is_kernels = data_is_kernels),
                     f_args)
      clust_i <- do.call(multi_omic_clustering, temp_args)
      
      # Clustering metrics 
      silh_i <- list()
      if (data_is_kernels) {
        # While the unnormalized linear kernel could be used to compute 
        # silhouette in the original space, other kernels cannot. 
        silh_i <- NULL
      } else if (!is.null(silhouette_dissimilarities)) {
        for (j in 1:length(silhouette_dissimilarities)) {
          silh_i[[j]] <- clustering_metrics(clust_i, 
                                            dat = NULL, 
                                            by = by,
                                            clustering_dissimilarity = silhouette_dissimilarities[[j]], 
                                            cluster_size_table = FALSE, 
                                            silhouette_min_cluster_size = 0.0,
                                            distance_metric = "euclidean")$metrics
          silh_i[[j]]$metric[silh_i[[j]]$metric == "Silhouette"] <- paste0(names(silhouette_dissimilarities)[j], "_Silhouette")
        }
        silh_i <- Reduce(rbind, silh_i)
      } else {
        for (j in 1:length(dat_i)) {
          silh_i[[j]] <- clustering_metrics(clust_i, 
                                            dat = as.data.frame(dat_i[[j]]), 
                                            by = by,
                                            clustering_dissimilarity = NULL, 
                                            cluster_size_table = FALSE, 
                                            silhouette_min_cluster_size = 0.0,
                                            distance_metric = "euclidean")$metrics
          silh_i[[j]]$metric[silh_i[[j]]$metric == "Silhouette"] <- paste0(names(dat_list)[j], "_Silhouette")
          #colnames(silh_i[[j]])[colnames(silh_i[[j]]) == "Silhouette"] <- paste0(names(dat_list)[j], "_Silhouette")
        }
        silh_i <- Reduce(rbind, silh_i)
      }
      
      # Survival evaluation
      if (!is.null(survival_data)) {
        temp_args <- c(list(event_data = survival_data, 
                            clusters = clust_i, 
                            parallel = 1,
                            by = by),
                       f_args)
        survival_i <- do.call(survival_evaluation, temp_args)
      } else {
        survival_i <- NULL
      }
      
      # Clustering association analysis
      if (!is.null(association_data)) {
        temp_args <- c(list(clusters = clust_i, 
                            association_data = association_data, 
                            by = by,
                            parallel = 1),
                       f_args)
        association_i <- do.call(association_analysis_cv, temp_args)
      } else {
        association_i <- NULL
      }
      # Return
      out_i <- list()
      out_i$clusters <- clust_i
      out_i$internal_metrics <- silh_i
      out_i$survival <- survival_i
      out_i$association <- association_i
      if (!is.null(attributes(clust_i)$extra_output)) {
        if (!is.null(attributes(clust_i)$extra_output$mkkm_mr_weights)) {
          out_i$mkkm_mr_weights <- attributes(clust_i)$extra_output$mkkm_mr_weights
        }
      }
      out_i
    }, finally = close_parallel_cluster(parallel_clust))
    out$clusters <- data.table::setDT(out$clusters)
    
    # Clustering stability evaluation
    out$stability <- stability_eval(out$clusters, 
                                    by = by[!by %in% "fold"], 
                                    parallel = parallel)
    return(out)
  } else {
    stop("Not implemented.")
  }
}

multi_view_cv_fold <- function(dat_list, nfolds = 5, nruns = 2, ...) {
  if (!all(Reduce("&", lapply(dat_list[-1], function(x) x$id == dat_list[[1]]$id)))) {
    warning("Colnames in all views do not match.")
    stop("Cross-validation for missing sample views not implemented.")
  }
  return(cv_fold(dat_list = dat_list[1], nfolds = nfolds, nruns = nruns, ...))
}

subset_cv_data <- function(dat_list, cv_index, data_is_kernels = FALSE) {
  dat_i <- list()
  non_data_cols <- list()
  if(data_is_kernels) {
    for (j in 1:length(dat_list)) {
      if (sum(grepl("^dim[0-9]+$", colnames(dat_list[[j]]))) > nrow(dat_list[[j]])) {
        stop("Input kernels are not square!")
      }
      ij_ind <- match(cv_index$id, dat_list[[j]]$id)
      dat_i[[j]] <- as.matrix(as.data.frame(dat_list[[j]])[ij_ind, paste0("dim", ij_ind)])
      
      temp <- merge(cv_index, dat_list[[j]], by = "id")
      sel <- grep("^dim[0-9]+$", colnames(temp))
      if ("data.table" %in% class(temp)) {
        non_data_cols[[j]] <- temp[,-..sel]
      } else {
        non_data_cols[[j]] <- temp[,-sel]
      }
    }
  } else {
    for (j in 1:length(dat_list)) {
      dat_i[[j]] <- merge(cv_index, dat_list[[j]], by = "id")
      sel <- grep("^dim[0-9]+$", colnames(dat_i[[j]]))
      if ("data.table" %in% class(dat_i[[j]])) {
        non_data_cols[[j]] <- dat_i[[j]][,-..sel]
        dat_i[[j]] <- as.matrix(dat_i[[j]][,..sel])
      } else {
        non_data_cols[[j]] <- dat_i[[j]][,-sel]
        dat_i[[j]] <- as.matrix(dat_i[[j]][,sel])
      }
    }
  }
  names(dat_i) <- names(dat_list)
  names(non_data_cols) <- names(dat_list)
  return(list(dat_i = dat_i, non_data_cols = non_data_cols))
}

embarrassingly_parallel_pipeline <- function(dat_list, 
                                             cv_index,
                                             fold = 6, 
                                             run = 1,
                                             survival_data = NULL,
                                             association_data = NULL, 
                                             multi_omic_methods = NULL, 
                                             parallel = 1, 
                                             data_is_kernels = FALSE, 
                                             silhouette_dissimilarities = NULL,
                                             ...) {
  cvi <- cv_index[cv_index$fold == fold & cv_index$run == run,]
  if (length(multi_omic_methods) == 1) {
    dat_i <- subset_cv_data(dat_list, cvi, data_is_kernels)
    by <- c("run", "fold", "m", "k", "mkkm_mr_lambda")
    
    # multi-omic clustering
    # 1) multi-view clustering
    # 2) multi-view integration/embedding + clustering
    clust_i <- multi_omic_clustering(dat_list_clust = dat_i$dat_i,
                                     non_data_cols = dat_i$non_data_cols, 
                                     multi_view_methods = multi_omic_methods,
                                     data_is_kernels = data_is_kernels,
                                     ...)
    
    # Clustering metrics 
    silh_i <- list()
    if (data_is_kernels) {
      # While the unnormalized linear kernel could be used to compute 
      # silhouette in the original space, other kernels cannot. 
      silh_i <- NULL
    } else if (!is.null(silhouette_dissimilarities)) {
      for (j in 1:length(silhouette_dissimilarities)) {
        silh_i[[j]] <- clustering_metrics(clust_i, 
                                          dat = NULL, 
                                          by = by,
                                          clustering_dissimilarity = silhouette_dissimilarities[[j]], 
                                          cluster_size_table = FALSE, 
                                          silhouette_min_cluster_size = 0.0,
                                          distance_metric = "euclidean")$metrics
        silh_i[[j]]$metric[silh_i[[j]]$metric == "Silhouette"] <- paste0(names(silhouette_dissimilarities)[j], "_Silhouette")
      }
      silh_i <- Reduce(rbind, silh_i)
    } else {
      for (j in 1:length(dat_i$dat_i)) {
        silh_i[[j]] <- clustering_metrics(clust_i, 
                                          dat = as.data.frame(dat_i$dat_i[[j]]), 
                                          by = by,
                                          clustering_dissimilarity = NULL, 
                                          cluster_size_table = FALSE, 
                                          silhouette_min_cluster_size = 0.0,
                                          distance_metric = "euclidean")$metrics
        silh_i[[j]]$metric[silh_i[[j]]$metric == "Silhouette"] <- paste0(names(dat_list)[j], "_Silhouette")
      }
      silh_i <- Reduce(rbind, silh_i)
    }
    
    # Survival evaluation
    if (!is.null(survival_data)) {
      survival_i <- survival_evaluation(event_data = survival_data, 
                                        clusters = clust_i, 
                                        parallel = 1,
                                        by = by,
                                        ...)
    } else {
      survival_i <- NULL
    }
    
    # Clustering association analysis
    if (!is.null(association_data)) {
      association_i <- association_analysis_cv(clusters = clust_i, 
                                               association_data = association_data, 
                                               by = by,
                                               parallel = 1,
                                               ...)
    } else {
      association_i <- NULL
    }
    # Return
    out <- list()
    out$clusters <- clust_i
    out$internal_metrics <- silh_i
    out$survival <- survival_i
    out$association <- association_i
    if (!is.null(attributes(clust_i)$extra_output)) {
      if (!is.null(attributes(clust_i)$extra_output$mkkm_mr_weights)) {
        out$mkkm_mr_weights <- attributes(clust_i)$extra_output$mkkm_mr_weights
      }
    }
    out$clusters <- data.table::setDT(out$clusters)
    return(out)
  } else {
    stop("Not implemented.")
  }
}

#' Scoring of dimensionality reduction and clustering pipeline output
#'
#' Computes averages of metrics from pipeline output and also returns the
#' best combination based on a weighted sum of metrics.
#'
#' Metrics are renamed for convenience: 
#' \itemize{
#'   \item [Train/Test]Stability[Jaccard/ARI/NMI]
#'   \item [NMI/ARI/ChisqRR].<batch>
#'   \item [NMI/ARI].<subtype>
#'   \item ...
#' }
#'
#' @param res \code{\link{COPS}} output
#' @param by character vector containing column names to group analysis by
#' @param wsum an expression that indicates how a combined score is computed 
#' @param significance_level p-value cutoff for computing rejection rates
#' @param summarise If FALSE, adds \code{"run"} and \code{"fold"} to \code{by}. 
#'   By default the metrics are averaged across runs and folds. 
#' @param format_names If TRUE, formats internally used method names etc. to more 
#'   user friendly names.
#'
#' @return Returns a \code{list} containing a \code{data.frame} \code{$all} of all scores and
#'         a single row \code{$best} with the best score according to \code{wsum}.
#' @export
#' 
#' @examples library(COPS)
#' library(parallel)
#' 
#' res <- COPS(ad_ge_micro_zscore, 
#' association_data = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap", "tsne"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' metric = "euclidean",
#' n_clusters = 2:4)
#' 
#' scores <- clusteval_scoring(res, wsum = Silhouette - GSE.nmi, summarise = TRUE)
#' 
#' best <- get_best_result(res, scores)
#' head(best$embedding)
#' head(best$clusters)
#' 
#' @importFrom plyr join ddply
#' @importFrom reshape2 melt dcast
#' @importFrom stats sd as.formula
clusteval_scoring <- function(res,
                              by = c("datname", "drname", "k", "m", "mkkm_mr_lambda"),
                              wsum = TrainStabilityJaccard + Silhouette,
                              significance_level = 0.05,
                              summarise = TRUE, 
                              format_names = TRUE) {
  uf <- unique(res$clusters[["fold"]])
  if (summarise == FALSE) {
    if (!"run" %in% by) {
      by <- c(by, "run")
    }
    if (!"fold" %in% by) {
      by <- c(by, "fold")
    }
    non_reference_fold <- uf
  } else {
    non_reference_fold <- uf[uf %in% unique(res$clusters[["cv_index"]])]
  }
  
  # Internal metrics
  if (!is.null(res$internal_metrics)) {
    by_internal <- by[by %in% colnames(res$internal_metrics)]
    mean_internals <- plyr::ddply(res$internal_metrics[res$internal_metrics[["fold"]] %in% non_reference_fold,], 
                                  c(by_internal, "metric"), 
                                  function(x) data.frame(mean = mean(x$value, na.rm = TRUE), 
                                                         sd = sd(x$value, na.rm = TRUE)))
    mean_internals <- reshape2::dcast(mean_internals, 
                                      as.formula(paste(paste(by_internal, collapse = "+"), "~ metric")), 
                                      value.var = "mean")
  } else {
    mean_internals <- NULL
  }
  
  # Batch label chi-squared tests (deprecated)
  if (!is.null(res$chisq_pval)) {
    # Average number of chi-squared test pvalues under threshold in all given labels
    by_chisq <- by[by %in% colnames(res$chisq_pval)]
    if (summarise) {
      chisq_f <- function(x) data.frame(chisqRR = mean(x$p < significance_level, na.rm = TRUE))
      chisq_rr <- plyr::ddply(res$chisq_pval[res$chisq_pval[["fold"]] %in% non_reference_fold,], 
                              c(by_chisq, "label"), chisq_f)
      chisq_rr$label <- paste0("ChisqRR.", chisq_rr$label)
      chisq_rr <- reshape2::dcast(chisq_rr, 
                                  as.formula(paste(paste(by_chisq, collapse = "+"), "~ label")), 
                                  value.var = "chisqRR")
    } else {
      chisq_f <- function(x) data.frame(chisq.p = mean(x$p, na.rm = TRUE))
      chisq_rr <- plyr::ddply(res$chisq_pval[res$chisq_pval[["fold"]] %in% non_reference_fold,], 
                              c(by_chisq, "label"), chisq_f)
      chisq_rr$label <- paste0("chisq.p.", chisq_rr$label)
      chisq_rr <- reshape2::dcast(chisq_rr, 
                                  as.formula(paste(paste(by_chisq, collapse = "+"), "~ label")), 
                                  value.var = "chisq.p")
    }
  } else {
    chisq_rr <- NULL
  }
  
  # Batch label associations (deprecated)
  if (!is.null(res$batch_association)) {
    if (nrow(res$batch_association) > 0) {
      by_bassoc <- by[by %in% colnames(res$batch_association)]
      bassoc <- plyr::ddply(res$batch_association[res$batch_association[["fold"]] %in% non_reference_fold,], 
                              c(by_bassoc, "label"), 
                              function(x) data.frame(NMI = mean(x$nmi, na.rm = TRUE), ARI = mean(x$ari, na.rm = TRUE)))
      bassoc_nmi <- reshape2::dcast(bassoc, 
                                    as.formula(paste(paste(by_bassoc, collapse = "+"), "~ label")), 
                                    value.var = "NMI")
      colnames(bassoc_nmi)[!colnames(bassoc_nmi) %in% by_bassoc] <- paste0("NMI.", colnames(bassoc_nmi)[!colnames(bassoc_nmi) %in% by_bassoc])
      bassoc_ari <- reshape2::dcast(bassoc, 
                                    as.formula(paste(paste(by_bassoc, collapse = "+"), "~ label")), 
                                    value.var = "ARI")
      colnames(bassoc_ari)[!colnames(bassoc_ari) %in% by_bassoc] <- paste0("ARI.", colnames(bassoc_ari)[!colnames(bassoc_ari) %in% by_bassoc])
    } else {
      bassoc_nmi <- NULL
      bassoc_ari <- NULL
    }
  } else {
    bassoc_nmi <- NULL
    bassoc_ari <- NULL
  }
  
  # Subtype label associations (deprecated)
  if (!is.null(res$subtype_association)) {
    if (nrow(res$subtype_association) > 0) {
      by_sassoc <- by[by %in% colnames(res$subtype_association)]
      sassoc <- plyr::ddply(res$subtype_association[res$subtype_association[["fold"]] %in% non_reference_fold,], 
                            c(by_sassoc, "label"), 
                            function(x) data.frame(NMI = mean(x$nmi, na.rm = TRUE), ARI = mean(x$ari, na.rm = TRUE)))
      sassoc_nmi <- reshape2::dcast(sassoc, 
                                    as.formula(paste(paste(by_sassoc, collapse = "+"), "~ label")), 
                                    value.var = "NMI")
      colnames(sassoc_nmi)[!colnames(sassoc_nmi) %in% by_sassoc] <- paste0("NMI.", colnames(sassoc_nmi)[!colnames(sassoc_nmi) %in% by_sassoc])
      sassoc_ari <- reshape2::dcast(sassoc, 
                                    as.formula(paste(paste(by_sassoc, collapse = "+"), "~ label")), 
                                    value.var = "ARI")
      colnames(sassoc_ari)[!colnames(sassoc_ari) %in% by_sassoc] <- paste0("ARI.", colnames(sassoc_ari)[!colnames(sassoc_ari) %in% by_sassoc])
    } else {
      sassoc_nmi <- NULL
      sassoc_ari <- NULL
    }
  } else {
    sassoc_nmi <- NULL
    sassoc_ari <- NULL
  }
  
  # Survival likelihood ratio test
  if (!is.null(res$survival)) {
    if (summarise) {
      # TODO: p-value like handling, c-index easier?
      warning("Currently summary of survival p-values is not implemented.")
      survival <- res$survival[res$survival[["fold"]] %in% non_reference_fold,]
      colnames(survival)[colnames(survival) == "cluster_significance"] <- "SurvivalPValue"
      
    }
    #by_survival <- by[by %in% colnames(res$survival)]
    #survival <- plyr::ddply(res$survival, 
    #                        by_survival, 
    #                        function(x) data.frame(SurvivalPValue = mean(x$cluster_significance, na.rm = TRUE)))
    survival <- res$survival[res$survival[["fold"]] %in% non_reference_fold,]
    colnames(survival)[colnames(survival) == "cluster_significance"] <- "SurvivalPValue"
  } else {
    survival <- NULL
  }
  
  # Module correlation scores
  if (!is.null(res$modules)) {
    by_modules <- by[by %in% colnames(res$modules)]
    modules <- plyr::ddply(res$modules[res$modules[["fold"]] %in% non_reference_fold,], 
                            by_modules, 
                            function(x) data.frame(Module_score = mean(x$Module_score, na.rm = TRUE)))
  } else {
    modules <- NULL
  }
  
  # Variable association tests
  if (!is.null(res$association)) {
    if (summarise) {
      assoc_string <- "\\.nmi$|\\.ari$"
    } else {
      assoc_string <- "\\.nmi$|\\.ari$|\\.p$"
    }
    by_association <- by[by %in% colnames(res$association)]
    association <- plyr::ddply(res$association[res$association[["fold"]] %in% non_reference_fold,], 
                               by_association, 
                               function(x) as.data.frame(lapply(x[grepl(assoc_string, colnames(x))], mean, na.rm = TRUE)))
    if (summarise) {
      association_rr <- plyr::ddply(res$association[res$association[["fold"]] %in% non_reference_fold,], 
                                    by_association, 
                                    function(x) as.data.frame(lapply(x[grepl("\\.p$", colnames(x))], function(x) mean(x < significance_level, na.rm = TRUE))))
      colnames(association_rr) <- gsub("\\.p$", ".rr", colnames(association_rr))
      association <- plyr::join(association, association_rr, by = by_association)
    }
  } else {
    association <- NULL
  }
  
  # Cluster sizes
  if (!is.null(res$cluster_sizes)) {
    by_cluster_sizes <- by[by %in% colnames(res$cluster_sizes)]
    cluster_sizes <- res$cluster_sizes[res$cluster_sizes[["fold"]] %in% non_reference_fold,]
    cluster_cols <- grep("X[0-9]+", colnames(cluster_sizes))
    csize_names <- paste0("Cluster_", gsub("^X", "", colnames(cluster_sizes))[cluster_cols], "_size")
    colnames(cluster_sizes)[cluster_cols] <- csize_names
    smallest_cluster_size <- apply(as.data.frame(cluster_sizes)[,csize_names], 1, min, na.rm = TRUE)
    cluster_sizes[["Smallest_cluster_size"]] <- smallest_cluster_size
    cluster_sizes <- plyr::ddply(cluster_sizes, 
                                 by_cluster_sizes, 
                                 function(x) as.data.frame(lapply(x[c(csize_names, "Smallest_cluster_size")], mean, na.rm = TRUE)))
  } else {
    cluster_sizes <- NULL
  }
  
  if (!is.null(res$stability)) {
    # Stability
    by_stability <- by[by %in% colnames(res$stability)]
    stability <- res$stability[res$stability[["fold"]] %in% non_reference_fold,]
    stab_col_ind <- match(c("train_jsc", "train_nmi", "train_ari", "test_jsc", "test_nmi", "test_ari"), colnames(stability))
    colnames(stability)[stab_col_ind[!is.na(stab_col_ind)]] <- c("TrainStabilityJaccard", "TrainStabilityNMI", "TrainStabilityARI",
                                           "TestStabilityJaccard", "TestStabilityNMI", "TestStabilityARI")[!is.na(stab_col_ind)]
    stability <- plyr::ddply(stability, 
                             by_stability, 
                             function(x) as.data.frame(lapply(x[grepl("Stability", colnames(x))], mean, na.rm = TRUE)))
  } else {
    stability <- NULL
  }
  
  # Combine all metrics
  out <- list(mean_internals, 
              chisq_rr, 
              bassoc_nmi, 
              bassoc_ari, 
              sassoc_nmi, 
              sassoc_ari, 
              survival, 
              modules,
              association,
              cluster_sizes,
              stability)
  if (!is.null(res$mkkm_mr_weights) & !summarise) {
    out$mkkm_mr_weights <- res$mkkm_mr_weights[res$mkkm_mr_weights[["fold"]] %in% non_reference_fold,]
  }
  #out <- Reduce(plyr::join, out[!sapply(out, is.null)])
  out <- Reduce(function(x,y) plyr::join(x, y, 
                                         by = intersect(by, 
                                                        intersect(colnames(x), 
                                                                  colnames(y))), 
                                         type = "full"), 
                out[!sapply(out, is.null)])
  if(format_names) {
    out <- format_scores(out)
  }
  
  # Scoring
  score <- substitute(wsum)
  
  out$wsum <- eval(score, out)
  best <- out[which.max(out$wsum),]
  
  if ("run" %in% colnames(out)) out[["run"]] <- factor(out[["run"]])
  if ("fold" %in% colnames(out)) out[["fold"]] <- factor(out[["fold"]])
  if ("k" %in% colnames(out)) out[["k"]] <- factor(out[["k"]])
  
  out <- list(all = out[order(out$wsum, decreasing = TRUE),], 
              best = best)
  
  return(out)
}

#aov(value ~ run * fold, method_evaluations$internal_metrics[method_evaluations$internal_metrics$k == 2 & method_evaluations$internal_metrics$m == "hierarchical" & method_evaluations$internal_metrics$drname == "pca4" & method_evaluations$internal_metrics$metric == "Silhouette",])

# Helper for pipeline verbosity
time_taken_string <- function(start) {
  time <- Sys.time() - start
  paste(sprintf("%.2f", time), attributes(time)$units)
}
