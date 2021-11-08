#' Pipeline combining dimensionality reduction, clustering, cross-validation and evaluation
#'
#' Combines \code{\link{cv_fold}}, \code{\link{cv_dimred}}, \code{\link{cv_clusteval}} and
#' \code{\link{stability_eval}}.
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
#' @examples library(parallel)
#' library(COPS)
#' 
#' # DR-CL
#' res <- dimred_clusteval_pipeline(ad_ge_micro_zscore, 
#' batch_label = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap", "tsne"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' metric = "euclidean", 
#' n_clusters = 2:4)
#' 
#' # CL
#' res <- dimred_clusteval_pipeline(ad_ge_micro_zscore, 
#' batch_label = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("none"), 
#' cluster_methods = c("hierarchical"), 
#' metric = "correlation", 
#' n_clusters = 2:4)
#' 
#' # BK-CL
#' res <- dimred_clusteval_pipeline(ad_ge_micro_zscore, 
#' batch_label = ad_studies, 
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
dimred_clusteval_pipeline <- function(dat, 
                                      nfolds, 
                                      nruns, 
                                      batch_label = NULL,
                                      subtype_label = NULL,
                                      clinical_data = NULL,
                                      survival_data = NULL,
                                      module_eigs = NULL,
                                      verbose = TRUE,
                                      parallel = 1,
                                      pathway_enrichment_method = "none",
                                      #dim_reduction = TRUE,
                                      #pre_clust_cv = FALSE,
                                      ...) {
  if ("list" %in% class(dat)) {
    dat_list <- dat
  } else {
    dat_list <- list(dat)
  }
  
  pipeline_start <- Sys.time()
  if(verbose) print("Processing data sets ..."); flush.console()
  
  #if (is.null(names(batch_label))) names(batch_label) <- colnames(dat_list[[1]])
  # Manage batch label column names if batch labels given
  batch_label_names <- NULL
  if (!is.null(batch_label)) {
    batch_label_names <- colnames(batch_label)
    if (is.null(batch_label_names)) batch_label_names <- "batch_label"
  }
  
  subtype_label_names <- NULL
  if (!is.null(subtype_label)) {
    subtype_label_names <- colnames(subtype_label)
    if (is.null(subtype_label_names)) subtype_label_names <- "subtype_label"
  }
  
  if (pathway_enrichment_method != "none") {
    # Collect gene names for later
    gene_id_list <- lapply(dat_list, rownames)
  }
  
  # Convert data to data.table to optimize memory usage
  for (i in 1:length(dat_list)) {
    
    id <- colnames(dat_list[[i]])
    
    dat_list[[i]] <- data.table::as.data.table(t(dat_list[[i]]))
    #data.table::setDT(dat_list[[i]])
    colnames(dat_list[[i]]) <- paste0("dim", 1:ncol(dat_list[[i]]))
    dat_list[[i]]$id <- id
    data.table::setkey(dat_list[[i]], id)
    
    # Add batch label(s) as separate column(s)
    if (!is.null(batch_label)) {
      if (is.null(dim(batch_label))) {
        # 1d array given as batch label
        if (!is.null(names(batch_label))) {
          batch_id <- names(batch_label)
        } else {
          # Assume data and batch labels are in the same order
          batch_id <- id
        }
        batch_table <- cbind(batch_label = as.character(batch_label), c())
      } else { 
        # 2d batch_label
        batch_id <- rownames(batch_label)
        batch_table <- batch_label
      }
      #data.table::setDT(batch_table)
      batch_table <- data.table::as.data.table(batch_table)
      if (!all(id %in% batch_id)) {
        stop("All sample IDs in data do not match with batch label sample IDs.")
      }
      batch_table$id <- batch_id
      dat_list[[i]] <- plyr::join(dat_list[[i]], batch_table, by = "id")
    }
    
    if (!is.null(subtype_label)) {
      if (is.null(dim(subtype_label))) {
        if (!is.null(names(subtype_label))) {
          subtype_id <- names(subtype_label)
        } else {
          # Assume data and labels are in the same order
          subtype_id <- id
        }
        subtype_table <- cbind(subtype_label = as.character(subtype_label), c())
      } else { # matrix batch_label
        subtype_id <- rownames(subtype_label)
        subtype_table <- subtype_label
      }
      #data.table::setDT(subtype_table)
      subtype_table <- data.table::as.data.table(subtype_table)
      if (!all(id %in% subtype_id)) {
        warning("All subtype label sample IDs do not match with data.")
      }
      subtype_table$id <- subtype_id
      dat_list[[i]] <- plyr::join(dat_list[[i]], subtype_table, by = "id")
    }
  }
  
  #if (vertical_parallelization)
  
  # Create cross validation folds
  cv_index <- cv_fold(dat_list = dat_list, 
                      nfolds = nfolds, 
                      nruns = nruns, 
                      ...)
  
  # Pathway enrichment
  if (pathway_enrichment_method != "none") {
    if (length(dat_list) > 1) {
      stop("Multiple data sets with pathway enrichment enabled is not implemented.")
      # TODO: make it compatible. Currently datname is being used as pwname. 
    }
    pw_start <- Sys.time()
    if(verbose) print("Starting pathway enrichment ..."); flush.console()
    dat_pw <- cv_pathway_enrichment(dat_list, 
                                    cv_index, 
                                    gene_id_list, 
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
  } else{
    # Dimensionality reduction
    dimred_start <- Sys.time()
    if(verbose) print("Starting dimensionality reduction ..."); flush.console()
    dat_embedded <- cv_dimred(dat_list, 
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
  
  if (is.null(clinical_data)) {
    if(!is.null(batch_label)) {
      clinical_data <- batch_table
    }
    if(!is.null(subtype_label)) {
      if (is.null(clinical_data)) {
        clinical_data <- subtype_table
      } else {
        clinical_data <- merge(clinical_data, subtype_table, by = "id")
      }
    }
    if (!is.null(clinical_data)) {
      clinical_data <- as.data.frame(clinical_data)
      rownames(clinical_data) <- clinical_data$id
      clinical_data$id <- NULL
    }
  }
  
  if (!is.null(clinical_data)) {
    clinical_analysis_start <- Sys.time()
    if(verbose) print("Starting clinical variable association analysis ..."); flush.console()
    dat_clinical_score <- clinical_evaluation(dat_clustered$clusters, 
                                              clinical_data, 
                                              parallel = parallel, 
                                              ...)
    if(verbose) print(paste("Finished clinical variable association analysis in",
                            time_taken_string(clinical_analysis_start))); flush.console()
  }
  
  # Return
  out <- list(embedding = dat_embedded, 
              clusters = dat_clustered$clusters, 
              internal_metrics = dat_clustered$metrics,
              stability = dat_stability)
  if (!is.null(survival_data)) out$survival <- dat_survival
  if (!is.null(module_eigs)) out$modules <- dat_gm_score
  if (!is.null(clinical_data)) out$clinical <- dat_clinical_score
  out$cluster_sizes <- dat_clustered$cluster_sizes
  
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}


#' @describeIn clusteval_scoring Retrieves best clustering from CV results based on scores. 
#' In practice retrieves reference fold result from first run matching the best results. 
#' 
#' @param res result from \code{\link{dimred_clusteval_pipeline}}
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

COPS <- function(dat, 
                 nfolds, 
                 nruns, 
                 clinical_data = NULL,
                 survival_data = NULL,
                 module_eigs = NULL,
                 verbose = TRUE,
                 parallel = 1,
                 pathway_enrichment_method = "none",
                 #dim_reduction = TRUE,
                 #pre_clust_cv = FALSE,
                 ...) {
  if ("list" %in% class(dat)) {
    dat_list <- dat
  } else {
    dat_list <- list(dat)
  }
  
  pipeline_start <- Sys.time()
  if(verbose) print("Processing data sets ..."); flush.console()
  
  if (pathway_enrichment_method != "none") {
    # Collect gene names for later
    gene_id_list <- lapply(dat_list, rownames)
  }
  
  # Convert data to data.table to optimize memory usage
  for (i in 1:length(dat_list)) {
    id <- colnames(dat_list[[i]])
    dat_list[[i]] <- data.table::as.data.table(t(dat_list[[i]]))
    #data.table::setDT(dat_list[[i]])
    colnames(dat_list[[i]]) <- paste0("dim", 1:ncol(dat_list[[i]]))
    dat_list[[i]]$id <- id
    data.table::setkey(dat_list[[i]], id)
  }
  
  #if (vertical_parallelization)
  
  # Create cross validation folds
  cv_index <- cv_fold(dat_list = dat_list, 
                      nfolds = nfolds, 
                      nruns = nruns, 
                      ...)
  
  # Pathway enrichment
  if (pathway_enrichment_method != "none") {
    if (length(dat_list) > 1) {
      stop("Multiple data sets with pathway enrichment enabled is not implemented.")
      # TODO: make it compatible. Currently datname is being used as pwname. 
    }
    pw_start <- Sys.time()
    if(verbose) print("Starting pathway enrichment ..."); flush.console()
    dat_pw <- cv_pathway_enrichment(dat_list, 
                                    cv_index, 
                                    gene_id_list, 
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
    dat_embedded <- cv_dimred(dat_list, 
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
  
  if (!is.null(clinical_data)) {
    clinical_analysis_start <- Sys.time()
    if(verbose) print("Starting clinical variable association analysis ..."); flush.console()
    dat_clinical_score <- clinical_evaluation(dat_clustered$clusters, 
                                              clinical_data, 
                                              parallel = parallel, 
                                              ...)
    if(verbose) print(paste("Finished clinical variable association analysis in",
                            time_taken_string(clinical_analysis_start))); flush.console()
  }
  
  # Return
  out <- list(embedding = dat_embedded, 
              clusters = dat_clustered$clusters, 
              internal_metrics = dat_clustered$metrics,
              stability = dat_stability)
  if (!is.null(survival_data)) out$survival <- dat_survival
  if (!is.null(module_eigs)) out$modules <- dat_gm_score
  if (!is.null(clinical_data)) out$clinical <- dat_clinical_score
  out$cluster_sizes <- dat_clustered$cluster_sizes
  
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}

vertical_pipeline <- function(dat_list, cv_ind, multi_omic = FALSE, ...) {
  parallel_clust <- setup_parallelization(parallel)
  #dat_list_fold <- list()
  #for (i in 1:length(dat_list) {
  #  dat_list_fold[[i]] <- merge(dat_list[[i]], cv_ind, by = "id")
  #}
  if (multi_omic) {
    if (!all(Reduce("&", lapply(dat_list[-1], function(x) x$id == dat_list[[1]]$id)))) {
      warning("Colnames in all views do not match.")
      stop("Cross-validation for missing sample views not implemented.")
    } else {
      cv_index <- cv_fold(dat_list = dat_list[1], nfolds = nfolds, nruns = nruns, ...)
      cv_index_split <- split(cv_index[[1]], by = c("run", "fold"))
      out <- foreach(i = 1:length(cv_index_split), 
                     .combine = c,
                     .export = c("dat_list", "cv_index"),
                     .packages = c("iClusterPlus")) %dopar% 
      {
        dat_i <- list()
        non_data_cols <- list()
        for (j in 1:length(dat_list)) {
          dat_i[[j]] <- merge(cv_index_split[[i]], dat_list[[j]])
          sel <- grep("^dim[0-9]+$", colnames(dat_i[[j]]))
          non_data_cols[[j]] <- dat_i[[j]][,-..sel]
          dat_i[[j]] <- as.matrix(dat_i[[j]][,..sel])
        }
        #multi
      }
    }
  } else {
    
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
#' @param res \code{\link{dimred_clusteval_pipeline}} output
#' @param by character vector containing column names to group analysis by
#' @param wsum an expression that indicates how a combined score is computed 
#' @param chisq_significance_level p-value cutoff for computing rejection rate of \code{chisq.test}
#' @param summarise If FALSE, adds \code{"run"} and \code{"fold"} to \code{by}. By default the metrics 
#'   are averaged across runs and folds. 
#'
#' @return Returns a \code{list} containing a \code{data.frame} \code{$all} of all scores and
#'         a single row \code{$best} with the best score according to \code{wsum}.
#' @export
#' 
#' @examples library(COPS)
#' library(parallel)
#' 
#' res <- dimred_clusteval_pipeline(ad_ge_micro_zscore, 
#' batch_label = ad_studies, 
#' parallel = 2, nruns = 2, nfolds = 5, 
#' dimred_methods = c("pca", "umap", "tsne"), 
#' cluster_methods = c("hierarchical", "kmeans"), 
#' metric = "euclidean",
#' n_clusters = 2:4)
#' 
#' scores <- clusteval_scoring(res, wsum = Silhouette - batch_label.nmi, summarise = TRUE)
#' 
#' best <- get_best_result(res, scores)
#' head(best$embedding)
#' head(best$clusters)
#' 
#' @importFrom plyr join ddply
#' @importFrom reshape2 melt dcast
#' @importFrom stats sd as.formula
clusteval_scoring <- function(res,
                              by = c("datname", "drname", "k", "m"),
                              wsum = TrainStabilityJaccard + Silhouette,
                              chisq_significance_level = 0.05,
                              summarise = TRUE) {
  if (summarise == FALSE) {
    by <- c(by, "run", "fold")
  }
  # Internal metrics
  if (!is.null(res$internal_metrics)) {
    by_internal <- by[by %in% colnames(res$internal_metrics)]
    mean_internals <- plyr::ddply(res$internal_metrics, 
                                  c(by_internal, "metric"), 
                                  function(x) data.frame(mean = mean(x$value, na.rm = TRUE), 
                                                         sd = sd(x$value, na.rm = TRUE)))
    mean_internals <- reshape2::dcast(mean_internals, 
                                      as.formula(paste(paste(by_internal, collapse = "+"), "~ metric")), 
                                      value.var = "mean")
  } else {
    mean_internals <- NULL
  }
  
  if (!is.null(res$chisq_pval)) {
    # Average number of chi-squared test pvalues under threshold in all given labels
    by_chisq <- by[by %in% colnames(res$chisq_pval)]
    if (summarise) {
      chisq_f <- function(x) data.frame(chisqRR = mean(x$p < chisq_significance_level, na.rm = TRUE))
      chisq_rr <- plyr::ddply(res$chisq_pval, c(by_chisq, "label"), chisq_f)
      chisq_rr$label <- paste0("ChisqRR.", chisq_rr$label)
      chisq_rr <- reshape2::dcast(chisq_rr, 
                                  as.formula(paste(paste(by_chisq, collapse = "+"), "~ label")), 
                                  value.var = "chisqRR")
    } else {
      chisq_f <- function(x) data.frame(chisq.p = mean(x$p, na.rm = TRUE))
      chisq_rr <- plyr::ddply(res$chisq_pval, c(by_chisq, "label"), chisq_f)
      chisq_rr$label <- paste0("chisq.p.", chisq_rr$label)
      chisq_rr <- reshape2::dcast(chisq_rr, 
                                  as.formula(paste(paste(by_chisq, collapse = "+"), "~ label")), 
                                  value.var = "chisq.p")
    }
  } else {
    chisq_rr <- NULL
  }
  
  # Batch label associations
  if (!is.null(res$batch_association)) {
    if (nrow(res$batch_association) > 0) {
      by_bassoc <- by[by %in% colnames(res$batch_association)]
      bassoc <- plyr::ddply(res$batch_association, 
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
  
  # Subtype label associations
  if (!is.null(res$subtype_association)) {
    if (nrow(res$subtype_association) > 0) {
      by_sassoc <- by[by %in% colnames(res$subtype_association)]
      sassoc <- plyr::ddply(res$subtype_association, 
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
  
  if (!is.null(res$stability)) {
    # Stability
    #by_stability <- by[by %in% colnames(res$stability)]
    stability <- res$stability
    stab_col_ind <- match(c("train_jsc", "train_nmi", "train_ari", "test_jsc", "test_nmi", "test_ari"), colnames(stability))
    colnames(stability)[stab_col_ind] <- c("TrainStabilityJaccard", "TrainStabilityNMI", "TrainStabilityARI",
                                           "TestStabilityJaccard", "TestStabilityNMI", "TestStabilityARI")
  } else {
    stability <- NULL
  }
  
  # Survival likelihood ratio test
  if (!is.null(res$survival)) {
    if (summarise) {
      warning("Currently summary of survival p-values is not implemented.")
    }
    #by_survival <- by[by %in% colnames(res$survival)]
    #survival <- plyr::ddply(res$survival, 
    #                        by_survival, 
    #                        function(x) data.frame(SurvivalPValue = mean(x$cluster_significance, na.rm = TRUE)))
    survival <- res$survival
    colnames(survival)[colnames(survival) == "cluster_significance"] <- "SurvivalPValue"
  } else {
    survival <- NULL
  }
  
  # Module correlation scores
  if (!is.null(res$modules)) {
    by_modules <- by[by %in% colnames(res$modules)]
    modules <- plyr::ddply(res$modules, 
                            by_modules, 
                            function(x) data.frame(Module_score = mean(x$Module_score, na.rm = TRUE)))
  } else {
    modules <- NULL
  }
  
  # Clinical variable association tests
  if (!is.null(res$clinical)) {
    if (summarise) {
      warning("Currently summary of clinical p-values is not implemented.")
    }
    clinical <- res$clinical
  } else {
    clinical <- NULL
  }
  
  if (!is.null(res$cluster_sizes)) {
    cluster_sizes <- res$cluster_sizes
    cluster_cols <- which(!is.na(suppressWarnings(as.numeric(colnames(cluster_sizes)))))
    csize_names <- paste0("Cluster_", colnames(cluster_sizes)[cluster_cols], "_size")
    colnames(cluster_sizes)[cluster_cols] <- csize_names
  } else {
    cluster_sizes <- NULL
  }
  
  # Combine all metrics
  out <- list(mean_internals, 
              chisq_rr, 
              bassoc_nmi, 
              bassoc_ari, 
              sassoc_nmi, 
              sassoc_ari, 
              stability, 
              survival, 
              modules,
              clinical,
              cluster_sizes)
  #out <- Reduce(plyr::join, out[!sapply(out, is.null)])
  out <- Reduce(function(x,y) plyr::join(x, y, by = intersect(by, intersect(colnames(x), colnames(y)))), 
                out[!sapply(out, is.null)])
  
  # Scoring
  score <- substitute(wsum)
  
  out$wsum <- eval(score, out)
  best <- out[which.max(out$wsum),]
  
  return(list(all = out[order(out$wsum, decreasing = TRUE),], 
              best = best))
}

#aov(value ~ run * fold, method_evaluations$internal_metrics[method_evaluations$internal_metrics$k == 2 & method_evaluations$internal_metrics$m == "hierarchical" & method_evaluations$internal_metrics$drname == "pca4" & method_evaluations$internal_metrics$metric == "Silhouette",])

# Helper for pipeline verbosity
time_taken_string <- function(start) {
  time <- Sys.time() - start
  paste(sprintf("%.2f", time), attributes(time)$units)
}
