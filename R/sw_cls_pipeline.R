#' Pipeline combining dimensionality reduction, clustering, cross-validation and evaluation
#'
#' Combines \code{\link{cv_fold}}, \code{\link{cv_dimred}}, \code{\link{cv_clusteval}} and
#' \code{\link{stability_eval}}.
#'
#' @param dat_list list of pre-processed data sets to run pipeline for
#' @param nfolds number of cross-validation folds
#' @param nruns number of cross-validation replicates
#' @param batch_label vector or matrix with categorical variables on columns
#' @param verbose if \code{TRUE}, prints progress notifications
#' @param parallel sets up and registers \code{parallel} number of threads for supported operations
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
dimred_clusteval_pipeline <- function(dat_list, 
                                      nfolds, 
                                      nruns, 
                                      batch_label = NULL,
                                      verbose = TRUE,
                                      parallel = 1,
                                      #dim_reduction = TRUE,
                                      #pre_clust_cv = FALSE,
                                      ...) {
  pipeline_start <- Sys.time()
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    # Avoid warnings related to %dopar%
    foreach::registerDoSEQ()
  }
  if(verbose) print("Processing data sets ..."); flush.console()
  
  #if (is.null(names(batch_label))) names(batch_label) <- colnames(dat_list[[1]])
  # Manage batch label column names if batch labels given
  batch_label_names <- NULL
  if (!is.null(batch_label)) {
    batch_label_names <- colnames(batch_label)
    if (is.null(batch_label_names)) batch_label_names <- "batch_label"
  }
  
  # Convert data to data.table to optimize memory usage
  for (i in 1:length(dat_list)) {
    id <- colnames(dat_list[[i]])
    
    dat_list[[i]] <- as.data.table(t(dat_list[[i]]))
    #data.table::setDT(dat_list[[i]])
    colnames(dat_list[[i]]) <- paste0("dim", 1:ncol(dat_list[[i]]))
    dat_list[[i]]$id <- id
    setkey(dat_list[[i]], id)
    
    # Add batch label(s) as separate column(s)
    if (!is.null(batch_label)) {
      if (is.null(dim(batch_label))) {
        if (!is.null(names(batch_label))) {
          batch_id <- names(batch_label)
        } else {
          # Assume data and batch labels are in the same order
          batch_id <- id
        }
        batch_table <- cbind(batch_label = as.character(batch_label), c())
        #rownames(batch_label) <- batch_id
      } else { # matrix batch_label
        batch_id <- rownames(batch_label)
        batch_table <- batch_label
      }
      #data.table::setDT(batch_table)
      batch_table <- as.data.table(batch_table)
      if (!all(id %in% batch_id)) {
        stop("All batch label sample IDs do not match with data.")
      }
      batch_table$id <- batch_id
      dat_list[[i]] <- merge(dat_list[[i]], batch_table, by = "id")
    }
  }
  
  # Create cross validation folds
  cv_index <- cv_fold(dat_list = dat_list, nfolds = nfolds, nruns = nruns, ..., batch_label = batch_label)
  
  # Dimensionality reduction
  dimred_start <- Sys.time()
  if(verbose) print("Starting dimensionality reduction ..."); flush.console()
  dat_embedded <- cv_dimred(dat_list, cv_index, ...)
  if(verbose) print(paste("Finished dimensionality reduction in",
                          time_taken_string(dimred_start))); flush.console()
  
  # Clustering evaluation
  clusteval_start <- Sys.time()
  if(verbose) print("Starting clustering analysis ..."); flush.console()
  dat_clustered <- cv_clusteval(dat_embedded, ..., batch_label_names = batch_label_names)
  if(verbose) print(paste("Finished clustering analysis in",
                          time_taken_string(clusteval_start))); flush.console()
  
  # Clustering stability evaluation
  stability_test_start <- Sys.time()
  if(verbose) print("Starting clustering stability analysis ..."); flush.console()
  dat_stability <- stability_eval(dat_clustered$clusters, ...)
  if(verbose) print(paste("Finished clustering stability analysis in",
                          time_taken_string(stability_test_start))); flush.console()
  
  # Return
  out <- list(embedding = dat_embedded, 
              clusters = dat_clustered$clusters, 
              internal_metrics = dat_clustered$metrics,
              chisq_pval = dat_clustered$chisq_pval,
              stability = dat_stability)
  
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}

#' Cross-validation fold permutation
#' 
#' Creates cross-validation folds of data for downstream analysis. 
#'
#' @param dat_list list of data matrices with samples on columns
#' @param nfolds number of cross-validation folds
#' @param nruns number of cross-validation replicates
#' @param batch_label batch_label vector or \code{data.frame} with column \code{"batch_label"}
#' @param stratified_cv if \code{TRUE}, try to maximize separation of batch labels within folds
#' @param mixed_cv if \code{TRUE}, try to minimize separation of batch labels within folds
#' @param ... extra arguments are ignored
#'
#' @return list of data.frames with added columns "fold", "run" and "cv_index" as well as 
#'         duplicated rows of the original data corresponding to different folds.
#' @export
#'
#' @importFrom plyr join
cv_fold <- function(dat_list, 
                    nfolds = 5, 
                    nruns = 2, 
                    #batch_label = NULL, 
                    stratified_cv = FALSE, 
                    mixed_cv = FALSE,
                    ...) {
  out <- list()
  for (i in 1:length(dat_list)) {
    folded <- list()
    for (j in 1:nruns) {
      if (!is.null(dat_list[[i]]$batch_label) & (stratified_cv | mixed_cv)) {
        a_ind <- lapply(table(dat_list[[i]]$batch_label), function(x) sample(1:x, x))
        b_ind <- sample(1:length(unique(dat_list[[i]]$batch_label)), length(unique(dat_list[[i]]$batch_label)))
        c_ind <- cumsum(table(dat_list[[i]]$batch_label)[unique(dat_list[[i]]$batch_label)[b_ind]])
        cv_index <- c()
        for (u in 1:length(b_ind)) {
          un <- unique(dat_list[[i]]$batch_label)[b_ind[u]]
          cv_index[dat_list[[i]]$batch_label == un] <- a_ind[[un]] + ifelse(u > 1, c_ind[u-1], 0)
        }
        if (stratified_cv) {
          # Stratified cv folds such that holdout set labels mostly do not match to rest of data
          cv_index <- cv_index %/% -(length(dat_list[[i]]$batch_label) %/% -nfolds) + 1
        } else {
          # Mixed cv folds such that labels are evenly distributed within folds
          cv_index <- cv_index %% nfolds + 1
        }
      } else {
        # Completely random folds
        cv_index <- sample(1:nrow(dat_list[[i]])) %% nfolds + 1
      }
      # Got index, create folds +1 extra "fold" with whole data
      # TODO: reference is the same accross all runs, maybe include it only once? 
      #       Possible incompatibility with external methods.
      folded[[j]] <- list()
      for (f in 1:(nfolds+1)) {
        # TODO: fix downstream support so that test set can be included too
        #tempfold <- dat_list[[i]][cv_index != f, ]
        #tempfold$fold <- f
        #tempfold$run <- j
        #tempfold$cv_index <- cv_index[cv_index != f]
        folded[[j]][[f]] <- data.table(fold = f, run = j, 
                                       cv_index = cv_index[cv_index != f], 
                                       id = dat_list[[i]]$id[cv_index != f])
      }
      folded[[j]] <- data.table::rbindlist(folded[[j]])
    }
    out[[i]] <- data.table::rbindlist(folded)
  }
  names(out) <- names(dat_list)
  return(out)
}

#' Reduce dimensionality and cluster input data
#'
#' Convenience function that runs a specific setting of the dimensionality reduction
#' and clustering analysis pipeline. Accepts \code{\link{clusteval_scoring}} output
#' \code{$best} as input.
#'
#' @param dat_list \code{list} of pre-processed data sets, original input to pipeline or
#'                 containing element \code{data_id} corresponding to the preferred input
#' @param best a single row \code{data.frame} which specifies the preferred methods
#' @param data_id name or index of preferred element of \code{dat_list}
#' @param dimred_method name of dimensionality reduction method for
#'                      \code{\link{dim_reduction_suite}}
#' @param dimred_dim number of dimensions parameter for \code{\link{dim_reduction_suite}}
#' @param dimred_perp perplexity parameter for \code{\link{dim_reduction_suite}}
#' @param cluster_method clustering method name
#' @param k number of clusters parameter
#' @param metric distance metric parameter for clustering
#'
#' @return Returns a \code{list} with \code{$embedding} corresponding to the embedding
#'         generated using \code{\link{dim_reduction_suite}} and \code{$clustering}
#'         corresponding to clustering on \code{$embedding} using
#'         \code{\link[clValid]{clValid}}.
#' @export
#' @importFrom clValid clValid clusters
#' @importFrom stats cutree
dimred_cluster <- function(dat_list,
                            best = NULL,
                            data_id = 1,
                            dimred_method = "pca",
                            dimred_dim = 2,
                            dimred_perp = 30,
                            cluster_method = "kmeans",
                            k = 2,
                            metric = "euclidean") {
  if (!is.null(best)) {
    dat <- dat_list[[best$datname]]
    
    # Get dimensionality reduction method from name
    temp <- regexpr("^[A-Za-z]*", best$drname)
    dimred_method <- substr(best$drname, temp,
                            temp + attributes(temp)$match.length - 1)
    # Get number of dimensions from name
    temp <- regexpr("[0-9]*$", best$drname)
    dimred_dim <- substr(best$drname, temp,
                         temp + attributes(temp)$match.length - 1)
    dimred_dim <- as.numeric(dimred_dim)
    dimred_perp <- dimred_dim
    
    cluster_method <- as.character(best$m)
    k <- best$k
  } else {
    if (is.list(dat_list)) {
      dat <- dat_list[[data_id]]
    } else if (is.matrix(dat_list)|is.data.frame(dat_list)) {
      dat <- dat_list
    } else {
      stop("Unsupported data format")
    }
  }
  if (dimred_method == "original") {
    embed <- dat
  } else {
    embed <- dim_reduction_suite(t(dat),
                                 dimred_methods = dimred_method,
                                 output_dimensions = dimred_dim,
                                 tsne_perplexities = dimred_perp,
                                 include_original = FALSE)[[1]]
  }
  
  out <- suppressWarnings(clValid::clValid(embed,
                                           nClust = k,
                                           clMethods = cluster_method,
                                           metric = metric,
                                           validation=c("internal"),
                                           maxitems = Inf))
  
  # TODO: create generic unpacker for clValid output or use clustering methods directly
  temp <- clValid::clusters(out, cluster_method)
  
  if (cluster_method %in% c("hierarchical", "diana", "agnes")) {
    cluster <- cutree(temp, k)
  } else if (cluster_method == "pam") {
    cluster <- temp[[1]]$clustering
  } else if (cluster_method == "kmeans") {
    cluster <- temp[[1]]$cluster
  } else if (cluster_method == "sota") {
    cluster <- temp[[1]]$clust
  } else {
    stop(paste("Unsupported clustering method:", as.character(cluster_method)))
  }
  
  return(list(embedding = embed, clustering = cluster))
}

#' Scoring of dimensionality reduction and clustering pipeline output
#'
#' Computes averages of metrics given in pipeline output and also returns the
#' best combination based on a weighted sum of metrics.
#'
#' Default weights have
#' been chosen to emphasize stability and  presence of batch effect. The
#' internal metrics have wildly varying scales and thus the weights will
#' need to be guessed manually almost always.
#'
#' @param input \code{\link{dimred_clusteval_pipeline}} output
#' @param by character vector containing column names to group analysis by
#' @param chisq_significance_level p-value cutoff for computing rejection rate of
#'                                 \code{chisq.test}
#' @param w_connectivity weight for \code{\link[clValid]{connectivity}} score
#' @param w_dunn weight for \code{\link[clValid]{dunn}} score
#' @param w_silhouette weight for \code{silhouette} score
#' @param w_stability weight for \code{\link{stability_eval}} score
#' @param w_batch_effect weight for \code{chisq.test} rejection rate,
#'                       ignored unless batch_label was supplied as vector or
#'                       contained column named "batch_label"
#'
#' @return Returns a \code{list} containing a \code{data.frame} \code{$all} of all scores and
#'         a single row \code{$best} with the best score
#' @export
#' @importFrom plyr join ddply
#' @importFrom reshape2 melt
#' @importFrom stats sd
clusteval_scoring <- function(input,
                              by = c("datname", "drname", "k", "m"),
                              chisq_significance_level = 0.05,
                              w_connectivity = -0.1,
                              w_dunn = 0.1,
                              w_silhouette = 0.5,
                              w_stability = 1.0,
                              w_batch_effect = -1.0) {
  out <- list()
  
  by_internal <- by[by %in% colnames(input$internal_metrics)]
  out[[1]] <- plyr::ddply(input$internal_metrics, 
                          c(by_internal, "metric"), 
                          function(x) data.frame(mean = mean(x$value), 
                                                 sd = sd(x$value)))
  
  # Average number of chi-squared test pvalues under threshold in all given labels
  by_chisq <- by[by %in% colnames(input$chisq_pval)]
  out[[2]] <- plyr::ddply(input$chisq_pval, 
                          c(by_chisq), 
                          function(x) data.frame(metric = "ChiSqRR_batch_label",
                                                 mean = mean(x$p < chisq_significance_level),
                                                 sd = NA))
  
  by_stability <- by[by %in% colnames(input$stability)]
  out[[3]] <- input$stability[by_stability]
  out[[3]]$metric <- "JaccardDist"
  out[[3]]$mean <- input$stability$jdist_train
  out[[3]]$sd <- NA
  
  out <- Reduce(plyr::rbind.fill, out)
  
  ## Rank based on weighted sum of metrics
  # Connectivity [0,Inf] minimize
  # Dunn [0,Inf] maximize
  # Silhouette [-1,1] maximize
  # JaccardDist [0,1] minimize
  # ChiSqRR [0,1] minimize
  
  wsumtable <- data.frame(metric = c("Connectivity", "Dunn", "Silhouette",
                                     "JaccardDist", "ChiSqRR_batch_label"),
                          weight = c(w_connectivity, w_dunn, w_silhouette,
                                     -1*w_stability, w_batch_effect))
  
  wsumfun <- function(x) {
    temp <- invisible(plyr::join(x, wsumtable, by = "metric"))
    return(data.frame(wsum = sum(temp$mean * temp$weight), 
                      Connectivity = temp$mean[grep("Connectivity", temp$metric)], 
                      Dunn = temp$mean[grep("Dunn", temp$metric)], 
                      Silhouette = temp$mean[grep("Silhouette", temp$metric)], 
                      JaccardDist = temp$mean[grep("JaccardDist", temp$metric)],
                      ChiSqRR_batch_label = temp$mean[grep("ChiSqRR_batch_label", temp$metric)]))
  }
  
  by_all <- by[by %in% colnames(out)]
  scores <- plyr::ddply(out, by_all, wsumfun)
  
  best <- scores[which.max(scores$wsum),]
  
  return(list(all = scores[order(scores$wsum, decreasing = TRUE),], 
              best = best))
}

#aov(value ~ run * fold, method_evaluations$internal_metrics[method_evaluations$internal_metrics$k == 2 & method_evaluations$internal_metrics$m == "hierarchical" & method_evaluations$internal_metrics$drname == "pca4" & method_evaluations$internal_metrics$metric == "Silhouette",])

# Helper for pipeline verbosity
time_taken_string <- function(start) {
  time <- Sys.time() - start
  paste(sprintf("%.2f", time), attributes(time)$units)
}
