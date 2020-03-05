#' Pipeline combining dimensionality reduction, clustering, cross-validation and evaluation
#'
#' Combines \code{\link{dim_reduction_suite}}, \code{\link{clusteval_cv}} and
#' \code{\link{stability_eval}}.
#'
#' @param dat_list list of pre-processed data sets to run pipeline for
#' @param batch_label vector or matrix with categorical variables on columns
#' @param verbose if \code{TRUE}, prints progress notifications
#' @param parallel sets up and registers \code{parallel} number of threads for supported operations
#' @param dim_reduction set to \code{FALSE} to skip the dimensionality reduction
##' @param pre_clust_cv set to \code{TRUE} if data is already organized into folds from previous step
#' @param ... extra arguments are passed to pipeline components where appropriate
#'
#' @return Returns a \code{list} of pipeline component outputs for given settings and input data sets
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom utils flush.console
dimred_clusteval_pipeline <- function(dat_list,
                                      batch_label = NULL,
                                      verbose = TRUE,
                                      parallel = 1,
                                      dim_reduction = TRUE,
                                      pre_clust_cv = FALSE,
                                      ...) {
  pipeline_start <- Sys.time()
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    # Avoid warnings related to %dopar%
    foreach::registerDoSEQ()
  }
  out <- list()
  for (i in 1:length(dat_list)) {
    if(verbose) print(paste0("Starting pipeline for data set #", i)); flush.console()
    out_i <- list()
    if (dim_reduction) {
      dimred_start <- Sys.time()
      if(verbose) print(paste("Starting dimensionality reduction")); flush.console()
      out_i$embed_set <- dim_reduction_suite(dat_list[[i]], ...)
      if(verbose) print(paste("Finished dimensionality reduction in",
                              time_taken_string(dimred_start))); flush.console()
    } else {
      # INCOMPLETE
      # Pre computed embeddings for a specific setting?
      out_i$embed_set <- dat_list[[i]]
    }
    out_i$evals <-list()
    for (j in 1:length(out_i$embed_set)) {
      if(verbose) print(paste("Starting clustering and evaluation for embedding:",
                              names(out_i$embed_set)[j]));flush.console()
      out_j <- list()

      # Clustering + internal evaluation
      clusteval_start <- Sys.time()
      if (!pre_clust_cv) {
        out_j$embed_eval <- clusteval_cv(out_i$embed_set[[j]], batch_label, ...)
      } else {
        stop("pre_clust_cv is not implemented yet")
        # TODO: Run clustering for each predefined fold
      }
      if(verbose) print(paste("Finished clustering analysis in",
                              time_taken_string(clusteval_start)));flush.console()
      # Evaluate stability
      stability_test_start <- Sys.time()
      out_j$embed_stab <- stability_eval(out_j$embed_eval$clust, parallel = parallel > 1, ...)
      if(verbose) print(paste("Finished stability analysis in",
                              time_taken_string(stability_test_start)));flush.console()

      # Batch effect estimation
      if (!is.null(batch_label)) {
        batch_effect_start <- Sys.time()
        out_j$embed_baef <- class_associations(out_i$embed_set[[j]], batch_label, ...)
        if(verbose) print(paste("Finished batch effect estimation in",
                                time_taken_string(batch_effect_start)));flush.console()
      }
      if (is.null(names(out_i$embed_set)))  {
        out_i$evals[[j]] <- out_j
      } else {
        out_i$evals[[names(out_i$embed_set)[j]]] <- out_j
      }
    }
    if (is.null(names(dat_list)))  {
      out[[i]] <- out_i
    } else {
      out[[names(dat_list)[i]]] <- out_i
    }
  }
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}

dimred_clusteval_pipeline2 <- function(dat_list, 
                                      nfolds, 
                                      nruns, 
                                      batch_label = NULL,
                                      verbose = TRUE,
                                      parallel = 1,
                                      dim_reduction = TRUE,
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
  
  # Create cross validation folds
  dat_folded <- cv_fold(dat_list, nfolds, nruns, ..., batch_label = batch_label)
  
  # Dimensionality reduction
  dimred_start <- Sys.time()
  if(verbose) print("Starting dimensionality reduction ..."); flush.console()
  dat_folded <- cv_dimred(dat_folded, ...)
  if(verbose) print(paste("Finished dimensionality reduction in",
                          time_taken_string(dimred_start))); flush.console()
  
  # Clustering evaluation
  clusteval_start <- Sys.time()
  if(verbose) print("Starting clustering analysis ..."); flush.console()
  dat_clustered <- cv_clusteval(dat_folded, ..., batch_label_names = batch_label_names)
  if(verbose) print(paste("Finished clustering analysis in",
                          time_taken_string(clusteval_start))); flush.console()
  
  # Clustering stability evaluation
  stability_test_start <- Sys.time()
  if(verbose) print("Starting clustering stability analysis ..."); flush.console()
  dat_stability <- stability_eval2(dat_clustered$clusters, ...)
  if(verbose) print(paste("Finished clustering stability analysis in",
                          time_taken_string(stability_test_start))); flush.console()
  
  # Return
  out <- list(embedding = dat_folded, 
              clusters = dat_clustered$clusters, 
              internal_metrics = dat_clustered$metrics,
              chisq_pval = dat_clustered$chisq_pval,
              stability = dat_stability)
  
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  if(verbose) print(paste("Finished pipeline in",
                          time_taken_string(pipeline_start)));flush.console()
  return(out)
}

cv_fold <- function(dat_list, 
                    nfolds = 5, 
                    nruns = 2, 
                    batch_label = NULL, 
                    stratified_cv = FALSE, 
                    mixed_cv = FALSE,
                    ...) {
  out <- list()
  for (i in 1:length(dat_list)) {
    temp <- as.data.frame(t(dat_list[[i]]))
    colnames(temp) <- paste0("dim", 1:ncol(temp))
    temp$id <- rownames(temp)
    # Add batch label(s) as separate column(s)
    if (!is.null(batch_label)) {
      if (is.null(dim(batch_label))) {
        if (!is.null(names(batch_label))) {
          batch_id <- names(batch_label)
        } else {
          # Assume data and batch labels are in the same order
          batch_id <- colnames(dat_list[[i]])
        }
        batch_label <- cbind(batch_label = as.character(batch_label), c())
        rownames(batch_label) <- batch_id
      }
      batch_label <- as.data.frame(batch_label)
      if (any(!(rownames(batch_label) %in% colnames(dat_list[[i]])))) {
        # Assume that mismatches are due to processing ...
        # Assume data and batch labels are in the same order
        rownames(batch_label) <- colnames(dat_list[[i]])
      }
      batch_label$id <- rownames(batch_label)
      temp <- plyr::join(temp, batch_label, by = "id")
    }
    folded <- data.frame()
    for (j in 1:nruns) {
      if (!is.null(batch_label) & (stratified_cv | mixed_cv)) {
        a_ind <- lapply(table(batch_label), function(x) sample(1:x, x))
        b_ind <- sample(1:length(unique(batch_label)), length(unique(batch_label)))
        c_ind <- cumsum(table(batch_label)[unique(batch_label)[b_ind]])
        cv_index <- c()
        for (u in 1:length(b_ind)) {
          un <- unique(batch_label)[b_ind[u]]
          cv_index[batch_label == un] <- a_ind[[un]] + ifelse(u > 1, c_ind[u-1], 0)
        }
        if (stratified_cv) {
          # Stratified cv folds such that holdout set labels mostly do not match to rest of data
          cv_index <- cv_index %/% -(length(batch_label) %/% -nfolds) + 1
        } else {
          # Mixed cv folds such that labels are evenly distributed within folds
          cv_index <- cv_index %% nfolds + 1
        }
      } else {
        # Completely random folds
        cv_index <- sample(1:nrow(temp)) %% nfolds + 1
      }
      # Got index, create folds +1 extra "fold" with whole data
      # TODO: reference is the same accross all runs, maybe include it only once?
      for (f in 1:(nfolds+1)) {
        # TODO: fix downstream support so that test set can be included too
        tempfold <- temp[cv_index != f, ]
        tempfold$fold <- f
        tempfold$run <- j
        tempfold$cv_index <- cv_index[cv_index != f]
        folded <- rbind(folded, tempfold)
      }
    }
    out[[i]] <- folded
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
    dat <- dat_list[[best$data_set]]

    # Get dimensionality reduction method from name
    temp <- regexpr("^[A-Za-z]*", best$dimred)
    dimred_method <- substr(best$dimred, temp,
                            temp + attributes(temp)$match.length - 1)
    # Get number of dimensions from name
    temp <- regexpr("[0-9]*$", best$dimred)
    dimred_dim <- substr(best$dimred, temp,
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
    embed <- dim_reduction_suite(dat,
                                 dimred_methods = dimred_method,
                                 output_dimensions = dimred_dim,
                                 tsne_perplexities = dimred_perp,
                                 include_original = FALSE)[[1]]
  }

  out <- suppressWarnings(clValid::clValid(t(embed),
                                           nClust = k,
                                           clMethods = cluster_method,
                                           metric = metric,
                                           validation=c()))

  # TODO: create generic unpacker for clValid output or use clustering methods directly
  temp <- clValid::clusters(out, as.character(best$m))

  if (best$m %in% c("hierarchical", "diana", "agnes")) {
    cluster <- cutree(temp, k)
  } else if (best$m == "pam") {
    cluster <- temp[[1]]$clustering
  } else if (best$m == "kmeans") {
    cluster <- temp[[1]]$cluster
  } else if (best$m == "sota") {
    cluster <- temp[[1]]$clust
  } else {
    stop(paste("Unsupported clustering method:", as.character(cluster_method)))
  }

  return(list(embedding = t(embed), clustering = cluster))
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
                              chisq_significance_level = 0.05,
                              w_connectivity = -0.1,
                              w_dunn = 0.1,
                              w_silhouette = 0.5,
                              w_stability = 1.0,
                              w_batch_effect = -1.0) {
  out <- list()
  # For each data set index i
  for (i in 1:length(input)) {
    out[[i]] <- list()

    # For each dimension reduction method index j
    for (j in 1:length(input[[i]]$evals)) {
      out_j <- plyr::ddply(input[[i]]$evals[[j]]$embed_eval$metrics,
                           c("k", "m", "metric"),
                           function(x) data.frame(mean = mean(x$value), sd = sd(x$value)))

      temp <- input[[i]]$evals[[j]]$embed_stab
      colnames(temp)[colnames(temp) == "jdist"] <- "mean"
      temp$sd <- NA
      temp$metric <- "JaccardDist"
      out_j <- rbind(out_j, temp[,c(1,2,5,3,4)])

      # Compute chisq test rejection rate if batch_label was given
      if (!is.null(input[[i]]$evals[[j]]$embed_eval$chisq_pval))  {
        temp <- reshape2::melt(input[[i]]$evals[[j]]$embed_eval$chisq_pval,
                               id.vars = c("run", "fold", "k", "m"),
                               variable.name = "metric",
                               value.name = "p")
        temp$metric <- gsub("^p", "ChiSqRR_", temp$metric)
        chisq_rr <- plyr::ddply(temp,
                                c("k", "m", "metric"),
                                function(x) data.frame(mean = mean(x$p < chisq_significance_level),
                                                       sd = NA))
      }

      out[[i]][[j]] <- rbind(out_j, chisq_rr)
    }
    if(!is.null(names(input[[i]]$evals))) names(out[[i]]) <- names(input[[i]]$evals)

    # TODO: add class_associations output as a separate thing?
    #       Or report best combination and run separately.

    if (!is.null(names(out[[i]]))) {
      dimred_method <- rep(names(out[[i]]), sapply(out[[i]], nrow))
    } else {
      dimred_method <- rep(1:length(out[[i]]), sapply(out[[i]], nrow))
    }
    out[[i]] <- Reduce(rbind, out[[i]])
    out[[i]]$dimred <- dimred_method
  }
  if(!is.null(names(input))) names(out) <- names(input)

  # Combine everything into one data.frame
  if (!is.null(names(out))) {
    data_set_name <- rep(names(out), sapply(out, nrow))
  } else {
    data_set_name <- rep(1:length(out), sapply(out, nrow))
  }
  out <- Reduce(rbind, out)
  out$data_set <- data_set_name

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

  scores <- plyr::ddply(out, c("data_set", "dimred", "k", "m"), wsumfun)

  best <- scores[which.max(scores$wsum),]

  return(list(all = scores[order(scores$wsum, decreasing = TRUE),], 
              best = best))
}

# Helper for pipeline verbosity
time_taken_string <- function(start) {
  time <- Sys.time() - start
  paste(sprintf("%.2f", time), attributes(time)$units)
}
