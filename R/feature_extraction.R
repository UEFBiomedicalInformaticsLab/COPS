#' Dimensionality reduction suite
#'
#' Applies multiple dimensionality reduction techniques to input data.
#'
#' Method options
#' \itemize{
#' \item "none" - Returns original data as is
#' \item "pca" - Principal Component Analysis
#' \item "tsne" - t-distributed stochastic neighbor embedding
#' \item "umap" - Uniform Manifold Approximation and Projection for Dimension Reduction
#' }
#' By default also appends original data to outputs.
#'
#' @param dat Data matrix, features on columns and samples on rows.
#' @param dimred_methods Vector of method names, see details for options.
#' @param output_dimensions Vector of dimensionalities to compute using each applicable method.
#' @param pca_dims PCA specific output dimensions.
#' @param umap_dims UMAP specific output dimensions.
#' @param tsne_perplexities Vector of t-SNE perplexity settings to generate embeddings with.
#' @param tsne_pca Whether to apply PCA before t-SNE, which massively boosts performance.
#' @param umap_neighbors UMAP parameter, affects manifold computation.
#' @param initial_dims Number of principal components used in t-SNE and UMAP.
#' @param ... Extra arguments are ignored.
#'
#' @return Returns a \code{list} of embeddings, elements are named based on methods used
#' @export
#' @importFrom FactoMineR PCA
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
dim_reduction_suite <- function(dat,
                                dimred_methods = c("pca", "umap"), 
                                output_dimensions = NULL, 
                                pca_dims = c(2),
                                umap_dims = c(2),
                                tsne_perplexities = c(45),
                                tsne_pca = TRUE, 
                                umap_neighbors = 20,
                                initial_dims = 50,
                                ...) {
  # TODO: redirect extra arguments
  # eargs <- list(...)

  out <- list()
  if ("none" %in% dimred_methods) {
    out$original <- dat
    dimred_methods <- dimred_methods[dimred_methods != "none"]
  }
  
  for (m in dimred_methods) {
    if (m %in% c("pca", "umap")) {
      if (is.null(output_dimensions)) {
        if (m == "pca") {
          dims <- pca_dims
        } else if (m == "umap") {
          dims <- umap_dims
        }
      } else {
        dims <- output_dimensions
      }
    } else if (m == "tsne") {
      # Rtsne throws error if perplexity is too high compared to data size 
      # We can remove the offending values (possibly yielding 0 t-SNE based reductions)
      dims <- tsne_perplexities[3 * tsne_perplexities < nrow(dat) - 1]
    } else {
      warning(paste("Unsupported method:", m))
      dims <- c()
    }
    if (m == "pca") pca_temp <- FactoMineR::PCA(dat,
                                                scale.unit = FALSE,
                                                ncp = max(dims),
                                                graph = FALSE)
    # Handle everything in the same loop
    for (d in dims) {
      if (m == "pca") {
        temp <- pca_temp$ind$coord[,1:d]
        colnames(temp) <- paste0("dim", 1:d)
      } else if (m == "tsne") {
        if (3 * d > dim(dat)[1] - 1) stop("t-SNE perplexity is too high.")
        if (tsne_pca) {
          tsne_pca_temp <- FactoMineR::PCA(dat, scale.unit = FALSE, ncp = min(c(initial_dims, dim(dat))), graph = FALSE)
          tsne_input <- tsne_pca_temp$ind$coord
        } else {
          tsne_input <- dat
        }
        temp <- Rtsne::Rtsne(tsne_input,
                             dims = 2,
                             perplexity = d,
                             #initial_dims = min(100, dim(dat)[2]),
                             check_duplicates = FALSE,
                             pca = FALSE,
                             partial_pca = FALSE,
                             verbose = FALSE)$Y
        colnames(temp) <- paste0("dim", 1:2)
        #rownames(temp) <- colnames(dat)
      } else if (m == "umap") {
        temp <- uwot::umap(dat,
                           n_neighbors = umap_neighbors,
                           n_components = d,
                           pca = min(c(initial_dims, dim(dat))),
                           verbose = FALSE,
                           init = "normlaplacian")
        colnames(temp) <- paste0("dim", 1:d)
        #rownames(temp) <- colnames(dat)
      } else {
        # never run
        temp <- NA
      }
      out[[paste0(m,d)]] <- temp
    }
  }
  return(out)
}

#' Dimensionality reduction on cross-validated data sets
#'
#' @param dat_list A list of data.tables.
#' @param cv_index A data.frame indicating cv folds and runs such as returned by \code{\link{cv_fold}}.
#' @param cv_split_data Can be set to FALSE if \code{dat_list} elements already contain the columns \code{"run"} and \code{"fold"}.
#' @param ... Extra arguments are passed to \code{\link{dim_reduction_suite}}.
#'
#' @return list of data sets
#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table data.table
cv_dimred <- function(dat_list, 
                      cv_index, 
                      cv_split_data = TRUE, 
                      parallel = 1, 
                      by = c("run", "fold"), 
                      ...) {
  temp_list <- list()
  if (cv_split_data) {
    for (i in 1:length(cv_index)) {
      temp <- cv_index[[i]]
      datname <- names(cv_index)[i]
      if (is.null(datname)) datname <- i
      temp$datname <- datname
      temp <- split_by_safe(temp, by)
      temp <- lapply(temp, function(x) as.data.frame(merge(dat_list[[datname]], x, by = "id")))
      temp_list <- c(temp_list, temp)
    }
  } else {
    for (i in 1:length(dat_list)) {
      temp <- split_by_safe(dat_list[[i]], by)
      temp_list <- c(temp_list, temp)
    }
  }
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(foreach(i = temp_list, 
                 .combine = c,
                 .export = c("dim_reduction_suite"), #"dat_list"),
                 .packages = c("FactoMineR", "Rtsne", "uwot", "plyr")) %dopar% {
    sel <- grep("^dim[0-9]+$", colnames(i))
    dr_temp <- dim_reduction_suite(i[,sel], ...)
    dr_temp <- lapply(dr_temp, function(x) cbind(i[,-sel], as.data.frame(x)))
    dr_temp
  }, finally = close_parallel_cluster(parallel_clust))
  return(out)
}
