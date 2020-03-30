#' Dimensionality reduction suite
#'
#' Applies multiple dimensionality reduction techniques to input data.
#'
#' Method options
#' \itemize{
#' \item "pca" - Principal Component Analysis
#' \item "tsne" - t-distributed stochastic neighbor embedding
#' \item "umap" - Uniform Manifold Approximation and Projection for Dimension Reduction
#' }
#' By default also appends original data to outputs.
#'
#' @param dat data matrix, features on columns and samples on rows
#' @param dimred_methods vector of method names, see details for options
#' @param output_dimensions vector of dimensionalities to compute using each applicable method
#' @param pca_dims PCA specific output dimensions
#' @param umap_dims UMAP specific output dimensions
#' @param tsne_perplexities vector of t-SNE perplexity settings to generate embeddings with
#' @param umap_neighbors scalar UMAP parameter, affects manifold resolution
#' @param include_original if \code{TRUE}, includes original data in output list
#' @param ... extra arguments are ignored currently
#'
#' @return Returns a \code{list} of embeddings, elements are named based on methods used
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
dim_reduction_suite <- function(dat,
                                dimred_methods = c("pca", "umap"), # TODO: fix "tsne" with data.table?
                                output_dimensions = NULL, #c(2:4),
                                pca_dims = 2:6,
                                umap_dims = 2:10,
                                tsne_perplexities = c(5,30,50),
                                umap_neighbors = 20,
                                include_original = FALSE,
                                ...) {
  # TODO: redirect extra arguments
  # eargs <- list(...)

  out <- list()
  if (include_original) {
    out$original <- dat
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
        temp <- Rtsne::Rtsne(dat,
                             dims = 2,
                             perplexity = d,
                             initial_dims = min(50, dim(dat)[2]),
                             check_duplicates = FALSE,
                             verbose = FALSE)$Y
        colnames(temp) <- paste0("dim", 1:2)
        #rownames(temp) <- colnames(dat)
      } else if (m == "umap") {
        temp <- uwot::umap(dat,
                           n_neighbors = umap_neighbors,
                           n_components = d,
                           pca = min(50, dim(dat)[2]),
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
#' @param dat_list list of data sets
#' @param ... arguments passed to \code{\link{dim_reduction_suite}}
#'
#' @return list of data sets
#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom data.table data.table
cv_dimred <- function(dat_list, cv_index, ...) {
  temp_list <- list()
  for (i in 1:length(cv_index)) {
    temp <- cv_index[[i]]
    datname <- names(cv_index)[i]
    temp$datname <- datname
    if (is.null(temp$datname)) temp$datname <- i
    temp <- split(temp, by = c("run", "fold"))
    temp <- lapply(temp, function(x) as.data.frame(merge(dat_list[[datname]], x, by = "id")))
    temp_list <- c(temp_list, temp)
  }
  
  out <- foreach(i = temp_list, 
                 .combine = c,
                 .export = c("dim_reduction_suite"), #"dat_list"),
                 .packages = c("FactoMineR", "Rtsne", "uwot", "plyr")) %dopar% {
    sel <- grep("^dim[0-9]+$", colnames(i))
    dr_temp <- dim_reduction_suite(i[,sel], ...)
    dr_temp <- lapply(dr_temp, function(x) cbind(i[,-sel], as.data.frame(x)))
    dr_temp
  }
  return(out)
}
