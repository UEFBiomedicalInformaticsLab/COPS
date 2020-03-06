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
#' @param dat data matrix, features on rows and samples on columns
#' @param dimred_methods vector of method names, see details for options
#' @param output_dimensions vector of dimensionalities to compute using each applicable method
#' @param tsne_perplexities vector of t-SNE perplexity settings to generate embeddings with
#' @param umap_neighbors scalar UMAP parameter, affects manifold resolution
#' @param include_original if \code{TRUE}, includes original data in output list
#' @param parallel NOT IMPLEMENTED
#' @param ... extra arguments are ignored currently
#'
#' @return Returns a \code{list} of embeddings, elements are named based on methods used
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
dim_reduction_suite <- function(dat,
                                dimred_methods = c("pca", "tsne", "umap"),
                                output_dimensions = c(2:4),
                                tsne_perplexities = c(5,30,50),
                                umap_neighbors = 20,
                                include_original = TRUE,
                                parallel = FALSE,
                                ...) {
  # TODO: redirect extra arguments
  # eargs <- list(...)

  # TODO: parallelize with foreach
  #out <- foreach(i = prod(length(dimred_methods, dimensions)))

  out <- list()
  if (include_original) {
    temp <- dat
    rownames(temp) <- paste0("dim", 1:nrow(temp))
    out$original <- temp
  }
  
  for (m in dimred_methods) {
    if (m %in% c("pca", "umap")) {
      dims <- output_dimensions
    } else if (m == "tsne") {
      # Rtsne throws error if perplexity is too high compared to data size
      dims <- tsne_perplexities[3 * tsne_perplexities < ncol(dat) - 1]
    } else {
      warning(paste("Unsupported method:", m))
      dims <- c()
    }
    if (m == "pca") pca_temp <- FactoMineR::PCA(t(dat),
                                                scale.unit = FALSE,
                                                ncp = max(dims),
                                                graph = FALSE)
    # Handle everything in the same loop
    for (d in dims) {
      if (m == "pca") {
        temp <- pca_temp$ind$coord[,1:d]
        colnames(temp) <- paste0("dim", 1:d)
      } else if (m == "tsne") {
        temp <- Rtsne::Rtsne(t(dat),
                             dims = 2,
                             perplexity = d,
                             initial_dims = min(50, dim(dat)[1]),
                             check_duplicates = FALSE,
                             verbose = FALSE)$Y
        colnames(temp) <- paste0("dim", 1:2)
        rownames(temp) <- colnames(dat)
      } else if (m == "umap") {
        temp <- uwot::umap(t(dat),
                           n_neighbors = umap_neighbors,
                           n_components = d,
                           pca = min(50, dim(dat)[1]),
                           verbose = FALSE,
                           init = "normlaplacian")
        colnames(temp) <- paste0("dim", 1:d)
        rownames(temp) <- colnames(dat)
      } else {
        # never run
        temp <- NA
      }
      out[[paste0(m,d)]] <- t(temp) # keep samples on columns
    }
  }
  return(out)
}

# Input has to be a list of data.frames with columns corresponding to genes and also run and fold ids
# Output has to be a list of data.frames because embeddings may have different dimensions
cv_dimred <- function(dat_list, ...) {
  temp_list <- list()
  for (i in 1:length(dat_list)) {
    temp <- dat_list[[i]]
    temp$datname <- names(dat_list)[i]
    if (is.null(temp$datname)) temp$datname <- i
    #temp <- plyr::dlply(temp, c("run", "fold"), function(x) x)
    temp <- split(temp, temp[c("run", "fold")])
    temp_list <- c(temp_list, temp)
  }
  
  #keep_columns <- c("id", "run", "fold", "cv_index", "datname")
  out <- foreach(i = 1:length(temp_list),
                      .combine = c,
                      .export = c("dim_reduction_suite"),
                      .packages = c("FactoMineR", "Rtsne", "uwot", "plyr")) %dopar% {
    temp <- temp_list[[i]]
    rownames(temp) <- temp$id
    #temp <- temp[!colnames(temp) %in% c("id", "run", "fold", "cv_index", "datname")]
    temp <- dim_reduction_suite(t(temp[grepl("^dim[0-9]+$", colnames(temp))]), ...)
    for (j in 1:length(temp)) {
      temp[[j]] <- as.data.frame(t(temp[[j]]))
      temp[[j]]$id <- rownames(temp[[j]])
      temp[[j]] <- plyr::join(temp[[j]], 
                              #temp_list[[i]][keep_columns[keep_columns %in% colnames(temp_list[[i]])]], 
                              temp_list[[i]][!grepl("^dim[0-9]+$", colnames(temp_list[[i]]))],
                              by = "id")
    }
    temp
  }
  return(out)
}
