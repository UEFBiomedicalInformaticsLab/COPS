#' Multi-omic clustering via multi-view clustering or integration
#'
#' @param dat_list_clust 
#' @param non_data_cols 
#' @param multi_view_methods 
#' @param n_clusters 
#'
#' @return
#' @export
#' @importFrom iClusterPlus iClusterPlus
#' @importFrom IntNMF nmf.mnnals
#' @importFrom MOFA2 create_mofa get_default_data_options get_default_model_options get_default_training_options prepare_mofa run_mofa
#' @importFrom reticulate use_python
multi_omic_clustering <- function(dat_list_clust, 
                                  non_data_cols,
                                  multi_view_methods = "iClusterPlus",
                                  n_clusters = 2, 
                                  distance_metric = "euclidean", 
                                  correlation_method = "spearman",
                                  nmf_maxiter = 200,
                                  nmf_st.count = 20,
                                  nmf_n.ini = 30,
                                  nmf_ini.nndsvd = TRUE,
                                  mofa_scale_views = FALSE,
                                  mofa_convergence_mode = "medium",
                                  mofa_maxiter = 1000,
                                  mofa_environment = NULL,
                                  mofa_threads = 1,
                                  anf_neighbors = 20,
                                  kkmeans_maxiter = 100,
                                  kkmeans_n_init = 100,
                                  kkmeans_parallel = 1, 
                                  mkkm_mr_lambda = 1, 
                                  mkkm_mr_tolerance = 1e-8, 
                                  mkkm_mr_parallel = 1, 
                                  data_is_kernels = FALSE, 
                                  ...) {
  res <- list()
  if("iClusterPlus" %in% multi_view_methods) {
    if (length(dat_list_clust) > 4) stop("iClusterPlus only supports up to four views.")
    if (length(dat_list_clust) >= 1) dt1 <- dat_list_clust[[1]] else dt1 <- NULL
    if (length(dat_list_clust) >= 2) dt2 <- dat_list_clust[[2]] else dt2 <- NULL
    if (length(dat_list_clust) >= 3) dt3 <- dat_list_clust[[3]] else dt3 <- NULL
    if (length(dat_list_clust) == 4) dt4 <- dat_list_clust[[4]] else dt4 <- NULL
    
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- iClusterPlus::iClusterPlus(dt1, dt2, dt3, dt4, 
                                               type = rep("gaussian", length(dat_list_clust)),
                                               K = k-1)
        k_res <- data.frame(m = "iClusterPlus", 
                            k = k,
                            cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], k_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("IntNMF" %in% multi_view_methods) {
    for (k in n_clusters) {
      k_res <- tryCatch({
        nmf_view_weights <- sapply(dat_list_clust, function(x) mean(sqrt(apply(x^2, 1, sum))))
        temp_res <- nmf.mnnals(dat_list_clust, 
                               k = k, 
                               maxiter = nmf_maxiter,
                               st.count = nmf_st.count,
                               n.ini = nmf_n.ini,
                               ini.nndsvd = nmf_ini.nndsvd,
                               wt = nmf_view_weights)
        k_res <- data.frame(m = "IntNMF", 
                            k = k,
                            cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], k_res)
        }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("MOFA2" %in% multi_view_methods) {
    temp_res <- tryCatch({
      Sys.setenv(OMP_NUM_THREADS=mofa_threads)
      Sys.setenv(MKL_NUM_THREADS=mofa_threads)
      if (!is.null(mofa_environment)) {
        reticulate::use_python(mofa_environment)
      }
      
      mofa_obj <- MOFA2::create_mofa(lapply(dat_list_clust, t))
      data_opts <- MOFA2::get_default_data_options(mofa_obj)
      data_opts$scale_views <- mofa_scale_views
      model_opts <- MOFA2::get_default_model_options(mofa_obj)
      train_opts <- MOFA2::get_default_training_options(mofa_obj)
      train_opts$convergence_mode <- mofa_convergence_mode
      train_opts$maxiter <- mofa_maxiter
      train_opts$verbose <- FALSE
      
      mofa_obj <- MOFA2::prepare_mofa(
        object = mofa_obj,
        data_options = data_opts,
        model_options = model_opts,
        training_options = train_opts
      )
      
      mofa_obj <- MOFA2::run_mofa(mofa_obj, outfile = NULL, save_data = FALSE)
      
      mofa_embedding <- mofa_obj@expectations$Z$group1
      colnames(mofa_embedding) <- paste0("dim", 1:ncol(mofa_embedding))
      mofa_embedding <- as.data.frame(mofa_embedding)
      #mofa_embedding$id <- rownames(mofa_embedding)
      mofa_embedding$drname <- "MOFA2"
      
      mofa_diss <- clustering_dissimilarity(mofa_embedding, distance_metric, correlation_method)
      mofa_cops_clust <- clustering_analysis(cbind(mofa_embedding, non_data_cols[[1]]), 
                                             n_clusters = n_clusters,
                                             clustering_dissimilarity = mofa_diss,
                                             ...)
      mofa_cops_clust
    }, error = function(e) return(NULL))
    if(!is.null(temp_res)) if(nrow(temp_res) > 1) res <- c(res, list(temp_res))
  }
  if ("ANF" %in% multi_view_methods) {
    temp_res <- ANF::ANF(K = anf_neighbors)
  }
  if ("kkmeanspp" %in% multi_view_methods) {
    # Kernel k-means++
    if (data_is_kernels) {
      multi_omic_kernels_linear <- dat_list_clust
    } else {
      # Just linear for now
      multi_omic_kernels_linear <- lapply(dat_list_clust, function(x) (x) %*% t(x))
      multi_omic_kernels_linear <- lapply(multi_omic_kernels_linear, center_kernel)
      multi_omic_kernels_linear <- lapply(multi_omic_kernels_linear, normalize_kernel)
    }
    
    # Average kernel
    multi_omic_kernels_linear <- Reduce('+', multi_omic_kernels_linear) / length(multi_omic_kernels_linear)
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- kernel_kmeans(multi_omic_kernels_linear, k, n_initializations = kkmeans_n_init, maxiter = kkmeans_maxiter)
        temp_res <- data.frame(m = "kkmeanspp", k = k, cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("mkkm_mr" %in% multi_view_methods) {
    if (data_is_kernels) {
      kernels <- dat_list_clust
    } else {
      # Just linear for now
      kernels <- lapply(dat_list_clust, function(x) (x) %*% t(x))
      kernels <- lapply(kernels, center_kernel)
      kernels <- lapply(kernels, normalize_kernel)
    }
    for (k in n_clusters) {
      k_res <- tryCatch({
        # Optimize combined kernel
        optimal_kernel <- mkkm_mr(kernels, 
                                  k = k, 
                                  lambda = mkkm_mr_lambda, 
                                  tolerance = mkkm_mr_tolerance, 
                                  parallel = mkkm_mr_parallel)
        # Run k-means++
        temp_res <- kernel_kmeans(optimal_kernel$K, 
                                  n_k = k, 
                                  n_initializations = kkmeans_n_init, 
                                  maxiter = kkmeans_maxiter,
                                  parallel = kkmeans_parallel)
        temp_res <- data.frame(m = "mkkm_mr", k = k, cluster = temp_res$clusters, kernel_mix = paste(optimal_kernel$mu, collapse = ";"))
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  return(plyr::rbind.fill(res))
}