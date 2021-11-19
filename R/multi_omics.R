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
      mofa_obj <- MOFA2::create_mofa(dat_list_clust)
      data_opts <- MOFA2::get_default_data_options(mofa_obj)
      data_opts$scale_views <- mofa_scale_views
      model_opts <- MOFA2::get_default_model_options(mofa_obj)
      train_opts <- MOFA2::get_default_training_options(mofa_obj)
      train_opts$convergence_mode <- mofa_convergence_mode
      train_opts$maxiter <- mofa_maxiter
      
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
      mofa_embedding$id <- rownames(mofa_embedding)
      mofa_embedding$drname <- "MOFA2"
      
      mofa_diss <- clustering_dissimilarity(temp, distance_metric, correlation_method)
      mofa_cops_clust <- COPS::clustering_analysis(mofa_embedding, 
                                                   n_clusters = n_clusters,
                                                   clustering_dissimilarity = mofa_diss,
                                                   ...)
      cbind(non_data_cols[[1]], mofa_cops_clust)
    }, error = function(e) return(NULL))
    if(!is.null(temp_res)) if(nrow(temp_res) > 1) res <- c(res, list(temp_res))
  }
  return(plyr::rbind.fill(res))
}