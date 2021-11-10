#' Multi-omic clustering via multi-view clustering or integration
#'
#' @param dat_list_clust 
#' @param non_data_cols 
#' @param multi_view_methods 
#' @param n_clusters 
#'
#' @return
#' @export
multi_omic_clustering <- function(dat_list_clust, 
                                  non_data_cols,
                                  multi_view_methods = "iClusterPlus",
                                  n_clusters = 2, 
                                  nmf_maxiter = 200,
                                  nmf_st.count = 20,
                                  nmf_n.ini = 30,
                                  nmf_ini.nndsvd = TRUE,
                                  ...) {
  res <- list()
  if("iClusterPlus" %in% multi_view_methods) {
    if (length(dat_list_clust) > 4) stop("iClusterPlus only supports up to four views.")
    if (length(dat_list_clust) >= 1) dt1 <- dat_list_clust[[1]] else dt1 <- NULL
    if (length(dat_list_clust) >= 2) dt2 <- dat_list_clust[[2]] else dt2 <- NULL
    if (length(dat_list_clust) >= 3) dt3 <- dat_list_clust[[3]] else dt3 <- NULL
    if (length(dat_list_clust) == 4) dt4 <- dat_list_clust[[4]] else dt4 <- NULL
    
    for (k in n_clusters) {
      temp_res <- iClusterPlus::iClusterPlus(dt1, dt2, dt3, dt4, 
                                             type = rep("gaussian", length(dat_list_clust)),
                                             K = k-1)
      k_res <- data.frame(m = "iClusterPlus", 
                          k = k,
                          cluster = temp_res$clusters)
      k_res <- cbind(non_data_cols[[1]], k_res)
      res <- c(res, list(k_res))
    }
  }
  if ("IntNMF" %in% multi_view_methods) {
    for (k in n_clusters) {
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
      k_res <- cbind(non_data_cols[[1]], k_res)
      res <- c(res, list(k_res))
    }
  }
  return(plyr::rbind.fill(res))
}