#' Multi-omic clustering via multi-view clustering or integration
#'
#' @param dat_list List of input \code{data.frame}s for input.
#' @param meta_data A single \code{data.frame} or a list that includes meta data for 
#'   each view. If a list is provided, at the moment only the first element is 
#'   used (by appending to clustering output). 
#' @param multi_omic_methods Vector of algorithm names to be applied. See details. 
#' @param n_clusters Integer vector of number of clusters to output. 
#' @param distance_metric Distance metric for clustering factorized data 
#'   (only for MOFA).
#' @param correlation_method Correlation method for \code{distance_metric}, 
#'   if applicable.
#' @param standardize_data If set, standardizes data before clustering. 
#' @param non_negativity_transform Vector of transformation names for IntNMF. 
#'   See details below. 
#' @param view_distributions A vector specifying the distribution to use for each view. 
#'   Used by iCluster+, iClusterBayes and MOFA2. 
#'   Options are "gaussian", "bernoulli" and "poisson". 
#' @param icp_lambda iCluster+ L1 penalty for each view. 
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#' @param icp_burnin iCluster+ number of MCMC burn in samples for approximating 
#'   joint distribution of latent variables. 
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#' @param icp_draw iCluster+ number of MCMC samples to draw after burn in for 
#'   approximating joint distribution of latent variables.  
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#' @param icp_maxiter iCluster+ maximum number of Newton-Rhapson (EM) iterations. 
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#'   \code{\link[iClusterPlus]{iClusterBayes}}. 
#' @param icp_sdev iCluster+ MCMC random walk standard deviation. 
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#' @param icp_eps iCluster+ algorithm convergence threshold. 
#'   See \code{\link[iClusterPlus]{iClusterPlus}}.
#' @param icb_burnin iClusteBayes number of samples for MCMC burn in. 
#'   See \code{\link[iClusterPlus]{iClusterBayes}}.
#' @param icb_draw iClusteBayes number of MCMC samples to draw after burn in. 
#'   See \code{\link[iClusterPlus]{iClusterBayes}}.
#' @param icb_sdev iClusteBayes MCMC random walk standard deviation. 
#'   See \code{\link[iClusterPlus]{iClusterBayes}}.
#' @param icb_thin iClusteBayes MCMC thinning, only one sample in every icb_thin 
#'   samples will be used. 
#'   See \code{\link[iClusterPlus]{iClusterBayes}}.
#' @param nmf_maxiter Maxiter for IntNMF. See 
#'   \code{\link[IntNMF]{nmf.mnnals}}.
#' @param nmf_st.count Count stability for IntNMF. 
#'   See \code{\link[IntNMF]{nmf.mnnals}}.
#' @param nmf_n.ini Number of initializations for IntNMF. 
#'   See \code{\link[IntNMF]{nmf.mnnals}}.
#' @param nmf_ini.nndsvd If set, IntNMF uses NNDSVD for initialization. 
#'   See \code{\link[IntNMF]{nmf.mnnals}}.
#' @param nmf_scaling Omic weights that are used for scaling. Defaults to the 
#'   Frobenius norm ratio similarly to Chalise et al. 2017.
#' @param mofa_convergence_mode MOFA convergence threshold. 
#'   See \code{\link[MOFA2]{get_default_training_options}}.
#' @param mofa_maxiter MOFA maximum iterations. 
#'   See \code{\link[MOFA2]{get_default_training_options}}.
#' @param mofa_environment If set, uses the specified Python environment 
#'   (with mofapy). Defaults to basilisk. 
#' @param mofa_lib_path Path to libpython. May be required if using non-default 
#'   \code{mofa_environment}. 
#' @param anf_neighbors Number of neighbours to use in knn-graph. 
#' @param kkmeans_algorithm See \code{\link{kernel_kmeans}}.
#' @param kkmeans_refine See \code{\link{kernel_kmeans}}.
#' @param kkmeans_tol See \code{\link{kernel_kmeans}}.
#' @param kkmeans_maxiter See \code{\link{kernel_kmeans}}.
#' @param kkmeans_n_init See \code{\link{kernel_kmeans}}.
#' @param mkkm_mr_lambda Regularization parameter for \code{\link{mkkm_mr}}.
#' @param mkkm_mr_tolerance Convergence threshold for \code{\link{mkkm_mr}}.
#' @param mkkm_mr_mosek If set, uses \code{Rmosek} for convex optimization 
#'   instead of \code{CVXR} for \code{\link{mkkm_mr}}.
#' @param ecmc_a Regularization parameter for \code{\link{ECMC}}.
#' @param ecmc_b Regularization parameter for \code{\link{ECMC}}.
#' @param ecmc_eps Convergence threshold for \code{\link{ECMC}}.
#' @param ecmc_maxiter Maximum number of iterations for \code{\link{ECMC}}.
#' @param ecmc_mkkm_mr If set, uses \code{\link{mkkm_mr}} on consensus kernels 
#'   obtained from \code{\link{ECMC}}. Otherwise uses the average kernel and 
#'     kernel k-means. 
#' @param data_is_kernels If \code{TRUE}, input data is assumed to be kernel matrices. 
#'   Otherwise kernels are computed based on input data and the 
#'   \code{kernels} parameter of \code{\link{get_multi_omic_kernels}}. 
#' @param zero_var_removal If set, removes all zero variance features 
#'   from the data. It is called fold-wise, because this is assumed to be run 
#'   inside CV. 
#' @param mvc_threads Number of threads to use for supported operations. 
#' @param gene_id_list List of gene/feature names for each view. If set, matches 
#'   pipeline standardized feature names ("dim1", "dim2", ...) to names on the list. 
#'   Required for pathway kernels. 
#' @param preprocess_data If the input data has already been processed by the 
#'   \code{\link{COPS}}-pipeline, this should be disabled. 
#' @param ... Arguments are passed to \code{\link{clustering_analysis}} 
#'   when using MOFA and \code{\link{get_multi_omic_kernels}} when using kernel 
#'   methods. 
#'
#' @return \code{data.frame} of clustering results
#' 
#' @details
#' Supported methods:
#' \itemize{
#'   \item "ANF" - Affinity Network Fusion \code{\link[ANF]{ANF}}
#'   \item "iClusterPlus" or "iCluster+" - \code{\link[iClusterPlus]{iClusterPlus}}. Supports only up to 4 views. 
#'   \item "iClusterBayes" -  code{\link[iClusterPlus]{iClusterBayes}}. Supports only up to 6 views
#'   \item "IntNMF" - Integrative Non-negative Matrix Factorization 
#'     \code{\link[IntNMF]{nmf.mnnals}}.
#'   \item "average_kernel" - kernel k-means with average kernel. 
#'   \item "mkkm_mr" - Multiple Kernel K-Means with Matrix-induced Regularization 
#'     \code{\link{mkkm_mr}}.
#'   \item "ECMC" - Enhanced Consensus Multi-view Clustering \code{\link{ECMC}}.
#'   \item "MOFA2" - Multi-Omics Factor Analysis. 
#'     See \code{vignette("getting_started_R", "MOFA2")}. 
#'     Resulting factorization is clustered with single-view algorithms by using 
#'     \code{\link{clustering_analysis}}.
#' }
#' 
#' For supported kernels see \code{\link{get_multi_omic_kernels}}: 
#' 
#' NMF non-negativity transform may be necessary if non-negativity was not 
#' considered while pre-processing the data. There are a few convenience 
#' functions included to transform the data as needed: 
#' \itemize{
#'   \item "logistic" - \code{1/(1 + exp(-x))}, maps input from (-Inf,Inf) to [0,1]. 
#'     Used for e.g. microarray data or methylation M-values. 
#'   \item "rank" - ranks values and divides by length, maps input from 
#'     (-Inf,Inf) to [0,1].
#'   \item "offset2" - adds 2 to input. Useful for e.g. copy number alterations 
#'     (assuming no alterations lower than -2). 
#' }
#' 
#' @export
#' @importFrom ANF ANF
#' @importFrom iClusterPlus iClusterPlus iClusterBayes
#' @importFrom IntNMF nmf.mnnals
#' @importFrom MOFA2 create_mofa get_default_data_options get_default_model_options get_default_training_options prepare_mofa run_mofa
#' @importFrom reticulate use_python
multi_omic_clustering <- function(
    dat_list, 
    meta_data = NULL,
    multi_omic_methods = "ANF",
    n_clusters = 2, 
    distance_metric = "euclidean", 
    correlation_method = "spearman",
    standardize_data = FALSE,
    non_negativity_transform = rep_len("none", length(dat_list)),
    view_distributions = rep_len("gaussian", length(dat_list)),
    icp_lambda = rep(0.03, length(dat_list)),
    icp_burnin = 100,
    icp_draw = 200,
    icp_maxiter = 20, 
    icp_sdev = 0.05,
    icp_eps = 1e-4,
    icb_burnin = 1000,
    icb_draw = 1200,
    icb_sdev = 0.5, 
    icb_thin = 1,
    nmf_maxiter = 200,
    nmf_st.count = 20,
    nmf_n.ini = 30,
    nmf_ini.nndsvd = TRUE,
    nmf_scaling = "F-ratio",
    mofa_convergence_mode = "medium",
    mofa_maxiter = 1000,
    mofa_environment = NULL,
    mofa_lib_path = NULL,
    anf_neighbors = 20,
    kkmeans_algorithm = "spectral_qr", 
    kkmeans_refine = TRUE, 
    kkmeans_maxiter = 100, 
    kkmeans_n_init = 100, 
    kkmeans_tol = 1e-8, 
    mkkm_mr_lambda = 1, 
    mkkm_mr_tolerance = 1e-8, 
    mkkm_mr_mosek = FALSE, 
    mkkm_mr_mosek_verbosity = 0L, 
    ecmc_a = 1, 
    ecmc_b = 1, 
    ecmc_eps = 1e-6,
    ecmc_maxiter = 100,
    ecmc_mkkm_mr = TRUE, 
    data_is_kernels = FALSE, 
    zero_var_removal = TRUE,
    mvc_threads = 1,
    gene_id_list = NULL,
    preprocess_data = TRUE, 
    ...
) {
  if (preprocess_data) {
    dat_processed <- data_preprocess(dat_list)
    dat_list <- dat_processed[["dat_list"]]
    gene_id_list <- dat_processed[["gene_id_list"]]
    
    meta_data <- list()
    for (j in 1:length(dat_list)) {
      sel <- grep("^dim[0-9]+$", colnames(dat_list[[j]]))
      
      if (data_is_kernels & length(sel) > nrow(dat_list[[j]])) {
        stop("Input kernels are not square!")
      }
      if ("data.table" %in% class(dat_list[[j]])) {
        meta_data[[j]] <- dat_list[[j]][,-..sel]
      } else {
        meta_data[[j]] <- dat_list[[j]][,-sel]
      }
      dat_list[[j]] <- as.matrix(as.data.frame(dat_list[[j]])[,sel])
    }
  } else {
    if (is.null(meta_data)) meta_data <- data.frame(id = rownames(dat_list[[1]]))
    if ("data.frame" %in% class(meta_data)) meta_data <- list(meta_data)
  }
  extra_output <- NULL # For returning things like view weights
  if (zero_var_removal & !data_is_kernels) {
    dat_list <- lapply(dat_list, function(x) x[, apply(x, 2, var) > 0, drop = FALSE])
  }
  if (standardize_data) {
    dat_list <- lapply(dat_list, scale)
  }
  if (any(multi_omic_methods %in% c("kkmeans", "kkmeanspp", "mkkm_mr", "ECMC"))) {
    multi_omic_kernels <- get_multi_omic_kernels(
      dat_list, 
      data_is_kernels = data_is_kernels, 
      gene_id_list = gene_id_list, 
      preprocess_data = FALSE,
      mvc_threads = mvc_threads, 
      ...
    )
  }
  res <- list()
  if(any(c("iClusterPlus", "iClusterBayes") %in% multi_omic_methods)) {
    icp_view_types <- view_distributions
    icp_view_types[icp_view_types == "bernoulli"] <- "binomial"
  }
  
  if(any(c("iCluster+", "iClusterPlus") %in% multi_omic_methods)) {
    if (length(dat_list) > 4) stop("iClusterPlus only supports up to four views.")
    if (length(dat_list) >= 1) dt1 <- dat_list[[1]] else dt1 <- NULL
    if (length(dat_list) >= 2) dt2 <- dat_list[[2]] else dt2 <- NULL
    if (length(dat_list) >= 3) dt3 <- dat_list[[3]] else dt3 <- NULL
    if (length(dat_list) == 4) dt4 <- dat_list[[4]] else dt4 <- NULL
    
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- iClusterPlus::iClusterPlus(
          dt1, 
          dt2, 
          dt3, 
          dt4, 
          type = icp_view_types,
          K = k-1,
          lambda = icp_lambda,
          n.burnin = icp_burnin,
          n.draw = icp_draw,
          maxiter = icp_maxiter, 
          sdev = icp_sdev,
          eps = icp_eps)
        k_res <- data.frame(
          m = "iClusterPlus", 
          k = k,
          cluster = temp_res$clusters)
        if (ncol(meta_data[[1]]) > 0) k_res <- cbind(meta_data[[1]], k_res)
        k_res
      }, error = function(e) {warning(e); return(NULL)})
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if("iClusterBayes" %in% multi_omic_methods) {
    if (length(dat_list) > 6) stop("iClusterPlus only supports up to six views.")
    if (length(dat_list) >= 1) dt1 <- dat_list[[1]] else dt1 <- NULL
    if (length(dat_list) >= 2) dt2 <- dat_list[[2]] else dt2 <- NULL
    if (length(dat_list) >= 3) dt3 <- dat_list[[3]] else dt3 <- NULL
    if (length(dat_list) >= 4) dt4 <- dat_list[[4]] else dt4 <- NULL
    if (length(dat_list) >= 5) dt5 <- dat_list[[5]] else dt5 <- NULL
    if (length(dat_list) == 6) dt6 <- dat_list[[6]] else dt6 <- NULL
    
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- iClusterPlus::iClusterBayes(
          dt1, 
          dt2, 
          dt3, 
          dt4, 
          dt5, 
          dt6, 
          type = icp_view_types,
          K = k-1,
          n.burnin = icb_burnin,
          n.draw = icb_draw,
          sdev = icb_sdev, 
          thin = icb_thin)
        k_res <- data.frame(
          m = "iClusterBayes", 
          k = k,
          cluster = temp_res$clusters)
        if (ncol(meta_data[[1]]) > 0) k_res <- cbind(meta_data[[1]], k_res)
        k_res
      }, error = function(e) {warning(e); return(NULL)})
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("IntNMF" %in% multi_omic_methods) {
    dat_list_nmf <- dat_list
    for (i in 1:length(dat_list_nmf)) {
      if (non_negativity_transform[i] == "logistic") {
        dat_list_nmf[[i]] <- 1/(1 + exp(-dat_list_nmf[[i]]))
      }
      if (non_negativity_transform[i] == "rank") {
        dat_list_nmf[[i]] <- apply(
          dat_list_nmf[[i]], 
          2, 
          function(x) rank(x) / length(x))
      }
      if (non_negativity_transform[i] == "offset2") {
        dat_list_nmf[[i]] <- dat_list_nmf[[i]] + 2
      }
    }
    for (k in n_clusters) {
      k_res <- tryCatch({
        if (nmf_scaling == "F-ratio") {
          # First minmax scale each matrix to [0,1]
          for (i in 1:length(dat_list_nmf)) {
            dat_list_nmf[[i]] <- dat_list_nmf[[i]] - min(dat_list_nmf[[i]])
            dat_list_nmf[[i]] <- dat_list_nmf[[i]] / max(dat_list_nmf[[i]])
          }
          nmf_view_weights <- sapply(dat_list_nmf, Matrix::norm, type = "F")
          nmf_view_weights <- max(nmf_view_weights) / nmf_view_weights
        } else {
          nmf_view_weights <- 1 / sapply(dat_list_nmf, max)
        }
        
        temp_res <- nmf.mnnals(
          dat_list_nmf, 
          k = k, 
          maxiter = nmf_maxiter,
          st.count = nmf_st.count,
          n.ini = nmf_n.ini,
          ini.nndsvd = nmf_ini.nndsvd,
          wt = nmf_view_weights)
        k_res <- data.frame(
          m = "IntNMF", 
          k = k,
          cluster = temp_res$clusters)
        if (ncol(meta_data[[1]]) > 0) k_res <- cbind(meta_data[[1]], k_res)
        k_res
        }, error = function(e) {warning(e); return(NULL)})
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("MOFA2" %in% multi_omic_methods) {
    temp_res <- tryCatch(
      {
        Sys.setenv(OMP_NUM_THREADS=mvc_threads)
        Sys.setenv(MKL_NUM_THREADS=mvc_threads)
        if (!is.null(mofa_lib_path)) {
          Sys.setenv(
            LD_LIBRARY_PATH = paste(
              mofa_lib_path, 
              Sys.getenv("LD_LIBRARY_PATH"), 
              sep = ":"
            )
          )
        }
        if (!is.null(mofa_environment)) {
          reticulate::use_virtualenv(mofa_environment)
        }
        
        mofa_obj <- MOFA2::create_mofa(lapply(dat_list, t))
        data_opts <- MOFA2::get_default_data_options(mofa_obj)
        model_opts <- MOFA2::get_default_model_options(mofa_obj)
        mofa_view_names <- names(model_opts$likelihoods)
        model_opts$likelihoods <- view_distributions
        names(model_opts$likelihoods) <- mofa_view_names
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
        
        mofa_diss <- clustering_dissimilarity_from_data(
          mofa_embedding, 
          distance_metric, 
          correlation_method)
        if (ncol(meta_data[[1]]) > 0) {
          mofa_embedding <- cbind(mofa_embedding, meta_data[[1]])
        }
        mofa_cops_clust <- clustering_analysis(
          mofa_embedding, 
          n_clusters = n_clusters,
          clustering_dissimilarity = mofa_diss,
          ...)
        mofa_cops_clust
      }, 
      error = function(e) {warning(e); return(NULL)})
    if(!is.null(temp_res)) if(nrow(temp_res) > 1) res <- c(res, list(temp_res))
  }
  if ("ANF" %in% multi_omic_methods) {
    aff_dist <- lapply(dat_list, dist)
    aff_dist <- lapply(aff_dist, as.matrix)
    aff_list <- lapply(aff_dist, ANF::affinity_matrix, k = anf_neighbors)
    aff_mat <- ANF::ANF(aff_list, K = anf_neighbors)
    for (k in n_clusters) {
      temp_res <- ANF::spectral_clustering(aff_mat, k)
      temp_res <- data.frame(m = "ANF", k = k, cluster = temp_res)
      if (ncol(meta_data[[1]]) > 0) k_res <- cbind(meta_data[[1]], temp_res)
      res <- c(res, list(k_res))
    }
  }
  if ("average_kernel" %in% multi_omic_methods) {
    multi_omic_kernels <- Reduce('+', multi_omic_kernels) / length(multi_omic_kernels)
    if (kkmeans_algorithm %in% c("spectral", "spectral_qr")) {
      eigs <- eigen(multi_omic_kernels, symmetric = TRUE)$vectors
    } else {
      eigs <- NULL
    }
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- kernel_kmeans(
          K = multi_omic_kernels, 
          n_k = k, 
          algorithm = kkmeans_algorithm, 
          spectral_qr_refine = kkmeans_refine, 
          kernel_eigen_vectors = eigs, 
          max_iter = kkmeans_maxiter, 
          num_init = kkmeans_n_init, 
          tol = kkmeans_tol, 
          parallel = mvc_threads
        )
        temp_res <- data.frame(
          m = "kkmeans", 
          k = k, 
          cluster = temp_res$clusters
        )
        if (ncol(meta_data[[1]]) > 0) temp_res <- cbind(meta_data[[1]], temp_res)
        temp_res
      }, error = function(e) {warning(e); return(NULL)})
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("mkkm_mr" %in% multi_omic_methods) {
    if (is.null(extra_output)) extra_output <- list()
    if (is.null(extra_output$mkkm_mr_weights)) {
      extra_output$mkkm_mr_weights <- data.frame()
    }
    for (k in n_clusters) {
      for (lambda_i in mkkm_mr_lambda) {
        k_res <- tryCatch({
          # Optimize combined kernel
          optimal_kernel <- mkkm_mr(
            multi_omic_kernels, 
            k = k, 
            lambda = lambda_i, 
            tolerance = mkkm_mr_tolerance, 
            parallel = mvc_threads,
            use_mosek = mkkm_mr_mosek, 
            mosek_verbosity = mkkm_mr_mosek_verbosity
          )
          temp_res <- kernel_kmeans(
            K = optimal_kernel$K, 
            n_k = k, 
            algorithm = kkmeans_algorithm, 
            spectral_qr_refine = kkmeans_refine, 
            kernel_eigen_vectors = optimal_kernel$H, 
            max_iter = kkmeans_maxiter, 
            num_init = kkmeans_n_init, 
            tol = kkmeans_tol, 
            parallel = mvc_threads
          )
          temp_res <- data.frame(
            m = "mkkm_mr", 
            k = k, 
            mkkm_mr_lambda = lambda_i, 
            cluster = temp_res$clusters)
          extra_output$mkkm_mr_weights <- rbind(
            extra_output$mkkm_mr_weights, 
            data.frame(
              m = "mkkm_mr", 
              k = k, 
              mkkm_mr_lambda = lambda_i, 
              kernel_mix = paste0(
                names(multi_omic_kernels), 
                ":", 
                optimal_kernel$mu, 
                collapse = ";"
              )
            )
          )
          if (ncol(meta_data[[1]]) > 0) {
            temp_res <- cbind(meta_data[[1]], temp_res)
          }
          temp_res
        }, error = function(e) {warning(e); return(NULL)})
        if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
      }
    }
  }
  if ("ECMC" %in% multi_omic_methods) {
    for (k in n_clusters) {
      k_res <- tryCatch({
        # Find consensus kernels
        consensus_res <- ECMC(
          multi_omic_kernels, 
          a = ecmc_a, 
          b = ecmc_b,
          eps = ecmc_eps,
          maxiter = ecmc_maxiter,
          solver = ecmc_solver,
          parallel = mvc_threads)
        if (ecmc_mkkm_mr) {
          if (is.null(extra_output)) extra_output <- list()
          if (is.null(extra_output$mkkm_mr_weights)) {
            extra_output$mkkm_mr_weights <- data.frame()
          }
          # Optimize combined kernel
          optimal_kernel <- mkkm_mr(
            consensus_res$C, 
            k = k, 
            lambda = mkkm_mr_lambda, 
            tolerance = mkkm_mr_tolerance, 
            parallel = mvc_threads,
            use_mosek = mkkm_mr_mosek, 
            mosek_verbosity = mkkm_mr_mosek_verbosity
          )
          extra_output$mkkm_mr_weights <- rbind(
            extra_output$mkkm_mr_weights, 
            data.frame(
              m = "mkkm_mr", 
              k = k, 
              kernel_mix = paste0(
                names(multi_omic_kernels), 
                ":", 
                optimal_kernel$mu, 
                collapse = ";"
              )
            )
          )
          Kk <- optimal_kernel$K
          eigs_k <- optimal_kernel$H
        } else {
          Kk <- consensus_res$C_sum
          eigs_k <- NULL
        }
        temp_res <- kernel_kmeans(
          K = Kk, 
          n_k = k, 
          algorithm = kkmeans_algorithm, 
          spectral_qr_refine = kkmeans_refine, 
          kernel_eigen_vectors = eigs_k, 
          max_iter = kkmeans_maxiter, 
          num_init = kkmeans_n_init, 
          tol = kkmeans_tol, 
          parallel = mvc_threads
        )
        temp_res <- data.frame(
          m = "ecmc",
          k = k, 
          cluster = temp_res$clusters, 
          kernel_mix = paste(optimal_kernel$mu, collapse = ";"))
        if (ncol(meta_data[[1]]) > 0) {
          temp_res <- cbind(meta_data[[1]], temp_res)
        }
        temp_res
      }, error = function(e) {warning(e); return(NULL)})
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  out <- plyr::rbind.fill(res)
  if (!is.null(extra_output)) attributes(out)$extra_output <- extra_output
  attributes(out)$multi_omic = TRUE # used for formatting
  return(out)
}