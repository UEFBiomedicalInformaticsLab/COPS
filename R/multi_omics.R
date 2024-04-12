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
#' @param kernels Character vector of kernel names to use for different views. 
#'   See details. 
#' @param kernels_center Logical vector specifying which kernels should be 
#'   centered. Repeated for each view if length 1.
#' @param kernels_normalize Logical vector specifying which kernels should be 
#'   normalized Repeated for each view if length 1.
#' @param kernels_scale_norm Logical vector specifying which kernels should be 
#'   scaled to unit F-norm. Repeated for each view if length 1.
#' @param kernel_gammas Numeric vector specifying gamma for the gaussian kernel. 
#' @param pathway_networks List of \code{igraph} objects containing pathway 
#'   networks. Required for pathway kernels. 
#' @param pathway_node_betweenness_endpoints whether to include shortest path 
#'   endpoints in betweenness. Including it results in more non-zero weights 
#'   in BWK and PAMOGK. 
#' @param pamogk_restart Restart probability for PAMOGK RWR. 
#' @param pamogk_seeds Seed selection strategy for PAMOGK RWR, one of: 
#'   "discrete", "continuous", "threshold". See details below. 
#' @param pamogk_seed_under_threshold z-score threshold for under-expressed. 
#' @param pamogk_seed_over_threshold z-score threshold for over-expressed. 
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
#' @param data_is_kernels If set, input data is assumed to be kernel matrices. 
#'   Otherwise kernels are computed based on input data and the 
#'   \code{kernels} parameter. 
#' @param zero_var_removal If set, removes all zero variance features 
#'   from the data. It is called fold-wise, because this is assumed to be run 
#'   inside CV. 
#' @param mvc_threads Number of threads to use for supported operations. 
#' @param gene_id_list List of gene/feature names for each view. If set, matches 
#'   pipeline standardized feature names ("dim1", "dim2", ...) to names on the list. 
#'   Required for pathway kernels. 
#' @param preprocess_data If the input data has already been processed by the 
#'   \code{\link{COPS}}-pipeline, this should be disabled. 
#' @param ... Arguments are passed to \code{\link{clustering_analysis}} when using MOFA. 
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
#' Supported kernels: 
#' \itemize{
#'   \item "linear" - Linear kernel based on standard dot product. 
#'   \item "gaussian", "rbf" - Gaussian kernel, a.k.a. radial basis function.
#'   \item "jaccard" - Kernel based on Jaccard index. Used for binary features. 
#'   \item "tanimoto" - For now, this is identical to "jaccard".
#'   \item "BWK" - Betweenness Weighted Kernel. Uses pathway networks to compute 
#'     betweenness centrality which is used to weight features in a linear 
#'     pathway kernels. 
#'   \item "PAMOGK" - PAthway Multi-Omics Graph Kernel (Tepeli et al. 2021). 
#'     Uses z-scores, RWR and shortest paths in pathway networks to create 
#'     pathway kernels. 
#'   \item "PIK" - Pathway Induced Kernel (Manica et al. 2019). Uses pathway 
#'     network adjacency matrices (specifically normalized Laplacians) to define 
#'     pathway kernels. 
#' }
#' Please note that for pathway kernels, the input data must always be mapped to 
#' genes and that the names must match with the gene names in the pathways. 
#' The default set of pathways is KEGG molecular pathways with gene symbols. 
#' 
#' PAMOGK RWR seed weight options:
#' \itemize{
#'   \item "discrete" - 1 if |z| > t, 0 otherwise. 
#'   \item "continuous" - z
#'   \item "threshold" - z if |z| > t, 0 otherwise
#' }
#' Regardless of the option, the seeds are divided into two sets based on the 
#' sign of the z-score. Each set has a separate smoothing step and the end 
#' result is two different kernels per pathway per omic. This impacts the 
#' RWR label smoothing by changing the initial distribution. 
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
    kernels = rep_len("linear", length(dat_list)),
    kernels_center = TRUE,
    kernels_normalize = TRUE,
    kernels_scale_norm = FALSE,
    kernel_gammas = rep_len(0.5, length(dat_list)),
    pathway_networks = NULL,
    pathway_node_betweenness_endpoints = TRUE, 
    pathway_first_shortest_path = FALSE, 
    pamogk_restart = 0.7,
    pamogk_seeds = "discrete", 
    pamogk_seed_under_threshold = qnorm(0.025), 
    pamogk_seed_over_threshold = qnorm(0.975), 
    kkmeans_algorithm = "qr_decomposition", 
    kkmeans_refine = TRUE, 
    kkmeans_maxiter = 100, 
    kkmeans_n_init = 100, 
    kkmeans_tol = 1e-8, 
    mkkm_mr_lambda = 1, 
    mkkm_mr_tolerance = 1e-8, 
    mkkm_mr_mosek = FALSE, 
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
    # Rare binary features such as some somatic mutations could end up missing 
    # in some of the folds. They cause issues and should be removed. 
    dat_list <- lapply(dat_list, function(x) x[,apply(x, 2, var) > 0])
  }
  if (standardize_data) {
    dat_list <- lapply(dat_list, scale)
  }
  if (any(multi_omic_methods %in% c("kkmeans", "kkmeanspp", "mkkm_mr", "ECMC"))) {
    # In fold centering and normalization
    if (length(kernels_center) != length(dat_list)) {
      kernels_center <- rep_len(kernels_center, length(dat_list))
    }
    if (length(kernels_normalize) != length(dat_list)) {
      kernels_normalize <- rep_len(kernels_normalize, length(dat_list))
    }
    if (length(kernels_scale_norm) != length(dat_list)) {
      kernels_scale_norm <- rep_len(kernels_scale_norm, length(dat_list))
    }
  }
  if (data_is_kernels) {
    multi_omic_kernels <- dat_list
    multi_omic_kernels[kernels_center] <- lapply(
      multi_omic_kernels[kernels_center],
      center_kernel)
    multi_omic_kernels[kernels_normalize] <- lapply(
      multi_omic_kernels[kernels_normalize],
      normalize_kernel)
    multi_omic_kernels[kernels_scale_norm] <- lapply(
      multi_omic_kernels[kernels_scale_norm],
      scale_kernel_norm)
  } else if (any(multi_omic_methods %in% c("kkmeans", "kkmeanspp", "mkkm_mr", "ECMC"))) {
    # Pathway-based kernels need pathway networks
    if (any(kernels %in% c("PIK", "BWK", "PAMOGK")) &
        is.null(pathway_networks)) {
      w1 <- "No pathway networks specified for pathway kernel."
      w2 <- "Defaulting to KEGG pathways mapped to gene symbols."
      warning(paste(w1, w2))
      kegg_nets <- KEGG_networks()
      for (net_i in 1:length(kegg_nets)) {
        kegg_entrez <- gsub("hsa:", "", names(igraph::V(kegg_nets[[net_i]])))
        kegg_symbol <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db, 
          kegg_entrez, 
          "SYMBOL", 
          "ENTREZID")
        kegg_nets[[net_i]] <- igraph::set.vertex.attribute(
          kegg_nets[[net_i]], 
          "name", 
          value = kegg_symbol)
      }
      pathway_networks <- kegg_nets
    }
    if (any(kernels %in% c("BWK", "PAMOGK"))) {
      nw_weights <- node_betweenness_parallel(
        pathway_networks, 
        mvc_threads, 
        pathway_node_betweenness_endpoints = pathway_node_betweenness_endpoints, 
        pathway_first_shortest_path = pathway_first_shortest_path)
      nw_weights <- lapply(nw_weights, sqrt)
    }
    # Construct kernels
    multi_omic_kernels <- list()
    for (i in 1:length(dat_list)) {
      if (kernels[i] == "linear") {
        temp <- dat_list[[i]] %*% t(dat_list[[i]])
        if (kernels_center[i]) temp <- center_kernel(temp)
        if (kernels_normalize[i]) temp <- normalize_kernel(temp)
        if (kernels_scale_norm[i]) temp <- scale_kernel_norm(temp)
        temp <- list(temp)
        names(temp) <- names(dat_list)[i]
        multi_omic_kernels <- c(multi_omic_kernels, temp)
      } else if (kernels[i] %in% c("gaussian", "rbf")) {
        temp <- exp(- kernel_gammas[i] * as.matrix(dist(dat_list[[i]]))**2)
        temp <- list(temp)
        names(temp) <- names(dat_list)[i]
        multi_omic_kernels <- c(multi_omic_kernels, temp)
      } else if (kernels[i] %in% c("jaccard", "tanimoto")) {
        temp <- jaccard_matrix(t(dat_list[[i]]))
        temp[is.nan(temp)] <- 0
        diag(temp) <- 1
        if (kernels_center[i]) temp <- center_kernel(temp)
        if (kernels_normalize[i]) temp <- normalize_kernel(temp)
        if (kernels_scale_norm[i]) temp <- scale_kernel_norm(temp)
        temp <- list(temp)
        names(temp) <- names(dat_list)[i]
        multi_omic_kernels <- c(multi_omic_kernels, temp)
      } else if (kernels[i] %in% c("PIK", "BWK", "PAMOGK")) {
        gene_col_ind <- as.integer(gsub("^dim", "", colnames(dat_list[[i]])))
        temp <- dat_list[[i]]
        colnames(temp) <- gene_id_list[[i]][gene_col_ind]
        if (kernels[i] == "PIK") {
          temp <- scale(temp, scale = TRUE) # z-scores
          temp <- PIK_from_networks(temp, pathway_networks, parallel = mvc_threads)
          names(temp) <- paste0(names(dat_list)[i], "_", names(temp))
          temp <- lapply(temp, as.matrix)
          if (kernels_normalize[i]) {
            temp <- lapply(temp, normalize_kernel)
            temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
          }
          if (kernels_scale_norm[i]) temp <- lapply(temp, scale_kernel_norm)
          multi_omic_kernels <- c(multi_omic_kernels, temp)
        } else if (kernels[i] == "BWK") {
          temp <- t(temp)
          temp <- lapply(nw_weights, function(w) weighted_linear_kernel(temp, w))
          names(temp) <- paste0(names(dat_list)[i], "_", names(temp))
          temp <- temp[!sapply(temp, is.null)]
          temp <- lapply(temp, as.matrix)
          temp <- temp[which(sapply(temp, function(x) var(as.vector(x))) > 0)]
          if (kernels_center[i]) temp <- lapply(temp, center_kernel)
          if (kernels_normalize[i]) {
            temp <- lapply(temp, normalize_kernel)
            temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
          }
          if (kernels_scale_norm[i]) temp <- lapply(temp, scale_kernel_norm)
          multi_omic_kernels <- c(multi_omic_kernels, temp)
        } else if (kernels[i] == "PAMOGK") {
          temp <- scale(temp, scale = TRUE) # z-scores
          if (mvc_threads > 1) {
            rwr_threads <- mvc_threads
          } else {
            rwr_threads <- NULL
          }
          
          if (pamogk_seeds == "discrete") {
            seed_up <- t(temp) > pamogk_seed_over_threshold
            seed_dn <- t(temp) < pamogk_seed_under_threshold
          } else if (pamogk_seeds == "continuous") {
            seed_up <- t(temp)
            seed_dn <- -t(temp)
            seed_up[seed_up < 0] <- 0
            seed_dn[seed_dn < 0] <- 0
          } else if (pamogk_seeds == "threshold") {
            seed_up <- seed_dn <- t(temp)
            seed_up[seed_up < pamogk_seed_over_threshold] <- 0
            seed_dn[seed_dn > pamogk_seed_under_threshold] <- 0
            seed_dn <- -seed_dn
          } else {
            stop(paste0(
              "pamogk_seeds option '", 
              pamogk_seeds, 
              "' not recognized. "
            ))
          }
          up_gene_ind <- seed_up > 0
          dn_gene_ind <- seed_dn > 0
          for (j in 1:length(pathway_networks)) {
            pw_gene_ind <- rownames(seed_up) %in% names(igraph::V(pathway_networks[[j]]))
            if (!any(pw_gene_ind)) next
            any_up_gene <- apply(up_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
            # Skip pathways where only one sample has seeds
            if (sum(any_up_gene)>1) {
              k_up <- dnet::dRWR(
                pathway_networks[[j]], 
                normalise = "laplacian",
                setSeeds = seed_up,
                restart = pamogk_restart,
                normalise.affinity.matrix = "none",
                parallel = mvc_threads > 1,
                multicores = rwr_threads)
              rownames(k_up) <- names(igraph::V(pathway_networks[[j]]))
              colnames(k_up) <- rownames(temp)
              k_up <- weighted_linear_kernel(as.matrix(k_up), nw_weights[[j]])
              if (!is.null(k_up)) {
                if (var(as.vector(k_up)) > 0) {
                  if (kernels_center[i]) k_up <- center_kernel(k_up)
                  if (kernels_normalize[i]) {
                    k_up <- normalize_kernel(k_up)
                    k_up[is.na(k_up)] <- 0
                  }
                  if (kernels_scale_norm[i]) k_up <- scale_kernel_norm(k_up)
                  multi_omic_kernels <- c(multi_omic_kernels, list(k_up))
                }
              }
            }
            # Same check for down genes
            any_dn_gene <- apply(dn_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
            if (sum(any_dn_gene)>1) {
              k_dn <- dnet::dRWR(
                pathway_networks[[j]], 
                normalise = "laplacian",
                setSeeds = seed_dn,
                restart = pamogk_restart,
                normalise.affinity.matrix = "none",
                parallel = mvc_threads > 1,
                multicores = rwr_threads)
              rownames(k_dn) <- names(igraph::V(pathway_networks[[j]]))
              colnames(k_dn) <- rownames(temp)
              k_dn <- weighted_linear_kernel(as.matrix(k_dn), nw_weights[[j]])
              if (!is.null(k_dn)) {
                if (var(as.vector(k_dn)) > 0) {
                  if (kernels_center[i]) k_dn <- center_kernel(k_dn)
                  if (kernels_normalize[i]) {
                    k_dn <- normalize_kernel(k_dn)
                    k_dn[is.na(k_dn)] <- 0
                  }
                  if (kernels_scale_norm[i]) k_dn <- scale_kernel_norm(k_dn)
                  multi_omic_kernels <- c(multi_omic_kernels, list(k_dn))
                }
              }
            }
          }
        } 
      } else {
        stop(paste0("Kernel \"", kernels[i], "\" is not supported."))
      }
    }
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
            use_mosek = mkkm_mr_mosek)
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
            use_mosek = mkkm_mr_mosek)
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