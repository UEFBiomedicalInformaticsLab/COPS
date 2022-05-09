#' Multi-omic clustering via multi-view clustering or integration
#'
#' @param dat_list_clust 
#' @param non_data_cols 
#' @param multi_view_methods 
#' @param n_clusters 
#'
#' @return
#' @export
#' @importFrom ANF ANF
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
                                  standardize_data = FALSE,
                                  non_negativity_transform = rep_len("none", length(dat_list_clust)),
                                  icp_view_types = rep_len("gaussian", length(dat_list_clust)),
                                  icp_bayes_burnin = 1000,
                                  icp_bayes_draw = 1200,
                                  nmf_maxiter = 200,
                                  nmf_st.count = 20,
                                  nmf_n.ini = 30,
                                  nmf_ini.nndsvd = TRUE,
                                  mofa_scale_views = FALSE,
                                  mofa_likelihoods = rep_len("gaussian", length(dat_list_clust)), 
                                  mofa_convergence_mode = "medium",
                                  mofa_maxiter = 1000,
                                  mofa_environment = NULL,
                                  mofa_lib_path = NULL,
                                  anf_neighbors = 20,
                                  kernels = rep_len("linear", length(dat_list_clust)),
                                  kernels_center = TRUE,
                                  kernels_normalize = TRUE,
                                  kernel_gammas = rep_len(0.5, length(dat_list_clust)),
                                  pathway_networks = NULL,
                                  pamogk_restart = 0.7,
                                  kkmeans_maxiter = 100,
                                  kkmeans_n_init = 100,
                                  mkkm_mr_lambda = 1, 
                                  mkkm_mr_tolerance = 1e-8, 
                                  mkkm_mr_initialization = TRUE, 
                                  mkkm_mr_mosek = FALSE,
                                  ecmc_a = 1, 
                                  ecmc_b = 1, 
                                  ecmc_eps = 1e-6,
                                  ecmc_maxiter = 100,
                                  ecmc_solver = "MOSEK",
                                  ecmc_mkkm_mr = TRUE, 
                                  data_is_kernels = FALSE, 
                                  foldwise_zero_var_removal = TRUE,
                                  mvc_threads = 1,
                                  gene_id_list = NULL,
                                  ...) {
  if (foldwise_zero_var_removal & !data_is_kernels) {
    # Rare binary features such as some somatic mutations could end up missing 
    # in some of the folds. They cause issues and should be removed. 
    dat_list_clust <- lapply(dat_list_clust, function(x) x[,apply(x, 2, var) > 0])
  }
  if (standardize_data) {
    dat_list_clust <- lapply(dat_list_clust, scale)
  }
  if (any(multi_view_methods %in% c("kkmeans", "kkmeanspp", "mkkm_mr", "ECMC"))) {
    # In fold centering and normalization
    if (length(kernels_center) != length(dat_list_clust)) {
      kernels_center <- rep_len(kernels_center, length(dat_list_clust))
    }
    if (length(kernels_normalize) != length(dat_list_clust)) {
      kernels_normalize <- rep_len(kernels_normalize, length(dat_list_clust))
    }
  }
  if (data_is_kernels) {
    multi_omic_kernels <- dat_list_clust
    multi_omic_kernels[kernels_center] <- lapply(multi_omic_kernels[kernels_center],
                                                 center_kernel)
    multi_omic_kernels[kernels_normalize] <- lapply(multi_omic_kernels[kernels_normalize],
                                                    normalize_kernel)
  } else if (any(multi_view_methods %in% c("kkmeans", "kkmeanspp", "mkkm_mr", "ECMC"))) {
    # Pathway-based kernels need pathway networks
    if (any(kernels %in% c("PIK", "BWK", "PAMOGK")) &
        is.null(pathway_networks)) {
      w1 <- "No pathway networks specified for pathway kernel."
      w2 <- "Defaulting to KEGG pathways mapped to gene symbols."
      warning(paste(w1, w2))
      kegg_nets <- KEGG_networks()
      for (net_i in 1:length(kegg_nets)) {
        kegg_entrez <- gsub("hsa:", "", names(igraph::V(kegg_nets[[net_i]])))
        kegg_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, kegg_entrez, 
                                             "SYMBOL", "ENTREZID")
        kegg_nets[[net_i]] <- igraph::set.vertex.attribute(kegg_nets[[net_i]] , 
                                                           "name", 
                                                           value = kegg_symbol)
      }
      pathway_networks <- kegg_nets
    }
    if (any(kernels %in% c("BWK", "PAMOGK"))) {
      nw_weights <- node_betweenness_parallel(pathway_networks, mvc_threads)
      nw_weights <- lapply(nw_weights, sqrt)
    }
    # Construct kernels
    multi_omic_kernels <- list()
    for (i in 1:length(dat_list_clust)) {
      if (kernels[i] == "linear") {
        temp <- dat_list_clust[[i]] %*% t(dat_list_clust[[i]])
        if (kernels_center[i]) temp <- center_kernel(temp)
        if (kernels_normalize[i]) temp <- normalize_kernel(temp)
        multi_omic_kernels <- c(multi_omic_kernels, list(temp))
      } else if (kernels[i] == "gaussian") {
        temp <- exp(- kernel_gammas[i] * as.matrix(dist(dat_list_clust[[i]]))**2)
        multi_omic_kernels <- c(multi_omic_kernels, list(temp))
      } else if (kernels[i] %in% c("jaccard", "tanimoto")) {
        temp <- jaccard_matrix(t(dat_list_clust[[i]]))
        temp[is.nan(temp)] <- 0
        diag(temp) <- 1
        if (kernels_center[i]) temp <- center_kernel(temp)
        if (kernels_normalize[i]) temp <- normalize_kernel(temp)
        multi_omic_kernels <- c(multi_omic_kernels, list(temp))
      } else if (kernels[i] %in% c("PIK", "BWK", "PAMOGK")) {
        gene_col_ind <- as.integer(gsub("^dim", "", colnames(dat_list_clust[[i]])))
        colnames(dat_list_clust[[i]]) <- gene_id_list[[i]][gene_col_ind]
        if (kernels[i] == "PIK") {
          temp <- dat_list_clust[[i]]
          temp <- scale(temp, scale = TRUE) # z-scores
          temp <- PIK_from_networks(temp, pathway_networks, parallel = mvc_threads)
          names(temp) <- paste0(names(dat_list_clust)[i], "_", names(temp))
          temp <- lapply(temp, as.matrix)
          if (kernels_normalize[i]) {
            temp <- lapply(temp, normalize_kernel)
            temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
          }
          multi_omic_kernels <- c(multi_omic_kernels, temp)
        } else if (kernels[i] == "BWK") {
          temp <- t(dat_list_clust[[i]])
          temp <- lapply(nw_weights, function(w) weighted_linear_kernel(temp, w))
          names(temp) <- paste0(names(dat_list_clust)[i], "_", names(temp))
          temp <- temp[!sapply(temp, is.null)]
          temp <- lapply(temp, as.matrix)
          temp <- temp[which(sapply(temp, function(x) var(as.vector(x))) > 0)]
          if (kernels_center[i]) temp <- lapply(temp, center_kernel)
          if (kernels_normalize[i]) {
            temp <- lapply(temp, normalize_kernel)
            temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
          }
          multi_omic_kernels <- c(multi_omic_kernels, temp)
        } else if (kernels[i] == "PAMOGK") {
          temp <- dat_list_clust[[i]]
          temp <- scale(temp, scale = TRUE) # z-scores
          if (mvc_threads > 1) {
            rwr_threads <- mvc_threads
          } else {
            rwr_threads <- NULL
          }
          for (j in 1:length(pathway_networks)) {
            k_up <- dnet::dRWR(pathway_networks[[j]], 
                               normalise = "laplacian",
                               setSeeds = t(temp) > qnorm(0.975),
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
                multi_omic_kernels <- c(multi_omic_kernels, list(k_up))
              }
            }
            k_dn <- dnet::dRWR(pathway_networks[[j]], 
                               normalise = "laplacian",
                               setSeeds = t(temp) < qnorm(0.025),
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
                multi_omic_kernels <- c(multi_omic_kernels, list(k_dn))
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
  if("iClusterPlus" %in% multi_view_methods) {
    if (length(dat_list_clust) > 4) stop("iClusterPlus only supports up to four views.")
    if (length(dat_list_clust) >= 1) dt1 <- dat_list_clust[[1]] else dt1 <- NULL
    if (length(dat_list_clust) >= 2) dt2 <- dat_list_clust[[2]] else dt2 <- NULL
    if (length(dat_list_clust) >= 3) dt3 <- dat_list_clust[[3]] else dt3 <- NULL
    if (length(dat_list_clust) == 4) dt4 <- dat_list_clust[[4]] else dt4 <- NULL
    
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- iClusterPlus::iClusterPlus(dt1, dt2, dt3, dt4, 
                                               type = icp_view_types,
                                               K = k-1)
        k_res <- data.frame(m = "iClusterPlus", 
                            k = k,
                            cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], k_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if("iClusterBayes" %in% multi_view_methods) {
    if (length(dat_list_clust) > 6) stop("iClusterPlus only supports up to six views.")
    if (length(dat_list_clust) >= 1) dt1 <- dat_list_clust[[1]] else dt1 <- NULL
    if (length(dat_list_clust) >= 2) dt2 <- dat_list_clust[[2]] else dt2 <- NULL
    if (length(dat_list_clust) >= 3) dt3 <- dat_list_clust[[3]] else dt3 <- NULL
    if (length(dat_list_clust) >= 4) dt4 <- dat_list_clust[[4]] else dt4 <- NULL
    if (length(dat_list_clust) >= 5) dt5 <- dat_list_clust[[5]] else dt5 <- NULL
    if (length(dat_list_clust) == 6) dt6 <- dat_list_clust[[6]] else dt6 <- NULL
    
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- iClusterPlus::iClusterBayes(dt1, dt2, dt3, dt4, dt5, dt6, 
                                                type = icp_view_types,
                                                K = k-1,
                                                n.burnin = icp_bayes_burnin,
                                                n.draw = icp_bayes_draw)
        k_res <- data.frame(m = "iClusterBayes", 
                            k = k,
                            cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], k_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("IntNMF" %in% multi_view_methods) {
    for (i in 1:length(dat_list_clust)) {
      if (non_negativity_transform[i] == "logistic") {
        dat_list_clust[[i]] <- 1/(1 + exp(-dat_list_clust[[i]]))
      }
      if (non_negativity_transform[i] == "rank") {
        dat_list_clust[[i]] <- apply(dat_list_clust[[i]], 2, function(x) rank(x) / length(x))
      }
      if (non_negativity_transform[i] == "offset2") {
        dat_list_clust[[i]] <- dat_list_clust[[i]] + 2
      }
    }
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
      Sys.setenv(OMP_NUM_THREADS=mvc_threads)
      Sys.setenv(MKL_NUM_THREADS=mvc_threads)
      if (!is.null(mofa_lib_path)) {
        Sys.setenv(LD_LIBRARY_PATH = paste(mofa_lib_path, 
                                           Sys.getenv("LD_LIBRARY_PATH"), 
                                           sep = ":"))
      }
      if (!is.null(mofa_environment)) {
        reticulate::use_virtualenv(mofa_environment)
      }
      
      mofa_obj <- MOFA2::create_mofa(lapply(dat_list_clust, t))
      data_opts <- MOFA2::get_default_data_options(mofa_obj)
      data_opts$scale_views <- mofa_scale_views
      model_opts <- MOFA2::get_default_model_options(mofa_obj)
      mofa_view_names <- names(model_opts$likelihoods)
      model_opts$likelihoods <- mofa_likelihoods
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
      
      mofa_diss <- clustering_dissimilarity_from_data(mofa_embedding, distance_metric, correlation_method)
      mofa_cops_clust <- clustering_analysis(cbind(mofa_embedding, non_data_cols[[1]]), 
                                             n_clusters = n_clusters,
                                             clustering_dissimilarity = mofa_diss,
                                             ...)
      mofa_cops_clust
    }, error = function(e) return(NULL))
    if(!is.null(temp_res)) if(nrow(temp_res) > 1) res <- c(res, list(temp_res))
  }
  if ("ANF" %in% multi_view_methods) {
    #aff_list <- lapply(dat_list_clust, knn_g, k = ann_knn, jaccard_kernel = ann_jk)
    aff_dist <- lapply(dat_list_clust, dist)
    aff_dist <- lapply(aff_dist, as.matrix)
    aff_list <- lapply(aff_dist, ANF::affinity_matrix, k = anf_neighbors)
    aff_mat <- ANF::ANF(aff_list, K = anf_neighbors)
    for (k in n_clusters) {
      temp_res <- ANF::spectral_clustering(aff_mat, k)
      temp_res <- data.frame(m = "ANF", k = k, cluster = temp_res)
      k_res <- cbind(non_data_cols[[1]], temp_res)
      res <- c(res, list(k_res))
    }
  }
  if ("kkmeans" %in% multi_view_methods) {
    # Average kernel
    multi_omic_kernels <- Reduce('+', multi_omic_kernels) / length(multi_omic_kernels)
    approximation <- eigen(multi_omic_kernels, symmetric = TRUE)$vectors
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- kernel_kmeans_algorithm(multi_omic_kernels, 
                                            n_k = k, 
                                            init = apply(approximation[,1:k], 1, which.max), 
                                            maxiter = kkmeans_maxiter)
        temp_res <- data.frame(m = "kkmeans", k = k, cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("kkmeanspp" %in% multi_view_methods) {
    # Average kernel
    multi_omic_kernels <- Reduce('+', multi_omic_kernels) / length(multi_omic_kernels)
    for (k in n_clusters) {
      k_res <- tryCatch({
        temp_res <- kernel_kmeans(multi_omic_kernels, 
                                  k, 
                                  n_initializations = kkmeans_n_init, 
                                  maxiter = kkmeans_maxiter)
        temp_res <- data.frame(m = "kkmeanspp", k = k, cluster = temp_res$clusters)
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("mkkm_mr" %in% multi_view_methods) {
    for (k in n_clusters) {
      k_res <- tryCatch({
        # Optimize combined kernel
        optimal_kernel <- mkkm_mr(multi_omic_kernels, 
                                  k = k, 
                                  lambda = mkkm_mr_lambda, 
                                  tolerance = mkkm_mr_tolerance, 
                                  parallel = mvc_threads,
                                  use_mosek = mkkm_mr_mosek)
        if (mkkm_mr_initialization) {
          temp_res <- kernel_kmeans_algorithm(optimal_kernel$K, 
                                              n_k = k, 
                                              init = apply(optimal_kernel$H, 1, which.max), 
                                              maxiter = kkmeans_maxiter)
        } else {
          # Run k-means++ (random initialization)
          temp_res <- kernel_kmeans(optimal_kernel$K, 
                                    n_k = k, 
                                    n_initializations = kkmeans_n_init, 
                                    maxiter = kkmeans_maxiter,
                                    parallel = mvc_threads)
        }
        temp_res <- data.frame(m = "mkkm_mr", k = k, cluster = temp_res$clusters, 
                               kernel_mix = paste(optimal_kernel$mu, collapse = ";"))
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  if ("ECMC" %in% multi_view_methods) {
    for (k in n_clusters) {
      k_res <- tryCatch({
        # Find consensus kernels
        consensus_res <- ECMC(multi_omic_kernels, 
                              a = ecmc_a, 
                              b = ecmc_b,
                              eps = ecmc_eps,
                              maxiter = ecmc_maxiter,
                              solver = ecmc_solver,
                              parallel = mvc_threads)
        if (ecmc_mkkm_mr) {
          # Optimize combined kernel
          optimal_kernel <- mkkm_mr(consensus_res$C, 
                                    k = k, 
                                    lambda = mkkm_mr_lambda, 
                                    tolerance = mkkm_mr_tolerance, 
                                    parallel = mvc_threads,
                                    use_mosek = mkkm_mr_mosek)
          if (mkkm_mr_initialization) {
            temp_res <- kernel_kmeans_algorithm(optimal_kernel$K, 
                                                n_k = k, 
                                                init = apply(optimal_kernel$H, 1, which.max), 
                                                maxiter = kkmeans_maxiter)
          } else {
            # Run k-means++ (random initialization)
            temp_res <- kernel_kmeans(optimal_kernel$K, 
                                      n_k = k, 
                                      n_initializations = kkmeans_n_init, 
                                      maxiter = kkmeans_maxiter,
                                      parallel = mvc_threads)
          }
        } else {
          approximation <- eigen(consensus_res$C_sum, symmetric = TRUE)$vectors
          temp_res <- kernel_kmeans_algorithm(consensus_res$C_sum, 
                                              n_k = k, 
                                              init = apply(approximation[,1:k], 1, which.max), 
                                              maxiter = kkmeans_maxiter)
        }
        temp_res <- data.frame(m = "ecmc", k = k, cluster = temp_res$clusters, 
                               kernel_mix = paste(optimal_kernel$mu, collapse = ";"))
        cbind(non_data_cols[[1]], temp_res)
      }, error = function(e) return(NULL))
      if(!is.null(k_res)) if(nrow(k_res) > 1) res <- c(res, list(k_res))
    }
  }
  return(plyr::rbind.fill(res))
}