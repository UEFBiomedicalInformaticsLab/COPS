library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

dat_list <- list(expr = log2(tbrca_norm[combined_gene_filter,]+1))

batch_label <- brca_norm_batch
subtype_label <- brca_norm_subtypes_all

# Manage batch label column names if batch labels given
batch_label_names <- NULL
if (!is.null(batch_label)) {
  batch_label_names <- colnames(batch_label)
  if (is.null(batch_label_names)) batch_label_names <- "batch_label"
}

subtype_label_names <- NULL
if (!is.null(subtype_label)) {
  subtype_label_names <- colnames(subtype_label)
  if (is.null(subtype_label_names)) subtype_label_names <- "subtype_label"
}

# Convert data to data.table to optimize memory usage
for (i in 1:length(dat_list)) {
  id <- colnames(dat_list[[i]])
  
  dat_list[[i]] <- data.table::as.data.table(t(dat_list[[i]]))
  #data.table::setDT(dat_list[[i]])
  colnames(dat_list[[i]]) <- paste0("dim", 1:ncol(dat_list[[i]]))
  dat_list[[i]]$id <- id
  data.table::setkey(dat_list[[i]], id)
  
  # Add batch label(s) as separate column(s)
  if (!is.null(batch_label)) {
    if (is.null(dim(batch_label))) {
      if (!is.null(names(batch_label))) {
        batch_id <- names(batch_label)
      } else {
        # Assume data and batch labels are in the same order
        batch_id <- id
      }
      batch_table <- cbind(batch_label = as.character(batch_label), c())
    } else { # matrix batch_label
      batch_id <- rownames(batch_label)
      batch_table <- batch_label
    }
    #data.table::setDT(batch_table)
    batch_table <- data.table::as.data.table(batch_table)
    if (!all(id %in% batch_id)) {
      stop("All sample IDs in data do not match with batch label sample IDs.")
    }
    batch_table$id <- batch_id
    dat_list[[i]] <- plyr::join(dat_list[[i]], batch_table, by = "id")
  }
  
  if (!is.null(subtype_label)) {
    if (is.null(dim(subtype_label))) {
      if (!is.null(names(subtype_label))) {
        subtype_id <- names(subtype_label)
      } else {
        # Assume data and labels are in the same order
        subtype_id <- id
      }
      subtype_table <- cbind(subtype_label = as.character(subtype_label), c())
    } else { # matrix batch_label
      subtype_id <- rownames(subtype_label)
      subtype_table <- subtype_label
    }
    #data.table::setDT(subtype_table)
    subtype_table <- data.table::as.data.table(subtype_table)
    if (!all(id %in% subtype_id)) {
      warning("All subtype label sample IDs do not match with data.")
    }
    subtype_table$id <- subtype_id
    dat_list[[i]] <- plyr::join(dat_list[[i]], subtype_table, by = "id")
  }
}
cv_index <- COPS::cv_fold(dat_list, nfolds = NFOLDS, nruns = NRUNS)

clusters <- data.frame()
metrics <- data.frame()
chisqp <- data.frame()
bassoc <- data.frame()
sassoc <- data.frame()
for (i in unique(cv_index$expr$run)) {
  for (j in unique(cv_index$expr$fold)) {
    dat <- as.data.frame(merge(dat_list$expr, cv_index$expr[run == i & fold == j,], by = "id"))
    
    temp <- dat[grepl("^dim[0-9]+$", colnames(dat))]
    temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
    rownames(temp) <- dat$id
    
    clusters_ij <- data.frame()
    metrics_ij <- data.frame()
    for (k in 1:length(NCLUSTERS)) {
      # SC3 only accepts input in the form of SingleCellExperiment 
      workaround <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(scale(temp, center = TRUE, scale = FALSE)))) # centering necessary?
      SummarizedExperiment::rowData(workaround)$feature_symbol <- colnames(temp)
      workaround <- SC3::sc3(workaround, ks = NCLUSTERS[k], gene_filter = FALSE, n_cores = NULL)
      clust_k <- cutree(workaround@metadata$sc3$consensus[[1]]$hc, NCLUSTERS[k])
      silh_k <- workaround@metadata$sc3$consensus[[1]]$silhouette
      metrics_ij <- rbind(metrics_ij, data.frame(m = "sc3", k = NCLUSTERS[k], 
                                           metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
      clusters_ij <- rbind(clusters_ij, data.frame(id = rownames(temp), m = "sc3", 
                                             k = NCLUSTERS[k], cluster = clust_k))
    }
    
    out_list_ij <- list()
    
    # Combine outputs with metadata
    out_list_ij$clusters <- plyr::join(dat[!grepl("^dim[0-9]+$", colnames(dat))], 
                                    clusters_ij, 
                                    by = "id")
    out_list_ij$clusters <- out_list_ij$clusters[!is.na(out_list_ij$clusters$cluster),] # potential issue with reference fold missing later
    
    out_list_ij$metrics <- metrics_ij
    
    # Insert meta data (if present)
    out_list_ij$metrics$run <- dat$run[1]
    out_list_ij$metrics$fold <- dat$fold[1]
    out_list_ij$metrics$datname <- dat$datname[1]
    out_list_ij$metrics$drname <- dat$drname[1]
    
    # Pearson's chi-squared test
    f1 <- function(x, c1, c2) {
      temp <- tryCatch(suppressWarnings(chisq.test(x[[c1]], x[[c2]])), error = function(e) NULL)
      if (!is.null(temp)) {
        temp <- data.frame(p = NA)
      } else {
        temp <- data.frame(p = temp$p.value)
      }
      temp$run <- x$run[1]
      temp$fold <- x$fold[1]
      temp$datname <- x$datname[1]
      temp$drname <- x$drname[1]
      temp$k <- x$k[1]
      temp$m <- x$m[1]
      temp$label <- c2
      return(temp)
    }
    # Other contingency table based metrics
    f2 <- function(x, c1, c2) {
      temp <- data.frame(nmi = aricode::NMI(x[[c1]], x[[c2]]),
                         ari = aricode::ARI(x[[c1]], x[[c2]]))
      temp$run <- x$run[1]
      temp$fold <- x$fold[1]
      temp$datname <- x$datname[1]
      temp$drname <- x$drname[1]
      temp$k <- x$k[1]
      temp$m <- x$m[1]
      temp$label <- c2
      return(temp)
    }
    
    if (!is.null(batch_label_names)) {
      batch_label_chisq <- list()
      batch_label_assoc <- list()
      for (a in 1:length(batch_label_names)) {
        batch_label_chisq[[a]] <- Reduce("rbind", 
                                         lapply(split(out_list_ij$clusters, 
                                                      out_list_ij$clusters[c("k", "m")]), 
                                                f1, c1 = "cluster", c2 = batch_label_names[a]))
        batch_label_assoc[[a]] <- Reduce("rbind", 
                                         lapply(split(out_list_ij$clusters, 
                                                      out_list_ij$clusters[c("k", "m")]), 
                                                f2, c1 = "cluster", c2 = batch_label_names[a]))
      }
      out_list_ij$chisq_pval <- Reduce("rbind", batch_label_chisq)
      out_list_ij$batch_association <- Reduce("rbind", batch_label_assoc)
    }
    
    if (!is.null(subtype_label_names)) {
      subtype_label_assoc <- list()
      for (b in 1:length(subtype_label_names)) {
        temp <- out_list_ij$clusters[!is.na(out_list_ij$clusters[[subtype_label_names[b]]]),]
        subtype_label_assoc[[b]] <- Reduce("rbind", 
                                           lapply(split(temp, 
                                                        temp[c("k", "m")]), 
                                                  f2, c1 = "cluster", c2 = subtype_label_names[b]))
      }
      out_list_ij$subtype_association <- Reduce("rbind", subtype_label_assoc)
    }
    clusters <- rbind(clusters, out_list_ij$clusters)
    metrics <- rbind(metrics, out_list_ij$metrics)
    chisqp <- rbind(chisqp, out_list_ij$chisq_pval)
    bassoc <- rbind(bassoc, out_list_ij$batch_association)
    sassoc <- rbind(sassoc, out_list_ij$subtype_association)
  }
}


foreach::registerDoSEQ()
# Clustering stability evaluation
clust_stability <- COPS::stability_eval(as.data.frame(clusters), by = c("run", "k", "m"), by2 = "fold")

# Survival evaluation
clust_survival <- COPS::survival_evaluation(dat_survival, clusters, by = c("run", "fold", "k", "m"))

# Gene module correlation evaluation
clust_gm_score <- COPS::module_evaluation(clusters, module_eigs = MEs, module_cor_threshold = 0.25,
                                          by = c("run", "fold", "k", "m"))

# Return
out <- list(clusters = clusters, 
            internal_metrics = metrics,
            chisq_pval = chisqp,
            batch_association = bassoc,
            subtype_association = sassoc,
            stability = clust_stability)
out$survival <- clust_survival
out$modules <- clust_gm_score

scores_sc3 <- COPS::clusteval_scoring(out, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_sc3$all, paste0(path_intermediate_results, "/brca/sc3/scores.csv"))
write.csv(clusters, gzfile(paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz")))







