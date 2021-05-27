library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

args <- commandArgs(trailingOnly = TRUE)

cv_run <- as.integer(args[1])
cv_fold <- as.integer(args[2])

path <- paste0(path_intermediate_results, "/brca/sc3/")

cv_index <- read.csv(paste0(path, "cv_index.csv"), row.names = 1, header = TRUE)

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

dat <- as.data.frame(merge(dat_list$expr, cv_index[cv_index$run == cv_run & cv_index$fold == cv_fold,], by = "id"))

temp <- dat[grepl("^dim[0-9]+$", colnames(dat))]
temp <- temp[sapply(temp, function(x) all(!is.na(x)))]
rownames(temp) <- dat$id

# Alternative
clusters <- data.frame()

# Metrics collected to data.frame
metrics <- data.frame()

for (j in 1:length(n_clusters)) {
  # SC3 only accepts input in the form of SingleCellExperiment 
  hack <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(temp)))
  SummarizedExperiment::rowData(hack)$feature_symbol <- colnames(temp)
  hack <- SC3::sc3(hack, ks = n_clusters[j], gene_filter = FALSE, n_cores = NULL)
  clust_k <- cutree(hack@metadata$sc3$consensus[[1]]$hc, n_clusters[j])
  silh_k <- hack@metadata$sc3$consensus[[1]]$silhouette
  metrics <- rbind(metrics, data.frame(m = "sc3", k = n_clusters[j], 
                                       metric = "Silhouette", value = mean(silh_k[,"sil_width"])))
  clusters <- rbind(clusters, data.frame(id = rownames(temp), m = "sc3", 
                                         k = n_clusters[j], cluster = clust_k))
}

out_list <- list()

# Combine outputs with metadata
out_list$clusters <- plyr::join(dat[!grepl("^dim[0-9]+$", colnames(dat))], 
                                clusters, 
                                by = "id")
out_list$clusters <- out_list$clusters[!is.na(out_list$clusters$cluster),] # potential issue with reference fold missing later

out_list$metrics <- metrics

# Insert meta data (if present)
out_list$metrics$run <- dat$run[1]
out_list$metrics$fold <- dat$fold[1]
out_list$metrics$datname <- dat$datname[1]
out_list$metrics$drname <- dat$drname[1]

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
  for (i in 1:length(batch_label_names)) {
    batch_label_chisq[[i]] <- Reduce("rbind", 
                                     lapply(split(out_list$clusters, 
                                                  out_list$clusters[c("k", "m")]), 
                                            f1, c1 = "cluster", c2 = batch_label_names[i]))
    batch_label_assoc[[i]] <- Reduce("rbind", 
                                     lapply(split(out_list$clusters, 
                                                  out_list$clusters[c("k", "m")]), 
                                            f2, c1 = "cluster", c2 = batch_label_names[i]))
  }
  out_list$chisq_pval <- Reduce("rbind", batch_label_chisq)
  out_list$batch_association <- Reduce("rbind", batch_label_assoc)
}

if (!is.null(subtype_label_names)) {
  subtype_label_assoc <- list()
  for (i in 1:length(subtype_label_names)) {
    temp <- out_list$clusters[!is.na(out_list$clusters[[subtype_label_names[i]]]),]
    subtype_label_assoc[[i]] <- Reduce("rbind", 
                                       lapply(split(temp, 
                                                    temp[c("k", "m")]), 
                                              f2, c1 = "cluster", c2 = subtype_label_names[i]))
  }
  out_list$subtype_association <- Reduce("rbind", subtype_label_assoc)
}

if (!file.exists(paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz"))) {
  write.table(out_list$clusters, paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz"), append = FALSE, sep = ",", dec = ".", row.names = FALSE)
} else {
  write.table(out_list$clusters, paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}
if (!file.exists(paste0(path_intermediate_results, "/brca/sc3/metrics.csv.gz"))) {
  write.table(out_list$metrics, paste0(path_intermediate_results, "/brca/sc3/metrics.csv.gz"), append = FALSE, sep = ",", dec = ".", row.names = FALSE)
} else {
  write.table(out_list$metrics, paste0(path_intermediate_results, "/brca/sc3/metrics.csv.gz"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}
if (!file.exists(paste0(path_intermediate_results, "/brca/sc3/chisq_pval.csv.gz"))) {
  write.table(out_list$chisq_pval, paste0(path_intermediate_results, "/brca/sc3/chisq_pval.csv.gz"), append = FALSE, sep = ",", dec = ".", row.names = FALSE)
} else {
  write.table(out_list$chisq_pval, paste0(path_intermediate_results, "/brca/sc3/chisq_pval.csv.gz"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}
if (!file.exists(paste0(path_intermediate_results, "/brca/sc3/bassoc.csv.gz"))) {
  write.table(out_list$batch_association, paste0(path_intermediate_results, "/brca/sc3/bassoc.csv.gz"), append = FALSE, sep = ",", dec = ".", row.names = FALSE)
} else {
  write.table(out_list$batch_association, paste0(path_intermediate_results, "/brca/sc3/bassoc.csv.gz"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}
if (!file.exists(paste0(path_intermediate_results, "/brca/sc3/sassoc.csv.gz"))) {
  write.table(out_list$subtype_association, paste0(path_intermediate_results, "/brca/sc3/sassoc.csv.gz"), append = FALSE, sep = ",", dec = ".", row.names = FALSE)
} else {
  write.table(out_list$subtype_association, paste0(path_intermediate_results, "/brca/sc3/sassoc.csv.gz"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}










