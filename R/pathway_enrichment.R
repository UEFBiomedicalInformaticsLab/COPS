#' Transform a gene-level data matrix into path-level information
#' 
#' This is a utility function which wraps several methods to extract pathway-based features from gene feature data
#' 
#' @param x a numeric matrix representing gene features profiles. Genes on the rows and samples on the columns.
#' @param enrichment_method options: "DiffRank", "GSVA", "RWRFGSEA"
#' @param gene_set_list list of gene sets (character vectors) corresponding to pathway annotations or gene sets, 
#'   if \code{NULL} gene sets are retrieved with \code{\link[msigdbr]{msigdbr}}.
#' @param batch_label_pw batch labels for batch-wise enrichment, normal enrichment performed if NULL
#' @param min_size a numeric value indicating the minimum size of gene sets included
#' @param max_size a numeric value indicating the maximum size of gene sets included
#' @param parallel a numeric value indicating the number of processors to use when doing the calculations in parallel.
#' @param verbose controls verbosity
#' @param gene_key_x if \code{is.null(gene_set_list)} and \code{x} rownames are not gene symbols, this specifies the column name in 
#'   \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} to translate pathway gene symbols to. The default value results in gene symbol based gene 
#'   sets when \code{is.null(gene_set_list)}. 
#' @param gsva_kcdf distribution name for \code{\link[GSVA]{gsva}} empirical probability distribution kernel.
#' @param gs_subcats if \code{is.null(gene_set_list)}, this character vector indicates \code{\link[msigdbr]{msigdbr}}) gene set subcategory 
#'   names that are included in the analysis.
#' @param ... extra arguments are passed to \code{\link[COPS]{RWRFGSEA}}
#' 
#' @return a list of \code{data.frame}s corresponding to the transformed features 
#'   based on the selected gene sets (only supports GO, KEGG and REACTOME at the moment)
#' @export
#' 
#' @examples library(parallel)
#' library(COPS)
#' 
#' ## GSVA example
#' ad_gsva <- genes_to_pathways(
#'     as.matrix(ad_ge_micro_zscore), 
#'     "GSVA", 
#'     parallel = 2, 
#'     gene_key_x = "ENSEMBL", 
#'     gs_subcats = "CP:KEGG")
#' # batch-wise
#' ad_gsva <- genes_to_pathways(
#'     as.matrix(ad_ge_micro_zscore), 
#'     "GSVA", 
#'     batch_label_pw = ad_studies$GSE, 
#'     parallel = 2, 
#'     gene_key_x = "ENSEMBL", 
#'     gs_subcats = "CP:KEGG")
#' 
#' ## DiffRank example with batch-wise enrichment
#' ad_diffrank <- genes_to_pathways(
#'     ad_ge_micro_zscore, 
#'     "DiffRank", 
#'     parallel = 2, 
#'     gene_key_x = "ENSEMBL", 
#'     gs_subcats = "CP:KEGG")
#' # batch-wise
#' ad_diffrank <- genes_to_pathways(
#'     ad_ge_micro_zscore, 
#'     "DiffRank", 
#'     batch_label_pw = ad_studies$GSE, 
#'     parallel = 2, 
#'     gene_key_x = "ENSEMBL", 
#'     gs_subcats = "CP:KEGG")
#' 
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr filter
#' @importFrom msigdbr msigdbr
genes_to_pathways <- function(
    x, 
    enrichment_method = "GSVA",
    gene_set_list = NULL,
    batch_label_pw = NULL, 
    min_size = 5, 
    max_size = 200, 
    parallel = 1,
    verbose = FALSE,
    gene_key_x = "SYMBOL",
    gs_subcats = c("GO:BP", "GO:MF", "CP:KEGG", "CP:REACTOME"),
    gsva_kcdf = "Gaussian",
    ...
) {
  #x <- t(x)
  if (is.null(gene_set_list)) {
    # extract pathways information from msigdb (https://www.gsea-msigdb.org/)
    db_annots <- msigdbr::msigdbr(species = "Homo sapiens")
    db_annots <- dplyr::filter(
      db_annots, 
      grepl(
        paste0("^", paste(gs_subcats, collapse = "$|^"), "$"), 
        gs_subcat
      )
    )
    
    if (gene_key_x != "SYMBOL") {
      db_annots$gene_id <- suppressMessages(
        as.character(
          AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db, 
            db_annots$human_gene_symbol, 
            column = gene_key_x, 
            keytype = "SYMBOL"
          )
        )
      )
    } else {
      db_annots$gene_id <- db_annots$human_gene_symbol
    }
    
    gene_set_list <- lapply(
      split(db_annots, db_annots$gs_name), 
      function(a) a$gene_id)
  }
  
  gene_set_list <- gene_set_list[
    which(
      sapply(gene_set_list, length) <= max_size & 
        sapply(gene_set_list, length) >= min_size
    )
  ]
  
  if(!is.null(batch_label_pw)) {
    # Batch-wise enrichment analysis
    if (any(table(batch_label_pw) < 10)) {
      warning(paste(
        "Batch-wise pathway enrichment results may be unreliable due to low", 
        "number of samples in some batches."))
    }
    if (any(table(batch_label_pw) < 2)) {
      stop(paste(
        "Batch-wise pathway enrichment failed, all batches must have at least", 
        "2 samples."))
    }
    batch_dat_list <- lapply(
      unique(batch_label_pw), 
      function(a) x[, which(a == as.character(batch_label_pw)), drop = FALSE])
    
    # Apply the specified enrichment methods
    if (tolower(enrichment_method) == "gsva") {
      enriched_dat <- lapply(
        batch_dat_list, 
        gsva_wrapper, 
        gene_set_list = gene_set_list, 
        verbose = verbose, 
        parallel = parallel,
        kcdf = gsva_kcdf,
        min_size = min_size, 
        max_size = max_size)
    } else if (tolower(enrichment_method) == "diffrank") {
      enriched_dat <- lapply(
        batch_dat_list, 
        DiffRank, 
        gene_set_list = gene_set_list, 
        parallel = parallel)
    } else if (gsub("[-_]", "", tolower(enrichment_method)) == "rwrfgsea") {
      enriched_dat <- lapply(
        batch_dat_list, 
        RWRFGSEA, 
        gene_set_list = gene_set_list, 
        parallel = parallel, 
        verbose = verbose, 
        ...)
    } else {
      stop(paste("Unsupported pathway enrichment method:", enrichment_method))
    }
    enriched_dat <- t(Reduce(plyr::rbind.fill.matrix, lapply(enriched_dat, t)))
  } else {
    # Regular enrichment analysis
    if (tolower(enrichment_method) == "gsva") {
      enriched_dat <- suppressWarnings(
        gsva_wrapper(
          x = x, 
          gene_set_list = gene_set_list, 
          verbose = verbose, 
          parallel = parallel, 
          kcdf = gsva_kcdf,
          min_size = min_size, 
          max_size = max_size))
    } else if (tolower(enrichment_method) == "diffrank") {
      enriched_dat <- COPS::DiffRank(x, gene_set_list, parallel)
    } else if (gsub("[-_]", "", tolower(enrichment_method)) == "rwrfgsea") {
      enriched_dat <- RWRFGSEA(
        x, 
        gene_set_list = gene_set_list, 
        parallel = parallel, 
        verbose = verbose, 
        ...)
    } else {
      stop(paste("Unsupported pathway enrichment method:", enrichment_method))
    }
  }
  
  gs_category <- sapply(strsplit(rownames(enriched_dat), "_"), function(a) a[[1]])
  
  out <- lapply(unique(gs_category), function(a) enriched_dat[gs_category == a,])
  names(out) <- unique(gs_category) # converted to dat_name by subsample_dimred
  
  # Format output
  if (enrichment_method ==  "GSVA") {
    names(out) <- paste0(names(out), "_GSVA")
  } else if (enrichment_method ==  "DiffRank") {
    names(out) <- paste0(names(out), "_DiffRank")
  } else if (enrichment_method ==  "RWRFGSEA") {
    names(out) <- paste0(names(out), "_RWRFGSEA")
  }
  
  out <- out[sapply(out, nrow) > 0]
  
  return(out)
}

#' @importFrom GSVA gsva
#' @importFrom BiocParallel MulticoreParam SerialParam
gsva_wrapper <- function(
    x, 
    gene_set_list = NULL,
    min_size = 5, 
    max_size = 200, 
    parallel = 1,
    verbose = FALSE,
    kcdf = "Gaussian"
) {
  if (parallel > 1) {
    bioc_parallel <- BiocParallel::MulticoreParam(workers = parallel)
  } else {
    bioc_parallel <- BiocParallel::SerialParam()
  }
  # Check GSVA function arguments and adjust inputs
  gsva_api <- args(GSVA::gsva)
  if (names(as.list(gsva_api))[1] == "param") {
    gsva_args <- list(
      param = GSVA::gsvaParam(
        exprData = x, 
        geneSets = gene_set_list, 
        minSize = min_size, 
        maxSize = max_size, 
        maxDiff = TRUE, 
        kcdf = kcdf
      ), 
      verbose = verbose, 
      BPPARAM = bioc_parallel
    )
  } else {
    # Assume older version
    gsva_args <- list(
      expr = x, 
      gset.idx.list = gene_set_list, 
      min.sz = min_size, 
      max.sz = max_size, 
      mx.diff = TRUE, 
      kcdf = kcdf, 
      verbose = verbose, 
      BPPARAM = bioc_parallel
    )
  }
  out <- do.call(GSVA::gsva, args = gsva_args)
  return(out)
}


#' Cross-validation wrapped pathway enrichment
#'
#' This function wraps the pathway enrichment based data transformations for use in a 'cross-validated' clustering pipeline. 
#' Assumes a parallel backend has been registered for \code{\link[foreach]{foreach}} with \code{\link[doParallel]{registerDoParallel}}. 
#' Setting \code{parallel} does nothing. 
#'
#' @param dat_list list of data.tables corresponding gene expression matrices where the gene names have been replaced by dimension IDs (dim1, dim2, ...).
#' @param sub_index list of data.frames corresponding to cross-validation fold indicators as produced by \code{\link{subsampling}}
#' @param gene_id_list list of gene name vectors of the corresponding columns in dat_list
#' @param parallel ignored and set to 1 for spawned subprocesses
#' @param ... arguments passed on to \code{\link[COPS]{genes_to_pathways}}
#'
#' @return \code{list} of \code{data.frame}s with extracted pathway features
#' @export
subsample_pathway_enrichment <- function(
    dat_list, 
    sub_index, 
    gene_id_list, 
    parallel = 1, 
    ...
) {
  temp_list <- list()
  for (i in 1:length(sub_index)) {
    temp <- sub_index[[i]]
    datname <- names(sub_index)[i]
    if (is.null(datname)) datname <- i
    temp$datname <- datname
    temp <- merge(dat_list[[datname]], temp, by = "id")
    temp <- split_by_safe(temp, by = c("run", "fold"))
    temp <- lapply(temp, as.data.frame)
    temp <- lapply(
      temp, 
      function(x) list(dat = x, gene_ids = gene_id_list[[i]]))
    temp_list <- c(temp_list, temp)
  }
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(
    foreach(
      i = temp_list, 
      .combine = c,
      .export = c("genes_to_pathways"), #"dat_list"),
      .packages = c(
        "GSVA", 
        "fgsea", 
        "dnet", 
        "msigdbr", 
        "AnnotationDbi", 
        "org.Hs.eg.db")
      ) %dopar% {
        sel <- grep("^dim[0-9]+$", colnames(i$dat))
        temp <- i$dat[, sel]
        colnames(temp) <- i$gene_ids
        # Parallel = 1 for subprocesses
        pw_temp <- genes_to_pathways(t(temp), parallel = 1, ...)
        if(length(pw_temp) == 0) stop("Pathway enrichment failed, empty result.")
        pw_temp <- lapply(pw_temp, function(x) {
         temp2 <- as.data.frame(t(x))
         if (ncol(temp2)>0) {
           # replace pw names with dimX for compatibility
           colnames(temp2) <- paste0("dim", 1:ncol(temp2)) 
           return(cbind(i$dat[,-sel], temp2))
         } else {
           return(data.frame())
         }
        })
        # It is possible for the enrichment to return 0 enriched features.
        # Here we remove all cv-folds where this happens, which is
        # unsafe if reference fold is removed but other folds remain. 
        # TODO: investigate this issue further
        pw_temp <- pw_temp[sapply(pw_temp, ncol) > 0] 
        for (j in 1:length(pw_temp)) {
         pw_temp[[j]]$datname <- names(pw_temp)[j]
        }
        pw_temp
      }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(out)
}

#' @describeIn genes_to_pathways DiffRank by Wang et al. BMC Medical Genomics 2019
#' 
#' @param x gene feature matrix, samples on columns and genes on rows.
#' @param gene_set_list list of gene sets with gene names that correspond to rows in \strong{x}
#' @param parallel a numeric value indicating the number of processors to use when doing the calculations in parallel.
#' 
#' @return \code{data.frame} of pathway activities for each sample
#' @export
DiffRank <- function(
    x, 
    gene_set_list, 
    parallel = 1
) {
  ranks <- apply(x, 2, function(a) order(order(a)))
  out <- matrix(NA, length(gene_set_list), ncol(x))
  rownames(out) <- names(gene_set_list)
  colnames(out) <- colnames(x)
  
  parallel_clust <- setup_parallelization(parallel)
  
  out <- tryCatch(
    foreach(
      i = gene_set_list, 
      .combine = rbind,
      .export = c(),
      .multicombine = TRUE,
      .maxcombine = max(length(gene_set_list), 2)
    ) %dopar% {
      ind <- which(rownames(x) %in% i)
      # Compare mean ranks of pw genes vs non-pw genes
      pw_mean <- apply(ranks[ind,,drop=FALSE] - nrow(ranks)/2, 2, mean)
      npw_mean <- apply(ranks[-ind,,drop=FALSE] - nrow(ranks)/2, 2, mean)
      pw_mean - npw_mean
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  
  rownames(out) <- names(gene_set_list)
  out <- out[!apply(out, 1, function(a) all(is.na(a))),]
  
  return(out)
}

