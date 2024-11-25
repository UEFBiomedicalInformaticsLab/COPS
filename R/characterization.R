#' Heatmap annotation automation
#' 
#' Generates stacked \code{\link[ComplexHeatmap]{HeatmapAnnotation}}s from a 
#' list of annotations. 
#' 
#' Characters and factors are assigned colors from categorical palettes 
#' according to the number of categories. NOTE: if a palette does not have 
#' enough colors, the colors are assigned from the next palette in the list. 
#'
#' @param annotations list of variables
#' @param factor_color_sets list of color sets (not in order)
#' @param ... passed onto \code{\link[ComplexHeatmap]{HeatmapAnnotation}}
#'
#' @return \code{\link[ComplexHeatmap]{`HeatmapList-class`}}
#' @export 
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation %v%
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pals watlington kelly
heatmap_annotations <- function(
    annotations, 
    factor_color_sets = list(
      RColorBrewer::brewer.pal(8, "Dark2"), 
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(9, "Pastel1"),
      pals::watlington(n = 16), 
      pals::kelly(n = 22)
    ), 
    ...
) {
  a_types <- sapply(annotations, class)
  for (i in names(annotations)) {
    if (a_types[i] %in% c("character", "integer")) {
      annotations[[i]] <- factor(annotations[[i]])
    }
  }
  a_types <- sapply(annotations, class)
  
  a_sizes <- sapply(annotations, function(x) length(table(x)))
  
  factor_palette_size <- sapply(factor_color_sets, length)
  factor_palette_size_reserved <- rep(0, length(factor_palette_size))
  
  a_palette_ind <- c()
  a_palette_offset <- c()
  
  heat_annotations <- list()
  heat_annotation_colors <- list()
  
  cont_colors <- c("blue", "green", "orange", "red")
  cont_taken <- rep(FALSE, length(cont_colors))
  
  for (i in names(annotations)) {
    if (a_types[i] == "factor") {
      fpi <- which(factor_palette_size - factor_palette_size_reserved >= a_sizes[i])
      if (length(fpi) == 0) stop(paste0("No large enough palettes left for ", i, "!"))
      fpi <- fpi[1]
      a_palette_ind[i] <- fpi
      a_palette_offset[i] <- factor_palette_size_reserved[fpi]
      factor_palette_size_reserved[fpi] <- factor_palette_size_reserved[fpi] + a_sizes[i]
      
      col_inds <- (1+a_palette_offset[i]):(a_sizes[i]+a_palette_offset[i])
      heat_annotation_colors[[i]] <- factor_color_sets[[fpi]][col_inds]
      
      names(heat_annotation_colors[[i]]) <- names(table(annotations[[i]]))
    } else {
      cpi <- which(!cont_taken)
      if (length(cpi) == 0) stop("Too few continuous scales defined!")
      cpi <- cpi[1]
      cont_taken[cpi] <- TRUE
      heat_annotation_colors[[i]] <- circlize::colorRamp2(
        c(
          min(annotations[[i]], na.rm = TRUE), 
          max(annotations[[i]], na.rm = TRUE)
        ), 
        c("white", cont_colors[cpi]))
    }
    
    heat_anno_args <- list(
      annotations[[i]], 
      list(heat_annotation_colors[[i]])
    )
    names(heat_anno_args) <- c(i, "col")
    names(heat_anno_args[["col"]]) <- i
    
    heat_annotations[[i]] <- do.call(
      ComplexHeatmap::HeatmapAnnotation, 
      args = c(heat_anno_args, list(...))
    )
  }
  
  heat_annotations_combined <- Reduce(
    ComplexHeatmap::`%v%`, 
    heat_annotations)
  return(heat_annotations_combined)
}

#' Annotated heatmap
#' 
#' Generate a heatmap with annotations.
#'
#' @param dat feature matrix or list of matrices, samples on columns
#' @param variable_list list of column annotation variables
#' @param feature_names optional character vector matching to a subset of \code{dat} rownames
#' @param center whether to apply row-wise center, passed to \link[base]{scale} 
#' @param scale whether to apply row-wise scaling, passed to \link[base]{scale} 
#' @param show_column_names whether to include colnames (sample names)
#' @param show_row_names whether to include rownames (gene names)
#' @param show_column_dend whether to include column dendrogram
#' @param show_row_dend whether to include row dendrogram
#' @param row_names_side where row names should be located
#' @param legend_names names to use for \code{dat} color legends
#' @param color_breaks manual color breaks in form \code{c(min, middle, max)}
#' @param heatmap_color_function passed to \code{\link[ComplexHeatmap]{Heatmap}} as \code{col}
#' @param ... other arguments passed to \code{\link[ComplexHeatmap]{Heatmap}} 
#'
#' @return \code{\link[ComplexHeatmap]{`HeatmapList-class`}}
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap %v%
#' @importFrom circlize colorRamp2
heatmap_annotated <- function(
    dat, 
    variable_list = list(), 
    feature_names = NULL,
    center = TRUE, scale = FALSE, 
    show_column_names = FALSE, 
    show_row_names = TRUE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    row_names_side = "left", 
    legend_names = NULL, 
    color_breaks = NULL,
    heatmap_color_function = circlize::colorRamp2(color_breaks, c("blue", "white", "red")), 
    ...
) {
  if (!"list" %in% class(dat)) dat <- list(dat)
  if (length(variable_list) > 0) {
    col_annots <- heatmap_annotations(variable_list)
    heatmap_list <- list(col_annots)
  } else {
    heatmap_list <- list()
  }
  if (is.null(legend_names)) legend_names <- names(dat)
  if (is.null(legend_names)) legend_names <- paste0("input", 1:length(dat))
  for (i in 1:length(dat)) {
    if (!is.null(feature_names)) {
      dat[[i]] <- dat[[i]][intersect(feature_names, rownames(dat[[i]])),]
    }
    if (center | scale) {
      dat[[i]] <- t(scale(t(dat[[i]]), center = center, scale = scale))
    }
    if (is.null(color_breaks)) {
      color_breaks <- c(min(dat[[i]]), mean(dat[[i]]), max(dat[[i]]))
    }
    heatmap_list <- c(
      heatmap_list, 
      list(ComplexHeatmap::Heatmap(
        dat[[i]], 
        name = legend_names[i],
        show_column_names = show_column_names, 
        show_row_names = show_row_names,
        show_column_dend = show_column_dend,
        show_row_dend = show_row_dend,
        row_names_side = row_names_side,
        col = heatmap_color_function,
        ...
      ))
    )
  }
  heatmap_combined <- Reduce(ComplexHeatmap::`%v%`, heatmap_list)
  return(heatmap_combined)
}

#' Univariate group difference tests for feature selection
#' 
#' Computes various statistics for molecular differences between groups. 
#'
#' @param dat molecular data, samples on columns, rows must be named
#' @param group categorical variable matching \code{dat} columns
#' @param remove_zero_var whether to remove zero variance genes
#' @param parallel number of threads
#'
#' @return \code{list} of test results
#' @export
univariate_features <- function(
    dat, 
    group, 
    remove_zero_var = TRUE, 
    parallel = 1
) {
  dat <- t(dat)
  if (remove_zero_var) {
    nzv <- which(apply(dat, 2, var) > 0)
    dat <- dat[,nzv]
  }
  gene_names <- colnames(dat)
  if (is.null(gene_names)) stop("No feature names found.")
  # Combine threaded output
  cfun <- function(a,b) {
    out <- list()
    for (i in c("kw_test_p", "anova_p", "kw_test_p_adj", "anova_p_adj")) {
      out[[i]] <- c(a[[i]],b[[i]])
    }
    for (i in c("t_test_t", "w_test_w", "t_test_p", "w_test_p",
                "correlation", "median", "mean")) {
      out[[i]] <- rbind(a[[i]],b[[i]])
    }
    return(out)
  }
  parallel_clust <- setup_parallelization(parallel)
  out <- tryCatch(foreach(
    j = lapply(gene_names, function(x) dat[,x]), 
    .combine = cfun) %dopar% 
  {
    out_j <- list()
    out_j$anova_p <- summary(aov(j ~ group))[[1]][["Pr(>F)"]][1]
    out_j$kw_test_p <- kruskal.test(j ~ group)$p.value
    out_j$t_test_t <- rep(NA, length(unique(group)))
    names(out_j$t_test_t) <- unique(group)
    cor_j <- out_j$w_test_p <- out_j$w_test_w <- out_j$t_test_p <- out_j$t_test_t
    for (k in unique(group)) {
      k_indicator <- group == k
      if (sum(k_indicator) > 1) {
        t_res <- t.test(
          j[k_indicator], 
          j[!k_indicator], 
          alternative = "two.sided", 
          var.equal = FALSE)
        out_j$t_test_t[k] <- t_res$statistic
        out_j$t_test_p[k] <- t_res$p.value
        w_res <- wilcox.test(
          j[k_indicator], j[!k_indicator], 
          alternative = "two.sided")
        out_j$w_test_w[k] <- w_res$statistic
        out_j$w_test_p[k] <- w_res$p.value
      }
      cor_j[k] <- cor(k_indicator, j, method = "spearman")
    }
    out_j$correlation <- data.frame(t(cor_j))
    out_j$median <- data.frame(t(tapply(j, group, median, na.rm = TRUE)))
    out_j$mean <- data.frame(t(tapply(j, group, mean, na.rm = TRUE)))
    out_j
  }, finally = close_parallel_cluster(parallel_clust))
  # Vectorized computations for signal-to-noise ratio
  out$snr <- matrix(
    NA, 
    nrow = length(gene_names), 
    ncol = length(unique(group)))
  out$effect <- matrix(
    NA, 
    nrow = length(gene_names), 
    ncol = length(unique(group)))
  colnames(out$snr) <- colnames(out$effect) <- unique(group)
  rownames(out$snr) <- rownames(out$effect) <- gene_names
  
  for (k in unique(group)) {
    mu_k <- apply(dat[group == k, gene_names, drop = FALSE], 2, mean)
    mu_nk <- apply(dat[group != k, gene_names, drop = FALSE], 2, mean)
    sigma_k <- apply(dat[group == k, gene_names, drop = FALSE], 2, sd)
    sigma_nk <- apply(dat[group != k, gene_names, drop = FALSE], 2, sd)
    out$snr[,k] <- (mu_k - mu_nk) / (sigma_k + sigma_nk)
    out$effect[,k] <- (mu_k - mu_nk)
  }
  # Adjust p-values
  out$kw_test_p_adj <- p.adjust(out$kw_test_p, method = "BH")
  out$anova_p_adj <- p.adjust(out$anova_p, method = "BH")
  
  out$t_p_adj <- apply(out$t_test_p, 2, p.adjust, method = "BH")
  out$w_p_adj <- apply(out$w_test_p, 2, p.adjust, method = "BH")
  # Add gene names
  for (i in c("kw_test_p", "anova_p", "kw_test_p_adj", "anova_p_adj")) {
    names(out[[i]]) <- gene_names
  }
  for (i in c("t_test_t", "w_test_w", "t_test_p", "w_test_p",
              "correlation", "median", "mean", "t_p_adj", "w_p_adj")) {
    rownames(out[[i]]) <- gene_names
  }
  
  return(out)
}
