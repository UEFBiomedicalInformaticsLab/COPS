#' Pre-process event data for survival analysis
#' 
#' Use when survival data is not in time and event format e.g., TCGA clinical data.
#'
#' @param event_data A data.frame that contains survival times, event data and covariates.
#' @param event_time_name Name of the column that contains time to event. 
#' @param follow_up_time_name Name of the column that contains time to end of followup for a censored data point. 
#' @param event_field_name Name of the column that contains event indicators.
#' @param event_name Event column value that corresponds to the event of interest.
#' @param event_time_cutoff Upper cutoff value for time.
#' @param event_time_lower_cutoff Lower cutoff value for time.
#'
#' @return \code{data.frame} with "time" and "event" columns added to input
#' @export
survival_preprocess <- function(
    event_data,
    event_time_name = "days_to_death",
    follow_up_time_name = "days_to_last_followup",
    event_field_name = "vital_status",
    event_name = "Dead", 
    event_time_cutoff = Inf,
    event_time_lower_cutoff = 0
) {
  event_data[[event_time_name]] <- suppressMessages(
    as.numeric(event_data[[event_time_name]]))
  event_data[[follow_up_time_name]] <- suppressMessages(
    as.numeric(event_data[[follow_up_time_name]]))
  etime_na_ind <- is.na(event_data[[event_time_name]])
  ftime_na_ind <- is.na(event_data[[follow_up_time_name]])
  clinical_filter <- 
    (event_data[[event_time_name]] <= event_time_cutoff |  etime_na_ind) & 
    (event_data[[event_time_name]] > event_time_lower_cutoff | etime_na_ind) & 
    (event_data[[follow_up_time_name]] <= event_time_cutoff | ftime_na_ind) & 
    (event_data[[follow_up_time_name]] > event_time_lower_cutoff | ftime_na_ind)
  event_data <- event_data[clinical_filter,]
  
  event_data$time <- event_data[[follow_up_time_name]]
  event_data$time[!etime_na_ind[clinical_filter]] <- 
    event_data[[event_time_name]][!etime_na_ind[clinical_filter]]
  event_data$event <- event_data[[event_field_name]] == event_name
  
  return(event_data)
}

#' Likelihood ratio analysis of Cox PH survival models for clustering results
#'
#' @param event_data A data.frame that contains survival times, event data and covariates.
#' @param clusters A data.frame or data.table containing clustering information.
#' @param survival_time_col Name of the column in \code{event_data} that contains survival time. 
#' @param survival_event_col Name of the column in \code{event_data} that contains event indicators.
#' @param survival_covariate_names Names of covariate columns in \code{event_data}. 
#' @param row_id Name of column in \code{event_data} that matches sample IDs in \code{clusters}.
#' @param by Vector of column names that identify a single clustering result in \code{clusters}.
#' @param parallel number of parallel threads
#' @param ... Extra arguments are ignored.
#'
#' @return \code{data.table} of Cox PH LRT-test p-values and Harrell's concordance index values
#' @export
#' 
#' @importFrom survival Surv coxph
subsample_survival_evaluation <- function(
    event_data, 
    clusters, 
    survival_time_col = "time", 
    survival_event_col = "event", 
    survival_covariate_names = NULL,
    row_id = "ID", 
    by = c("run", "fold", "datname", "drname", "k", "m"), 
    parallel = 1, 
    ...
) {
  by <- by[by %in% colnames(clusters)]
  clust_list <- split_by_safe(clusters, by)
  parallel_clust <- setup_parallelization(parallel)
  out <- tryCatch(
    foreach(
      clust = clust_list,
      .combine = function(...) data.table::rbindlist(list(...)),
      .export = c("by"),
      .packages = c("survival"),
      .multicombine = TRUE,
      .maxcombine = max(length(clust_list), 2)
    ) %dopar% {
      survival_ind <- match(event_data[[row_id]], clust$id)
      temp <- event_data[!is.na(survival_ind),]
      temp$cluster <- NA
      temp$cluster <- clust$cluster[survival_ind[!is.na(survival_ind)]]
      temp$cluster <- factor(temp$cluster)
      if ("data.table" %in% class(clust)) {
       out_i <- data.frame(
         clust[1,..by], 
         cluster_significance = NA,
         concordance_index = NA)
      } else {
       out_i <- data.frame(
         clust[1,by], 
         cluster_significance = NA,
         concordance_index = NA)
      }
      cluster_sizes <- table(temp$cluster)
      if (length(cluster_sizes) > 1 & all(cluster_sizes > 1)) {
        covariates_in_temp <- sapply(
          survival_covariate_names, 
          function(x) length(table(temp[[x]])) > 1)
        
        model_formula <- paste0(
          "survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
          paste(
            c(survival_covariate_names[covariates_in_temp], "cluster"), 
            collapse = " + "
          )
        )
        out_i_vals <- tryCatch(
          {
            model <- survival::coxph(as.formula(model_formula), data = temp)
            model_formula0 <- paste0(
              "survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
              paste(
                survival_covariate_names[covariates_in_temp], 
                collapse = " + "
              )
            )
            if (all(!covariates_in_temp)) {
              model_formula0 <- paste(model_formula0, "+ 1")
            }
            p_val <- tryCatch(
              {
                model0 <- survival::coxph(as.formula(model_formula0), data = temp)
                res <- anova(model, model0, test="LRT")
                p_col <- which(colnames(res) %in% c("P(>|Chi|)", "Pr(>|Chi|)"))
                res[2, p_col]
              }, 
              error = function(e) {warning(e);return(NA)}
            )
            ci <- tryCatch(
              survival::concordance(model)[["concordance"]], 
              error = function(e) {warning(e);return(NA)})
            list(p_val = p_val, ci = ci)
          }, 
          error = function(e) {warning(e);return(list(p_val = NA, ci = NA))}
        )
        out_i$cluster_significance <- out_i_vals[["p_val"]]
        out_i$concordance_index <- out_i_vals[["ci"]]
      }
      out_i
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(out)
}
