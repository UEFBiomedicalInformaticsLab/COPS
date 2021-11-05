#' Pre-process event data for survival analysis
#'
#' @param event_data A data.frame that contains survival times, event data and covariates.
#' @param event_time_name Name of the column that contains time to event. 
#' @param follow_up_time_name Name of the column that contains time to end of followup for a censored data point. 
#' @param event_field_name Name of the column that contains event indicators.
#' @param event_name Event column value that corresponds to the event of interest.
#' @param event_time_cutoff Upper cutoff value for time.
#' @param event_time_lower_cutoff Lower cutoff value for time.
#'
#' @return
#' @export
survival_preprocess <- function(event_data,
                                event_time_name = "days_to_death",
                                follow_up_time_name = "days_to_last_followup",
                                event_field_name = "vital_status",
                                event_name = "Dead", 
                                event_time_cutoff = Inf,
                                event_time_lower_cutoff = 0) {
  event_data[[event_time_name]] <- as.numeric(event_data[[event_time_name]])
  event_data[[follow_up_time_name]] <- as.numeric(event_data[[follow_up_time_name]])
  clinical_filter <- (event_data[[event_time_name]] <= event_time_cutoff | is.na(event_data[[event_time_name]])) & 
                     (event_data[[event_time_name]] > event_time_lower_cutoff | is.na(event_data[[event_time_name]])) & 
                     (event_data[[follow_up_time_name]] <= event_time_cutoff | is.na(event_data[[follow_up_time_name]])) & 
                     (event_data[[follow_up_time_name]] > event_time_lower_cutoff | is.na(event_data[[follow_up_time_name]]))
  event_data <- event_data[clinical_filter,]
  
  event_data$time <- event_data[[follow_up_time_name]]
  event_data$time[!is.na(event_data[[event_time_name]])] <- event_data[[event_time_name]][!is.na(event_data[[event_time_name]])]
  event_data$event <- event_data[[event_field_name]] == event_name
  
  return(event_data)
}

#' Likelihood ratio analysis of Cox PH survival models for clustering results
#'
#' @param event_data A data.frame that contains survival times, event data and covariates.
#' @param clusters A data.frame or data.table containing clustering information.
#' @param survival_time_col Name of the column in \code{event_data} that contains survival time. 
#' @param survival_event_col Name of the column in \code{event_data} that contains event indicators.
#' @param survival_covariate_names Nme of covariate columns. 
#' @param row_id Name of column in \code{event_data} that matches sample IDs in \code{clusters}.
#' @param by Vector of column names that identify a single clustering result in \code{clusters}.
#' @param ... Extra arguments are ignored.
#'
#' @return
#' @export
#' 
#' @importFrom survival Surv coxph
survival_evaluation <- function(event_data, 
                                clusters, 
                                survival_time_col = "time", 
                                survival_event_col = "event", 
                                survival_covariate_names = c("age", "stage"),
                                row_id = "ID", 
                                by = c("run", "fold", "datname", "drname", "k", "m"), 
                                parallel = 1, 
                                ...) {
  parallel_clust <- setup_parallelization(parallel)
  
  if (data.table::is.data.table(clusters)) {
    clust_list <- split(clusters, by = by)
  } else {
    clust_list <- split(clusters, clusters[, by])
  }
  
  out <- tryCatch(foreach(clust = clust_list,
                 .combine = function(...) data.table::rbindlist(list(...)),
                 .export = c("by"),
                 .packages = c("survival"),
                 .multicombine = TRUE,
                 .maxcombine = length(clust_list)) %dopar% {
                   survival_ind <- match(event_data[[row_id]], clust$id)
                   temp <- event_data[!is.na(survival_ind),]
                   temp$cluster <- NA
                   temp$cluster <- clust$cluster[survival_ind[!is.na(survival_ind)]]
                   temp$cluster <- factor(temp$cluster)
                   
                   out_i <- data.frame(clust[1,..by], cluster_significance = NA)
                   out_i$id <- NULL
                   out_i$cluster <- NULL
                   if (length(table(temp$cluster)) > 1) {
                     covariates_in_temp <- sapply(survival_covariate_names, function(x) length(table(temp[[x]])) > 1)
                     
                     model_formula <- paste0("survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
                                             paste(c(survival_covariate_names[covariates_in_temp], "cluster"), collapse = " + "))
                     model <- survival::coxph(as.formula(model_formula), data = temp)
                     #coef(summary(model))[,"Pr(>|z|)"]
                     model_formula0 <- paste0("survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
                                              paste(survival_covariate_names[covariates_in_temp], collapse = " + "))
                     model0 <- survival::coxph(as.formula(model_formula0), data = temp)
                     res <- anova(model, model0, test="LRT")
                     out_i$cluster_significance <- res[["P(>|Chi|)"]][2]
                   }
                   
                   out_i
                 }, finally = if(parallel > 1) parallel::stopCluster(parallel_clust))
  return(out)
}
