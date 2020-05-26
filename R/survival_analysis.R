


survival_preprocess <- function(event_data,
                                event_time_name = "days_to_death",
                                follow_up_time_name = "days_to_last_followup",
                                event_field_name = "vital_status",
                                event_name = "Dead", 
                                event_time_cutoff = Inf) {
  event_data[[event_time_name]] <- as.numeric(event_data[[event_time_name]])
  event_data[[follow_up_time_name]] <- as.numeric(event_data[[follow_up_time_name]])
  clinical_filter <- (event_data[[event_time_name]] <= event_time_cutoff | is.na(event_data[[event_time_name]])) & 
                     (event_data[[follow_up_time_name]] <= event_time_cutoff | is.na(event_data[[follow_up_time_name]]))
  event_data <- event_data[clinical_filter,]
  
  event_data$time <- event_data[[follow_up_time_name]]
  event_data$time[!is.na(event_data[[event_time_name]])] <- event_data[[event_time_name]][!is.na(event_data[[event_time_name]])]
  event_data$event <- event_data[[event_field_name]] == event_name
  
  return(event_data)
}

survival_evaluation <- function(event_data, 
                                clusters, 
                                survival_time_col = "time", 
                                survival_event_col = "event", 
                                survival_covariate_names = c("age", "stage"),
                                patient_id = "patient", 
                                ...) {
  clust_list <- split(clusters, by = c("run", "fold", "datname", "drname", "k", "m"))
  
  out <- foreach(clust = clust_list,
                 .combine = function(...) data.table::rbindlist(list(...)),
                 .export = c(),
                 .packages = c("survival"),
                 .multicombine = TRUE,
                 .maxcombine = length(clust_list)) %dopar% {
                   survival_ind <- match(event_data[[patient_id]], substr(clust$id, 1, 12))
                   temp <- event_data[!is.na(survival_ind),]
                   temp$cluster <- NA
                   temp$cluster <- clust$cluster[survival_ind[!is.na(survival_ind)]]
                   temp$cluster <- factor(temp$cluster)
                   
                   covariates_in_temp <- sapply(survival_covariate_names, function(x) length(table(temp[[x]])) > 1)
                   
                   model_formula <- paste0("survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
                                           paste(c(survival_covariate_names[covariates_in_temp], "cluster"), collapse = " + "))
                   model <- survival::coxph(as.formula(model_formula), data = temp)
                   #coef(summary(model))[,"Pr(>|z|)"]
                   model_formula0 <- paste0("survival::Surv(", survival_time_col, ", ",  survival_event_col, ") ~ ", 
                                           paste(survival_covariate_names[covariates_in_temp], collapse = " + "))
                   model0 <- survival::coxph(as.formula(model_formula0), data = temp)
                   res <- anova(model, model0, test="LRT")
                   res <- data.frame(clust[1,], cluster_significance = res[["P(>|Chi|)"]][2])
                   res$id <- NULL
                   res$cluster <- NULL
                   res
                 }
  return(out)
}
