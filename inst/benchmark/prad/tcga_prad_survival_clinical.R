## Load data
source("load_config.R")
source("prad/prad_default_parameters.R")
source("prad/tcga_prad_mrna_data.R") 

# Model selection
few_na <- which(sapply(prad_norm@colData, function(x) mean(is.na(x))) < 0.1)
few_na_not1 <- which(sapply(prad_norm@colData[few_na], function(x) length(table(x))) > 1)
table(prad_norm@colData$patient.histological_type, useNA = "always")
table(prad_norm@colData$patient.icd_10, useNA = "always")
table(prad_norm@colData$patient.icd_o_3_histology, useNA = "always")
table(prad_norm@colData$patient.icd_o_3_site, useNA = "always")
table(prad_norm@colData$patient.residual_tumor, useNA = "always")
table(prad_norm@colData$patient.prior_dx, useNA = "always")
table(prad_norm@colData$pathology_N_stage, useNA = "always")
table(prad_norm@colData$pathology_T_stage, useNA = "always")

few_na <- which(sapply(prad_subtype_all, function(x) mean(is.na(x) | x == "NA")) < 0.1)
few_na_not1 <- which(sapply(prad_subtype_all[few_na], function(x) length(table(x))) > 1)

which(sapply(prad_survival, function(x) mean(is.na(x) | x == "NA")) < 0.2)
which(sapply(prad_subtype_all, function(x) mean(is.na(x) | x == "NA")) < 0.2)


prad_surv_temp <- prad_survival
prad_surv_temp$PSA <- prad_surv_temp$patient.clinical_cqcf.psa_result_preop
prad_surv_temp$age <- prad_surv_temp$years_to_birth
prad_surv_temp$GC <- prad_surv_temp$Gleason_category
prad_surv_temp$T <- prad_surv_temp$pathology_T_stage
prad_surv_temp$N <- prad_surv_temp$pathology_N_stage

model1 <- survival::coxph(survival::Surv(time, event) ~ PSA, data = prad_surv_temp)
summary(model1) # *
model2 <- survival::coxph(survival::Surv(time, event) ~ age, data = prad_surv_temp)
summary(model2) # .
model3 <- survival::coxph(survival::Surv(time, event) ~ GC, data = prad_surv_temp)
summary(model3) #
model4 <- survival::coxph(survival::Surv(time, event) ~ T, data = prad_surv_temp)
summary(model4) #
model5 <- survival::coxph(survival::Surv(time, event) ~ N, data = prad_surv_temp)
summary(model5) # *

model12 <- survival::coxph(survival::Surv(time, event) ~ PSA + age, data = prad_surv_temp)
model15 <- survival::coxph(survival::Surv(time, event) ~ PSA + N, data = prad_surv_temp)
model25 <- survival::coxph(survival::Surv(time, event) ~ age + N, data = prad_surv_temp)

model125 <- survival::coxph(survival::Surv(time, event) ~ PSA + age + N, data = prad_surv_temp)

prad_survival_models <- list(model1, model2, model3, model4, model5, model12, model15, model25, model125)
prad_survival_aic <- sapply(prad_survival_models, stats::AIC)
prad_survival_bic <- sapply(prad_survival_models, stats::BIC)
prad_survival_pvalue_hi <- sapply(prad_survival_models, function(x) min(summary(x)$coef[,"Pr(>|z|)"]))

model0 <- survival::coxph(survival::Surv(time, event) ~ patient.clinical_cqcf.psa_result_preop, 
                          data = prad_surv_temp[!is.na(prad_surv_temp$pathology_N_stage),])
anova(model0, model15, test="LRT") # .

prad_model_variables <- lapply(prad_survival_models, function(x) colnames(attributes(x$terms)$factors))
prad_model_variables <- sapply(prad_model_variables, paste, collapse = " + ")

xtable::print.xtable(xtable::xtable(data.frame(Variables = prad_model_variables, 
                                               AIC = prad_survival_aic, 
                                               BIC = prad_survival_bic,
                                               pvalue = prad_survival_pvalue_hi),
                                    digits = c(NA, NA, 2, 2, 4)),
                     include.rownames = FALSE)



