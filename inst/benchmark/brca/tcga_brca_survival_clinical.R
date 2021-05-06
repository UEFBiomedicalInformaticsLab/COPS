## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
#source("load_config.R")
source("brca/brca_default_parameters.R")
source("brca/tcga_brca_mrna_data.R") 

few_na <- which(sapply(brca_norm@colData, function(x) mean(is.na(x) | x == "NA")) < 0.1)
few_na_not1 <- which(sapply(brca_norm@colData[few_na], function(x) length(table(x))) > 1)
which(sapply(dat_survival, function(x) mean(is.na(x) | x == "NA")) < 0.25)

dat_survival$pathology <- dat_survival$BRCA_Pathology

model1 <- survival::coxph(survival::Surv(time, event) ~ age, 
                          data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])
model2 <- survival::coxph(survival::Surv(time, event) ~ stage, 
                          data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])
model3 <- survival::coxph(survival::Surv(time, event) ~ pathology, 
                          data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])

model12 <- survival::coxph(survival::Surv(time, event) ~ age + stage, 
                           data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])
model13 <- survival::coxph(survival::Surv(time, event) ~ age + pathology, 
                           data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])
model23 <- survival::coxph(survival::Surv(time, event) ~ stage + pathology, 
                           data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])

model123 <- survival::coxph(survival::Surv(time, event) ~ age + stage + pathology, 
                            data = dat_survival[!dat_survival$BRCA_Subtype_PAM50 %in% c("NA", "Normal"),])

brca_survival_models <- list(model1, model2, model3, model12, model13, model23, model123)
brca_survival_aic <- sapply(brca_survival_models, stats::AIC)
brca_survival_bic <- sapply(brca_survival_models, stats::BIC)
brca_survival_pvalue_hi <- sapply(brca_survival_models, function(x) min(summary(x)$coef[,"Pr(>|z|)"]))

brca_model_variables <- lapply(brca_survival_models, function(x) colnames(attributes(x$terms)$factors))
brca_model_variables <- sapply(brca_model_variables, paste, collapse = " + ")

xtable::print.xtable(xtable::xtable(data.frame(Variables = brca_model_variables, 
                                               AIC = brca_survival_aic, 
                                               BIC = brca_survival_bic,
                                               pvalue = brca_survival_pvalue_hi),
                                    digits = c(NA, NA, 2, 2, 4)),
                     include.rownames = FALSE)