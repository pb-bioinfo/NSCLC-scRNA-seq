library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(glmmTMB)


obj_macro_T <- qread("obj_macro_T.qs")
#####
#glmm SPP1, with batch

covariates =  c("histology_v4","sex_v2","smoking_v4","Grade_v4")
primary_var = "Race"
random_batch = TRUE
seurat_obj = obj_macro_T
outcome_target= "SPP1"

spp1_data <- FetchData(seurat_obj, vars = "SPP1", layer = "counts")
spp1_counts <- spp1_data$SPP1


model_df <- data.frame(
  Counts  = as.numeric(spp1_counts),
  Patient = as.factor(seurat_obj$Sample_ID), # Random intercept
  Depth   = seurat_obj$nCount_RNA,           # Sequencing depth offset
  stringsAsFactors = FALSE
)

model_df[[primary_var]] <- seurat_obj[[primary_var, drop = TRUE]]

cov_with_batch = c(covariates,"batch")
covariate_df = seurat_obj@meta.data[,cov_with_batch]

total_df = cbind(model_df,covariate_df)

cov_total   <- c(primary_var,covariates)

spp1_total_glmm_spp1_high


#  cov_total   <- c("Race", "histology_v3", "sex_v2", "smoking_v2", "Grade_v3")
#  primary_var <- "Race"
other_covs  <- setdiff(cov_total, primary_var)

# Model 1: Race only
# Models 2a-2d: Race + 1 covariate
# Model 3: All covariates
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in other_covs) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
  #    model_specs[[paste0("Model_2_Plus_", cv)]] <- c(cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- cov_total

all_terms_list <- list()


for (m_name in names(model_specs)) {
  current_covs <- model_specs[[m_name]]
  
  if (random_batch == FALSE){
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient) + offset(log(Depth))")
  }else{
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient)+ (1 | batch) + offset(log(Depth))")
  }
  glmm_formula <- as.formula(fml_str)
  
  
  # Run pseudobulk linear model with error trapping
  fit_result <- tryCatch(
    
    mixed_mod <- glmmTMB(glmm_formula, 
                         data = total_df, 
                         family = nbinom2),
    error = function(e) {
      warning("Failed fitting ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    
    summ = summary(fit_result)
    
    res_df <- as.data.frame(summ$coefficients$cond)
    #original is    Estimate Std. Error     z value     Pr(>|z|)
    
    cond_df <- res_df %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_target,
        Formula    = fml_str,
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, 
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    
    
    
    all_terms_list[[m_name]] <- cond_df
    
    
  }
}

df_all_terms    <- do.call(rbind, all_terms_list)

rownames(df_all_terms)    <- NULL

spp1_total_glmm_spp1_high <- df_all_terms

write.csv(spp1_total_glmm_spp1_high,"spp1_total_glmm_batchv2.csv")



#####
#glmm spp1 without batch
covariates =  c("histology_v4","sex_v2","smoking_v4","Grade_v4")
primary_var = "Race"
random_batch = FALSE
seurat_obj = obj_macro_T
outcome_target= "SPP1"

spp1_data <- FetchData(seurat_obj, vars = "SPP1", layer = "counts")
spp1_counts <- spp1_data$SPP1


model_df <- data.frame(
  Counts  = as.numeric(spp1_counts),
  Patient = as.factor(seurat_obj$Sample_ID), # Random intercept
  Depth   = seurat_obj$nCount_RNA,           # Sequencing depth offset
  stringsAsFactors = FALSE
)

model_df[[primary_var]] <- seurat_obj[[primary_var, drop = TRUE]]

cov_with_batch = c(covariates,"batch")
covariate_df = seurat_obj@meta.data[,cov_with_batch]

total_df = cbind(model_df,covariate_df)

cov_total   <- c(primary_var,covariates)

spp1_total_glmm_spp1_high


#  cov_total   <- c("Race", "histology_v3", "sex_v2", "smoking_v2", "Grade_v3")
#  primary_var <- "Race"
other_covs  <- setdiff(cov_total, primary_var)

# Model 1: Race only
# Models 2a-2d: Race + 1 covariate
# Model 3: All covariates
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in other_covs) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
  #    model_specs[[paste0("Model_2_Plus_", cv)]] <- c(cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- cov_total

all_terms_list <- list()


for (m_name in names(model_specs)) {
  current_covs <- model_specs[[m_name]]
  
  if (random_batch == FALSE){
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient) + offset(log(Depth))")
  }else{
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient)+ (1 | batch) + offset(log(Depth))")
  }
  glmm_formula <- as.formula(fml_str)
  
  
  # Run pseudobulk linear model with error trapping
  fit_result <- tryCatch(
    
    mixed_mod <- glmmTMB(glmm_formula, 
                         data = total_df, 
                         family = nbinom2),
    error = function(e) {
      warning("Failed fitting ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    
    summ = summary(fit_result)
    
    res_df <- as.data.frame(summ$coefficients$cond)
    #original is    Estimate Std. Error     z value     Pr(>|z|)
    
    cond_df <- res_df %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_target,
        Formula    = fml_str,
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, 
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    
    
    
    all_terms_list[[m_name]] <- cond_df
    
    
  }
}

df_all_terms    <- do.call(rbind, all_terms_list)

rownames(df_all_terms)    <- NULL

spp1_total_glmm_spp1_high <- df_all_terms

write.csv(spp1_total_glmm_spp1_high,"spp1_total_glmm_no_batchv4.csv")



#####
#glmm SPP1 LUAD only
covariates =  c("histology_v4","sex_v2","smoking_v4","Grade_v4")
primary_var = "Race"
random_batch = FALSE
seurat_obj = subset(obj_macro_T, histology_v4 == "LUAD")
outcome_target= "SPP1"

spp1_data <- FetchData(seurat_obj, vars = "SPP1", layer = "counts")
spp1_counts <- spp1_data$SPP1


model_df <- data.frame(
  Counts  = as.numeric(spp1_counts),
  Patient = as.factor(seurat_obj$Sample_ID), # Random intercept
  Depth   = seurat_obj$nCount_RNA,           # Sequencing depth offset
  stringsAsFactors = FALSE
)

model_df[[primary_var]] <- seurat_obj[[primary_var, drop = TRUE]]

cov_with_batch = c(covariates,"batch")
covariate_df = seurat_obj@meta.data[,cov_with_batch]

total_df = cbind(model_df,covariate_df)

cov_total   <- c(primary_var,covariates)

spp1_total_glmm_spp1_high


#  cov_total   <- c("Race", "histology_v3", "sex_v2", "smoking_v2", "Grade_v3")
#  primary_var <- "Race"
other_covs  <- setdiff(cov_total, primary_var)

# Model 1: Race only
# Models 2a-2d: Race + 1 covariate
# Model 3: All covariates
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in other_covs) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
  #    model_specs[[paste0("Model_2_Plus_", cv)]] <- c(cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- cov_total

all_terms_list <- list()


for (m_name in names(model_specs)) {
  current_covs <- model_specs[[m_name]]
  
  if (random_batch == FALSE){
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient) + offset(log(Depth))")
  }else{
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient)+ (1 | batch) + offset(log(Depth))")
  }
  glmm_formula <- as.formula(fml_str)
  
  
  # Run pseudobulk linear model with error trapping
  fit_result <- tryCatch(
    
    mixed_mod <- glmmTMB(glmm_formula, 
                         data = total_df, 
                         family = nbinom2),
    error = function(e) {
      warning("Failed fitting ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    
    summ = summary(fit_result)
    
    res_df <- as.data.frame(summ$coefficients$cond)
    #original is    Estimate Std. Error     z value     Pr(>|z|)
    
    cond_df <- res_df %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_target,
        Formula    = fml_str,
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, 
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    
    
    
    all_terms_list[[m_name]] <- cond_df
    
    
  }
}

df_all_terms    <- do.call(rbind, all_terms_list)

rownames(df_all_terms)    <- NULL

spp1_total_glmm_spp1_high <- df_all_terms

write.csv(spp1_total_glmm_spp1_high,"spp1_glmm_luad_v2.csv")



#####

#glmm SPP1, race array, using Pop_CEU as primary variable
obj_macro_array <- subset(obj_macro_T,Sample_ID %in% race_array_meta_short$Sample_ID)

lookup_vec_array <- setNames(race_array_meta_short$`Pop3(CEU)`, race_array_meta_short$Sample_ID)
obj_macro_array@meta.data$Pop_CEU <- lookup_vec_array[as.character(obj_macro_array@meta.data$Sample_ID)]


covariates =  c("histology_v4","sex_v2","smoking_v4","Grade_v4")
primary_var = "Pop_CEU"
random_batch = FALSE
seurat_obj =obj_macro_array
outcome_target= "SPP1"

spp1_data <- FetchData(seurat_obj, vars = "SPP1", layer = "counts")
spp1_counts <- spp1_data$SPP1


model_df <- data.frame(
  Counts  = as.numeric(spp1_counts),
  Patient = as.factor(seurat_obj$Sample_ID), # Random intercept
  Depth   = seurat_obj$nCount_RNA,           # Sequencing depth offset
  stringsAsFactors = FALSE
)

model_df[[primary_var]] <- seurat_obj[[primary_var, drop = TRUE]]

cov_with_batch = c(covariates,"batch")
covariate_df = seurat_obj@meta.data[,cov_with_batch]

total_df = cbind(model_df,covariate_df)

cov_total   <- c(primary_var,covariates)

spp1_total_glmm_spp1_high


#  cov_total   <- c("Race", "histology_v3", "sex_v2", "smoking_v2", "Grade_v3")
#  primary_var <- "Race"
other_covs  <- setdiff(cov_total, primary_var)

# Model 1: Race only
# Models 2a-2d: Race + 1 covariate
# Model 3: All covariates
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in other_covs) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
  #    model_specs[[paste0("Model_2_Plus_", cv)]] <- c(cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- cov_total

all_terms_list <- list()


for (m_name in names(model_specs)) {
  current_covs <- model_specs[[m_name]]
  
  if (random_batch == FALSE){
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient) + offset(log(Depth))")
  }else{
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient)+ (1 | batch) + offset(log(Depth))")
  }
  glmm_formula <- as.formula(fml_str)
  
  
  # Run pseudobulk linear model with error trapping
  fit_result <- tryCatch(
    
    mixed_mod <- glmmTMB(glmm_formula, 
                         data = total_df, 
                         family = nbinom2),
    error = function(e) {
      warning("Failed fitting ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    
    summ = summary(fit_result)
    
    res_df <- as.data.frame(summ$coefficients$cond)
    #original is    Estimate Std. Error     z value     Pr(>|z|)
    
    cond_df <- res_df %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_target,
        Formula    = fml_str,
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, 
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    
    
    
    all_terms_list[[m_name]] <- cond_df
    
    
  }
}

df_all_terms    <- do.call(rbind, all_terms_list)

rownames(df_all_terms)    <- NULL

spp1_total_glmm_spp1_high <- df_all_terms

write.csv(spp1_total_glmm_spp1_high,"spp1_glmm_array.csv")



#####
#glmm SPP1 using self-reported race, on the subset with race array data
covariates =  c("histology_v4","sex_v2","smoking_v4","Grade_v4")
primary_var = "Race"
random_batch = FALSE
seurat_obj =obj_macro_array
outcome_target= "SPP1"

spp1_data <- FetchData(seurat_obj, vars = "SPP1", layer = "counts")
spp1_counts <- spp1_data$SPP1


model_df <- data.frame(
  Counts  = as.numeric(spp1_counts),
  Patient = as.factor(seurat_obj$Sample_ID), # Random intercept
  Depth   = seurat_obj$nCount_RNA,           # Sequencing depth offset
  stringsAsFactors = FALSE
)

model_df[[primary_var]] <- seurat_obj[[primary_var, drop = TRUE]]

cov_with_batch = c(covariates,"batch")
covariate_df = seurat_obj@meta.data[,cov_with_batch]

total_df = cbind(model_df,covariate_df)

cov_total   <- c(primary_var,covariates)

spp1_total_glmm_spp1_high


#  cov_total   <- c("Race", "histology_v3", "sex_v2", "smoking_v2", "Grade_v3")
#  primary_var <- "Race"
other_covs  <- setdiff(cov_total, primary_var)

# Model 1: Race only
# Models 2a-2d: Race + 1 covariate
# Model 3: All covariates
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in other_covs) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
  #    model_specs[[paste0("Model_2_Plus_", cv)]] <- c(cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- cov_total

all_terms_list <- list()


for (m_name in names(model_specs)) {
  current_covs <- model_specs[[m_name]]
  
  if (random_batch == FALSE){
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient) + offset(log(Depth))")
  }else{
    fml_str    <- paste("Counts ~", paste(current_covs, collapse = " + "), "+ (1 | Patient)+ (1 | batch) + offset(log(Depth))")
  }
  glmm_formula <- as.formula(fml_str)
  
  
  # Run pseudobulk linear model with error trapping
  fit_result <- tryCatch(
    
    mixed_mod <- glmmTMB(glmm_formula, 
                         data = total_df, 
                         family = nbinom2),
    error = function(e) {
      warning("Failed fitting ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    
    summ = summary(fit_result)
    
    res_df <- as.data.frame(summ$coefficients$cond)
    #original is    Estimate Std. Error     z value     Pr(>|z|)
    
    cond_df <- res_df %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_target,
        Formula    = fml_str,
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, 
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    
    
    
    all_terms_list[[m_name]] <- cond_df
    
    
  }
}

df_all_terms    <- do.call(rbind, all_terms_list)

rownames(df_all_terms)    <- NULL

spp1_total_glmm_spp1_high <- df_all_terms

write.csv(spp1_total_glmm_spp1_high,"spp1_glmm_array_self_reported.csv")

