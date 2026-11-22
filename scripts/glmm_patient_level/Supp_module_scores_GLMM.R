library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(glmmTMB)

sce_score <- qread("sce_macro.qs")

#####
#M2 score glmm

outcome_col     = "M2_score1"
seurat_obj = sce_score

primary_var     = "Race"
covariates      = c("histology_v4", "sex_v2", "smoking_v2", "Grade_v4")
sample_col      = "Sample_ID"
depth_col       = "nCount_RNA"
use_depth_offset = FALSE
family          = gaussian()
random_batch = FALSE

meta <- seurat_obj@meta.data



# 1. Verify required columns exist
req_cols <- c(outcome_col, primary_var, sample_col, covariates)
if (use_depth_offset) req_cols <- c(req_cols, depth_col)

missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 2. Build base modeling dataframe
base_df <- data.frame(
  Outcome = as.numeric(meta[[outcome_col]]),
  Primary = meta[[primary_var]],
  Patient = as.factor(meta[[sample_col]]),
  stringsAsFactors = FALSE
)
colnames(base_df)[1:2] <- c(outcome_col, primary_var)

# Append covariates
for (cv in covariates) {
  base_df[[cv]] <- meta[[cv]]
}

if (use_depth_offset) {
  base_df$Depth_Offset <- log(as.numeric(meta[[depth_col]]))
}

if (random_batch){
  base_df$batch <- meta$batch
}

# 3. Define the hierarchical model specifications
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in covariates) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- c(primary_var, covariates)

# 4. Iterate across models
all_terms_list <- list()

for (m_name in names(model_specs)) {
  current_vars <- model_specs[[m_name]]
  
  # Subset to complete cases for the active formula
  active_cols <- c(outcome_col, "Patient", current_vars)
  if (use_depth_offset) active_cols <- c(active_cols, "Depth_Offset")
  
  fit_df <- base_df[complete.cases(base_df[, active_cols, drop = FALSE]), , drop = FALSE]
  
  # Construct formula
  fixed_effects <- paste(current_vars, collapse = " + ")
  offset_term   <- if (use_depth_offset) " + offset(Depth_Offset)" else ""
  
  
  if (random_batch){
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient) + (1 | batch)")
  }else{
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient)")
  }
  
  
  model_formula <- as.formula(formula_str)
  
  message("Fitting ", m_name, " (", nrow(fit_df), " cells across ", 
          length(unique(fit_df$Patient)), " patients)...")
  print(colnames(fit_df))
  # Fit model with error trapping
  fit_result <- tryCatch(
    glmmTMB(
      formula = model_formula,
      data    = fit_df,
      family  = family
    ),
    error = function(e) {
      warning("Fitting failed for ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    sum_fit <- summary(fit_result)
    cond_coefs <- as.data.frame(sum_fit$coefficients$cond)
    
    # Extract diagnostics and tidy dataframe
    coef_df <- cond_coefs %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_col,
        Formula    = formula_str,
        N_Cells    = nrow(fit_df),
        N_Patients = length(unique(fit_df$Patient)),
        AIC        = AIC(fit_result),
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, N_Cells, N_Patients, AIC,
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    all_terms_list[[m_name]] <- coef_df
  }
}

final_df <- do.call(rbind, all_terms_list)
rownames(final_df) <- NULL

write.csv(final_df,"glmm_m2_resultsv3.csv")


#####
#M1 score glmm


outcome_col     = "M1_score1"
seurat_obj = sce_score

primary_var     = "Race"
covariates      = c("histology_v4", "sex_v2", "smoking_v2", "Grade_v4")
sample_col      = "Sample_ID"
depth_col       = "nCount_RNA"
use_depth_offset = FALSE
family          = gaussian()
random_batch = FALSE

meta <- seurat_obj@meta.data



# 1. Verify required columns exist
req_cols <- c(outcome_col, primary_var, sample_col, covariates)
if (use_depth_offset) req_cols <- c(req_cols, depth_col)

missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 2. Build base modeling dataframe
base_df <- data.frame(
  Outcome = as.numeric(meta[[outcome_col]]),
  Primary = meta[[primary_var]],
  Patient = as.factor(meta[[sample_col]]),
  stringsAsFactors = FALSE
)
colnames(base_df)[1:2] <- c(outcome_col, primary_var)

# Append covariates
for (cv in covariates) {
  base_df[[cv]] <- meta[[cv]]
}

if (use_depth_offset) {
  base_df$Depth_Offset <- log(as.numeric(meta[[depth_col]]))
}

if (random_batch){
  base_df$batch <- meta$batch
}

# 3. Define the hierarchical model specifications
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in covariates) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- c(primary_var, covariates)

# 4. Iterate across models
all_terms_list <- list()

for (m_name in names(model_specs)) {
  current_vars <- model_specs[[m_name]]
  
  # Subset to complete cases for the active formula
  active_cols <- c(outcome_col, "Patient", current_vars)
  if (use_depth_offset) active_cols <- c(active_cols, "Depth_Offset")
  
  fit_df <- base_df[complete.cases(base_df[, active_cols, drop = FALSE]), , drop = FALSE]
  
  # Construct formula
  fixed_effects <- paste(current_vars, collapse = " + ")
  offset_term   <- if (use_depth_offset) " + offset(Depth_Offset)" else ""
  
  
  if (random_batch){
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient) + (1 | batch)")
  }else{
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient)")
  }
  
  
  model_formula <- as.formula(formula_str)
  
  message("Fitting ", m_name, " (", nrow(fit_df), " cells across ", 
          length(unique(fit_df$Patient)), " patients)...")
  print(colnames(fit_df))
  # Fit model with error trapping
  fit_result <- tryCatch(
    glmmTMB(
      formula = model_formula,
      data    = fit_df,
      family  = family
    ),
    error = function(e) {
      warning("Fitting failed for ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    sum_fit <- summary(fit_result)
    cond_coefs <- as.data.frame(sum_fit$coefficients$cond)
    
    # Extract diagnostics and tidy dataframe
    coef_df <- cond_coefs %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_col,
        Formula    = formula_str,
        N_Cells    = nrow(fit_df),
        N_Patients = length(unique(fit_df$Patient)),
        AIC        = AIC(fit_result),
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, N_Cells, N_Patients, AIC,
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    all_terms_list[[m_name]] <- coef_df
  }
}

final_df <- do.call(rbind, all_terms_list)
rownames(final_df) <- NULL

write.csv(final_df,"glmm_m1_resultsv3.csv")

#####
#glmm on cytotoxic score
sce_cd8 <- qread("sce_cd8.qs")


qsave(sce_cd8,"sce_cd8.qs")

outcome_col     = "cytotoxic_score1"
seurat_obj = sce_cd8

primary_var     = "Race"
covariates      = c("histology_v4", "sex_v2", "smoking_v2", "Grade_v4")
sample_col      = "Sample_ID"
depth_col       = "nCount_RNA"
use_depth_offset = FALSE
family          = gaussian()
random_batch = FALSE

meta <- seurat_obj@meta.data



# 1. Verify required columns exist
req_cols <- c(outcome_col, primary_var, sample_col, covariates)
if (use_depth_offset) req_cols <- c(req_cols, depth_col)

missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 2. Build base modeling dataframe
base_df <- data.frame(
  Outcome = as.numeric(meta[[outcome_col]]),
  Primary = meta[[primary_var]],
  Patient = as.factor(meta[[sample_col]]),
  stringsAsFactors = FALSE
)
colnames(base_df)[1:2] <- c(outcome_col, primary_var)

# Append covariates
for (cv in covariates) {
  base_df[[cv]] <- meta[[cv]]
}

if (use_depth_offset) {
  base_df$Depth_Offset <- log(as.numeric(meta[[depth_col]]))
}

if (random_batch){
  base_df$batch <- meta$batch
}

# 3. Define the hierarchical model specifications
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in covariates) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- c(primary_var, covariates)

# 4. Iterate across models
all_terms_list <- list()

for (m_name in names(model_specs)) {
  current_vars <- model_specs[[m_name]]
  
  # Subset to complete cases for the active formula
  active_cols <- c(outcome_col, "Patient", current_vars)
  if (use_depth_offset) active_cols <- c(active_cols, "Depth_Offset")
  
  fit_df <- base_df[complete.cases(base_df[, active_cols, drop = FALSE]), , drop = FALSE]
  
  # Construct formula
  fixed_effects <- paste(current_vars, collapse = " + ")
  offset_term   <- if (use_depth_offset) " + offset(Depth_Offset)" else ""
  
  
  if (random_batch){
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient) + (1 | batch)")
  }else{
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient)")
  }
  
  
  model_formula <- as.formula(formula_str)
  
  message("Fitting ", m_name, " (", nrow(fit_df), " cells across ", 
          length(unique(fit_df$Patient)), " patients)...")
  print(colnames(fit_df))
  # Fit model with error trapping
  fit_result <- tryCatch(
    glmmTMB(
      formula = model_formula,
      data    = fit_df,
      family  = family
    ),
    error = function(e) {
      warning("Fitting failed for ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    sum_fit <- summary(fit_result)
    cond_coefs <- as.data.frame(sum_fit$coefficients$cond)
    
    # Extract diagnostics and tidy dataframe
    coef_df <- cond_coefs %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_col,
        Formula    = formula_str,
        N_Cells    = nrow(fit_df),
        N_Patients = length(unique(fit_df$Patient)),
        AIC        = AIC(fit_result),
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, N_Cells, N_Patients, AIC,
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    all_terms_list[[m_name]] <- coef_df
  }
}

final_df <- do.call(rbind, all_terms_list)
rownames(final_df) <- NULL

write.csv(final_df,file = "glmm_cyto_resultsv3.csv")



#####

#exhausted score glmm

outcome_col     = "exhausted_score1"
seurat_obj = sce_cd8

primary_var     = "Race"
covariates      = c("histology_v4", "sex_v2", "smoking_v2", "Grade_v4")
sample_col      = "Sample_ID"
depth_col       = "nCount_RNA"
use_depth_offset = FALSE
family          = gaussian()
random_batch = FALSE

meta <- seurat_obj@meta.data



# 1. Verify required columns exist
req_cols <- c(outcome_col, primary_var, sample_col, covariates)
if (use_depth_offset) req_cols <- c(req_cols, depth_col)

missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 2. Build base modeling dataframe
base_df <- data.frame(
  Outcome = as.numeric(meta[[outcome_col]]),
  Primary = meta[[primary_var]],
  Patient = as.factor(meta[[sample_col]]),
  stringsAsFactors = FALSE
)
colnames(base_df)[1:2] <- c(outcome_col, primary_var)

# Append covariates
for (cv in covariates) {
  base_df[[cv]] <- meta[[cv]]
}

if (use_depth_offset) {
  base_df$Depth_Offset <- log(as.numeric(meta[[depth_col]]))
}

if (random_batch){
  base_df$batch <- meta$batch
}

# 3. Define the hierarchical model specifications
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in covariates) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- c(primary_var, covariates)

# 4. Iterate across models
all_terms_list <- list()

for (m_name in names(model_specs)) {
  current_vars <- model_specs[[m_name]]
  
  # Subset to complete cases for the active formula
  active_cols <- c(outcome_col, "Patient", current_vars)
  if (use_depth_offset) active_cols <- c(active_cols, "Depth_Offset")
  
  fit_df <- base_df[complete.cases(base_df[, active_cols, drop = FALSE]), , drop = FALSE]
  
  # Construct formula
  fixed_effects <- paste(current_vars, collapse = " + ")
  offset_term   <- if (use_depth_offset) " + offset(Depth_Offset)" else ""
  
  
  if (random_batch){
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient) + (1 | batch)")
  }else{
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient)")
  }
  
  
  model_formula <- as.formula(formula_str)
  
  message("Fitting ", m_name, " (", nrow(fit_df), " cells across ", 
          length(unique(fit_df$Patient)), " patients)...")
  print(colnames(fit_df))
  # Fit model with error trapping
  fit_result <- tryCatch(
    glmmTMB(
      formula = model_formula,
      data    = fit_df,
      family  = family
    ),
    error = function(e) {
      warning("Fitting failed for ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    sum_fit <- summary(fit_result)
    cond_coefs <- as.data.frame(sum_fit$coefficients$cond)
    
    # Extract diagnostics and tidy dataframe
    coef_df <- cond_coefs %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_col,
        Formula    = formula_str,
        N_Cells    = nrow(fit_df),
        N_Patients = length(unique(fit_df$Patient)),
        AIC        = AIC(fit_result),
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, N_Cells, N_Patients, AIC,
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    all_terms_list[[m_name]] <- coef_df
  }
}

final_df <- do.call(rbind, all_terms_list)
rownames(final_df) <- NULL

write.csv(final_df,file = "glmm_exhausted_resultsv3.csv")

#####

seurat_obj  <- qread("sce_T.qs")
outcome_col     = "score_TLS1"

primary_var     = "Race"
covariates      = c("histology_v4", "sex_v2", "smoking_v2", "Grade_v4")
sample_col      = "Sample_ID"
depth_col       = "nCount_RNA"
use_depth_offset = FALSE
family          = gaussian()
random_batch = FALSE

meta <- seurat_obj@meta.data



# 1. Verify required columns exist
req_cols <- c(outcome_col, primary_var, sample_col, covariates)
if (use_depth_offset) req_cols <- c(req_cols, depth_col)

missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 2. Build base modeling dataframe
base_df <- data.frame(
  Outcome = as.numeric(meta[[outcome_col]]),
  Primary = meta[[primary_var]],
  Patient = as.factor(meta[[sample_col]]),
  stringsAsFactors = FALSE
)
colnames(base_df)[1:2] <- c(outcome_col, primary_var)

# Append covariates
for (cv in covariates) {
  base_df[[cv]] <- meta[[cv]]
}

if (use_depth_offset) {
  base_df$Depth_Offset <- log(as.numeric(meta[[depth_col]]))
}

if (random_batch){
  base_df$batch <- meta$batch
}

# 3. Define the hierarchical model specifications
model_specs <- list(
  "Model_1_Unadjusted" = primary_var
)

for (cv in covariates) {
  model_specs[[paste0("Model_2_Plus_", cv)]] <- c(primary_var, cv)
}

model_specs[["Model_3_Fully_Adjusted"]] <- c(primary_var, covariates)

# 4. Iterate across models
all_terms_list <- list()

for (m_name in names(model_specs)) {
  current_vars <- model_specs[[m_name]]
  
  # Subset to complete cases for the active formula
  active_cols <- c(outcome_col, "Patient", current_vars)
  if (use_depth_offset) active_cols <- c(active_cols, "Depth_Offset")
  
  fit_df <- base_df[complete.cases(base_df[, active_cols, drop = FALSE]), , drop = FALSE]
  
  # Construct formula
  fixed_effects <- paste(current_vars, collapse = " + ")
  offset_term   <- if (use_depth_offset) " + offset(Depth_Offset)" else ""
  
  
  if (random_batch){
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient) + (1 | batch)")
  }else{
    formula_str   <- paste0(outcome_col, " ~ ", fixed_effects, offset_term, " + (1 | Patient)")
  }
  
  
  model_formula <- as.formula(formula_str)
  
  message("Fitting ", m_name, " (", nrow(fit_df), " cells across ", 
          length(unique(fit_df$Patient)), " patients)...")
  print(colnames(fit_df))
  # Fit model with error trapping
  fit_result <- tryCatch(
    glmmTMB(
      formula = model_formula,
      data    = fit_df,
      family  = family
    ),
    error = function(e) {
      warning("Fitting failed for ", m_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(fit_result)) {
    sum_fit <- summary(fit_result)
    cond_coefs <- as.data.frame(sum_fit$coefficients$cond)
    
    # Extract diagnostics and tidy dataframe
    coef_df <- cond_coefs %>%
      tibble::rownames_to_column(var = "Term") %>%
      dplyr::rename(
        Estimate  = Estimate,
        Std_Error = `Std. Error`,
        z_value   = `z value`,
        p_value   = `Pr(>|z|)`
      ) %>%
      dplyr::mutate(
        Model_Name = m_name,
        Outcome    = outcome_col,
        Formula    = formula_str,
        N_Cells    = nrow(fit_df),
        N_Patients = length(unique(fit_df$Patient)),
        AIC        = AIC(fit_result),
        CI_Lower   = Estimate - 1.96 * Std_Error,
        CI_Upper   = Estimate + 1.96 * Std_Error
      ) %>%
      dplyr::select(
        Model_Name, Outcome, Formula, N_Cells, N_Patients, AIC,
        Term, Estimate, Std_Error, CI_Lower, CI_Upper, z_value, p_value
      )
    
    all_terms_list[[m_name]] <- coef_df
  }
}

final_df <- do.call(rbind, all_terms_list)
rownames(final_df) <- NULL

write.csv(final_df,file = "glmm_tls_resultsv3.csv")
