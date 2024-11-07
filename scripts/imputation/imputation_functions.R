############################################################
# HLA-COVID Analysis: Imputation and Filtering Functions
# Author: David García Valentín-Fernández
# Part 2: Core analysis functions
############################################################

#' HLA Imputation Function
#' Performs HLA imputation and generates quality metrics
#' @param hla_id HLA locus identifier
#' @param genetic_data Genetic data from hlaBED2Geno
#' @param clinical_data Clinical data frame
#' @param model_list HLA model list
#' @return List containing imputation results and metrics

impute_hla <- function(hla_id, genetic_data, clinical_data, model_list) {
  tryCatch({
    # Get appropriate model
    model <- get_hla_model(hla_id, model_list)
    if(is.null(model)) return(NULL)
    
    # Prepare data for imputation
    prepared_data <- prepare_imputation_data(
      genetic_data, 
      clinical_data
    )
    
    # Perform imputation
    prediction <- predict_hla(model, genetic_data, prepared_data$covariates)
    if(is.null(prediction)) return(NULL)
    
    # Calculate and return metrics
    metrics <- calculate_imputation_metrics(prediction, hla_id)
    
    return(list(
      genotypes = prediction$value,
      probabilities = prediction$postprob,
      metrics = metrics
    ))
  }, error = function(e) {
    log_error("HLA imputation", hla_id, e)
    return(NULL)
  })
}

#' Data Filtering Functions
#' A collection of functions to filter and clean the data

#' Filter missing values
#' @param data Data frame
#' @param columns Columns to check for missing values
#' @return Filtered data frame
filter_missing_values <- function(data, columns) {
  complete_cases <- complete.cases(data[, columns])
  return(data[complete_cases, ])
}

#' Filter by allele frequency
#' @param data Data frame
#' @param allele_col Column containing allele information
#' @param min_freq Minimum frequency threshold
#' @return Filtered data frame
filter_by_allele_frequency <- function(data, allele_col, min_freq = 10) {
  freq_table <- table(data[[allele_col]])
  valid_alleles <- names(freq_table)[freq_table >= min_freq]
  return(data[data[[allele_col]] %in% valid_alleles, ])
}

#' Process and scale covariates
#' @param data Data frame containing covariates
#' @return Processed data frame with scaled covariates
process_covariates <- function(data) {
  # Process sex
  data$sex_numeric <- as.numeric(factor(data$sexo_redcap, 
                                        levels = c("mujer", "hombre")))
  data$sex_numeric[is.na(data$sex_numeric)] <- get_mode(data$sex_numeric)
  
  # Process age
  data$age_scaled <- scale_variable(data$REDAD)
  
  # Process PCs
  pc_columns <- grep("^PC", names(data), value = TRUE)
  data[pc_columns] <- scale_pc_matrix(data[pc_columns])
  
  return(data)
}

#' Comprehensive data filtering pipeline
#' @param data Data frame to filter
#' @param allele_col Column containing allele information
#' @param required_vars Required variables to check
#' @param min_freq Minimum frequency threshold
#' @return List containing filtered data and filtering statistics
filter_data <- function(data, 
                        allele_col, 
                        required_vars = c("R_SEV4", "R_HOSP", "ASINT"),
                        min_freq = 10) {
  
  initial_count <- nrow(data)
  
  # Process covariates
  data_processed <- process_covariates(data)
  
  # Filter missing values
  data_complete <- filter_missing_values(data_processed, required_vars)
  after_missing_count <- nrow(data_complete)
  
  # Filter by frequency
  data_filtered <- filter_by_allele_frequency(data_complete, allele_col, min_freq)
  final_count <- nrow(data_filtered)
  
  # Calculate statistics
  stats <- calculate_filtering_stats(
    initial_count,
    after_missing_count,
    final_count
  )
  
  return(list(
    filtered_data = data_filtered,
    stats = stats
  ))
}

# Helper functions

#' Get HLA model
#' @param hla_id HLA locus identifier
#' @param model_list List of models
#' @return HLA model object
get_hla_model <- function(hla_id, model_list) {
  if(hla_id %in% names(model_list)) {
    return(hlaModelFromObj(model_list[[hla_id]]))
  }
  return(NULL)
}

#' Prepare data for imputation
#' @param genetic_data Genetic data
#' @param clinical_data Clinical data
#' @return List of prepared data
prepare_imputation_data <- function(genetic_data, clinical_data) {
  sample_ids <- genetic_data$sample.id
  matched_data <- clinical_data[match(sample_ids, clinical_data$id), ]
  processed_data <- process_covariates(matched_data)
  
  return(list(
    matched_data = matched_data,
    covariates = create_covariate_matrix(processed_data)
  ))
}

#' Create covariate matrix
#' @param data Processed data frame
#' @return Matrix of covariates
create_covariate_matrix <- function(data) {
  pc_cols <- grep("^PC.*_scaled$", names(data), value = TRUE)
  covariates <- cbind(
    sex = data$sex_numeric,
    age = data$age_scaled,
    as.matrix(data[, pc_cols])
  )
  return(covariates)
}

#' Calculate imputation metrics
#' @param prediction Prediction results
#' @param hla_id HLA locus identifier
#' @return List of metrics
calculate_imputation_metrics <- function(prediction, hla_id) {
  max_probs <- apply(prediction$postprob, 2, max)
  
  metrics <- list(
    low_conf = mean(max_probs < 0.5) * 100,
    med_conf = mean(max_probs >= 0.5 & max_probs < 0.75) * 100,
    high_conf = mean(max_probs >= 0.75) * 100
  )
  
  save_probability_plot(max_probs, hla_id)
  
  return(metrics)
}

#' Save probability distribution plot
#' @param probs Vector of probabilities
#' @param hla_id HLA locus identifier
save_probability_plot <- function(probs, hla_id) {
  pdf(file.path("plots", "severity", 
                paste0("HLA_", hla_id, "_probability_dist.pdf")))
  hist(probs, 
       main = paste("Posterior Probability Distribution - HLA-", hla_id),
       xlab = "Posterior Probability",
       ylab = "Frequency",
       breaks = 20,
       col = "lightblue",
       border = "darkblue")
  add_probability_thresholds()
  dev.off()
}

#' Calculate filtering statistics
#' @param initial Initial count
#' @param after_missing Count after missing value filtering
#' @param final Final count
#' @return List of statistics
calculate_filtering_stats <- function(initial, after_missing, final) {
  list(
    initial_count = initial,
    missing_filtered = initial - after_missing,
    frequency_filtered = after_missing - final,
    final_count = final,
    total_filtered_percent = ((initial - final) / initial) * 100
  )
}

#' Error logging function
#' @param process Process name
#' @param id Identifier
#' @param error Error object
log_error <- function(process, id, error) {
  message(sprintf("Error in %s for %s: %s", process, id, conditionMessage(error)))
}

# Utility functions
get_mode <- function(x) {
  as.numeric(names(which.max(table(x, useNA = "no"))))
}

scale_variable <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(scale(x))
}

scale_pc_matrix <- function(pc_data) {
  pc_matrix <- as.matrix(sapply(pc_data, as.numeric))
  return(scale(pc_matrix))
}

# Example usage:
# results <- lapply(hla_loci, function(locus) {
#     impute_hla(locus, genetic_data, clinical_data, model_list)
# })