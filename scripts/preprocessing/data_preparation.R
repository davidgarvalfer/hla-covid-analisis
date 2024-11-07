############################################################
# HLA-COVID Analysis
# Author: David García Valentín-Fernández
# Description: Analysis of HLA alleles and COVID-19 severity
############################################################

#' Setup Project Environment
#' Loads required packages and creates necessary directories

# 1. Required Packages
required_packages <- c(
  "HIBAG",
  "data.table",
  "ggplot2",
  "corrplot",
  "pROC",
  "qvalue",
  "dplyr",
  "tidyr",
  "gridExtra",
  "ggrepel",
  "reshape2",
  "caret",
  "factoextra",
  "car",
  "glmnet",
  "MASS",
  "brglm2"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

# 2. Project Structure
project_dirs <- c(
  "data/raw",              # Raw input data
  "data/processed",        # Processed data
  "results/severity",      # Severity analysis results
  "results/hospitalization", # Hospitalization analysis results
  "results/asymptomatic",   # Asymptomatic analysis results
  "results/multivariable",  # Multivariable analysis results
  "plots/severity",         # Severity plots
  "plots/hospitalization",  # Hospitalization plots
  "plots/asymptomatic",     # Asymptomatic plots
  "plots/multivariable"     # Multivariable plots
)

# Create directories if they don't exist
invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# 3. Data Loading and Initial Processing
#' Load and preprocess the clinical and genetic data
#' @param clinical_file Path to clinical data file
#' @param model_file Path to HLA model file
#' @param bed_prefix Prefix for BED files

load_and_process_data <- function(clinical_file, model_file, bed_prefix) {
  # Load clinical data
  clinical_data <- read.table(clinical_file, header = TRUE, sep = "\t")
  
  # Load HLA model
  load(model_file)
  
  # Load genetic data
  geno_data <- hlaBED2Geno(
    bed.fn = paste0(bed_prefix, ".bed"),
    fam.fn = paste0(bed_prefix, ".fam"),
    bim.fn = paste0(bed_prefix, ".bim")
  )
  
  # Clean numeric data
  clinical_data <- clean_numeric_data(clinical_data)
  
  # Define HLA loci
  hla_loci <- c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1")
  
  return(list(
    clinical = clinical_data,
    genetic = geno_data,
    loci = hla_loci,
    model = if(exists("model.list")) model.list else NULL
  ))
}

# 4. Helper Functions
#' Clean numeric data in clinical dataset
#' @param data Clinical data frame
#' @return Cleaned data frame

clean_numeric_data <- function(data) {
  # Clean age
  data$REDAD <- as.numeric(gsub(",", ".", as.character(data$REDAD)))
  
  # Clean PCs
  pc_cols <- grep("^PC", names(data), value = TRUE)
  for(col in pc_cols) {
    data[[col]] <- as.numeric(gsub(",", ".", as.character(data[[col]])))
  }
  
  return(data)
}

# 5. Data Validation
#' Validate required variables in clinical data
#' @param data Clinical data frame
#' @return Warning messages if variables are missing

validate_clinical_data <- function(data) {
  required_vars <- c("id", "R_SEV4", "R_HOSP", "ASINT")
  missing_vars <- required_vars[!required_vars %in% names(data)]
  
  if(length(missing_vars) > 0) {
    warning("Missing required variables: ", 
            paste(missing_vars, collapse = ", "))
  }
  
  cat("Data structure summary:\n")
  cat("- Number of samples:", nrow(data), "\n")
  cat("- Available variables:", paste(names(data), collapse = ", "), "\n")
}

# Example usage:
# data <- load_and_process_data(
#     clinical_file = "path/to/clinical.txt",
#     model_file = "path/to/model.RData",
#     bed_prefix = "path/to/genetic_data"
# )
# validate_clinical_data(data$clinical)