############################################################
# HLA-COVID Analysis: Analysis and Visualization Functions
# Author: David García Valentín-Fernández
# Part 3: Statistical analysis and visualization
############################################################

#' Statistical Analysis Functions
#' Collection of functions for statistical analysis of HLA data

#' Calculate adjusted odds ratio
#' @param data Data frame
#' @param allele Allele of interest
#' @param allele_col Column containing allele information
#' @param outcome_var Outcome variable
#' @return List with OR and confidence intervals
calculate_adjusted_or <- function(data, allele, allele_col, outcome_var = "R_SEV4") {
  tryCatch({
    # Create formula with covariates
    covariates <- paste(c("sex_numeric", "age_scaled", 
                          grep("^PC.*_scaled$", names(data), value = TRUE)), 
                        collapse = " + ")
    formula_str <- paste(outcome_var, "~", paste(allele_col, covariates, sep = " + "))
    model <- glm(as.formula(formula_str), data = data, family = binomial)
    
    # Extract coefficient
    coef_idx <- grep(paste0("^", allele_col), names(coef(model)))
    if(length(coef_idx) > 0) {
      coef <- coef(model)[coef_idx]
      se <- sqrt(diag(vcov(model)))[coef_idx]
      
      return(list(
        OR = exp(coef),
        CI_lower = exp(coef - 1.96*se),
        CI_upper = exp(coef + 1.96*se),
        P_value = summary(model)$coefficients[coef_idx, 4]
      ))
    }
    return(NULL)
  }, error = function(e) {
    log_error("OR calculation", allele, e)
    return(NULL)
  })
}

#' Perform association test
#' @param data Data frame
#' @param allele_col Allele column
#' @param outcome_var Outcome variable
#' @return Association test results
perform_association_test <- function(data, allele_col, outcome_var = "R_SEV4") {
  covariates <- paste(c("sex_numeric", "age_scaled", 
                        grep("^PC.*_scaled$", names(data), value = TRUE)), 
                      collapse = " + ")
  
  formula_null <- as.formula(paste(outcome_var, "~", covariates))
  formula_full <- as.formula(paste(outcome_var, "~", 
                                   paste(allele_col, covariates, sep = " + ")))
  
  model_null <- glm(formula_null, data = data, family = binomial)
  model_full <- glm(formula_full, data = data, family = binomial)
  
  lrt <- anova(model_null, model_full, test = "LRT")
  
  return(list(
    test_statistic = lrt$Deviance[2],
    p_value = lrt$`Pr(>Chi)`[2],
    df = lrt$Df[2]
  ))
}