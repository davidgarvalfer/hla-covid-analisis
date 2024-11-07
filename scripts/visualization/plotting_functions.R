#' Visualization Functions
#' Collection of functions for creating visualizations

#' Create forest plot
#' @param results Results data frame
#' @param locus HLA locus
#' @return ggplot object
create_forest_plot <- function(results, locus) {
  p <- ggplot(results, aes(x = reorder(Allele, OR))) +
    geom_point(aes(y = OR, size = Frequency, color = Effect)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = Effect), 
                  width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
    coord_flip() +
    scale_color_manual(values = c("Risk" = "red", 
                                  "Protective" = "blue", 
                                  "Non-significant" = "gray")) +
    scale_y_log10(limits = c(0.1, 10),
                  breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    theme_minimal() +
    labs(title = paste("Forest Plot -", locus),
         subtitle = "Adjusted odds ratios with 95% CI",
         x = "Allele",
         y = "Odds Ratio (log scale)",
         size = "Frequency",
         color = "Effect")
  
  return(p)
}

#' Create heatmap
#' @param data Data frame
#' @param locus HLA locus
#' @return ggplot object
create_heatmap <- function(data, locus) {
  p <- ggplot(data, aes(x = Outcome, y = Allele, fill = Residual)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0,
                         limits = c(-max(abs(data$Residual)),
                                    max(abs(data$Residual)))) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    labs(title = paste("Association Heatmap -", locus),
         subtitle = "Standardized residuals",
         x = "Clinical Outcome",
         y = "Allele",
         fill = "Residual")
  
  return(p)
}

#' Create correlation plot
#' @param data Data frame
#' @param locus HLA locus
#' @return ggplot object
create_correlation_plot <- function(data, locus) {
  p <- ggplot(data, 
              aes(x = Frequency, y = OR)) +
    geom_point(aes(color = P_value < 0.05)) +
    geom_text_repel(aes(label = Allele),
                    size = 3,
                    max.overlaps = Inf,
                    box.padding = 0.5) +
    geom_hline(yintercept = 1, 
               linetype = "dashed", 
               color = "red") +
    scale_y_log10() +
    scale_color_manual(values = c("grey", "red"),
                       labels = c("Non-significant", "Significant")) +
    theme_minimal() +
    labs(title = paste("Frequency vs Effect Size -", locus),
         subtitle = "Odds ratios by allele frequency",
         x = "Allele Frequency",
         y = "Odds Ratio (log scale)",
         color = "Significance")
  
  return(p)
}

#' Create summary plots
#' @param results Results list
#' @param locus HLA locus
#' @return List of plots
create_summary_plots <- function(results, locus) {
  plots <- list()
  
  # Forest plot
  plots$forest <- create_forest_plot(results$or_results, locus)
  
  # Heatmap
  plots$heatmap <- create_heatmap(results$residuals, locus)
  
  # Correlation plot
  plots$correlation <- create_correlation_plot(results$or_results, locus)
  
  return(plots)
}

#' Save plots
#' @param plots List of plots
#' @param locus HLA locus
save_plots <- function(plots, locus) {
  plot_types <- names(plots)
  
  for(type in plot_types) {
    filename <- file.path("plots", "severity",
                          paste0("HLA_", locus, "_", type, ".pdf"))
    
    ggsave(filename,
           plots[[type]],
           width = 10,
           height = 8,
           dpi = 300)
  }
}

#' Generate comprehensive analysis report
#' @param results Analysis results
#' @param locus HLA locus
generate_report <- function(results, locus) {
  report_file <- file.path("results", "severity",
                           paste0("HLA_", locus, "_report.txt"))
  
  sink(report_file)
  
  cat("ANALYSIS REPORT FOR HLA-", locus, "\n")
  cat("===============================\n\n")
  
  # Sample statistics
  cat("Sample Statistics:\n")
  cat("-----------------\n")
  print(results$sample_stats)
  cat("\n")
  
  # Association results
  cat("Association Results:\n")
  cat("-------------------\n")
  print(results$association_test)
  cat("\n")
  
  # Top findings
  cat("Significant Findings:\n")
  cat("-------------------\n")
  significant_results <- subset(results$or_results, P_value < 0.05)
  print(significant_results)
  
  sink()
}

# Example usage:
# results <- analyze_hla_data(data, "A")
# plots <- create_summary_plots(results, "A")
# save_plots(plots, "A")
# generate_report(results, "A")
