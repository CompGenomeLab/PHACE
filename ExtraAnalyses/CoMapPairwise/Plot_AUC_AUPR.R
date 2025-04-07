# Load necessary libraries
library(ggplot2)
library(ggsignif)
library(cowplot)

# Function to load and preprocess data
load_and_preprocess <- function(file_path) {
  load(file_path)
  vect_main <- as.data.frame(all)
  vect_main$AUC[which(vect_main$AUC == "SameVal")] <- 0
  vect_main$AUPR[which(vect_main$AUPR == "SameVal")] <- 0
  vect_main <- vect_main[-which(vect_main$Method == "MIp"),]
  return(vect_main)
}

# Function to generate and save plots
generate_and_save_plots <- function(vect_main, measure, file_name_suffix) {
  thrs <- c(0, 5)
  closes <- 8
  
  plot_list <- list()
  for (close in closes) {
    for (thr in thrs) {
      vect <- vect_main[which(vect_main$THR == thr),]
      vect[which(vect[, 1] == "ZMIp"), 1] <- "MIp"
      
      # Subset data by method and convert to numeric
      method_list <- c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap")
      values_list <- lapply(method_list, function(method) as.numeric(vect[which(vect[, 1] == method), measure]))
      
      # Create dataframe ensuring that all algorithms appear, even if they have no data
      df <- data.frame(
        Algorithm = rep(method_list, each = length(values_list[[1]])),
        Value = unlist(values_list)
      )
      
      # Setting factor levels to ensure the order of algorithms
      df$Algorithm <- factor(df$Algorithm, levels = method_list)
      
      p <- ggplot(df, aes(x = Algorithm, y = Value)) +
        geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5, color = "grey40") +
        geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.3, color = "black", fill = "lightblue") +
        theme_minimal() +
        labs(title = "", x = "", y = measure) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)
        ) +
        geom_signif(comparisons = list(c("PHACE", "DCA"), c("PHACE", "GaussDCA")), 
                    textsize = 8, vjust = -0.5, step_increase = 0.05, 
                    map_signif_level = TRUE, color = "black") +
        coord_cartesian(ylim = c(0, 1.1))
      
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
  
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  ggsave(sprintf("All_ROC_Paper_%s.pdf", file_name_suffix), combined_plot, width = 21, height = 10, units = "cm")
}

# Main execution block
file_path <- "/Users/nurdankuru/Desktop/PHACE_Rev1_Analyses/AllToolsROCMeasures_Revision1_All.RData"
vect_main <- load_and_preprocess(file_path)

# Generate and save plots for AUC
generate_and_save_plots(vect_main, "AUC", "Revision1_AUC")

# Generate and save plots for AUPR
generate_and_save_plots(vect_main, "AUPR", "Revision1_AUPR")
