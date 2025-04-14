rm(list=ls())
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

# Function to generate plots and return list
generate_plots <- function(vect_main, measure) {
  thrs <- c(0, 5)
  closes <- 8
  plot_list <- list()
  
  for (close in closes) {
    for (thr in thrs) {
      vect <- vect_main[which(vect_main$THR == thr),]
      vect[which(vect[, 1] == "ZMIp"), 1] <- "MIp"
      method_list <- c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap")
      values_list <- lapply(method_list, function(method) as.numeric(vect[which(vect[, 1] == method), measure]))
      
      df <- data.frame(
        Algorithm = rep(method_list, each = length(values_list[[1]])),
        Value = unlist(values_list)
      )
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
        geom_signif(
          comparisons = list(c("PHACE", "DCA"), c("PHACE", "GaussDCA")),
          y_position = c(1.03, 1.08),  # Adjust these values to control bar height
          textsize = 4,
          vjust = 0.75, 
          map_signif_level = TRUE,
          step_increase = 0,
          color = "black"
        ) +
        coord_cartesian(ylim = c(0, 1.1))
      
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
  return(plot_list)
}

# Main execution block
file_path <- "AllToolsROCMeasures_Revision1_162.RData"
vect_main <- load_and_preprocess(file_path)

# Get individual plots
auc_plots <- generate_plots(vect_main, "AUC")
aupr_plots <- generate_plots(vect_main, "AUPR")

# Combine into one figure (2 rows, 2 columns)
combined_plot <- plot_grid(
  plotlist = c(auc_plots, aupr_plots),
  labels = c("A", "", "B", ""),
  ncol = 2, nrow = 2
)

print(combined_plot)

# Save the figure
ggsave("All_ROC_Paper_Combined.png", combined_plot, width = 21, height = 20, units = "cm", dpi = 300)

