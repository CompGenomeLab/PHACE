# Clear the environment
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
  vect_main <- vect_main[-which(vect_main$Method == "CAPS_Bootstrap"),]
  rm(all)
  return(vect_main)
}

# Function to generate plots for each group
generate_plots <- function(data, ids, measure, plot_title) {
  inds <- unlist(lapply(ids, function(id) which(data$ID == id)))
  vect_main <- data[inds,]
  
  # Iterate through thresholds and closeness values
  plot_list <- list()
  for (thr in c(0, 5)) {
    vect <- vect_main[which(vect_main$THR == thr),]
    vect[which(vect[, 1] == "ZMIp"), 1] <- "MIp"
    
    # Collect values and match lengths correctly
    algorithms <- unique(vect[, 1])
    values_list <- lapply(algorithms, function(alg) as.numeric(vect[vect[, 1] == alg, measure]))
    
    # Ensure that all entries in values_list have values, otherwise, fill with NA
    max_length <- max(sapply(values_list, length))
    values_list <- lapply(values_list, function(x) {
      if (length(x) < max_length) {
        c(x, rep(NA, max_length - length(x)))
      } else {
        x
      }
    })
    
    # Create data frame
    df <- data.frame(
      Algorithm = rep(algorithms, each = max_length),
      Value = unlist(values_list)
    )
    
    p <- ggplot(df, aes(x = Algorithm, y = Value)) +
      geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5, color = "grey40") +
      geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.3, color = "black", fill = "lightblue") +
      geom_signif(comparisons = list(c("PHACE", "DCA"), c("PHACE", "GaussDCA")),
                  textsize = 8, vjust = -0.5, step_increase = 0.05, 
                  map_signif_level = TRUE, color = "black") +
      coord_cartesian(ylim = c(0, 1.1)) +
      labs(title = plot_title, x = "", y = measure) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14))
    plot_list[[length(plot_list) + 1]] <- p
  }
  return(plot_list)
}

# Load and preprocess data
vect_main_major <- load_and_preprocess("/Users/nurdankuru/Desktop/PHACE_Rev1_Analyses/AllToolsROCMeasures_Revision1_All.RData")

# Define groups
load("/Users/nurdankuru/Desktop/PHACE_Rev1_Analyses/Info_Revision1.RData")
group_ids <- list(
  "Num Seq Small" = data$ID[which(data$NumSeq <= 250)],
  "Num Seq Large" = data$ID[which(data$NumSeq > 250)],
  "Num Pos Small" = data$ID[which(data$NumPos <= 500)],
  "Num Pos Large" = data$ID[which(data$NumPos > 500)],
  "Tree Len Small" = data$ID[which(data$TreeLength <= 50)],
  "Tree Len Large" = data$ID[which(data$TreeLength > 50)]
)

# Generate and save plots
for (name_group in names(group_ids)) {
  plots <- generate_plots(vect_main_major, group_ids[[name_group]], "AUC", name_group)
  combined_plot <- plot_grid(plotlist = plots, ncol = 2)
  combined_plot
  ggsave(sprintf("All_ROC_Paper_%s.png", name_group), combined_plot, width = 10, height = 5, units = "in")
}

# Repeat for AUPR or other measures as needed
for (name_group in names(group_ids)) {
  plots <- generate_plots(vect_main_major, group_ids[[name_group]], "AUPR", name_group)
  combined_plot <- plot_grid(plotlist = plots, ncol = 2)
  combined_plot
  ggsave(sprintf("All_ROC_Paper_%s.png", name_group), combined_plot, width = 10, height = 5, units = "in")
}
