rm(list=ls())

library(ggplot2)
library(ggsignif)
library(cowplot)

# Load and preprocess function
load_and_preprocess <- function(file_path) {
  load(file_path)
  vect_main <- as.data.frame(all)
  vect_main$AUC[which(vect_main$AUC == "SameVal")] <- 0
  vect_main$AUPR[which(vect_main$AUPR == "SameVal")] <- 0
  vect_main <- vect_main[!(vect_main$Method %in% c("MIp", "CAPS_Bootstrap")), ]
  vect_main$Method[vect_main$Method == "ZMIp"] <- "MIp"
  rm(all)
  return(vect_main)
}

# Load data
vect_main_major <- load_and_preprocess("AllToolsROCMeasures.RData")
load("Info_Revision1_AllProt.RData")
ids <- unique(vect_main_major$ID)
data <- data[match(ids, data$ID),]

# Clean and categorize
data$NumSeq <- as.numeric(data$NumSeq)
data$NumPos <- as.numeric(data$NumPos)
data$TreeLength <- as.numeric(data$TreeLength)

# Thresholds
num_seq <- 250
num_pos <- 500
tree_len <- 50

# Group definitions
group_sets <- list(
  "NumSeq" = list(
    "Small" = data$ID[data$NumSeq <= num_seq[1]],
    "Large" = data$ID[data$NumSeq > num_seq[1]]
  ),
  "NumPos" = list(
    "Small" = data$ID[data$NumPos <= num_pos[1]],
    "Large" = data$ID[data$NumPos > num_pos[1]]
  ),
  "TreeLen" = list(
    "Small" = data$ID[data$TreeLength <= tree_len[1]],
    "Large" = data$ID[data$TreeLength > tree_len[1]]
  )
)

# Final fixed method order
method_order <- c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap")

# Updated generate_plot function
generate_plot <- function(data, ids, measure, thr_val) {
  inds <- which(data$ID %in% ids & data$THR == thr_val)
  vect <- data[inds,]
  vect$Method <- factor(vect$Method, levels = method_order)
  
  ggplot(vect, aes(x = Method, y = as.numeric(vect[[measure]]))) +
    geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5, color = "grey40") +
    geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.3, color = "black", fill = "lightblue") +
    geom_signif(
      comparisons = list(c("PHACE", "DCA"), c("PHACE", "GaussDCA")),
      y_position = c(1.03, 1.08),  # Adjust these values to control bar height
      textsize = 4,
      vjust = 0.75, 
      map_signif_level = TRUE,
      step_increase = 0,
      color = "black"
    ) +
    coord_cartesian(ylim = c(0, 1.1)) +
    labs(x = "", y = measure) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12)
    )
}

# Loop through attributes and create 3x4 plots with labels A, B, C
for (attribute in names(group_sets)) {
  plot_list <- list()
  panel_labels <- c("A", "B", "C")  # Only one label per row
  
  i <- 1
  for (size in c("Small", "Large")) {
    ids <- group_sets[[attribute]][[size]]
    
    # Generate 4 plots for this size category
    auc_all <- generate_plot(vect_main_major, ids, "AUC", 0)
    auc_5aa <- generate_plot(vect_main_major, ids, "AUC", 5)
    aupr_all <- generate_plot(vect_main_major, ids, "AUPR", 0)
    aupr_5aa <- generate_plot(vect_main_major, ids, "AUPR", 5)
    
    plot_list <- c(plot_list, list(auc_all, auc_5aa, aupr_all, aupr_5aa))
    i <- i + 1
  }
  
  # Combine into a 3x4 grid with row labels A–C only
  final_plot <- plot_grid(
    # First row with label "A"
    plot_grid(plotlist = plot_list[1:4], ncol = 4),
    # Second row with label "B"
    plot_grid(plotlist = plot_list[5:8], ncol = 4),
    ncol = 1,
    labels = c("A", "B"),
    label_size = 16,
    label_x = 0, label_y = 1
  )
  
  
  ggsave(sprintf("SuppFigure_%s_Combined.png", attribute), final_plot, width = 16, height = 10, units = "in")
}
