# Load required packages
library(ggplot2)
library(gridExtra)

# Load the data
load("Info_Revision1.RData")

# Convert variables to numeric
data$NumSeq <- as.numeric(data$NumSeq)
data$NumPos <- as.numeric(data$NumPos)
data$TreeLength <- as.numeric(data$TreeLength)

# Create professional histograms
p1 <- ggplot(data, aes(x = NumSeq)) +
  geom_histogram(bins = 30, fill = "#2c3e50", color = "white", alpha = 0.8) +
  labs(title = "Distribution of Number of Sequences",
       x = "Number of Sequences",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

p2 <- ggplot(data, aes(x = NumPos)) +
  geom_histogram(bins = 30, fill = "#3498db", color = "white", alpha = 0.8) +
  labs(title = "Distribution of Number of Positions",
       x = "Number of Positions",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

p3 <- ggplot(data, aes(x = TreeLength)) +
  geom_histogram(bins = 30, fill = "#e74c3c", color = "white", alpha = 0.8) +
  labs(title = "Distribution of Tree Length",
       x = "Tree Length",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Arrange plots in a grid
combined_plot <- grid.arrange(p1, p2, p3, ncol = 1)

# Save the combined plot as a high-quality PDF
ggsave("histograms_combined_162Prot.png", 
       plot = combined_plot,
       width = 8.5,  # inches
       height = 11,  # inches
       dpi = 300,   # high resolution
       device = "png")
