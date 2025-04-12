rm(list=ls())

vect_main <- read.csv("ROCPRMeasures_Revision1.csv")
vect_main <- as.data.frame(vect_main)

thrs <- c(0, 5)
closes <- 8

for (close in closes){
  for (thr in thrs){
    
    vect <- vect_main[which(vect_main$THR==thr),]
    
    w <- which(vect[,1]=="ZMIp")
    vect[w,1] <- "MIp"
    
    P1 <- vect[which(vect[,1]=="PHACE"),]
    P2 <- vect[which(vect[,1]=="DCA"),]
    P3 <- vect[which(vect[,1]=="GaussDCA"),]
    P4 <- vect[which(vect[,1]=="PSICOV"),]
    P5 <- vect[which(vect[,1]=="MIp"),]
    P6 <- vect[which(vect[,1]=="CAPS"),]
    P7 <- vect[which(vect[,1]=="CoMap"),]
    
    val1 <- as.numeric(P1$AUPR)
    val2 <- as.numeric(P2$AUPR)
    val3 <- as.numeric(P3$AUPR)
    val4 <- as.numeric(P4$AUPR)
    val5 <- as.numeric(P5$AUPR)
    val6 <- as.numeric(P6$AUPR)
    val7 <- as.numeric(P7$AUPR)
    
    library(ggplot2)
    library(ggsignif)
    
    df <- data.frame(
      Algorithm = rep(c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap"), each = length(val1)),
      Value = c(val1, val2, val3, val4, val5, val6, val7)
    )
    
    print(t.test(val1, val2))
    df$Algorithm <- factor(df$Algorithm, levels = c("PHACE", "DCA", "GaussDCA", "PSICOV", "CAPS", "MIp", "CoMap"))
    
    p <- ggplot(df, aes(x = Algorithm, y = Value)) +
      geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5, color = "grey40") +
      geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.3, color = "black", fill = "lightblue") +
      theme_minimal() +
      labs(title = "",
           x = "", y = "AUPR") +
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
    
    if (thr == 0) {
      p1 <- p
    } else if (thr == 5) {
      p2 <- p
    }
    
    
  }
  
}

library(cowplot)

combined_plot <- plot_grid(p1, p2, ncol = 2)
combined_plot 
ggsave(sprintf("AllTools_Revision1_AUPR.pdf"), combined_plot, width = ( 210/25.4), height = 4)


