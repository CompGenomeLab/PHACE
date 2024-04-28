library(cowplot)

load("Paper_ROC_Measures_AllTools.RData")
  
vect_main <- vect_all
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
      
      val1 <- as.numeric(P1$AUC)
      val2 <- as.numeric(P2$AUC)
      val3 <- as.numeric(P3$AUC)
      val4 <- as.numeric(P4$AUC)
      val5 <- as.numeric(P5$AUC)
      val6 <- as.numeric(P6$AUC)
      val7 <- as.numeric(P7$AUC)

      library(ggplot2)
      library(ggsignif)
      
      df <- data.frame(
        Algorithm = rep(c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap"), each = length(val1)),
        Value = c(val1, val2, val3, val4, val5, val6, val7)
      )
      
      df$Algorithm <- factor(df$Algorithm, levels = c("PHACE", "DCA", "GaussDCA", "PSICOV", "MIp", "CAPS", "CoMap"))

      p <- ggplot(df, aes(x = Algorithm, y = Value, color = Algorithm)) +
          geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5) +
          geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.2, color = "black") +
          theme_minimal() +
          labs(title = "",
               x = "", y = "AUC") +
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),   
                axis.text.y = element_text(size = 12),  
                axis.title.x = element_text(size = 14), 
                axis.title.y = element_text(size = 14)) + 
          theme(legend.position = "topright") +
          coord_cartesian(ylim = c(0, 1.09))
        
        if (thr == 0) {
          p1 <- p
        } else if (thr == 5) {
          p2 <- p
        }
      
      
    }
    
}
  
combined_plot <- plot_grid(p1, p2, ncol = 2)
combined_plot <- combined_plot +
  draw_label("a", x = 0, y = 1, size = 20, vjust = 1, hjust = 0) +
  draw_label("b", x = 0.5, y = 1, size = 20, vjust = 1, hjust = 0) 
combined_plot 

ggsave(sprintf("All_ROC_Paper.pdf"), combined_plot, width = (210/25.4), height = 4)
