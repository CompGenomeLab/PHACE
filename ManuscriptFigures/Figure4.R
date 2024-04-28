library(ggplot2)
library(ggsignif)

algs <- c("CAPS", "CoMap")
cols <- c("#F8766D", "#A58AFF", "#FB61D7")
measure <- "MCC"
# measure <- "F1"

plot_list <- list()
for (ind in 1:length(algs)){
  alg <- algs[ind]
  load(sprintf("Paper_Results_%s.RData", alg))
  
  vect_main <- vect
  vect_main <- as.data.frame(vect_main)
  rm(vect)
  
  phace_color <- cols[1]
  alg_color <- cols[(ind+1)]
  
  color_vector <- c(phace_color, alg_color)
  
  thrs <- c(0, 5)
  closes <- 8
  
  for (close in closes){
    for (thr in thrs){
      
      vect <- vect_main[which(vect_main$Thr==thr),]
      
      P1 <- vect[which(vect[,1]=="PHACE"),]
      P2 <- vect[which(vect[,1]==alg),]
      
      loc <- which(colnames(P1)==measure)
      val1 <- as.numeric(P1[,loc])
      val2 <- as.numeric(P2[,loc])
      
      val1[which(is.na(val1)==1)] <- 0
      val2[which(is.na(val2)==1)] <- 0
      
      df <- data.frame(
        Algorithm = rep(c("PHACE",alg), each = length(val1)),
        Value = c(val1, val2)
      )
      
      t_test_result <- t.test(val1, val2)
      p_value <- t_test_result$p.value
      
      df$Algorithm <- factor(df$Algorithm, levels = c("PHACE", alg))
      
      p <- ggplot(df, aes(x = Algorithm, y = Value, color = Algorithm)) +
        geom_jitter(position = position_jitter(width = 0.3), size = 1, alpha = 0.5) +
        geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.2, color = "black") +
        theme_minimal() +
        labs(title = "",
             x = "", y = measure) +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
              axis.text.y = element_text(size = 12),  
              axis.title.x = element_text(size = 14), 
              axis.title.y = element_text(size = 14)) + 
        scale_color_manual(values = color_vector) +
        geom_signif(comparisons = list(c("PHACE", alg)), textsize = 8, vjust = -0.5, step_increase = 0.05, map_signif_level = TRUE, color = "black") +
        theme(legend.position = "topright") +
        coord_cartesian(ylim = c(-0.2, 1))
      
      plot_list[[sprintf("%s_Thr%s", alg, thr)]] <- p
    }
  }
}


library(cowplot)

plot_list_Thr0 <- list()
plot_list_Thr5 <- list()

for (name in names(plot_list)) {
  if (grepl("_Thr0", name)) {
    plot_list_Thr0[[name]] <- plot_list[[name]]
  } else if (grepl("_Thr5", name)) {
    plot_list_Thr5[[name]] <- plot_list[[name]]
  }
}

combined_plot_Thr0 <- plot_grid(plotlist = plot_list_Thr0, ncol = 2, labels = "")
combined_plot_Thr5 <- plot_grid(plotlist = plot_list_Thr5, ncol = 2, labels = "")

combined_final <- plot_grid(combined_plot_Thr0, combined_plot_Thr5, nrow = 1, labels = c("A", "B"))
print(combined_final)

ggsave(sprintf("Paper_Pairwise_%s_%s.pdf", alg, measure), combined_final, width = ( 210/25.4), height = 4)


