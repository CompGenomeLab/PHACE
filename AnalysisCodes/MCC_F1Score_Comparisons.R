library(mltools)

ind <- as.numeric(args[1])
alg <- args[2]

ids <- list.files(sprintf("Pairwise%s", alg))
ids <- unlist(strsplit(ids, "_AllRes_PDB.RData"))
ids <- ids[ind]

close_vals <- 8
far_vals <- 16
thrs <- c(0, 5) 
vect <- c()
for (far_val in far_vals){
  for (close_val in close_vals){
    for (thr in thrs){
      for (id in ids){
        print(c(id, length(ids)))
        load(sprintf("Pairwise%s/%s_AllRes_PDB.RData", alg, id))
        
        ind <- c(which(colnames(mat)=="Pos1"), which(colnames(mat)=="Pos2"), which(colnames(mat)=="Dist"),
                 which(colnames(mat)=="Score2"), which(colnames(mat)=="Score3"),
                 which(colnames(mat)==alg))
        mat <- mat[,ind]
        
        val <- as.numeric(mat$Score3)
        val[which(as.numeric(mat$Score3)<=-1)] <- -1
        
        mat2 <- cbind(as.numeric(mat$Pos1), as.numeric(mat$Pos2), val)
        interaction_data <- as.data.frame(mat2)
        colnames(interaction_data) <- c("Pos1", "Pos2", "Score")
        
        # Calculate mean interaction score
        mean_score <- mean(interaction_data$Score)
        mean_score2 <- mean(c(interaction_data$Score, matrix(1, max(as.numeric(mat$Pos2)), 1)))
        interaction_data$Corrected_Score <- interaction_data$Score - mean_score
        interaction_data$Corrected_Score2 <- interaction_data$Score - mean_score2
        
        mat <- as.data.frame(mat)
        mat$PHACE <- interaction_data$Corrected_Score2
        
        nas <- which(is.na(mat$Dist)==1)
        if (length(nas)>0){
          mat <- mat[-nas,]
        }
        
        difs <- which(abs(as.numeric(mat$Pos1)-as.numeric(mat$Pos2))>=thr)
        mat <- mat[difs,]
        
        co <- which(as.numeric(mat$Dist)<=close_val)
        co <- setdiff(co, which(mat$Score2==-2))
        sub1 <- mat[co,]
        
        nonco <- which(as.numeric(mat$Dist)>=far_val)
        sub2 <- mat[-co,]
        
        mat <- rbind(sub1, sub2)
        
        if (length(sub1$Pos1)>=5 && length(sub2$Pos1)>=5){
          
          co <- which(as.numeric(mat$Dist)<=close_val)          
          nonco <- which(as.numeric(mat$Dist)>=far_val)
          
          co_pred_phace <- which(as.numeric(mat$ScoreUpd)>=0.25)
          nonco_pred_phace <- setdiff(1:length(mat$Pos1), co_pred_phace)
          
          tp_phace <- length(intersect(co, co_pred_phace))
          tn_phace <- length(intersect(nonco, nonco_pred_phace))
          fp_phace <- length(intersect(nonco, co_pred_phace))
          fn_phace <- length(intersect(co, nonco_pred_phace))
          mcc_phace <- mcc(TP = tp_phace, FP=fp_phace, TN=tn_phace, FN=fn_phace)
          
          pres_phace <- tp_phace/(tp_phace+fp_phace)
          rec_phace <- tp_phace/(tp_phace+fn_phace)
          f1_phace <- (2*pres_phace*rec_phace) / (pres_phace + rec_phace)
          
          ppv_phace <- tp_phace/(tp_phace+fp_phace)
          npv_phace <- tn_phace/(tn_phace+fn_phace)
          
          co_pred_alg <- which(is.na(mat$CoMap)!=1)
          nonco_pred_alg <- setdiff(1:length(mat$Pos1), co_pred_alg)
          
          tp_alg <- length(intersect(co, co_pred_alg))
          tn_alg <- length(intersect(nonco, nonco_pred_alg))
          fp_alg <- length(intersect(nonco, co_pred_alg))
          fn_alg <- length(intersect(co, nonco_pred_alg))
          mcc_alg <- mcc(TP = tp_alg, FP=fp_alg, TN=tn_alg, FN=fn_alg)
          
          pres_alg <- tp_alg/(tp_alg+fp_alg)
          rec_alg <- tp_alg/(tp_alg+fn_alg)
          f1_alg <- (2*pres_alg*rec_alg) / (pres_alg + rec_alg)
          
          ppv_alg <- tp_alg/(tp_alg+fp_alg)
          npv_alg <- tn_alg/(tn_alg+fn_alg)
          
          vect <- rbind(vect, c("PHACE", thr, 0.25,  close_val, far_val, id, tp_phace, fp_phace, tn_phace, fn_phace, pres_phace, rec_phace, mcc_phace, f1_phace, ppv_phace, npv_phace))
          vect <- rbind(vect, c(alg, thr, "Predefined", close_val, far_val, id, tp_alg, fp_alg, tn_alg, fn_alg, pres_alg, rec_alg, mcc_alg, f1_alg, ppv_alg, npv_alg))
          colnames(vect) <- c("Method", "Thr", "Thr2", "Close", "Far", "ID", "TP", "FP", "TN", "FN", "PRES", "REC", "MCC", "F1", "PPV", "NPV")
          
        }
      }
    }
  }
}


save(vect, file = sprintf("%s_MCC_F1.RData", alg))


