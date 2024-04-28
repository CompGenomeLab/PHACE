library(mltools)
library(AUC)
library(PRROC)
args = commandArgs(trailingOnly=TRUE)

ind <- as.numeric(args[1])
alg <- args[2]

ids <- list.files(sprintf("Pairwise%s", alg))
ids <- unlist(strsplit(ids, "_AllRes_PDB.RData"))
ids <- ids[ind]

vect <- c()
thrs <- c(0, 5)
close_vals <- 8
far_val <- 16
for (thr in thrs){
  for (close_val in close_vals){
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
      nas <- unique(c(nas, which(is.na(mat[,6])==1)))
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
      
      if (length(sub2$Pos1)>length(sub1$Pos1)){
        sx <- sort(as.numeric(sub2$Dist), decreasing = T, index.return=T)
        sub2 <- sub2[sx$ix, ]
        
        pos <- unique(c(sub2$Pos1, sub2$Pos2))
        for (max in c(5, 10, 15, 20, 25)){
          sub22 <- sub2
          chosen_pairs <- c()
          freq <- matrix(0, 1, length(pos))
          names(freq) <- pos
          kx <- 1
          while (length(chosen_pairs[,1])<length(sub1$Pos1) && kx <= length(sub2$Pos1)) {
            kx <- kx + 1
            p1 <- as.numeric(sub22[1,1])
            p2 <- as.numeric(sub22[1,2])
            
            if (freq[which(names(freq)==p1)]<max && freq[which(names(freq)==p2)]<max){
              freq[which(names(freq)==p1)] <- freq[which(names(freq)==p1)] + 1
              freq[which(names(freq)==p2)] <- freq[which(names(freq)==p2)] + 1
              chosen_pairs <- rbind(chosen_pairs, unlist(sub22[1,1:2]))
              sub22 <- sub22[-1,]
            } else {
              sub22 <- sub22[-1,]
            }
          }
          
          if (length(chosen_pairs[,1])==length(sub1$Pos1)){
            break
          }
        }
        
        print(max)
        chosen_pairs_var <- paste(chosen_pairs[,1], chosen_pairs[,2], sep = "-")
        sub2$Var <- paste(sub2$Pos1, sub2$Pos2, sep = "-")
        sub2 <- sub2[match(chosen_pairs_var, sub2$Var),]
        
      } else if (length(sub2$Pos1)<length(sub1$Pos1)) {
        s <- sort(as.numeric(sub1$Dist), decreasing = F, index.return = T)
        sub1 <- sub1[s$ix,]
        sub1 <- sub1[1:length(sub2$Pos1),]
      }
      
      mat <- rbind(sub1, sub2[,-which(colnames(sub2)=="Var")])
      save(mat, file=sprintf("Paper_ROC_Data_%s/%s_Thr_%s.RData", alg, id, thr))
      
      if (length(sub1$Pos1)>=10 && length(sub2$Pos1)>=10){
        print(id)
        y <- c(matrix(1, 1, length(sub1$Pos1)), matrix(-1, 1, length(sub2$Pos1)))
        
        pos <- length(which(y==1))
        neg <- length(which(y==-1))
        
        ord <- c(which(colnames(mat)=="PHACE"), which(colnames(mat)==alg))
        for (i in ord){
          nms <- colnames(mat)[i]
          val <- as.numeric(mat[,i])
          pp2 <- roc(val, as.factor(1 * (y == 1)))
          pr2 <- pr.curve(scores.class0 = val, weights.class0 = (1 * (y == 1)), curve = TRUE)
          vect <- rbind(vect, c(nms, id, thr, close_val, pos, neg, auc(pp2), pr2$auc.integral))
        }
        
        colnames(vect) <- c("Method", "ID", "THR", "CLOSE", "POS", "NEG", "AUC", "AUPR")
        save(vect, file = sprintf("Paper_ROC_Measures_%s/%s_ROC.RData", alg, id))
      }
    }
  }
}
