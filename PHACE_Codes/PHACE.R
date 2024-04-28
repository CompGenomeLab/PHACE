library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(mltools)
library("irr")

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}


wccc <- function(Y, X, W) {
  cov_w <- sum(W * (Y - weighted.mean(Y, W)) * (X - weighted.mean(X, W)))
  sd_Y_w <- sqrt(sum(W * (Y - weighted.mean(Y, W))^2))
  sd_X_w <- sqrt(sum(W * (X - weighted.mean(X, W))^2))
  
  wcc <- 2 * cov_w / (sd_Y_w^2 + sd_X_w^2 + (weighted.mean(Y, W) - weighted.mean(X, W))^2)
  return(wcc)
}

file_fasta <- sprintf("%s_MSA1.fasta", id)
file_fasta_org <- sprintf("%s_MaskedMSA.fasta", id)
file_nwk <- sprintf("%s_MSA1.treefile", id)

tr_org <- read.tree(file_nwk)
tree_info <- as.data.frame(as_tibble(tr_org))

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

fasta_org <- read.fasta(file = file_fasta_org)
msa_org <- fasta_org$ali

# Names of leaves
names_all <- tr_org[["tip.label"]]
msa <- msa[names_all, ]
# Number of total leaves&nodes
num_leaves <- length(tr_org[["tip.label"]])
num_nodes <- tr_org[["Nnode"]]

# Total number of positions from ancestralProbs file
total_pos <- length(msa[1,])

num_branch <- length(tree_info$branch.length)-1

load(sprintf("%s_TotalChange.RData", id))
mat_total_change <- mat_of_changes$total_change
mat_aa_change <- mat_of_changes$aa_change
mat_aff_branch_num <- mat_of_changes$aff_branch_num
mat_aff_branches <- mat_of_changes$aff_branches
mat_info <- mat_of_changes$info

mat_of_changes_combined_num <- mat_total_change
mat_of_changes_combined_num <-  matrix(as.numeric(mat_of_changes_combined_num), length(mat_of_changes_combined_num[,1]), length(mat_of_changes_combined_num[1,]))
row_means <- rowMeans(mat_of_changes_combined_num)
general_mean <- mean(mat_of_changes_combined_num)
row_means[which(row_means<general_mean)] <- general_mean
row_means <- row_means/general_mean
weight_branch <- 1/row_means

mat <- c() 
for (i1 in 1:(total_pos-1)) {
  for (i2 in (total_pos+1):total_pos) {
    weight_branch <- 1/row_means    
    com_gap <- length(unique(c(names(msa[which(msa[,i1]=="-"),i1]), names(msa[which(msa[,i2]=="-"),i2]))))
    cons1 <- length(unique(msa_org[,i1]))
    cons2 <- length(unique(msa_org[,i2]))
    
    data <- mat_info
    data <- as.data.frame(data)
    data$dif1 <- mat_total_change[,i1]
    data$dif2 <- mat_total_change[,i2]
    data$num_eff1 <- mat_aff_branch_num[,i1]
    data$num_eff2 <- mat_aff_branch_num[,i2]
    data$pos1_aff_br <- mat_aff_branches[,i1]
    data$pos2_aff_br <- mat_aff_branches[,i2]
    data$dif_aa1 <- mat_aa_change[,i1]
    data$dif_aa2 <- mat_aa_change[,i2]
    data <- data[,-c(1:2)]
    data$node <- mat_info[,2]
    data$branch_len <-  mat_info[,1]
    
    w1 <- which(data$dif_aa1=="X-")
    w2 <- which(data$dif_aa2=="X-")
    
    br1 <- setdiff(unlist(strsplit(data$pos1_aff_br[w1],"-")), "0")
    br2 <- setdiff(unlist(strsplit(data$pos2_aff_br[w2],"-")), "0")
    ww <- length(unique(c(br1, br2)))
    com_ww <- length(intersect(br1, br2))
    
    weight_inc1 <- max(1 - (ww-com_ww)/num_branch,0)

    dif1 <- data$dif1
    dif2 <- data$dif2
    
    if (min(cons1, cons2) == 1){
      score3 <- -1
    } else {
      print(c(i1, i2))
      part2 <- data
      
      dim_part2 <- length(part2$dif1)
      branch_lens <- as.numeric(part2$branch_len)
      subs <- sort(as.numeric(data$branch_len), decreasing = T)
      
      data <- as.data.frame(data)
      
      pos1_totalscore <- as.numeric(data$dif1)
      pos2_totalscore <- as.numeric(data$dif2)
      
      #####################################################
      ###### Score over total change ######################
      #####################################################
      
      dif_tot <- as.numeric(pos1_totalscore)-as.numeric(pos2_totalscore)
      loc1 <- which(dif_tot>=0.5)
      loc2 <- which(dif_tot<=-0.5)
      
      if (length(loc1)>0){
        part2_loc1 <- part2[loc1, ]
      }
      if (length(loc2)>0){
        part2_loc2 <- part2[loc2, ]
      }
      
      pos1_all_br <- setdiff(unlist(strsplit(part2$pos1_aff_br, "-")), "0")
      pos2_all_br <- setdiff(unlist(strsplit(part2$pos2_aff_br, "-")), "0")
      
      if (length(loc1)==1){
        lx <- loc1
        lx_orig <- loc1
        aff_lf <- as.numeric(part2_loc1$num_eff1)
        aff_br <- part2_loc1$pos1_aff_br
        aff_br <- setdiff(unlist(strsplit(aff_br, "-")), "0")
        if (aff_lf==1){
          alt <- msa_org[which(row.names(msa_org)==aff_br),i1]
          obs_org <- length(which(msa_org[,i1]==alt))
          if (obs_org==1){
            upgr <- upgr + 1
            pos1_totalscore[lx_orig] <- 0
            pos2_totalscore[lx_orig] <- 0
          }
        }
        
      } else if (length(loc1)>1){
        weight_ignore <- length(loc1)
        for (lx in 1:length(loc1)){
          lx_orig <- loc1[lx]
          aa_obs <- part2_loc1$dif_aa1[lx]
          sec_aa <- substr(aa_obs, 2, 2)
          obs <- length(which(msa[,i1]==sec_aa))
          aff_lf <- as.numeric(part2_loc1$num_eff1[lx])
          aff_br <- part2_loc1$pos1_aff_br[lx]
          aff_br <- setdiff(unlist(strsplit(aff_br, "-")), "0")
          if (aff_lf==1){
            alt <- msa_org[which(row.names(msa_org)==aff_br),i1]
            obs_org <- length(which(msa_org[,i1]==alt))
            if (obs_org==1){
              upgr <- upgr + 1
              pos1_totalscore[lx_orig] <- 0
              pos2_totalscore[lx_orig] <- 0
            }
          }
          
        }
      }
      
      
      if (length(loc2)==1){
        lx <- loc2
        lx_orig <- loc2
        aa_obs <- part2_loc2$dif_aa2
        sec_aa <- substr(aa_obs, 2, 2)
        obs <- length(which(msa[,i2]==sec_aa))
        aff_lf <- as.numeric(part2_loc2$num_eff2)
        aff_br <- part2_loc2$pos2_aff_br
        aff_br <- setdiff(unlist(strsplit(aff_br, "-")), "0")
        if (aff_lf==1){
          alt <- msa_org[which(row.names(msa_org)==aff_br),i2]
          obs_org <- length(which(msa_org[,i2]==alt))
          if (obs_org==1){
            upgr <- upgr + 1
            pos1_totalscore[lx_orig] <- 0
            pos2_totalscore[lx_orig] <- 0
          }
        }
        
      } else if (length(loc2)>1){
        weight_ignore <- length(loc2)
        for (lx in 1:length(loc2)){
          lx_orig <- loc2[lx]
          aa_obs <- part2_loc2$dif_aa2[lx]
          sec_aa <- substr(aa_obs, 2, 2)
          obs <- length(which(msa[,i2]==sec_aa))
          aff_lf <- as.numeric(part2_loc2$num_eff2[lx])
          aff_br <- part2_loc2$pos2_aff_br[lx]
          aff_br <- setdiff(unlist(strsplit(aff_br, "-")), "0")
          if (aff_lf==1){
            alt <- msa_org[which(row.names(msa_org)==aff_br),i2]
            obs_org <- length(which(msa_org[,i2]==alt))
            if (obs_org==1){
              upgr <- upgr + 1
              pos1_totalscore[lx_orig] <- 0
              pos2_totalscore[lx_orig] <- 0
            }
          }
          
        }
      }
      
      data$dif1 <- pos1_totalscore
      data$dif2 <- pos2_totalscore
      
      inc_gap1 <- intersect(intersect(which(data$dif_aa1=="X-"), which(data$dif_aa2!="--")),
                            which(data$dif_aa2!="X-"))
      inc_gap2 <- intersect(intersect(which(data$dif_aa2=="X-"), which(data$dif_aa1!="--")),
                            which(data$dif_aa1!="X-"))
      
      if (length(inc_gap1)>0){
        data[inc_gap1, 1] <- 0
      }
      if (length(inc_gap2)>0){
        data[inc_gap2, 2] <- 0
      }
      
      pos1_totalscore <- as.numeric(data$dif1)
      pos2_totalscore <- as.numeric(data$dif2)
      
      #####################################################
      ###### Score Computation and Assignment #############
      #####################################################
      
      if (sum(pos1_totalscore)>=1 && sum(pos2_totalscore)>=1 && max(pos1_totalscore)>=0.5 && max(pos2_totalscore)>=0.5){
        df <- abs(pos1_totalscore-pos2_totalscore)
        plc <- which(df>=0.5)
        weight_branch[plc] <- 1

        cc <- cbind(pos1_totalscore, pos2_totalscore)
        cc <- apply(cc, 1, max)
        cc2 <- cc
        cc[which(cc==0)] <- 1
        weight_branch2 <- sqrt(weight_branch*cc)
        score3 <- wccc(pos1_totalscore, pos2_totalscore, weight_branch2)*weight_inc1

      } else {
        score3 <- -1
      }
    }
    
    mat <- rbind(mat, c(i1, i2, score3))
  }
}

mat <- as.data.frame(mat)
colnames(mat) <- c("Pos1", "Pos2", "PHACE_Score")
save(mat, file= sprintf("%s_PHACE.RData",id))

