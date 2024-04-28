library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(mltools)
library("irr")

args = commandArgs(trailingOnly=TRUE)

id <- args[1]
print(id)

load(sprintf("Part1_AC/%s_MSA1.RData", id))
main <- data
rm(data)
load(sprintf("Part1_Gap/%s_MSA2.RData", id))
gap <- data
rm(data)

mat_of_all_changes_main <- main$mat_of_all_changes
mat_of_all_changes_main <- mat_of_all_changes_main[-1,]
mat_of_aa_changes_main <- main$mat_of_aa_changes
mat_of_aa_changes_main <- mat_of_aa_changes_main[-1,]
mat_of_aff_leaves_main <- main$mat_of_aff_leaves
mat_of_aff_leaves_main <- mat_of_aff_leaves_main[-1,]
mat_of_aff_branches_main <- main$mat_of_aff_branches
mat_of_aff_branches_main <- mat_of_aff_branches_main[-1,]

mat_of_all_changes_gap <- gap$mat_of_all_changes
mat_of_all_changes_gap <- mat_of_all_changes_gap[-1,]
mat_of_aa_changes_gap <- gap$mat_of_aa_changes
mat_of_aa_changes_gap <- mat_of_aa_changes_gap[-1,]
mat_of_aff_leaves_gap <- gap$mat_of_aff_leaves
mat_of_aff_leaves_gap <- mat_of_aff_leaves_gap[-1,]
mat_of_aff_branches_gap <- gap$mat_of_aff_branches
mat_of_aff_branches_gap <- mat_of_aff_branches_gap[-1,]

total_pos <- length(mat_of_aa_changes_main[1,])-2
print(total_pos)

mat_of_changes_combined <- matrix(0, length(mat_of_aa_changes_main[,1]), total_pos)
mat_of_aa_changes_combined <- matrix(0, length(mat_of_aa_changes_main[,1]), total_pos)
mat_of_num_aff_branches_combined <- matrix(0, length(mat_of_aa_changes_main[,1]), total_pos)
mat_of_aff_branches_combined <- matrix(0, length(mat_of_aa_changes_main[,1]), total_pos)
for (i1 in 1:total_pos) {

    print(c(i1, total_pos))
    dif1_main <- as.numeric(mat_of_all_changes_main[,(i1+1)])
    dif1_gap <- as.numeric(mat_of_all_changes_gap[,(i1+1)])

    dif_aa1_main <- mat_of_aa_changes_main[,(i1+1)]
    num_eff1_main <- as.numeric(mat_of_aff_leaves_main[,(i1+1)])
    pos1_aff_br_main <- mat_of_aff_branches_main[,(i1+1)]
  
    branch_len <- as.numeric(mat_of_aa_changes_main[,length(mat_of_aa_changes_main[1,])])
    branch_info <- mat_of_aa_changes_main[,1]
    
    dif_aa1_gap <- mat_of_aa_changes_gap[,(i1+1)]
    num_eff1_gap <- as.numeric(mat_of_aff_leaves_gap[,(i1+1)])
    pos1_aff_br_gap <- mat_of_aff_branches_gap[,(i1+1)]
    
    data_main <- cbind(dif1_main, dif_aa1_main, num_eff1_main, branch_len, pos1_aff_br_main, branch_info)

    
    data_gap <- cbind(dif1_gap, dif_aa1_gap, num_eff1_gap, branch_len, pos1_aff_br_gap, branch_info)
  
    data_main <- as.data.frame(data_main)
    data_gap <- as.data.frame(data_gap)
    
    gap_plc_i1 <- c(which(data_main$dif_aa1_main =="C-"), which(data_main$dif_aa1_main=="A-"))
    data_main$dif_aa1_main[gap_plc_i1] <- "--"
    
    loc_ind_gap1 <- data_gap$branch_info[which(data_gap$dif_aa1_gap=="CG")]
    
    data_main$dif1_main[match(loc_ind_gap1, data_main$branch_info)] <- data_gap$dif1_gap[which(data_gap$dif_aa1_gap=="CG")]
    data_main$num_eff1_main[match(loc_ind_gap1, data_main$branch_info)] <- data_gap$num_eff1_gap[which(data_gap$dif_aa1_gap=="CG")]
    data_main$pos1_aff_br_main[match(loc_ind_gap1, data_main$branch_info)] <- data_gap$pos1_aff_br_gap[which(data_gap$dif_aa1_gap=="CG")]
    data_main$dif_aa1_main[match(loc_ind_gap1, data_main$branch_info)] <- "X-"
    
    mat_of_changes_combined[,i1] <- data_main$dif1_main
    mat_of_aa_changes_combined[,i1] <- data_main$dif_aa1_main
    mat_of_num_aff_branches_combined[,i1] <- data_main$num_eff1_main
    mat_of_aff_branches_combined[,i1] <- data_main$pos1_aff_br_main
    if (i1==1){
      mat_of_info <- cbind(data_main$branch_len, data_main$branch_info)
    }
}

mat_of_changes <- list()
mat_of_changes$total_change <- mat_of_changes_combined
mat_of_changes$aa_change <-  mat_of_aa_changes_combined
mat_of_changes$aff_branch_num <- mat_of_num_aff_branches_combined
mat_of_changes$aff_branches <- mat_of_aff_branches_combined
mat_of_changes$info <- mat_of_info
save(mat_of_changes, file=sprintf("%s_TotalChange.RData", id))
