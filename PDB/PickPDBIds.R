load("INFO_PDBStr.RData")

all <- as.data.frame(all)
colnames(all) <- c("ID", "PDB_ID", "Length_Mapped", "Length_MSA", "DifferentPos_Length",
                   "DifferentPos_List", "LongerPosNum")

ratio <- as.numeric(all$Length_Mapped)/as.numeric(all$Length_MSA) # the ratio between mapped positions to the total sequence length
all <- all[which(ratio>=0.5),]
all <- all[-which(as.numeric(all$DifferentPos_Length)>=10),] # at most 10 different positison btw structure sequence and our query sequence

ids <- unique(all$ID)
best_pdb_all <- c()
for (id in ids){
  sub <- all[which(all$ID==id),]
  
  sub_df <- as.data.frame(sub)
  sub_df$Length_Mapped <- as.numeric(as.character(sub_df$Length_Mapped))
  sub_df$Length_MSA <- as.numeric(as.character(sub_df$Length_MSA))
  sub_df$DifferentPos_Length <- as.numeric(as.character(sub_df$DifferentPos_Length))
  sub_df$LongerPosNum <- as.numeric(as.character(sub_df$LongerPosNum))
  
  best_match <- sub_df[which(sub_df$Length_Mapped == max(sub_df$Length_Mapped)), ]
  best_match <- best_match[which(best_match$DifferentPos_Length == min(best_match$DifferentPos_Length)), ]
  if (nrow(best_match) > 1) {
    best_match <- best_match[which(best_match$LongerPosNum == min(best_match$LongerPosNum)), ]
  }
  
  best_pdb <- best_match[1,]
  best_pdb_all <- rbind(best_pdb_all, best_pdb)
}

save(best_pdb_all, file = "ID_PDB_List.RData")



