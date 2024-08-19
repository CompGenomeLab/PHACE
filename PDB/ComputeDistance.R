library(bio3d)
library(spaa)
library(dplyr)
library(tidyr)
library(broom)
library(stringr)
library(caret)

args = commandArgs(trailingOnly=TRUE)

load("ID_PDB_List.RData")

ids <- best_pdb_all$ID
ids2 <- best_pdb_all$PDB_ID

for (k in 1:length(ids)){
  id <- ids[k]
  id2 <- ids2[k]
  el_pos <- best_pdb_all$DifferentPos_List[k]
  el_pos <- as.numeric(unlist(strsplit(el_pos, "-")))
  
  file_fasta <- sprintf("/cta/groups/adebali/PHACT_Dataset/%s/1_psiblast/%s.fasta", id, id)
  file_str <- sprintf("/cta/users/nkuru/COEVOLUTION_DATA/Ref_and_Others2/Part2_Gap/PDB/PDBs/%s.pdb", id2)
  
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  ii <- grep(id, rownames(msa))
  lc <- which(msa[ii,]=="G")
  
  a <- read.pdb(file_str)
  data <- a[["atom"]]
  
  if (as.numeric(best_pdb_all$LongerPosNum)[k]>0){
    data <- data[-which(as.numeric(data$resno)>as.numeric(best_pdb_all$Length_MSA)[k]),]
  }
  
  data <- data[which(data$elety=="CB"),]
  CB_atom <- data
  
  if (length(unique(CB_atom$chain))>1){
    chs <- names(which(table(CB_atom$chain)==max(table(CB_atom$chain))))[1]
    CB_atom <- CB_atom[which(CB_atom$chain==chs),]
  }
  
  el_pos <- intersect(el_pos, CB_atom$resno)
  els <- c()
  if (length(el_pos)>0){
    el_pos <- as.numeric(el_pos)
    el_pos2 <- match(el_pos, CB_atom$resno)
    els <- c(els, el_pos2)
  }
  negs <- which(as.numeric(CB_atom$resno)<0)
  els <- c(els, negs)
  els <- unique(els)
  if (length(els)>0){
    els <- as.numeric(els)
    CB_atom <- CB_atom[-els,]
  }
  
  pos <- unique(as.numeric(CB_atom$resno))
  dist <- c()
  dist2 <- c()
  for (i in pos){
    for (j in pos){
      k1 <- which(CB_atom$resno==i)
      k2 <- which(CB_atom$resno==j)
      if (k1<k2){
        d <- sqrt((CB_atom$x[k1]-CB_atom$x[k2])^2 + (CB_atom$y[k1]-CB_atom$y[k2])^2 + (CB_atom$z[k1]-CB_atom$z[k2])^2)
        dist <- rbind(dist, c(i, j, d))
        dist2 <- rbind(rbind(dist2, c(i,j,d)), c(j,i,d))
      }
    }
  }
  
  a <- read.pdb(file_str)
  data <- a[["atom"]]
  if (as.numeric(best_pdb_all$LongerPosNum)[k]>0){
    data <- data[-which(as.numeric(data$resno)>as.numeric(best_pdb_all$Length_MSA)[k]),]
  }
  
  data2 <- data[which(data$elety=="CA"),]
  CB_atom2 <- data2
  
  if (length(unique(CB_atom2$chain))>1){
    chs <- names(which(table(CB_atom2$chain)==max(table(CB_atom2$chain))))[1]
    CB_atom2 <- CB_atom2[which(CB_atom2$chain==chs),]
  }
  
  pos2 <- unique(as.numeric(CB_atom2$resno))
  dist3 <- c()
  dist4 <- c()
  for (i in pos2){
    for (j in pos2){
      k1 <- which(CB_atom2$resno==i)
      k2 <- which(CB_atom2$resno==j)
      if (k1<k2){
        d <- sqrt((CB_atom2$x[k1]-CB_atom2$x[k2])^2 + (CB_atom2$y[k1]-CB_atom2$y[k2])^2 + (CB_atom2$z[k1]-CB_atom2$z[k2])^2)
        dist3 <- rbind(dist3, c(i, j, d))
        dist4 <- rbind(rbind(dist4, c(i,j,d)), c(j,i,d))
      }
    }
  }
  
  for (l in lc){
    # print(l)
    i1 <- which(dist3[,1]==l)
    i2 <- which(dist3[,2]==l)
    
    if (length(i1)>0 && length(i2)>0){ 
      dist <- rbind(dist, dist3[c(i1, i2),])
    }
    
    i1 <- which(dist4[,1]==l)
    i2 <- which(dist4[,2]==l)
    if (length(i1)>0 && length(i2)>0){
      dist2 <- rbind(dist2, dist4[c(i1, i2),])
    }
  }
  rm(dist3, dist4)
  
  save(dist, file=sprintf("Dists3/Dist_PDB_%s.RData", id))
}
  

