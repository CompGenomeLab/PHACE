library(bio3d)
library(spaa)
library(dplyr)
library(tidyr)
library(broom)
library(stringr)
library(AUC)
library(PRROC)
library(RColorBrewer)
library(caret)
library(mltools)

args = commandArgs(trailingOnly=TRUE)

convert_aa <- function(three_letter_codes) {
  three_letter_codes <- toupper(three_letter_codes)
  aa_mapping <- c("ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D",
                  "CYS" = "C", "GLU" = "E", "GLN" = "Q", "GLY" = "G",
                  "HIS" = "H", "ILE" = "I", "LEU" = "L", "LYS" = "K",
                  "MET" = "M", "PHE" = "F", "PRO" = "P", "SER" = "S",
                  "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V")
  
  one_letter_codes <- aa_mapping[three_letter_codes]
  
  invalid_codes <- three_letter_codes[is.na(one_letter_codes)]
  if (length(invalid_codes) > 0) {
    stop(paste("Invalid three-letter amino acid codes:", paste(invalid_codes, collapse = ", ")))
  }
  
  return(one_letter_codes)
}

load("Map_Id_PDB.RData")

files <- list.files("PDBs")
#alls <- unlist(strsplit(alls, ".pdb"))

#xx <- list.files("INFO3")
#xx <- unlist(strsplit(xx, "_Info.RData"))

ind <- as.numeric(args[1])

alls <- info_id_pdb$Id
#alls <- setdiff(alls, xx)

print(length(alls))

id <- alls[ind]

id2 <- info_id_pdb$PDB_Id[which(info_id_pdb$Id==id)]
id2 <- unlist(strsplit(id2, ";"))

file_fasta <- sprintf("%s.fasta", id)

fasta <- read.fasta(file = file_fasta)

msa <- fasta$ali

vect <- c()
for (ii in id2){
 print(is.element(sprintf("%s.pdb",ii), files))
 print(sprintf("%s.pdb",ii))
 if (is.element(sprintf("%s.pdb",ii), files)){
  file_str <- sprintf("PDBs/%s.pdb", ii)
  a <- read.pdb(file_str)
  data <- a[["atom"]]
  data <- data[which(data$elety=="CB"),]
  
  allowed_amino_acids <- c("ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D",
                           "CYS" = "C", "GLU" = "E", "GLN" = "Q", "GLY" = "G",
                           "HIS" = "H", "ILE" = "I", "LEU" = "L", "LYS" = "K",
                           "MET" = "M", "PHE" = "F", "PRO" = "P", "SER" = "S",
                           "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V")
  
  # Filter out unwanted amino acids
  data <- data[data$resid %in% names(allowed_amino_acids), ]
  
  if (length(which(data$resno<0))>=1){
    data <- data[-which(data$resno<0),]
  }
  
  CB_atom <- data
  
  el <- which(as.numeric(CB_atom$resno)>length(msa))
  if (length(el)>0){
    CB_atom <- CB_atom[-el,]
    
  }
  
  pos <- unique(as.numeric(CB_atom$resno))
  map <- convert_aa(CB_atom$resid[match(pos, CB_atom$resno)])
  
  len_main <- length(msa)
  
  if (max(pos)<=length(msa[1,])){
    df <- length(map) - sum(fasta$ali[pos]==map)
    
    w <- pos[which(fasta$ali[pos]!=map)]
    
    vect <- rbind(vect, c(id, ii, length(map), length(msa[1,]), df, paste(w, collapse = "-"), length(el)))
  } else {
    int <- intersect(pos, 1:length(msa))
    vect <- rbind(vect, c(id, ii, length(map), length(msa[1,]), "X", "LongerSeq", "X"))
  }
 }  
}

save(vect, file = sprintf("InfoFolder/%s_Info.RData", id))

