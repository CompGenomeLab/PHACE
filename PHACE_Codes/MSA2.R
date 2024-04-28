library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)

args = commandArgs(trailingOnly=TRUE)

id <- args[1]

file_fasta <- sprintf("%s_MaskedMSA.fasta", id)

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

for (i in 1:length(msa[1,])){
  gap <- which(msa[,i]=="-")
  oth <- which(msa[,i]!="-")
  
  msa[gap,i] <- "G"
  msa[oth,i] <- "C"
}

num_leaves <- length(msa[,1])
vect <- matrix(0, (num_leaves*2), 1)
k <- 0
for (i in seq(1, (num_leaves*2), 2)){
  k <- k + 1
  vect[i, 1] <- paste(">", row.names(msa)[k], sep = "")
}
k <- 0
for (i in seq(2, (num_leaves*2), 2)){
  k <- k + 1
  vect[i, 1] <- paste(msa[k,], collapse = "")
}

write.table(vect, quote = F, col.names = F, row.names = F, sprintf("MSA2/%s_MSA2.fasta", id))
