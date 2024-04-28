library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)

args = commandArgs(trailingOnly=TRUE)
id <- args[1]

res <- read.csv(sprintf("ToleranceScores/%s.csv", id))
eps <- 10^(-15)
res[,2:21] <- 1-log(res[,2:21] +eps)/log(eps)

file_fasta <- sprintf("%s_MaskedMSA.fasta", id)

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

for (i in 1:length(msa[1,])){
  sc <- res[i, 2:21]
  obs <- table(msa[which(msa[,i]!="-"),i])
  dom <- names(sort(obs, decreasing = T))[1]
  
  lc <- which(colnames(sc)==dom)
  dom_p2 <- names(sc)[which(as.numeric(sc)>as.numeric(sc[lc]))]
  
  dom_loc <- c()
  for (aa in c(dom, dom_p2)){
    dom_loc <- c(dom_loc, which(msa[,i]==aa))
  }
 
  msa[dom_loc,i] <- "XX"
  
  upd <- intersect(which(msa[,i]!="XX"), which(msa[,i]!="-"))
  msa[upd,i] <- "A"
  msa[dom_loc, i] <- "C"
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

write.table(vect, quote = F, col.names = F, row.names = F, sprintf("MSA1/%s_MSA1.fasta", id))
