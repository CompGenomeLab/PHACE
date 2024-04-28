#!/usr/bin/env Rscript
library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(Peptides)
source("./position_score.R")

args = commandArgs(trailingOnly=TRUE)

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

id <- args[1]

file_fasta <- sprintf("%s_MaskedMSA.fasta", id)
file_nwk <- sprintf("%s.treefile", id, id)
file_rst <- sprintf("%s.state", id, id)

output_name <- id

pos_chosen <- "all"

tr_org <- read.tree(file_nwk)
x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
x[,1] <- str_remove(x[,1], "Node")
colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )  

tree_info <- as.data.frame(as_tibble(tr_org))

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

# connections_1: Parent node, connections_2: connected node/leaf
connections_1 <- tree_info$parent
connections_2 <- tree_info$node

# Names of leaves
names_all <- tr_org[["tip.label"]]
msa <- msa[names_all, ]
# Number of total leaves&nodes
num_leaves <- length(tr_org[["tip.label"]])
num_nodes <- tr_org[["Nnode"]]

nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
nodes_raxml_prev <- nodes_raxml
num_nodes_prev <- num_nodes
num_leaves_prev <- num_leaves

# Total number of positions from ancestralProbs file
total_pos <- max(x$Site)
positions <- 1:total_pos
score_all <- matrix(0, total_pos, 21)

####################################################
####################################################
chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
leaf_names <- tree_info$label

chosen_nodes2 <- chosen_nodes

scores <- t(mapply(function(ps){position_score(ps, x, msa, trim_final, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes, nodes_raxml, num_leaves, total_pos, nodes_raxml)}, rep(positions)))

tolerance_scores <- matrix(unlist(scores), nrow = length(positions), ncol = 20, byrow = TRUE)
tolerance_scores <- cbind(positions, tolerance_scores)
colnames(tolerance_scores) <- c("Pos/AA", num_to_aa(1:20))

write.csv(tolerance_scores, quote = F, row.names = F, paste("ToleranceScores/", output_name, ".csv", sep = ""))

position_score <- function(ps, x, msa, trim_final, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes_prev, nodes_raxml_prev, num_leaves_prev, total_pos, nodes_raxml) {
  position <- ps
  # print(ps)
  
  b1 <- position + total_pos*(0:(num_nodes-1))
  TT <- x[b1,]
  
  node_info <- as.numeric(TT[,1])
  sort_node_info <- sort(node_info, decreasing = F, index.return=T)
  TT <- TT[sort_node_info$ix,]
  
  matrix_prob <- matrix(0, num_nodes, 20)
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  matrix_prob <- matrix_prob[nodes_raxml,]
  colnames(matrix_prob) <- names(sort(rr))
  
  el <- which(matrix_prob<0.01)
  matrix_prob[el] <- 0   
  
  msa_upd <- msa
  names_msa <- rownames(msa)
  
  tree_new <- tr_org
  tree_new_info <-  tree_info
  
  dd_node <- dist.nodes(tree_new)
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  
  root <- num_leaves + 1
  d_n <- dd_node[root, ]
  
  param <- mean(d_n)
  weights <- exp(-d_n^2/param^2)
  
  scores <- c()

  # Connections between leaves & nodes
  chosen_leaves <- tree_new_info[1:num_leaves,c("parent", "node")]
  # Connections between nodes & nodes
  chosen_nodes <- tree_new_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
  leaf_names <- tree_new_info$label
    
    # chosen_nodes2: ordered connections (for probability differences)
    chosen_nodes2 <- chosen_nodes
    
    position_vec <- msa_upd[, ps]
    
    position_num <- aa_to_num(position_vec)
    prob_leaves <- matrix(0, num_leaves, 20)
    prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
    
    gaps <- which(position_num == 21)
    
    root_node <- which(tree_new_info$parent==tree_new_info$node)-num_leaves
    
    diff_leaves <- matrix(0, num_leaves, 20)
    if (num_nodes == 1) {
      diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
      diff_nodes <- matrix(0, 1, 20)
      root_pr <- matrix_prob
    } else {
      diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
      diff_nodes <- matrix(0, num_nodes-1, 20)
      diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
      root_pr <- matrix_prob[root_node,]
    }
    

    ################## weights
    weight_leaf <- weights[1:num_leaves]
    weight_node <- weights[(num_leaves+1):length(weights)]
    
    score <- matrix(0,1,20)
    
    if (num_nodes != 1) {
      s1 <- sapply(1:20, function(ii){
        if (num_nodes == 2) {
          dif_pr <- diff_nodes[ii]
        } else {
          dif_pr <- diff_nodes[1:(num_nodes-1),ii]
        }
        dif_pr[dif_pr<0] <- 0
        sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
        score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
      })
    }
    
    ### NOVEL 29.03
    
    score_without_leaf <- score
    score_without_leaf <- score_without_leaf + root_pr
    score <- score + root_pr
    
    s2 <- sapply(1:20, function(ii){
      diff_lf <- diff_leaves[1:num_leaves,ii]
      diff_lf[gaps] <-  0
      diff_lf[diff_lf<0] <- 0
      
      
      s1 <- sum(weight_leaf * diff_lf)
      score[ii] <<- score[ii] + s1
    })
    
    scores <- (score)/(num_nodes+num_leaves)
    
    return(scores)
}
  
