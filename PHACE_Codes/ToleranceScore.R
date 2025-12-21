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
  amino_acids <- c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")
  aa <- toupper(aa)
  res <- match(aa, amino_acids)
  res[is.na(res)] <- 21
  return(res)
}

num_to_aa <- function(num) {
  amino_acids <- c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")
  aa <- ifelse(num == 21, "X", amino_acids[num])
  return(aa)
}

id <- args[1]
file_fasta <- args[2]
file_nwk <- args[3]
file_rst <- args[4]

if (!all(file.exists(file_fasta, file_nwk, file_rst))) {
  stop("One or more input files do not exist. Check your paths.")
}

output_name <- id

pos_chosen <- "all"

# Read tree file, NWK
tr_org <- read.tree(file_nwk)
tree_info <- as.data.frame(as_tibble(tr_org))

# Read state file
x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
x[,1] <- str_remove(x[,1], "Node")
colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )  

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

# ----------------------------
# Input validation checkpoints
# ----------------------------

# 1. Check that tree_info$label is not empty or all NA
if (all(is.na(tree_info$label)) || length(tree_info$label) == 0) {
  stop("tree_info$label is empty. Check that your .nwk file includes node labels (Node1, Node2, etc.).")
}

# 2. Check that the ASR probability table has expected shape
expected_cols <- 23  # 1: NodeID, 2: Site, 3: something like Posterior?, then 20 amino acids
if (ncol(x) != expected_cols) {
  stop("Unexpected number of columns in ASR file. Expected ", expected_cols,
       " but got ", ncol(x), ". Check your .state file format.")
}

# 3. Check that MSA dimensions match the tree
if (nrow(msa) != length(tr_org$tip.label)) {
  stop("MSA and tree tip label counts do not match: nrow(msa) = ", nrow(msa),
       ", but tree has ", length(tr_org$tip.label), " leaves.")
}

# 4. Optional informative messages for users
message("Tree info labels: ", length(tree_info$label), " entries")
message("ASR table size: ", paste(dim(x), collapse=" x "))
message("MSA size: ", paste(dim(msa), collapse=" x "))


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

#tolerance_scores <- matrix(unlist(scores), nrow = length(positions), ncol = 20, byrow = TRUE)
tolerance_scores <- cbind(positions, scores)
colnames(tolerance_scores) <- c("Pos/AA", num_to_aa(1:20))

write.csv(tolerance_scores, quote = F, row.names = F, paste("ToleranceScores/", output_name, ".csv", sep = ""))


