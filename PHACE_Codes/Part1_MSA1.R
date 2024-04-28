library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
source("./coev_diff_MSA1.R")
args = commandArgs(trailingOnly=TRUE)

id <- args[1]

aa_to_num <- function(aa) {
  amino_acids <- c("A", "C")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("A", "C")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

file_fasta <- sprintf("%s_MSA1.fasta", id)
file_nwk <- sprintf("%s_MSA1.treefile", id, id)
file_rst <- sprintf("%s_MSA1.state", id, id)

output_name <- id

tr_org <- read.tree(file_nwk)
x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
x[,1] <- str_remove(x[,1], "Node")
colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )
states <- x[, c(1:3)]

x <- x[,c(1:5)]

tree_info <- as.data.frame(as_tibble(tr_org))

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

print(c(id, dim(msa)))

fasta_node <- c()

# Names of leaves
names_all <- tr_org[["tip.label"]]
msa <- msa[names_all, ]
# Number of total leaves&nodes
num_leaves <- length(tr_org[["tip.label"]])
num_nodes <- tr_org[["Nnode"]]

nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]

# Total number of positions from ancestralProbs file
total_pos <- max(x$Site)

fasta_node <- matrix(0, num_nodes, total_pos)
for (ndd in 1:num_nodes){
  k1 <- states[which(states[,1]==ndd),3]
  fasta_node[ndd, ] <- k1
}
fasta_node2 <- matrix(0, num_nodes, total_pos)
for (i in 1:num_nodes){
  fasta_node2[i, ] <- fasta_node[nodes_raxml[i],]
}
fasta_node <- fasta_node2

num_branch <- length(tree_info$branch.length)-1
num_nodes_prev <- num_nodes
num_leaves_prev <- num_leaves

positions <- 1:total_pos

####################################################
names_msa <- row.names(msa)
trim_final <- matrix(0, num_leaves, total_pos)
trim_weight <- matrix(0, num_leaves, total_pos)
el <- c()
for (i in 1:total_pos){
  vec <- names_msa[which(msa[, i] == "-")]
  if (length(vec)>0){
    vec <- unique(vec)
    trim_final[1:length(vec),i] <- vec
    trim_weight[1:length(vec),i] <- vec
  }
  if (length(unique(msa[,i]))==1){
    el <- c(el, i)
  }
}
trim_msa <- trim_final

dd_node <- dist.nodes(tr_org)
root_dist <- dd_node[(num_leaves+1),]
root_dist <- root_dist[-(num_leaves+1)] 

root_node <- tree_info$label[num_leaves+1]


chosen_nodes2 <- tree_info[,c("parent", "node")]
chosen_nodes2 <- chosen_nodes2[-(num_leaves+1), ]
label_branch <- tree_info$label[-(num_leaves+1)]

conn_left <- chosen_nodes2[,2] #tree_new_info$parent
conn_right <- chosen_nodes2[,1] #tree_new_info$node

cc <- which(tree_info$parent==tree_info$node)
root_nodee <- tree_info$parent[cc]

levels <- matrix(0,num_branch,num_branch)
ii_all <- as.numeric(tree_info$node[(num_leaves+1)])
i<-1
mx<-0
while ((sum(ii_all>num_leaves)>=1)) {
  ndd <- ii_all
  if (sum(is.element(1:num_leaves, ndd))>=1){
    ndd <- ndd[-which(ndd<=num_leaves)]
  }
  mtc <- c()
  ii_all <- c()
  sapply(1:length(ndd), function(k) {
    ii <<- which(conn_right==ndd[k])
    jj <<- conn_left[ii]
    ii_all <<- c(ii_all, jj)
    sapply(jj, function(uu) {
      l1 <<- intersect(which(chosen_nodes2[,2]==ndd[k]),which(chosen_nodes2[,1]==uu))
      l2 <<- intersect(which(chosen_nodes2[,1]==ndd[k]),which(chosen_nodes2[,2]==uu))
      if (length(l2)>0){
        ff <<- paste(ndd[k], ",", uu, sep = "")
      } else{
        ff <<- paste(uu, ",", ndd[k], sep = "")
      }
      mtc <<- c(mtc, ff)
    })
    mx <<- max(mx, length(mtc))
  })
  levels[1:length(mtc),i] <- mtc
  i <- i+1
}
levels <- levels[,-c(i:num_branch)]
levels <- levels[-c((mx+1):num_branch),]
add <- matrix(0, length(levels[,1]),1)
add[1] <- c(sprintf("0, %s", root_nodee))
levels <- cbind(add, levels)




all_diffs <- c()
no_sig <- 0
mat_of_all_changes <- matrix(100, num_branch, total_pos)
mat_of_all_changes <- cbind(tree_info$label[-(num_leaves+1)], mat_of_all_changes)
mat_of_all_changes <- cbind(mat_of_all_changes, tree_info$branch.length[-(num_leaves+1)])
mat_of_aa_changes <- matrix(100, num_branch, total_pos)
mat_of_aa_changes <- cbind(tree_info$label[-(num_leaves+1)], mat_of_aa_changes)
mat_of_aa_changes <- cbind(mat_of_aa_changes, tree_info$branch.length[-(num_leaves+1)])
mat_of_aff_leaves <- matrix(100, num_branch, total_pos)
mat_of_aff_leaves <- cbind(tree_info$label[-(num_leaves+1)], mat_of_aff_leaves)
mat_of_aff_branches <- matrix(100, num_branch, total_pos)
mat_of_aff_branches <- cbind(tree_info$label[-(num_leaves+1)], mat_of_aff_leaves)
mat_of_all_changes_extra <- matrix(100, 1, total_pos)
mat_of_aa_changes_extra <- matrix(100, 1, total_pos)
mat_of_aff_leaves_extra <- matrix(100, 1, total_pos)
mat_of_aff_branches_extra <- matrix(100, 1, total_pos)
for (ps in 1:(total_pos)) {
    print(ps)
    position_vec1 <- msa[, ps]
    gap1 <- which(position_vec1=="-")
       
    trim_from_msa <- c() #unique(setdiff(trim_msa[,ps], "0"))
    trims <- c() #trim_from_msa
        
    if (length(trims)==num_leaves || length(trims)>=(num_leaves-3)){
          tree_new <- c()
          msa_upd <- c()
          no_sig <- no_sig + 1
          
    } else {
      if (length(trims)>0){
        tree_new <- drop.tip(tr_org, trims)
        tree_new_info <- as.data.frame(as_tibble(tree_new))
      } else {
          tree_new <- tr_org
          tree_new_info <- as.data.frame(as_tibble(tr_org))
      }
          
      msa_upd <- msa
      if (length(trims)>=1){
        for (y in 1:length(trims)){
          msa_upd <- msa_upd[-which(rownames(msa_upd)==trims[y]),]
          }
       }
    
      ord <- tree_info$label[num_leaves+1]
      num_nodes_prev <- num_nodes
      num_leaves_prev <- num_leaves
      nodes_raxml_prev <- nodes_raxml
        
      diff_all <- compute_difference(ps, ord, levels, tree_new, tree_new_info, num_nodes, num_leaves, x, total_pos, nodes_raxml, msa_upd, fasta_node)
      all_diffs[[sprintf("Pos_%s", ps)]] <- diff_all
      
      if (length(trims)==0){
        num_leaves_x <- as.numeric(diff_all$num_leaves)
        els <- setdiff(tree_info$label, tree_new_info$label)
        els2 <- match(els, mat_of_all_changes[,1])
        val <- diff_all$diff1
        val <- val[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]

        val2 <- diff_all$mat2
        val2 <- val2[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        val3 <- diff_all$eff_br
        val3 <- val3[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        val4 <- diff_all$aff_br
        val4 <- val4[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        mat_of_all_changes[, (ps+1)] <- val
        mat_of_aa_changes[, (ps+1)] <- val2
        mat_of_aff_leaves[, (ps+1)] <- val3
        mat_of_aff_branches[, (ps+1)] <- val4
      } else {
        num_leaves_x <- as.numeric(diff_all$num_leaves)
        els <- c(setdiff(tree_info$label, c(tree_new_info$label, tree_info$label[num_leaves+1])), tree_new_info$label[(num_leaves_x+1)])
        
        if (is.element(root_node, els)){
          els <- setdiff(els, root_node)
        }
        
        els2 <- match(els, mat_of_all_changes[,1])
        val <- diff_all$diff1
        val <- val[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        val2 <- diff_all$mat2
        val2 <- val2[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        val3 <- diff_all$eff_br
        val3 <- val3[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        val4 <- diff_all$aff_br
        val4 <- val4[match(tree_new_info$label[-(num_leaves_x+1)], intersect(tree_info$label[-(num_leaves+1)], tree_new_info$label[-(num_leaves_x+1)]))]
        mat_of_all_changes[-els2, (ps+1)] <- val
        mat_of_aa_changes[-els2, (ps+1)] <- val2
        mat_of_aff_leaves[-els2, (ps+1)] <- val3
        mat_of_aff_branches[-els2, (ps+1)] <- val4
      }
      
      mat_of_all_changes_extra[ps] <- diff_all$extra[1]
      mat_of_aa_changes_extra[ps] <- diff_all$extra[2]
      mat_of_aff_leaves_extra[ps] <- diff_all$extra[3]
      mat_of_aff_branches_extra[ps] <- diff_all$extra[4]
    }
}

mat_of_all_changes <- rbind(c("Node0", mat_of_all_changes_extra, 0), mat_of_all_changes)
mat_of_aa_changes <- rbind(c("Node0", mat_of_aa_changes_extra, 0), mat_of_aa_changes)
mat_of_aff_leaves <- rbind(c("Node0", mat_of_aff_leaves_extra), mat_of_aff_leaves)
mat_of_aff_branches <- rbind(c("Node0", mat_of_aff_branches_extra, 0), mat_of_aff_branches)

data <- list()
data$mat_of_all_changes <- mat_of_all_changes
data$mat_of_aa_changes <- mat_of_aa_changes
data$mat_of_aff_leaves <- mat_of_aff_leaves
data$mat_of_aff_branches <- mat_of_aff_branches

save(data, file=sprintf("Part1_AC/%s_MSA1.RData", id))

