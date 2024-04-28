compute_difference <- function(ps, ord, levels, tree_new, tree_new_info, num_nodes_prev, num_leaves_prev, x, total_pos, nodes_raxml_prev, msa_upd, fasta_node) {
  num_leaves <- length(tree_new[["tip.label"]])
  num_nodes <- tree_new[["Nnode"]]
  num_branch <- length(tree_new_info$parent)-1
  
  chosen_nodes2 <- tree_new_info[,c("parent", "node")]
  chosen_nodes2 <- chosen_nodes2[-(num_leaves+1), ]
  label_branch <- tree_new_info$label[-(num_leaves+1)]
  
  conn_left <- chosen_nodes2[,2] #tree_new_info$parent
  conn_right <- chosen_nodes2[,1] #tree_new_info$node
  
  cc <- which(tree_new_info$parent==tree_new_info$node)
  root_node <- tree_new_info$parent[cc]
  
  num_branch <- num_branch + 1
  
  diff <- matrix(0,num_branch,1)
  
  b1 <- ps + total_pos*(0:(num_nodes_prev-1))
  TT <- x[b1,]
  matrix_prob <- matrix(0, num_nodes, 2)
  
  node_info <- as.numeric(TT[,1])
  sort_node_info <- sort(node_info, decreasing = F, index.return=T)
  TT <- TT[sort_node_info$ix,]
  
  root <- TT[which(TT[,1]==substr(ord, 5, nchar(ord))), 4:5]
  nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_new_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
  nodes_raxml <- rbind(nodes_raxml, tree_new_info[num_leaves+1:num_nodes, "node"])
  
  nodes_raxml_org <- nodes_raxml
  
  elims <- setdiff(1:num_nodes_prev,nodes_raxml[1,])
  all_diff <- list()
  all_diff$elims <- elims
  elims2 <- c()
  if (length(elims)>0){
    elims2 <- c()
    for (i in elims){
      k <- as.numeric(names(nodes_raxml_prev)[which(nodes_raxml_prev==i)])
      elims2 <- c(elims2, (k-num_leaves_prev))
    }
    fasta_node_new <- fasta_node[-elims2,]
    for (w in 1:length(elims)){
      TT <- TT[-which(TT[,1]==elims[w]),]
    }
    nds <- TT[,1]
    for (k in 1:length(elims)){
      for (j in 1:length(nodes_raxml[1,])){
        if (nodes_raxml_org[1,j]>elims[k]){
          nodes_raxml[1,j] <- nodes_raxml[1,j]-1
        }
      }
    }
  } else {
    fasta_node_new <- fasta_node
  }
  all_diff$elims2 <- elims2
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  matrix_prob <- matrix_prob[nodes_raxml[1,],]
  
  root_prob <- matrix(0,1,2)
  root_prob[rr] <- root
  
  position_vec <- msa_upd[, ps]
  
  position_num <- aa_to_num(position_vec)
  prob_leaves <- matrix(0, num_leaves, 2)
  prob_leaves[cbind(which(position_num <= 2), position_num[which(position_num <= 2)])] <- 1
  
  all_prob <- rbind(prob_leaves, matrix_prob)
  
  chosen_nodes2 <- rbind(c(0, root_node), chosen_nodes2)
  label_branch <- c("Node0", label_branch)
  conn_left <- chosen_nodes2[,2] #tree_new_info$parent
  conn_right <- chosen_nodes2[,1] #tree_new_info$node
  
#  faa <- names(which(table(msa_upd[,ps])==max(table(msa_upd[,ps]))))
#  if (faa=="G"){
#    fa <- c(0, 1)
#  } else {
#    fa <- c(1, 0)
#  }
 
 
  faa <- names(which(table(msa_upd[,ps])==max(table(msa_upd[,ps]))))
  if (length(faa)==2){
    faa <- "C"
    fa <- c(1, 0)
  } else {
    if (faa=="G"){
      fa <- c(0, 1)
    } else {
      fa <- c(1, 0)
    }
  }
   
  nms <- num_to_aa(1:2)
  diff <- matrix(0,num_branch,1)
  mat2 <- matrix(0,num_branch,1)
  aff_brs <- matrix(0,num_branch,1)
  matrx <- matrix(0, num_branch, 5)
  diff_mat <- matrix(0, num_branch, 2)
  num_of_ch <- matrix(0, num_branch, 1)
  br_len <- tree_new_info$branch.length
  br_len <- br_len[-(num_leaves+1)]
  for (i in 1:num_branch) {
    lf <- chosen_nodes2[i,1]
    rg <- chosen_nodes2[i,2]
    
    if (lf==0){
      d_vec <- all_prob[rg, ] - fa
      first_aa <- faa
    } else {
      d_vec <- all_prob[rg, ] - all_prob[lf, ]
      first_aa <- fasta_node_new[(lf-num_leaves),ps]
    }
    
    if (rg>num_leaves){
      second_aa <- fasta_node_new[(rg-num_leaves),ps]
    } else if (rg<=num_leaves) {
      second_aa <- msa_upd[rg,ps]
    } 
    matrx[i, 1] <- first_aa
    matrx[i, 2] <- second_aa
    
    diff[i] <- sum(d_vec[d_vec>0])
    diff_mat[i, ] <- d_vec
    m1 <- first_aa
    m2 <- second_aa
    if (m1==m2){
      m1 <- "-"
      m2 <- "-"
    }
    mat2[i] <- paste(m1,m2, sep = "")
    if (m1==m2) {
      num_of_ch[i] <- 0
    } else {
      num_of_ch[i] <- 1
      aff_brs[i] <- label_branch[i]
    }
  }
  
  keep <- 0
  for (i in length(levels[1,]):1){
    xx <- levels[,i]
    if (sum(which(xx=="0"))>=1){
      xx <- xx[-which(xx=="0")]
    }
    ss <- c()
    for (j in 1:length(xx)) {
      x2 <- strsplit(xx[j],",")
      x2 <- as.numeric(unlist(x2))
      ss <- rbind(ss, x2)
    }
    for (p in 1:length(ss[,1])){
      lc <- intersect(which(conn_right==ss[p,1]),which(conn_left==ss[p,2]))
      aa_lc <- mat2[lc]
      rg <- ss[p,1]
      lf <- ss[p,2]
      if (aa_lc=="--"){
        keep <- num_of_ch[intersect(which(chosen_nodes2[,1]==rg), which(chosen_nodes2[,2]==lf))] + 1
        keep_lf <- paste(label_branch[intersect(which(chosen_nodes2[,1]==rg), which(chosen_nodes2[,2]==lf))], 
                           aff_brs[intersect(which(chosen_nodes2[,1]==rg), which(chosen_nodes2[,2]==lf))], sep="-")
        
        num_of_ch[intersect(which(chosen_nodes2[,1]==rg), which(chosen_nodes2[,2]==lf))] <- 0
        aff_brs[intersect(which(chosen_nodes2[,1]==rg), which(chosen_nodes2[,2]==lf))] <- 0
        prev <- which(chosen_nodes2[,2]==rg)
        aa_prev <- mat2[prev]
        if (length(prev)>0){
          iii<-1
          rg_prev <- chosen_nodes2[prev[iii],1]
          lf_prev <- rg
          num_of_ch[intersect(which(chosen_nodes2[,1]==rg_prev), which(chosen_nodes2[,2]==lf_prev))] <- 
            num_of_ch[intersect(which(chosen_nodes2[,1]==rg_prev), which(chosen_nodes2[,2]==lf_prev))] + keep
          aff_brs[intersect(which(chosen_nodes2[,1]==rg_prev), which(chosen_nodes2[,2]==lf_prev))] <- 
            paste(aff_brs[intersect(which(chosen_nodes2[,1]==rg_prev), which(chosen_nodes2[,2]==lf_prev))], keep_lf, sep = "-")
        }
      }
      
    }
  }

  dd_node <- dist.nodes(tree_new)
  root_dist <- dd_node[(num_leaves+1),]
  root_dist <- root_dist[-(num_leaves+1)] 

  all_diff$diff1 <- diff[2:length(diff)]
  all_diff$mat2 <- mat2[2:length(mat2)]
  all_diff$branch_len <- br_len
  all_diff$num_leaves <- num_leaves
  all_diff$num_nodes <- num_nodes
  all_diff$order <- tree_new_info$label
  all_diff$eff_br <- num_of_ch[2:length(num_of_ch)]
  all_diff$aff_br <- aff_brs[2:length(aff_brs)]
  all_diff$extra <- c(diff[1], mat2[1], num_of_ch[1], aff_brs[1])
  return(all_diff)
}

