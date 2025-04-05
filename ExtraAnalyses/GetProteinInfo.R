# GetInfo_AllProt.R
# This script extracts phylogenetic information from ancestral reconstruction log files. 
# It creates a summary .RData file containing phylogenetic information

# Process phylogenetic data from the log files
files <- list.files("PHACE_ProteinList")
data <- c()
kk <- 1

# Process files from index 1071 to 5024
for (id in files[1071:5024]){
  print(c(id, kk))
  kk <- kk + 1
  
  # Read and parse RAxML log file
  info <- readLines(sprintf("/PHACE_ProteinList/%s/%s.log", id, id))
  info <- unlist(info)
  
  # Extract tree length information
  len <- info[grep("Total tree length: ", info)]
  len <- unlist(strsplit(len, "Total tree length: "))
  len <- setdiff(len, "")
  
  # Extract alignment statistics
  alignment_info <- info[grep("Alignment has ", info)]
  matches <- regmatches(alignment_info, gregexpr("[0-9]+", alignment_info))
  numbers <- as.numeric(unlist(matches))
  
  num_seq <- numbers[1]  # Number of sequences
  num_pos <- numbers[2]  # Number of positions
  
  # Combine results
  data <- rbind(data, c(id, num_seq, num_pos, len))
}

# Convert to data frame and save results
data <- as.data.frame(data)
colnames(data) <- c("ID", "NumSeq", "NumPos", "TreeLength")
save(data, file = "Info_Revision1_AllProt.RData")










