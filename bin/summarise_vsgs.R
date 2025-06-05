pacman::p_load(dplyr, stringr)

args <- commandArgs(trailingOnly = TRUE)

fasta_files <- strsplit(gsub("\\[|\\]", "", args[1]), split = ",")[[1]]
fasta_files <- trimws(fasta_files)

# Initialize a vector to hold the counts
counts <- integer(length(fasta_files))

# Loop over the files
for (i in seq_along(fasta_files)) {
  # Read the file
  file_content <- readLines(fasta_files[i])

  # Count the occurrences of ">"
  counts[i] <- sum(str_count(file_content, ">"))
}

# Create a data frame with the results
vsg_count <- data.frame(sample = fasta_files, VSGs = counts)

# Simplify the file column to remove the path and the .fasta part
vsg_count$sample <- sub("\\_VSGs.fasta$", "", basename(vsg_count$sample))
vsg_count$sample <- sub("\\_ORF$", "", vsg_count$sample)

# Print the results
write.csv(vsg_count, './vsg_count.csv', row.names=FALSE)