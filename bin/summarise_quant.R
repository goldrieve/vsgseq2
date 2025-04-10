pacman::p_load(dplyr, stringr)

# Access the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define 'quant_dirs' from the command line arguments
# Remove the square brackets and split the argument into a list by a comma
quant_dirs <- strsplit(gsub("\\[|\\]", "", args[1]), split = ",")[[1]]

threshold <- as.numeric(args[2])

# Trim leading and trailing whitespaces from the file paths
quant_dirs <- trimws(quant_dirs)

# Append '/quant.sf' to the directory paths to get the file paths
my.files <- file.path(quant_dirs, "quant.sf")

# Read each 'quant.sf' file into a data frame and store them in a list
tpm <- lapply(my.files, function(file) {
  # Read the file
  read.table(file, header = TRUE)
})

# Extract TPM
tpm_extract <- function(DF) {
  DF <- DF %>% select(Name, TPM)
  return(DF)
}

# Apply the 'tpm_extract' function to each data frame in the list
tpm <- lapply(tpm, tpm_extract)

# Extract the last element from each file path and assign these names to the list of data frames
names <- sapply(strsplit(quant_dirs, split="\\/"), function(x) tail(x, n=1))
names(tpm) <- names

# Join all data frames in the list into a single wide data frame
wide_tpm <- dplyr:::reduce(tpm, full_join, by = "Name")
names <- append(names, 'vsg_id', after = 0)
names <- gsub("_quant", "", names)

# Remove '1.fq' from column names if present
names <- gsub("_1\\.fq", "", names)

colnames(wide_tpm) <- names

# Write the wide data frame to a CSV file
write.csv(wide_tpm, './tpm.csv', row.names=FALSE)

#Reads


# Read each 'quant.sf' file into a data frame and store them in a list
reads <- lapply(my.files, function(file) {
  # Read the file
  read.table(file, header = TRUE)
})

# Extract TPM
reads_extract <- function(DF) {
  DF <- DF %>% select(Name, NumReads)
  return(DF)
}

# Apply the 'tpm_extract' function to each data frame in the list
reads <- lapply(reads, reads_extract)

# Extract the last element from each file path and assign these names to the list of data frames
names <- sapply(strsplit(quant_dirs, split="\\/"), function(x) tail(x, n=1))
names(reads) <- names

# Join all data frames in the list into a single wide data frame
wide_reads <- dplyr:::reduce(reads, full_join, by = "Name")
names <- append(names, 'vsg_id', after = 0)
names <- gsub("_quant", "", names)

# Remove '1.fq' from column names if present
names <- gsub("_1\\.fq", "", names)

colnames(wide_reads) <- names

# Write the wide data frame to a CSV file
write.csv(wide_reads, './num_reads.csv', row.names=FALSE)

# Exclude the 'vsg_id' column and calculate the total sum of read counts for each sample
total_read_counts <- colSums(wide_reads[ , -1], na.rm = TRUE)
total_read_counts <- round(total_read_counts)

# Convert the result to a data frame for better readability
total_read_counts_df <- data.frame(Sample = names(total_read_counts), TotalReads = total_read_counts)

# Print a warning message for samples with total read counts below defined value
low_read_samples <- total_read_counts_df %>% filter(TotalReads < threshold)

# Write the result to a CSV file
write.csv(total_read_counts_df, './total_read_counts.csv', row.names = FALSE)

# Get the list of low-read samples
low_read_sample_names <- low_read_samples$Sample

# Filter out the low-read samples from the tpm_data
filtered_tpm_data <- wide_tpm %>%
  select(-all_of(low_read_sample_names))  # Remove columns corresponding to low-read samples

# Write the filtered data to a new CSV file
write.csv(filtered_tpm_data, './filtered_tpm.csv', row.names = FALSE)