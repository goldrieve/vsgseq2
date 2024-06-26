pacman::p_load(dplyr, stringr)

# Access the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define 'quant_dirs' from the command line arguments
# Remove the square brackets and split the argument into a list by a comma
quant_dirs <- strsplit(gsub("\\[|\\]", "", args[1]), split = ",")[[1]]

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
colnames(wide_tpm) <- names

# Write the wide data frame to a CSV file
write.csv(wide_tpm, './tpm.csv', row.names=FALSE)
