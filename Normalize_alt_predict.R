##

# GitHub: aagnolin (https://github.com/aagnolin/PauseMeter)

# Author:  Alberto Agnolin (alberto.agnolin.1@gmail.com)
# Name of script: Normalize_alt_predict
# Summary: this script normalizes ribosome profiling counts in alt_predict output files by multiplying the counts of each mapped ribosome profiling peak by the ratio 
# between the sum of counts of the file with the highest number of reads and the sum of counts of the one that is going to be normalized. In addition, it filters out 
# non coding RNAs that might still be present after mapping.

# NOTE: this script is meant to be used in conjunction with alt_predict (https://github.com/BiosystemsDataAnalysis/PausePredictionTools) and produces files that can
# be used with scripts PauseMeter.R and RiboScout.R.
# Read the README.md file for more information

##

# Load necessary libraries
library(dplyr)
library(stringr)

# Function to display help message
print_help <- function() {
  cat("Usage: Rscript normalize_riboseq.R <input_directory> <output_directory>\n\n")
  cat("Description:\n")
  cat("  This script normalizes all CSV files containing riboseq data in the input directory.\n")
  cat("  The normalized datasets are saved in the output directory.\n\n")
  cat("Arguments:\n")
  cat("  <input_directory>   The path to the directory containing input CSV files.\n")
  cat("  <output_directory>  The path to the directory where normalized CSV files will be saved.\n")
}

# Check if help option is provided
if (commandArgs(TRUE)[1] == "-h") {
  print_help()
  quit(save = "no", status = 0)
}

# Command-line arguments for input and output directories
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 2) {
  cat("Error: Incorrect number of arguments.\n\n")
  print_help()
  quit(save = "no", status = 1)
}

input_dir <- args[1]
output_dir <- args[2]

# Get a list of CSV files in the input directory
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Check if there are CSV files in the input directory
if (length(csv_files) == 0) {
  cat("Error: No CSV files found in the input directory.\n\n")
  print_help()
  quit(save = "no", status = 1)
}

# Load datasets, filter out rows with locus_tag starting with "BSU_", and calculate sum of counts for each
sum_counts <- numeric(length(csv_files))
datasets <- list()
for (i in seq_along(csv_files)) {
  dataset <- read.csv(csv_files[i])
  filtered_dataset <- dataset %>% filter(!grepl("^BSU_", locus_tag))
  sum_counts[i] <- sum(filtered_dataset$count)
  datasets[[i]] <- filtered_dataset
}

# Find the index of the dataset with the highest sum of counts
max_index <- which.max(sum_counts)

# Normalize datasets using the dataset with the highest sum of counts
normalized_datasets <- lapply(seq_along(datasets), function(i) {
  dataset <- datasets[[i]]
  sum_max <- sum(datasets[[max_index]]$count)
  Normalized_dataset <- mutate(dataset, Norm_count = dataset$count * (sum_counts[max_index] / sum(dataset$count)))
  # Get the filename of the input dataset
  input_filename <- basename(csv_files[i])
  
  # Construct the output filename with 'Normalized_' prefix
  output_filename <- paste0("Normalized_", input_filename)
  
  output_file <- file.path(output_dir, output_filename)
  write.csv(Normalized_dataset, output_file, row.names = FALSE)
  return(Normalized_dataset)
})