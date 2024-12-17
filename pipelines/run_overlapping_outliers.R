# Load required library
library(dplyr)

# Define file paths
file1 <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/out/pcadapt_bonferroni_outliers.csv"
file2 <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/out/outflank_outliers_info.csv"

# Read the files
outliers1 <- read.csv(file1)
outliers2 <- read.csv(file2)

# Ensure column names are consistent
colnames(outliers1) <- tolower(colnames(outliers1))
colnames(outliers2) <- tolower(colnames(outliers2))

# Select relevant columns for comparison
outliers1_subset <- outliers1 %>% select(locus, chromosome, position)
outliers2_subset <- outliers2 %>% select(locus)

# Find overlapping loci
overlapping_loci <- inner_join(outliers2_subset, outliers1_subset, by = "locus")

# Count the number of overlapping loci
num_overlaps <- nrow(overlapping_loci)

# Print results
cat("Number of overlapping loci:", num_overlaps, "\n")

# Save the overlapping loci with Chromosome and Position to a file
write.csv(overlapping_loci, 
          file = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/out/overlapping_outliers_with_positions.csv",
          row.names = FALSE)

cat("Overlapping loci with Chromosome and Position saved to 'overlapping_outliers_with_positions.csv'.\n")

