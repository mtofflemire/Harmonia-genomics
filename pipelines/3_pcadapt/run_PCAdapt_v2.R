# ==========================================================
# PCAdapt Analysis Pipeline for Harmonia Genomics SNP Dataset
# ==========================================================

# This script installs and loads required libraries, 
# processes genotype data, performs PCAdapt analysis,
# and generates a Manhattan plot highlighting Bonferroni outliers.

# ------------------------------
# 1. Install and Load Libraries
# ------------------------------
if (!requireNamespace("pcadapt", quietly = TRUE)) {
  install.packages("pcadapt")
}
if (!requireNamespace("qvalue", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("qvalue")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Load libraries
library(pcadapt)
library(qvalue)
library(ggplot2)
library(dplyr)

# ----------------------------
# 2. Define File Paths
# ----------------------------
# Input files
bed_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/filtered_no_outgroup_pruned.bed"
bim_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/filtered_no_outgroup_pruned.bim"

# Output files
bonferroni_output_plot <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/out/pcadapt_manhattan_plot_bonferroni.png"
bonferroni_outliers_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/out/pcadapt_bonferroni_outliers.txt"

# ----------------------------
# 3. Load Genotype Data
# ----------------------------
# Read genotype data in PLINK BED format
genotype_data <- read.pcadapt(bed_file, type = "bed")

# ---------------------------------
# 4. Perform PCAdapt Analysis
# ---------------------------------
# Run PCAdapt with 2 principal components (K = 2)
pcadapt_results <- pcadapt(genotype_data, K = 2)

# ---------------------------------
# 5. Visualize PCAdapt Results
# ---------------------------------
# Generate standard plots to inspect results
plot(pcadapt_results, option = "scores")      # PC scores
plot(pcadapt_results, option = "manhattan")  # Manhattan plot
plot(pcadapt_results, option = "qqplot")     # Q-Q plot
hist(pcadapt_results$pvalues,                # P-value distribution
     xlab = "p-values", 
     main = NULL, 
     breaks = 50, 
     col = "orange")
plot(pcadapt_results, option = "stat.distribution")  # Statistic distribution

# ---------------------------------
# 6. Bonferroni Correction
# ---------------------------------
# Adjust p-values using the Bonferroni method
padj <- p.adjust(pcadapt_results$pvalues, method = "bonferroni")
alpha <- 0.01  # Significance level
outliers <- which(padj < alpha)  # Identify significant outliers
cat("Number of Bonferroni outlier loci:", length(outliers), "\n")


# ----------------------------
# 7. Add Chromosome and SNP Information to Outliers
# This may be redundaant. Code right after pretty much does same thing. 
# ----------------------------
# Read the .bim file
bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim_data) <- c("CHR", "SNP", "CM", "BP", "REF", "ALT")

# Subset the .bim data for outlier loci
outlier_info <- bim_data[outliers, c("CHR", "SNP", "BP")]

# Save outliers with chromosome and SNP info to a file
write.table(outlier_info, bonferroni_outliers_file, 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Print out the first few outliers for verification
print(head(outlier_info))


# ==============================
# End of Analysis
# ==============================