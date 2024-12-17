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

# =======================================
# Plot publication-quality manhattan plot
# =======================================


# Create the Manhattan plot manually
manhattan_data <- data.frame(
  SNP = seq_along(pcadapt_results$pvalues),
  PValue = pcadapt_results$pvalues
)

# Plot the Manhattan plot using base R
plot(
  manhattan_data$SNP, 
  -log10(manhattan_data$PValue), 
  xlab = "Loci", 
  ylab = expression(-log[10](italic(p))),
  pch = 20, 
  col = "black",
  main = "Manhattan Plot with Highlighted Outliers"
)

# Add highlighted points for outliers
outliers <- which(padj < alpha)  # Outlier indices
points(
  outliers, 
  -log10(manhattan_data$PValue[outliers]), 
  col = "red", 
  pch = 19, 
  cex = 1
)



# ==========================================================
# Extract top 10 outlier loci from each Chromosome (first 8)
# ==========================================================

# Load libraries
library(dplyr)
library(purrr)

# File paths
file_path <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/pcadapt_bonferroni_outliers.csv"
output_dir <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results"

# Read the CSV file
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Add a new column for -log10(pvalue)
data <- data %>%
  mutate(neg_log_pvalue = -log10(pvalue))

# Extract the first 8 chromosomes and save as one combined file
top10_chr1_8 <- data %>%
  filter(Chromosome %in% 1:8) %>%              # Filter for Chromosomes 1 to 8
  group_by(Chromosome) %>%
  arrange(desc(neg_log_pvalue), .by_group = TRUE) %>%
  slice_head(n = 10) %>%                       # Top 10 per chromosome
  ungroup()

# Save the combined file for chromosomes 1-8
write.table(top10_chr1_8,
            file = file.path(output_dir, "top10_combined_chr1-8.csv"),
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

# Extract and save the top 10 loci for each remaining chromosome as separate files
data %>%
  filter(!Chromosome %in% 1:8) %>%             # Exclude chromosomes 1 to 8
  group_by(Chromosome) %>%
  arrange(desc(neg_log_pvalue), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  group_split() %>%                            # Split by chromosome
  walk2(
    group_keys(data %>% filter(!Chromosome %in% 1:8) %>% group_by(Chromosome))$Chromosome,
    ~ write.table(.x,
                  file = file.path(output_dir, paste0("top10_chr", .y, ".csv")),
                  sep = ",",
                  row.names = FALSE,
                  quote = FALSE)
  )

# Message for completion
cat("Top 10 loci for chromosomes 1-8 saved as 'top10_combined_chr1-8.csv'.\n")
cat("Top 10 loci for other chromosomes saved as individual files.\n")



# ==================================================================
# Plot manhattan plot with highlighting top 10 outlier loci (base R)
# ==================================================================


# Re-import the file with the correct delimiter
top10_outliers_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/top10_combined_chr1-8.csv"
top10_outliers <- read.table(top10_outliers_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Verify the structure
print(head(top10_outliers))
str(top10_outliers)

# Prepare Manhattan plot data
manhattan_data <- data.frame(
  SNP = 1:nrow(bim_data),                 # SNP indices as x-axis
  Position = bim_data$BP,                 # Positions from bim_data
  Chromosome = bim_data$CHR,              # Chromosome from bim_data
  PValue = pcadapt_results$pvalues        # P-values from pcadapt results
)

# Identify indices of outliers in the Manhattan data
outlier_indices <- which(
  manhattan_data$Position %in% top10_outliers$Position & 
    manhattan_data$Chromosome %in% top10_outliers$Chromosome
)

# Plot the Manhattan plot
plot(
  manhattan_data$SNP, 
  -log10(manhattan_data$PValue), 
  xlab = "Loci", 
  ylab = expression(-log[10](italic(p))),
  pch = 1, 
  col = "black",
  main = "Manhattan Plot with Highlighted Outliers"
)

# Highlight the outlier loci
if (length(outlier_indices) > 0) {
  points(
    outlier_indices, 
    -log10(manhattan_data$PValue[outlier_indices]), 
    col = "red", 
    pch = 1, 
    cex = 1
  )
  cat("Number of highlighted outliers:", length(outlier_indices), "\n")
} else {
  cat("No matching outliers found to highlight.\n")
}





# ========================================
# Save Manhattan Plot as High-Quality PNG
# ========================================

# Define the output file path
output_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/pcadapt_manhattan_plot_highres.png"

# Open a PNG device with high resolution
png(output_file, width = 5, height = 5, units = "in", res = 300)  # 300 dpi resolution

# Manhattan Plot
plot(
  manhattan_data$SNP, 
  -log10(manhattan_data$PValue), 
  xlab = "Loci", 
  ylab = expression(-log[10](italic(p))),
  pch = 1, 
  col = "black",
  main = "Manhattan Plot with Highlighted Outliers"
)

# Highlight the outlier loci
if (length(outlier_indices) > 0) {
  points(
    outlier_indices, 
    -log10(manhattan_data$PValue[outlier_indices]), 
    col = "red", 
    pch = 1, 
    cex = 1
  )
}

# Close the PNG device to save the file
dev.off()

cat("Plot saved as:", output_file, "\n")




# ==============================================
# Manhattan Plot with Centered Chromosome Labels
# ==============================================

# Define output file
output_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/pcadapt_manhattan_plot_highres_centered_labels.png"

# Open a PNG device for high-quality output
png(output_file, width = 5, height = 5, units = "in", res = 300)

# Prepare Manhattan data
manhattan_data <- data.frame(
  SNP = 1:nrow(bim_data),                 # SNP indices as x-axis
  Position = bim_data$BP,                 # Positions from bim_data
  Chromosome = as.factor(bim_data$CHR),   # Chromosomes as factors
  PValue = pcadapt_results$pvalues        # P-values from pcadapt results
)

# Assign alternating colors for chromosomes
chrom_colors <- rep(c("gray50", "gray20"), length.out = length(levels(manhattan_data$Chromosome)))
manhattan_data$color <- chrom_colors[as.numeric(manhattan_data$Chromosome)]

# Calculate chromosome midpoints for labels
chrom_midpoints <- tapply(manhattan_data$SNP, manhattan_data$Chromosome, function(x) (min(x) + max(x)) / 2)

# Plot Manhattan Plot with alternating colors
plot(
  manhattan_data$SNP, 
  -log10(manhattan_data$PValue),
  col = manhattan_data$color,    # Assign colors based on chromosome
  pch = 1,
  xlab = "Chromosomes",
  ylab = expression(-log[10](italic(p))),
  main = "Manhattan Plot with Centered Chromosome Labels",
  xaxt = "n"  # Remove default x-axis ticks
)

# Add chromosome labels centered under each chromosome
axis(1, at = chrom_midpoints, labels = levels(manhattan_data$Chromosome), las = 2)

# Highlight the outlier loci
if (length(outlier_indices) > 0) {
  points(
    outlier_indices, 
    -log10(manhattan_data$PValue[outlier_indices]), 
    col = "red", 
    pch = 1, 
    cex = 1
  )
}

# Close the PNG device
dev.off()

cat("Plot saved as:", output_file, "\n")









# ==============================
# End of Analysis
# ==============================





# ==============================
# End of Analysis
# ==============================