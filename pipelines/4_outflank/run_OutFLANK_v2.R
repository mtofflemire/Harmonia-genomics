# Load necessary libraries
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)

# === 1. Load Metadata and VCF ===
vcf_path <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/filtered_no_outgroup_pruned.vcf.gz"
meta_path <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/Harmonia_meta_3.csv"

# Load metadata
meta <- read.csv(meta_path)
head(meta)

# Confirm outgroup is removed in metadata
if ("out" %in% meta$SITES) {
  stop("Outgroup still present in metadata. Remove it and rerun.")
}

# Load VCF data
data <- read.vcfR(vcf_path)
geno <- extract.gt(data)  # Extract genotype information
dim(geno)  # Check dimensions

# === 2. Extract Chromosome and Position Information from VCF ===
vcf_metadata <- as.data.frame(data@fix[, c("CHROM", "POS")])  # Extract CHROM and POS columns
colnames(vcf_metadata) <- c("Chromosome", "Position")
vcf_metadata$Chromosome <- as.factor(vcf_metadata$Chromosome)
vcf_metadata$Position <- as.numeric(vcf_metadata$Position)
vcf_metadata$SNP <- 1:nrow(vcf_metadata)  # Assign sequential SNP identifiers

# === 3. Process Genotype Data ===
G <- geno
G[geno %in% c("0/0")] <- 0
G[geno %in% c("0/1")] <- 1
G[geno %in% c("1/1")] <- 2
G[is.na(G)] <- 9  # Set missing data to 9
tG <- t(G)  # Transpose: SNPs as columns, individuals as rows
dim(tG)  # Check dimensions

# Ensure metadata and genotype rows match
matching_order <- match(meta$ID, rownames(tG))
if (any(is.na(matching_order))) {
  stop("Mismatch between metadata IDs and genotype row names. Check alignment.")
}
tG_reordered <- tG[matching_order, , drop = FALSE]

# Confirm alignment
if (!identical(rownames(tG_reordered), as.character(meta$ID))) {
  stop("Row names of tG_reordered and meta$ID are not identical after reordering.")
}

# === 4. Subset Metadata and Genotypes for Selected Populations ===
# Define populations to include (excluding outgroup)
subpops <- unique(meta$SITES)  # Modify this if specific populations are required
submeta <- subset(meta, SITES %in% subpops)
subgen <- tG_reordered[meta$SITES %in% subpops, , drop = FALSE]

# Verify alignment after subsetting
if (!identical(rownames(subgen), as.character(submeta$ID))) {
  stop("Row names of subgen and submeta$ID are not identical after subsetting.")
}

# === 5. Calculate FST Matrix ===
fst <- MakeDiploidFSTMat(subgen, locusNames = 1:ncol(subgen), popNames = submeta$SITES)

# Add Chromosome and Position information to FST
fst$SNP <- as.character(1:nrow(fst))  # Ensure SNP column matches `vcf_metadata$SNP`
fst <- merge(fst, vcf_metadata, by = "SNP", all.x = TRUE)

# Verify data integrity
if (any(is.na(fst$Position))) {
  stop("Some SNPs do not have chromosome or position information.")
}

# === 6. Analyze FST Distribution and Identify Outliers ===
# Quick stats from OutFLANK analysis
hist(fst$FST, breaks = 50)
summary(fst$FST)  # Highest FST should be higher than the mean

# Run OutFLANK
OF <- OutFLANK(fst, LeftTrimFraction = 0.01, RightTrimFraction = 0.01,
               Hmin = 0.05, NumberOfSamples = length(unique(submeta$SITES)), qthreshold = 0.05)

# Plot OutFLANK results
OutFLANKResultsPlotter(OF, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005,
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Identify Outlier SNPs
P1 <- pOutlierFinderChiSqNoCorr(fst, Fstbar = OF$FSTNoCorrbar,
                                dfInferred = OF$dfInferred, qthreshold = 0.05, Hmin = 0.1)


outliers <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers)


plot(P1$LocusName,P1$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1$LocusName[outliers],P1$FST[outliers],col="magenta")

