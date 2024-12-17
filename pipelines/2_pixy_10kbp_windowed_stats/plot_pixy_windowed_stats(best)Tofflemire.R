install.packages("tidyverse")
install.packages("ggplot2")
library(tidyverse)

pixy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_pi.txt", header=TRUE, sep="\t")

pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)


pixy_df %>%
  filter(chromosome == 1) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic)) +
  geom_line(size = 0.25) +
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
  xlab("Position on Chromosome (Mb)") +
  ylab("Statistic Value") +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")







# Load libraries
library(tidyverse)

# Load your pixy output file
pixy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_pi.txt", header=TRUE, sep="\t")

# Plot nucleotide diversity (pi) for each population along chromosome 1
pixy_df %>%
  filter(chromosome == "chr1") %>%  # Adjust chromosome filter as needed
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_pi, color = pop)) +
  geom_line(size = 0.3) +
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(pi)) +  # Use expression to display the pi symbol
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")




# Load libraries
library(tidyverse)

# Load your pixy dxy output file
dxy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_dxy.txt", header=TRUE, sep="\t")

# Create a unique identifier for each population pair
dxy_df <- dxy_df %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "-"))

# Plot dxy for each population pair along chromosome 1
dxy_df %>%
  filter(chromosome == "chr3") %>%  # Adjust to the chromosome of interest
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_dxy, color = pop_pair)) +
  geom_line(size = 0.5) +
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(D[XY])) +  # Use expression to display Dxy symbol
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")




















# Load necessary libraries
library(tidyverse)

# Load your data (assuming they’re in the same structure as before)
pixy_pi <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_pi.txt", header=TRUE, sep="\t")
pixy_dxy <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_dxy.txt", header=TRUE, sep="\t")

# Combine pi and dxy data into a single data frame for plotting
pixy_df <- bind_rows(
  pixy_pi %>% mutate(statistic = "avg_pi", value = avg_pi),
  pixy_dxy %>% mutate(statistic = "avg_dxy", value = avg_dxy)
)

# Custom labeller for special characters in pi and dxy
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]"),
                             default = label_parsed)

# Plot summary statistics (pi and dxy) across all chromosomes
pixy_df %>%
  # Alternate coloring for chromosomes and order them properly
  mutate(chrom_color_group = case_when(
    as.numeric(chromosome) %% 2 != 0 ~ "even",
    chromosome == "X" ~ "even",
    chromosome == "Y" ~ "even",
    TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
  # Create plot
  ggplot(aes(x = (window_pos_1 + window_pos_2) / 2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller, value = label_value)) +
  xlab("Chromosome") +
  ylab("Statistic Value") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))










# Load necessary libraries
library(tidyverse)

# Load nucleotide diversity (pi) data
pixy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_pi.txt", header=TRUE, sep="\t")

# Plot nucleotide diversity (pi) for each population across all chromosomes
pixy_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_pi, color = pop)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 8) +  # Facet by chromosome, adjust ncol as needed
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(pi)) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")





# Load genetic divergence (dxy) data
dxy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_dxy.txt", header=TRUE, sep="\t")

# Create a unique identifier for each population pair
dxy_df <- dxy_df %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "-"))


# Plot dxy for each population pair across all chromosomes
dxy_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_dxy, color = pop_pair)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 4) +  # Facet by chromosome, adjust ncol as needed
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(D[XY])) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")








# Load necessary libraries
library(tidyverse)

# Load genetic divergence (dxy) data
dxy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_dxy.txt", header=TRUE, sep="\t")

# Filter out comparisons involving the outgroup
dxy_df <- dxy_df %>%
  filter(pop1 != "out" & pop2 != "out") %>%  # Exclude rows where either pop1 or pop2 is 'out'
  mutate(pop_pair = paste(pop1, pop2, sep = "-"))  # Create unique identifier for each population pair

# Plot dxy for each population pair across all chromosomes on a single row
dxy_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_dxy, color = pop_pair)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", nrow = 1) +  # Facet by chromosome in a single row
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(D[XY])) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")






# Create a custom color palette with more colors
custom_colors <- c("#377eb8", "green2", "red2")
                            
# Apply custom colors to the plot
dxy_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_dxy, color = pop_pair)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", nrow = 1) +
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(D[XY])) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = custom_colors)


# Load necessary libraries
library(tidyverse)

# Load nucleotide diversity (pi) data
pixy_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_pi.txt", header = TRUE, sep = "\t")

# Filter out the group labeled "out"
pixy_df <- pixy_df %>%
  filter(pop != "out")  # Exclude rows where pop is "out"

# Define a custom color palette
custom_colors <- c("#377eb8", "orangered", "purple4")  # Adjust as needed for the populations

# Plot nucleotide diversity (pi) for each population across all chromosomes
pixy_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000) %>%
  ggplot(aes(x = chr_position, y = avg_pi, color = pop)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 8) +  # Facet by chromosome
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(pi)) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +  # Remove legend for consistency
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = custom_colors)  # Apply the custom color palette




# Load necessary libraries
library(tidyverse)
library(patchwork)  # For combining plots

# Load and filter Tajima's D data
file_path_tajimas <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/output/ECN_individuals_tajimasD.Tajima.D"
tajimas_d_data <- read.table(file_path_tajimas, header = TRUE)

# Filter for the first 8 chromosomes for Tajima's D
selected_chromosomes_tajimas <- unique(tajimas_d_data$CHROM)[1:8]
subset_data_tajimas <- tajimas_d_data %>%
  filter(CHROM %in% selected_chromosomes_tajimas) %>%
  arrange(CHROM, BIN_START) %>%
  mutate(chr_position = (BIN_START + 50000) / 1000000)  # Midpoint in Mb assuming 100kb bins

# Plot Tajima's D with facet wrapping by chromosome and a single color
tajimas_d_plot <- ggplot(subset_data_tajimas, aes(x = chr_position, y = TajimaD, group = CHROM)) +
  geom_line(size = 0.4, color = "red") +  # Set single color for all lines
  facet_wrap(~ CHROM, scales = "free_x", nrow = 1) +  # Facet by chromosome in a single row
  labs(x = "Position on Chromosome (Mb)", y = "Tajima's D") +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(fill = "grey90", color = "black"),  # Box around each facet
        strip.text = element_text(face = "bold"),  # Bold facet labels
        panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
        panel.grid = element_blank(),  # Remove grid lines for clarity
        legend.position = "none")



# Display the combined plot
print(tajimas_d_plot)

# Save the combined plot as an image file
ggsave("combined_plot.png", combined_plot, width = 12, height = 10)












# Load necessary libraries
library(tidyverse)
library(patchwork)  # For combining plots

# Define file paths for each population's Tajima's D data
file_paths <- list(
  ECN = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/output/ECN_individuals_tajimasD.Tajima.D",
  WCN = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/output/WCN_individuals_tajimasD.Tajima.D",
  EUR_ = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/output/EUR_NAM_SAM_individuals_tajimasD.Tajima.D"
)

# Load and combine Tajima's D data for each population
tajimas_d_data <- map_dfr(file_paths, read.table, header = TRUE, .id = "Population")

# Filter for the first 8 chromosomes for Tajima's D
selected_chromosomes_tajimas <- unique(tajimas_d_data$CHROM)[1:8]
subset_data_tajimas <- tajimas_d_data %>%
  filter(CHROM %in% selected_chromosomes_tajimas) %>%
  arrange(Population, CHROM, BIN_START) %>%
  mutate(chr_position = (BIN_START + 50000) / 1000000)  # Midpoint in Mb assuming 100kb bins

# Plot Tajima's D with facet wrapping by chromosome and color by population
tajimas_d_plot <- ggplot(subset_data_tajimas, aes(x = chr_position, y = TajimaD, color = Population, group = interaction(Population, CHROM))) +
  geom_line(size = 0.5) +
  facet_wrap(~ CHROM, scales = "free_x", nrow = 1) +  # Facet by chromosome in a single row
  labs(x = "Position on Chromosome (Mb)", y = "Tajima's D") +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(fill = "grey90", color = "black"),  # Box around each facet
        strip.text = element_text(face = "bold"),  # Bold facet labels
        panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
        panel.grid = element_blank(),
        legend.position = "right") +  # Remove legend
  scale_color_manual(values = c("#377eb8", "orangered", "purple4"))  # Customize colors for populations

# Display the plot
print(tajimas_d_plot)

# Display the combined plot
print(tajimas_d_plot)






# Load necessary libraries
library(tidyverse)

# Load FST data
fst_df <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/data/pixy_fst.txt", 
                     header = TRUE, sep = "\t")

# Filter out comparisons involving "out" in pop2
fst_df <- fst_df %>%
  filter(pop2 != "out")  # Adjust as needed to exclude unwanted groups

# Define a population pair column for plotting
fst_df <- fst_df %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "-"))  # Combine pop1 and pop2 into one column

# Calculate midpoint of genomic windows for plotting
fst_df <- fst_df %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000)  # Convert to Mb

# Define a custom color palette for population pairs
fst_colors <- c("#377eb8", "green2", "red2","black")  # Adjust as needed

# Plot FST for each population pair across all chromosomes
fst_df %>%
  ggplot(aes(x = chr_position, y = avg_wc_fst, color = pop_pair)) +
  geom_line(size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 8) +  # Facet by chromosome
  xlab("Position on Chromosome (Mb)") +
  ylab(expression(F[ST])) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = fst_colors)  # Apply the custom color palette
