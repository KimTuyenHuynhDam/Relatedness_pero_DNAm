# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(seqinr)

# Load the dataset

nested_anova_data <- read_excel("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx")

# Filter the dataset
nested_anova_filtered <- nested_anova_data %>%
  filter(p_value < 0.05, FDR < 0.1) %>%
  filter(!is.na(SYMBOL))

# Load the chromosome sizes data
chromosome_sizes_df <- read_tsv("./data/sequence_report.tsv")

# Filter chromosome sizes and prepare for plotting
chromosome_sizes_df <- chromosome_sizes_df %>%
  filter(`Chromosome name` %in% c(as.character(1:22), "X")) %>%
  dplyr::select(`Chromosome name`, `Seq length`) %>%
  rename(seqnames = `Chromosome name`, Size = `Seq length`) %>%
  mutate(seqnames = factor(seqnames, levels = c(as.character(1:22), "X")))

# Merge chromosome sizes with nested_anova_filtered
combined_data <- nested_anova_filtered %>%
  left_join(chromosome_sizes_df, by = "seqnames")

# Ensure chromosome order: 1 to 22, X
chromosome_levels <- c(as.character(1:22), "X")

# Filter combined_data to keep only valid chromosomes
combined_data <- combined_data %>%
  filter(seqnames %in% chromosome_levels) %>%
  mutate(seqnames = factor(seqnames, levels = chromosome_levels)) # Set order

# Update chromosome sizes to match valid chromosomes
chromosome_sizes_df <- chromosome_sizes_df %>%
  filter(seqnames %in% chromosome_levels) %>%
  mutate(seqnames = factor(seqnames, levels = chromosome_levels))

# Custom color for the dataset
custom_color <- "blue"

# Plot
ggplot() +
  # Add chromosome lines to show lengths
  geom_segment(data = chromosome_sizes_df,
               aes(x = 0, xend = Size, y = seqnames, yend = seqnames),
               color = "black", size = 0.7) +
  
  # Add CGid positions with jitter
  geom_jitter(data = combined_data,
              aes(x = CGstart, y = seqnames),
              color = custom_color,
              height = 0.3, size = 1.5, alpha = 0.7) +
  
  # Add space between y-axis levels
  scale_y_discrete(limits = rev(levels(combined_data$seqnames))) + # Reverse levels, maintain order
  
  labs(
    title = "Position of top CpGs on Chromosomes",
    x = "Position (bp)",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(), # Remove gridlines
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save the plot with increased vertical size
ggsave("./plots/mapping_top_CpGs_BW_dataset.png", width = 12, height = 10)
