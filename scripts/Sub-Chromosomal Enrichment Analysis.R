# ==============================================================================
# Sub-Chromosomal Enrichment Analysis (Sliding Window)
# ==============================================================================

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
library(readr) # For read_tsv if needed

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------

# Load the nested ANOVA results (CpG statistics)
# Ensure this file contains: CGid, FDR, seqnames (Chromosome), CGstart (Position)
nested_anova_data <- read_excel("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx")


# Filter for valid chromosomes (1-23, X) and valid positions
valid_chromosomes <- c(as.character(1:23), "X")

clean_data <- nested_anova_data %>%
  filter(seqnames %in% valid_chromosomes) %>%
  filter(!is.na(CGstart)) %>%
  mutate(
    # Define "Significant" based on your study's threshold (e.g., FDR < 0.05)
    IsSignificant = ifelse(FDR < 0.05, TRUE, FALSE)
  )

# ------------------------------------------------------------------------------
# 2. Define Sliding Window Parameters
# ------------------------------------------------------------------------------
WINDOW_SIZE <- 2000000  # 2 Mb window size (Adjust to 1Mb=1000000 or 5Mb as needed)

# Calculate Global Background Statistics (The "Urn")
# M = Total number of significant CpGs in the entire genome (White balls)
M_total_sig <- sum(clean_data$IsSignificant)

# N = Total number of non-significant CpGs in the entire genome (Black balls)
N_total_nonsig <- nrow(clean_data) - M_total_sig

# Total CpGs
Total_CpGs <- nrow(clean_data)

print(glue("Global Background: {M_total_sig} significant CpGs out of {Total_CpGs} total CpGs."))

# ------------------------------------------------------------------------------
# 3. Perform Sliding Window Analysis (Hypergeometric Test)
# ------------------------------------------------------------------------------

window_stats <- clean_data %>%
  group_by(seqnames) %>%
  mutate(
    # Assign each CpG to a window bin
    WindowStart = floor(CGstart / WINDOW_SIZE) * WINDOW_SIZE,
    WindowEnd = WindowStart + WINDOW_SIZE
  ) %>%
  # Group by Chromosome AND Window
  group_by(seqnames, WindowStart, WindowEnd) %>%
  summarise(
    n_total_in_window = n(),              # k (Total draws)
    n_sig_in_window = sum(IsSignificant), # x (White balls drawn)
    .groups = "drop"
  ) %>%
  # Filter out sparse windows (e.g., windows with < 5 probes might act as noise)
  filter(n_total_in_window >= 5) %>%
  rowwise() %>%
  mutate(
    # Hypergeometric Test: phyper(q, m, n, k, lower.tail = FALSE)
    # q = number of significant CpGs in window - 1 (because lower.tail=FALSE is strictly >)
    # m = Total significant in genome
    # n = Total non-significant in genome
    # k = Total CpGs in this window
    p_value_enrichment = phyper(n_sig_in_window - 1, 
                                M_total_sig, 
                                N_total_nonsig, 
                                n_total_in_window, 
                                lower.tail = FALSE)
  ) %>%
  ungroup() %>%
  # Adjust for multiple testing across all windows
  mutate(
    FDR_Window = p.adjust(p_value_enrichment, method = "BH")
  )

# ------------------------------------------------------------------------------
# 4. Filter and Export Significant Windows
# ------------------------------------------------------------------------------

significant_windows <- window_stats %>%
  filter(FDR_Window < 0.05) %>%
  arrange(FDR_Window)

# Print summary
print(glue("Number of significant enriched windows found (FDR < 0.05): {nrow(significant_windows)}"))

if(nrow(significant_windows) > 0){
  print(head(significant_windows))
  # Save to CSV
  write.csv(significant_windows, "./tables/Significant_SubChromosomal_Clusters.csv", row.names = FALSE)
} else {
  print("No sub-chromosomal regions showed significant enrichment after FDR correction.")
}

# ------------------------------------------------------------------------------
# 5. Visualization (Enrichment Manhattan Plot)
# ------------------------------------------------------------------------------

# Create a cumulative position for continuous X-axis plotting
# (Simple plotting by facet is often clearer for windows)

plot_data <- window_stats %>%
  mutate(
    logP = -log10(p_value_enrichment),
    seqnames = factor(seqnames, levels = valid_chromosomes)
  )

enrichment_plot <- ggplot(plot_data, aes(x = WindowStart, y = logP, color = seqnames)) +
  geom_point(size = 1.5, alpha = 0.7) +
  
  # Threshold line for significance (approximate Bonferroni or FDR cutoff)
  # Here we just show a reference line, e.g., p < 0.001
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "grey50") +
  
  # Facet by chromosome to see the "Sub-chromosomal" structure
  facet_grid(~seqnames, scales = "free_x", space = "free_x") +
  
  scale_color_viridis_d(guide = "none") + # Or use custom colors
  
  labs(
    title = "Sub-chromosomal Enrichment of Relatedness-Associated CpGs",
    subtitle = paste0("Sliding Window Size: ", WINDOW_SIZE/1e6, " Mb. Y-axis: -log10(Enrichment P-value)"),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Hide messy x-axis labels in facet
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(size = 8, angle = 90) # Vertical chr labels
  )

# Save Plot
ggsave("./plots/SubChromosomal_Enrichment_Plot.png", enrichment_plot, width = 14, height = 6)

print(enrichment_plot)