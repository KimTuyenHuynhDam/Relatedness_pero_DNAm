# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)
library(openxlsx) # Required for writing Excel files

# ==============================================================================
# 1. DATA LOADING AND PREPARATION
# ==============================================================================

# Load the data

data <- read.xlsx("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx", sheet = "Sheet 1") %>%
  filter(seqnames %in% c(as.character(1:23), "X"))

# Calculate Global Totals (Background for Fisher's Test)
# We need these to compare each specific chromosome against the "rest of the genome"
total_CpGs_all   <- nrow(data)
total_sig_all    <- data %>% filter(FDR < 0.1) %>% nrow()
total_nonsig_all <- total_CpGs_all - total_sig_all
overall_ratio    <- total_sig_all / total_CpGs_all

# ==============================================================================
# 2. STATISTICAL ANALYSIS: FISHER'S EXACT TEST + FDR
# ==============================================================================

# Summarize data by chromosome
prepared_data <- data %>%
  distinct() %>%
  group_by(seqnames) %>% 
  summarize(
    total_CpGs = n(),
    significant_CpGs = sum(FDR < 0.1, na.rm = TRUE),
    non_significant_CpGs = total_CpGs - significant_CpGs,
    ratio = significant_CpGs / total_CpGs,
    
    # Expected counts (useful for visualization logic)
    expected_significant_CpGs = total_CpGs * overall_ratio,
    expected_non_significant_CpGs = total_CpGs - expected_significant_CpGs,
    deviation = significant_CpGs - expected_significant_CpGs,
    direction = ifelse(deviation > 0, "Higher", "Lower")
  ) %>%
  ungroup()

# Run Fisher's Exact Test row-by-row
fisher_results <- prepared_data %>%
  rowwise() %>%
  mutate(
    p_value_raw = {
      # Construct the 2x2 Contingency Table:
      #                 Significant        Non-Significant
      # In Chrom             A                   B
      # Rest of Genome       C                   D
      
      A <- significant_CpGs
      B <- non_significant_CpGs
      C <- total_sig_all - A
      D <- total_nonsig_all - B
      
      fisher_matrix <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
      
      # Run test
      test_result <- fisher.test(fisher_matrix)
      test_result$p.value
    }
  ) %>%
  ungroup() %>%
  # Correction for Multiple Testing (Benjamini-Hochberg)
  mutate(
    p_adj = p.adjust(p_value_raw, method = "BH"),
    is_significant_plot = p_adj < 0.05 # Significance flag based on FDR
  )

# Save the statistical results to Excel
# Ensure the directory exists
if(!dir.exists("./plots/Fisher-Exact-BW_dataset/")) dir.create("./plots/Fisher-Exact-BW_dataset/", recursive = TRUE)
write.xlsx(fisher_results, "./plots/Fisher-Exact-BW_dataset/fisher_exact_results.xlsx")

# ==============================================================================
# 3. VISUALIZATION
# ==============================================================================

# Prepare data for plotting
visual_data <- fisher_results %>%
  mutate(
    # Ensure standard chromosomal ordering (1..23, X)
    seqnames = factor(seqnames, levels = c(as.character(1:23), "X")),
    
    # Scale for bubble size: -log10 of the ADJUSTED p-value
    # Adding a tiny constant (e.g., 1e-300) avoids Inf if p_adj is 0, though rare in R
    neg_log10_p = -log10(p_adj+ 1e-300),
    # Create a label for significance: "*" if significant, empty if not
    sig_label = ifelse(p_adj < 0.05, "*", ""),
    # Create a column to control transparency (alpha)
    is_sig_alpha = ifelse(p_adj < 0.05, "Significant", "Not Significant")
  )

# Define Custom Theme
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

# --- Plot 1: Bar Plot (Observed vs Expected) ---
bar_plot<- ggplot(visual_data, aes(x = seqnames)) +
  
  # 1. THE BARS (Observed Count)
  # We map 'alpha' to significance so non-sig bars fade out
  geom_bar(aes(y = significant_CpGs, fill = direction, alpha = is_sig_alpha), 
           stat = "identity", width = 0.8) +
  
  # 2. THE POINTS (Expected Count)
  # We wrap 'fill' in aes() with a constant string to FORCE a legend entry
  geom_point(aes(y = expected_significant_CpGs, shape = "Expected Value"), 
             size = 3, color = "black", fill = "yellow", stroke = 1) +
  
  # 3. THE STARS (Significance Indicator)
  geom_text(aes(y = pmax(significant_CpGs, expected_significant_CpGs) + 5, label = sig_label), 
            size = 8, vjust = 0) +
  
  # 4. SCALES & LEGENDS
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue"), 
                    name = "Deviation Direction") +
  
  # Manual scale for shape to create the "Expected Value" legend item
  scale_shape_manual(name = "", values = c("Expected Value" = 21)) +
  
  # Manual scale for alpha (Solid for Sig, Faded for Non-Sig)
  scale_alpha_manual(values = c("Significant" = 1.0, "Not Significant" = 0.4), 
                     name = "Fisher's Test (FDR < 0.05)") +
  
  # 5. LABELS
  labs(
    title = "Observed vs. Expected Significant CpGs (FDR < 0.1)",
    subtitle = "Asterisk (*) denotes statistical significance (Fisher's Exact Test, FDR < 0.05)",
    x = "Chromosome",
    y = "Number of Significant CpGs"
  ) +
  
  # 6. LEGEND ORGANIZATION
  guides(
    fill = guide_legend(order = 1),
    alpha = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  ) +
  custom_theme

# Save and Print
ggsave("./plots/Fisher-Exact-BW_dataset/observed_vs_expected_barplot.png", bar_plot, width = 12, height = 7)
print(bar_plot)

# --- Plot 2: Scatter Plot (The main visualization) ---
scatter_plot <- ggplot(visual_data, aes(x = expected_significant_CpGs, y = significant_CpGs)) +
  # Reference line (Expected = Observed)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # Points: Shape = Significance (FDR < 0.05), Size = Magnitude of significance
  geom_point(aes(color = direction, 
                 size = neg_log10_p, 
                 shape = is_significant_plot), 
             alpha = 0.8) +
  
  # Labels
  geom_text_repel(aes(label = seqnames), size = 5, box.padding = 0.5, max.overlaps = Inf) +
  
  # Scales
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c("NS (FDR > 0.05)", "Sig (FDR < 0.05)")) + 
  scale_size_continuous(range = c(3, 10), name = "-Log10(FDR)") +
  
  # Labels and Guides
  labs(
    title = "Scatter Plot: Observed vs. Expected Significant CpGs",
    subtitle = "Fisher's Exact Test with Benjamini-Hochberg Correction",
    x = "Expected Significant CpGs",
    y = "Observed Significant CpGs",
    color = "Deviation",
    shape = "Significance"
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) +
  custom_theme

ggsave("./plots/Fisher-Exact-BW_dataset/scatter_plot_fisher.png", scatter_plot, width = 12, height = 8)
print(scatter_plot)

# --- Plot 3: Percentage Plot ---
percentage_plot <- ggplot(visual_data, aes(x = seqnames, y = ratio * 100, fill = direction)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(
    title = "Percentage of Significant CpGs per Chromosome",
    subtitle = "Relative to total CpGs surveyed on that chromosome",
    x = "Chromosome",
    y = "Percentage Significant CpGs (%)",
    fill = "Deviation"
  ) +
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/Fisher-Exact-BW_dataset/percentage_significant_CpGs_fisher.png", percentage_plot, width = 12, height = 6)
print(percentage_plot)

# --- Plot 4: Log P-value Plot ---
log_p_plot <- ggplot(visual_data, aes(x = seqnames, y = neg_log10_p)) +
  geom_point(aes(color = direction, size = ratio), alpha = 0.8) +
  # Threshold line for FDR = 0.05 (-log10(0.05) ≈ 1.3)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  annotate("text", x = 1, y = -log10(0.05) + 0.2, label = "FDR = 0.05", color = "darkgray", hjust = 0) +
  
  geom_text_repel(aes(label = seqnames), size = 5, box.padding = 0.5, max.overlaps = Inf) +
  
  labs(
    title = "-Log10(Adjusted P-value) by Chromosome",
    x = "Chromosome",
    y = "-Log10(FDR Adjusted P-value)",
    color = "Deviation",
    size = "Sig. Ratio"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/Fisher-Exact-BW_dataset/log_pvalue_plot_fisher.png", log_p_plot, width = 12, height = 6)
print(log_p_plot)