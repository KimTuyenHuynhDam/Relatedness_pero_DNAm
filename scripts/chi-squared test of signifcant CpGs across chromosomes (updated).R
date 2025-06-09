# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

# Load the data
data <- read.xlsx("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx", sheet = "Sheet 1") %>%
  filter(seqnames %in% c(as.character(1:22), "X"))

# Calculate the overall ratio of significant CpGs (FDR < 0.1)
total_CpGs <- nrow(data)  
significant_CpGs <- data %>%
  filter(FDR < 0.1) %>%
  nrow()

overall_ratio <- significant_CpGs / total_CpGs  

# Calculate chromosome-specific ratios
prepared_data <- data %>%
  distinct() %>%
  group_by(seqnames) %>% 
  summarize(
    total_CpGs = n(),
    significant_CpGs = sum(FDR < 0.1, na.rm = TRUE),
    non_significant_CpGs = total_CpGs - significant_CpGs,
    ratio = significant_CpGs / total_CpGs,
    expected_significant_CpGs = total_CpGs * overall_ratio,
    expected_non_significant_CpGs = total_CpGs - expected_significant_CpGs,
    deviation = significant_CpGs - expected_significant_CpGs,
    direction = ifelse(deviation > 0, "Higher", "Lower")
  ) %>%
  rowwise()

# Perform chi-squared test for each row
chi_squared_results <- prepared_data %>%
  rowwise() %>%
  mutate(
    p_value = {
      observed <- c(significant_CpGs, non_significant_CpGs)
      expected <- c(expected_significant_CpGs, expected_non_significant_CpGs)
      chisq_result <- chisq.test(x = observed, p = expected / sum(expected))
      chisq_result$p.value
    }
  ) %>%
  ungroup()

write.xlsx(chi_squared_results, "./plots/chi_squared-BW_dataset/chi_squared.xlsx")

# Prepare data for visualization
visual_data <- chi_squared_results %>%
  mutate(seqnames = factor(seqnames, levels = c(as.character(1:22), "X")))

# Define a theme with increased font size
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold")
  )

# Bar Plot
bar_plot <- ggplot(visual_data, aes(x = seqnames, fill = direction)) +
  geom_bar(aes(y = significant_CpGs), stat = "identity", alpha = 0.8, position = "dodge") +
  geom_point(aes(y = expected_significant_CpGs), shape = 21, size = 3, color = "black", fill = "yellow") +
  labs(
    title = "Observed vs. Expected Significant CpGs (FDR < 0.1)",
    x = "Chromosome",
    y = "Significant CpGs",
    fill = "Deviation Direction"
  ) +
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/chi_squared-BW_dataset/observed_vs_expected_barplot.png", bar_plot, width = 12, height = 6)
print(bar_plot)

# Scatter Plot
scatter_plot <- ggplot(visual_data, aes(x = expected_significant_CpGs, y = significant_CpGs)) +
  geom_point(aes(color = direction, size = -log10(p_value)), alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(aes(label = seqnames), size = 5) +
  labs(
    title = "Scatter Plot: Observed vs. Expected Significant CpGs",
    x = "Expected Significant CpGs",
    y = "Observed Significant CpGs",
    color = "Deviation Direction",
    size = "Significance (-log10 p-value)"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/chi_squared-BW_dataset/scatter_plot.png", scatter_plot, width = 12, height = 6)
print(scatter_plot)

# Percentage Plot
percentage_plot <- ggplot(visual_data, aes(x = seqnames, y = ratio * 100, fill = direction)) +
  geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
  labs(
    title = "Percentage of Significant CpGs (FDR < 0.1) per Chromosome",
    x = "Chromosome",
    y = "Percentage Significant CpGs (%)",
    fill = "Deviation Direction"
  ) +
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/chi_squared-BW_dataset/percentage_significant_CpGs.png", percentage_plot, width = 12, height = 6)
print(percentage_plot)

# -log10(P-value) Plot
log_p_plot <- ggplot(visual_data, aes(x = seqnames, y = -log10(p_value))) +
  geom_point(aes(color = direction, size = ratio), alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = seqnames), size = 5) +
  labs(
    title = "-log10(P-value) of Chromosome Deviations",
    x = "Chromosome",
    y = "-log10(P-value)",
    color = "Deviation Direction",
    size = "Significant CpGs Ratio"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  custom_theme

ggsave("./plots/chi_squared-BW_dataset/log_pvalue_plot.png", log_p_plot, width = 12, height = 6)
print(log_p_plot)

