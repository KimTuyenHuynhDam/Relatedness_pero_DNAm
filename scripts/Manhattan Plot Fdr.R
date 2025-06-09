# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)


# Step 1: Load the data
data <- read.xlsx("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx", sheet = "Sheet 1")

# Step 2: Data Preparation
manhattan_data <- data %>%
  dplyr:: select(CGid, FDR, estimate, seqnames, GeneRegionID) %>%   # Select relevant columns
  filter(!is.na(seqnames) & !is.na(FDR) & FDR > 0) %>%     # Remove rows with missing seqnames, FDR, or FDR <= 0
  mutate(
    logFDR = -log10(FDR),                                 # Transform FDR to -log10(FDR)
    color = ifelse(estimate > 0, "red", "blue"),        # Assign color based on estimate
    seqnames = factor(seqnames, levels = c(as.character(1:22), "X")), # Explicit chromosome order
    jittered_position = as.numeric(seqnames) + runif(n(), -0.25, 0.25) # Jitter x-axis positions
  ) %>%
  distinct(CGid, .keep_all = TRUE)                        # Keep unique values

# Step 3: Diagnostic Check for Missing Values
missing_values <- manhattan_data %>%
  filter(is.na(logFDR))
print("Rows excluded due to missing or out-of-range logFDR:")
print(missing_values)

# Retain all valid logFDR values within the range
manhattan_data <- manhattan_data %>%
  filter(!is.na(logFDR))

# Step 4: Annotate top 20 most significant points
top20 <- manhattan_data %>%
  arrange(FDR) %>%
  head(20)

manhattan_data_clean <- manhattan_data %>%
  filter(!is.na(seqnames) & !is.na(jittered_position))


# Step 5: Create Manhattan Plot
manhattan_plot <- ggplot(manhattan_data_clean, aes(x = jittered_position, y = logFDR, color = color)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggrepel::geom_text_repel(data = top20, aes(label = GeneRegionID),
                           size = 3, max.overlaps = 100, box.padding = 0.7, point.padding = 0.3) +
  labs(
    title = "Manhattan Plot of FDR-adjusted P-values",
    x = "Chromosome",
    y = "-log10(FDR)",
    caption = "Red: Positive correlation | Blue: Negative correlation"
  ) +
  scale_x_continuous(breaks = 1:23, labels = c(1:22, "X")) +
  scale_y_continuous(limits = c(0, 5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save and display the plot
ggsave("./plots/manhattan_plot_BW_dataset.png", manhattan_plot, width = 12, height = 8)
print(manhattan_plot)
