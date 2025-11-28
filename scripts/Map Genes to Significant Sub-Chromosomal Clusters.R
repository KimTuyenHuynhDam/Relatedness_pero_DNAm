# ==============================================================================
# Map Genes to Significant Sub-Chromosomal Clusters
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(glue)

# 1. Load the Data
# ----------------
# Load the significant windows you generated
sig_windows <- read_csv("./tables/Significant_SubChromosomal_Clusters.csv")

# Load the full CpG annotation data (contains gene symbols)
# Ensure column names match your file (e.g., SYMBOL, GeneRegionID)
nested_anova_data <- read_excel("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx")

# 2. Prepare the Annotation Data
# ------------------------------
# We need Chromosome, Position, Significance, and Gene Symbol
annotated_cpgs <- nested_anova_data %>%
  filter(!is.na(seqnames) & !is.na(CGstart)) %>%
  mutate(
    IsSignificant = ifelse(FDR < 0.05, TRUE, FALSE)
  )

# 3. Perform the Overlap (Map Windows -> CpGs -> Genes)
# ----------------------------------------------------
# We will loop through each significant window to find the genes inside

get_genes_in_window <- function(chr, start, end, cpg_data) {
  # Filter CpGs in this specific window
  window_cpgs <- cpg_data %>%
    filter(seqnames == as.character(chr)) %>%
    filter(CGstart >= start & CGstart < end)
  
  # 1. Get all genes present in this window (context)
  all_genes <- unique(window_cpgs$SYMBOL)
  all_genes <- all_genes[!is.na(all_genes)] # Remove NAs
  
  # 2. Get ONLY genes associated with SIGNIFICANT CpGs (the "drivers")
  sig_cpgs <- window_cpgs %>% filter(IsSignificant == TRUE)
  driver_genes <- unique(sig_cpgs$SYMBOL)
  driver_genes <- driver_genes[!is.na(driver_genes)]
  
  # Return a list or formatted string
  list(
    n_sig_cpgs_found = nrow(sig_cpgs),
    Driver_Genes = paste(driver_genes, collapse = "; "),
    All_Genes_In_Window = paste(all_genes, collapse = "; ")
  )
}

# Apply the function to each row of the significant windows dataframe
results <- sig_windows %>%
  rowwise() %>%
  mutate(gene_info = list(get_genes_in_window(seqnames, WindowStart, WindowEnd, annotated_cpgs))) %>%
  unnest_wider(gene_info) %>%
  ungroup()

# 4. View and Save Results
# ------------------------
# Sort by significance of the window
final_table <- results %>%
  arrange(FDR_Window) %>%
  dplyr:: select(seqnames, WindowStart, WindowEnd, FDR_Window, n_sig_in_window, Driver_Genes, All_Genes_In_Window)

# Print the top clusters and their driver genes
print(head(final_table))

# Save to CSV
write_csv(final_table, "./tables/Significant_Clusters_Gene_Map.csv")

print("Gene mapping complete. Results saved to './tables/Significant_Clusters_Gene_Map.csv'.")