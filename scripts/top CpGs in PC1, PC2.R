

# Load necessary libraries
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glmnet)
library(caret)
library(tidyr)
library(tibble)



# Load all necessary dataset

map = read.csv("./data/Peromyscus_maniculatus_bairdii.hu_pman_2.1.100.HorvathMammalMethylChip40.v1.csv") #annotation 
key = read.xlsx("./data/DNAm mice info -BW dataset.xlsx") # mice information
beta = read.csv("./data/sesame_data_BW.csv") #normalized beta sesame - DNAm information

# Ensure CGids are unique
rownames(beta) <- make.unique(rownames(beta))

# Step 1: Convert `beta` to long format
beta_long <- as.data.frame(beta) %>%
  rownames_to_column("RowID") %>%         # Retain row names in a temporary column
  mutate(CGid = X) %>%                    # Assign `X` values to `CGid`
  dplyr:: select(-X) %>%                          # Remove the original `X` column
  pivot_longer(-c(RowID, CGid),           # Keep `RowID` and `CGid` as fixed columns
               names_to = "SID",
               values_to = "Beta") %>%
  dplyr:: select(CGid, SID, Beta)                 # Reorder columns for clarity


# Step 2: Merge with `key` to add metadata
beta_with_metadata <- beta_long %>%
  inner_join(key, by = c("SID" = "SID")) %>% distinct()




# Ensure Sex is a factor and Age is numeric
complete_data <- beta_with_metadata %>%
  mutate(
    Sex = as.factor(Sex),
    Age = as.numeric(Age)
  )

# Perform nested model testing for each CpG
nested_model_results <- complete_data %>% 
  group_by(CGid) %>%
  nest() %>%
  mutate(
    # Fit reduced and full models
    reduced_model = map(data, ~ lm(Beta ~ Age + Sex, data = .)),
    full_model = map(data, ~ lm(Beta ~ Relatedness + Age + Sex, data = .)),
    # Perform ANOVA between models
    anova_result = map2(full_model, reduced_model, anova),
    # Extract p-value for Relatedness term
    p_value = map_dbl(anova_result, ~ .x$`Pr(>F)`[2])
  )

nested_model_results <- nested_model_results %>% 
  dplyr::select(CGid, p_value) %>% 
  as.data.frame() %>% 
  mutate(FDR = p.adjust(p_value, method = "BH")) 

# Filter significant CpGs
significant_cpgs <- nested_model_results %>% 
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  dplyr::select(CGid, FDR, p_value)

########################


# Create beta matrix: rows = samples, columns = CpGs
beta_wide <- complete_data %>%
  dplyr :: select(ExternalSampleName, CGid, Beta) %>%
  distinct() %>%
  pivot_wider(names_from = CGid, values_from = Beta) %>%
  arrange(ExternalSampleName)  # Ensure rows are sorted by Sample ID




#Convert to matrix
beta_mat <- as.matrix(beta_wide %>% dplyr :: select(-ExternalSampleName))


# Run PCA
pca_res <- prcomp(beta_mat, center = TRUE, scale. = TRUE)



# Loadings (contribution of each CpG to PCs)
loadings <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  rownames_to_column("CGid")

# Rank CpGs by absolute loading for PC1 and PC2
top_pc1 <- loadings %>% 
  arrange(desc(abs(PC1))) %>%
  slice_head(n = 369)   # top 200 CpGs for PC1 (adjust cutoff as needed)

top_pc2 <- loadings %>% 
  arrange(desc(abs(PC2))) %>%
  slice_head(n = 369)

# Save lists
write.xlsx(list(Top_PC1 = top_pc1, Top_PC2 = top_pc2), 
           "tables/Top_CpGs_PC1_PC2.xlsx", overwrite = TRUE)



# Annotate CpGs with your annotation map
annotated_pc1 <- map %>% merge(top_pc1, by ="CGid")
write.xlsx(annotated_pc1, "tables/Top_CpGs_PC1_annotated.xlsx", overwrite = TRUE )

annotated_pc2 <- map %>% merge(top_pc2, by ="CGid")
write.xlsx(annotated_pc2, "tables/Top_CpGs_PC2_annotated.xlsx", overwrite = TRUE )





# significant_cpgs comes from your nested model pipeline (FDR < 0.05)

# Intersections
pc1_related_overlap <- intersect(top_pc1$CpG, significant_cpgs$CGid)
pc2_related_overlap <- intersect(top_pc2$CpG, significant_cpgs$CGid)

# Summaries
summary_df <- tibble(
  Category = c("PC1 top CpGs", "PC2 top CpGs", "Relatedness CpGs",
               "PC1 ∩ Relatedness", "PC2 ∩ Relatedness"),
  Count = c(nrow(top_pc1), nrow(top_pc2), nrow(significant_cpgs),
            length(pc1_related_overlap), length(pc2_related_overlap))
)

write.xlsx(summary_df, "tables/CpG_overlap_summary.xlsx", overwrite = TRUE)


