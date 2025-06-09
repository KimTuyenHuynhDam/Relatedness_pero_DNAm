
# Load necessary libraries
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)
library(kinship2)
library(readxl)
library(dplyr)
library(ggplot2)



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




# Nested modeling, ANOVA, and extracting direction
nested_anova <- complete_data %>%
  group_by(CGid) %>%
  nest() %>%
  mutate(
    # Fit reduced and full models
    reduced_model = map(data, ~ lm(Beta ~ Sex + Age, data = .)),
    full_model = map(data, ~ lm(Beta ~ Relatedness + Sex + Age, data = .)),
    # Perform ANOVA to compare models
    anova_result = map2(full_model, reduced_model, anova),
    p_value = map_dbl(anova_result, ~ .x$`Pr(>F)`[2]),
    # Extract estimate (direction) of Relatedness from the full model
    full_model_tidy = map(full_model, broom::tidy),
    estimate = map_dbl(full_model_tidy, ~ .x$estimate[.x$term == "Relatedness"])
  )

# Finalize results with FDR adjustment and annotation
nested_anova_final <- nested_anova %>%
  dplyr:: select(CGid, p_value, estimate) %>% as.data.frame() %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  inner_join(map, by = "CGid")

# Export the final results to an Excel file
write.xlsx(nested_anova_final, "./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx")

