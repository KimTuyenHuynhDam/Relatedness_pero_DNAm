

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





# Select Top CpGs 
top_cpgs <- significant_cpgs %>%
  ungroup() %>%                # Remove grouping to avoid slicing within groups
  arrange(FDR) %>%             # Sort by FDR in ascending order
  pull(CGid)                   # Extract the CpG IDs



# Filter and reshape the data
ml_data <- complete_data %>%
  filter(CGid %in% top_cpgs) %>%
  dplyr::select(SID, CGid, Beta) %>%
  pivot_wider(names_from = CGid, values_from = Beta) %>%
  inner_join(
    complete_data %>%
      dplyr::select(SID, Age, Sex, Relatedness) %>% distinct(),
    by = "SID"
  )

# Convert Sex to numeric
ml_data <- ml_data %>%
  mutate(Sex = as.numeric(Sex))


set.seed(100)  # For reproducibility


# Perform feature selection using LASSO


cv_lasso <- cv.glmnet(as.matrix(ml_data %>%
                                  dplyr::select(all_of(top_cpgs), Age, Sex)), 
                      ml_data$Relatedness, 
                      alpha = 1, 
                      nfolds = 10)

best_lambda_lasso <- cv_lasso$lambda.min

lasso_model <- glmnet(as.matrix(ml_data %>%
                                  dplyr::select(all_of(top_cpgs), Age, Sex)), 
                      ml_data$Relatedness, 
                      alpha = 1, 
                      lambda = best_lambda_lasso)

# Extract non-zero coefficients
selected_features <- coef(lasso_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Feature") %>%
  filter(`s0` != 0) %>%
  pull(Feature)



selected_features <- setdiff(selected_features, c("(Intercept)", "Sex"))

# Prepare predictors (X) and target variable (y)
X <- ml_data %>%
  dplyr::select(all_of(selected_features)) %>%
  as.matrix()

# Standardize the predictors
X_scaled <- scale(X)
y <- ml_data$Relatedness


# Split data into 80% train and 20% test set

train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X_scaled[train_index, ]
y_train <- y[train_index]

X_test <- X_scaled[-train_index, ]
y_test <- y[-train_index]

# Train Elastic Net model
cv_model <- cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 10)
best_lambda <- cv_model$lambda.min
final_model <- glmnet(X_train, y_train, alpha = 0.5, lambda = best_lambda)

# Predict on test set
y_pred <- predict(final_model, newx = X_test, s = best_lambda) %>% as.vector()

# Compute performance metrics
test_rmse <- sqrt(mean((y_test - y_pred)^2))  # Root Mean Squared Error
test_r2 <- cor(y_test, y_pred)^2  # R-squared Score

# Print performance metrics
print(glue::glue("Test RMSE: {round(test_rmse, 3)}"))
print(glue::glue("Test R²: {round(test_r2, 3)}"))

# Plot observed vs predicted values
ggplot(data.frame(Observed = y_test, Predicted = y_pred), aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Observed vs Predicted Relatedness (Train-Test Split with Feature Selection)", x = "Observed Relatedness", y = "Predicted Relatedness") +
  theme_minimal() +
  annotate("text", x = min(y), y = max(y), label = glue::glue("RMSE: {round(test_rmse, 3)}\nR²: {round(test_r2, 3)}"), hjust = 0, vjust = 1, size = 5, color = "blue")


ggsave('./plots/Observed vs Predicted Relatedness (Train-Test Split with Feature Selection).jpeg')
