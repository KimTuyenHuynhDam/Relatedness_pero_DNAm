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



# Entropy + SD calculation
entropy_DNAm_BW <- complete_data %>%
  group_by(ExternalSampleName) %>%
  summarize(
    entropy = {
      values <- na.omit(Beta)
      bins <- hist(values, breaks=seq(0, 1, length.out=11), plot=FALSE)
      probs <- bins$counts / sum(bins$counts)
      entropy.empirical(probs, unit="log2")
    }
  ) %>%
  left_join(key, by = "ExternalSampleName") %>% distinct() 






# Simple scatterplot with regression line
p_age_entropy_BW <- ggplot(entropy_DNAm_BW, aes(x = entropy, y = Age)) +
  geom_point(color = "darkblue", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Relationship between Shannon Entropy and Chronological Age",
    x = "Shannon Entropy",
    y = "Chronological Age (years)"
  ) +
  theme_minimal(base_size = 14)
# Print and save
print(p_age_entropy_BW)

ggsave('./plots/entropy/Relationship between Shannon Entropy and  Chronological Age.jpeg')





entropy_DNAm_BW <- entropy_DNAm_BW %>%
  mutate(Sex = ifelse(Sex == 'F', 'Females', 'Males'))  %>%
  mutate(Sex = as.factor(Sex)) 


# Fit linear model
lm_age_entropy <- lm( Age ~ entropy , data = entropy_DNAm_BW)

# Predict 
entropy_DNAm_BW <- entropy_DNAm_BW %>%
  mutate(
    Age_predicted_from_entropy = predict(lm_age_entropy, newdata = .),
    Residual_Age_predicted_from_entropy = Age_predicted_from_entropy - Age
  )





write.xlsx(entropy_DNAm_BW, './DNAm_BW/BW mice tail - DNAm - residual entropy.xlsx')



ggplot(entropy_DNAm_BW, aes(x = Residual_Age_predicted_from_entropy, y = Relatedness)) +
  geom_point(color = "blue", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Relationship between Residual Age predicted from entropy and Relatedness",
    x = "Residual Age predicted from entropy",
    y = "Relatedness"
  ) +
  theme_minimal(base_size = 14)


ggsave('./plots/entropy/Relationship between Residual_Age_predicted_from_entropy and Relatedness.jpeg')



#####################################



plot_scatter <- function(
    data,
    xvar,
    yvar,
    xlab,
    ylab,
    facet_by = NULL,
    title = NULL,
    filename = NULL,
    width = 8,
    height = 6
) {
  library(ggplot2)
  library(ggpubr)
  x_sym <- rlang::sym(xvar)
  y_sym <- rlang::sym(yvar)
  
  p <- ggplot(data, aes(x = !!x_sym, y = !!y_sym, color = Sex)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, aes(group = Sex), linetype = "dashed") +
    stat_cor(aes(group = Sex), label.x.npc = "left", label.y.npc = "top", method = "pearson") +
    labs(
      title = title %||% paste("Correlation between", xvar, "and", yvar, "by Sex"),
      x = xlab,
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  # Add optional faceting
  if (!is.null(facet_by)) {
    facet_formula <- as.formula(paste("~", facet_by))
    p <- p + facet_wrap(facet_formula)
  }
  
  # Print in viewer
  print(p)
  
  # Save if filename given
  if (!is.null(filename)) {
    ggsave(
      filename = paste0("./plots/entropy/", filename, ".jpeg"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300
    )
  }
}




#Residual_Age_predicted_from_entropy

plot_scatter(
  data = entropy_DNAm_BW ,
  xvar = "Residual_Age_predicted_from_entropy",
  yvar = "Relatedness",
  xlab = 'Residual age',
  ylab = 'Relatedness',
  facet_by = NULL,
  title = "Residual age vs Relatedness",
  filename = "Residual_Age_predicted_from_entropy vs Relatedness"
)


###################
plot_boxplot_relatedness <- function(
    data,
    xvar,
    yvar,
    xlab,
    ylab,
    facet_by = NULL,
    title = NULL,
    filename = NULL,
    width = 8,
    height = 6,
    legend_title = NULL,
    show_legend = TRUE
) {
  library(ggplot2)
  library(ggpubr)
  library(rlang)
  
  x_sym <- sym(xvar)
  y_sym <- sym(yvar)
  
  # Check if xvar has exactly 2 levels
  x_levels <- unique(data[[xvar]])
  if (length(x_levels) != 2) {
    stop("xvar must have exactly two levels for t-test comparison with bracket.")
  }
  
  # Variance check for Welch vs Student's t-test
  f_test <- var.test(reformulate(xvar, yvar), data = data)
  var_equal <- f_test$p.value > 0.05
  
  # Prepare the p-value comparison bracket
  comp <- list(as.character(x_levels))
  
  # Plot
  p <- ggplot(data, aes(x = !!x_sym, y = !!y_sym, color = !!x_sym)) +
    geom_boxplot(size = 1.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.15), alpha = 0.4) +
    stat_compare_means(
      comparisons = comp,
      method = "t.test",
      method.args = list(var.equal = var_equal),
      label = "p.format"
    ) +
    labs(
      title = title %||% paste("Comparison of", yvar, "by", xvar),
      x = xlab,
      y = ylab,
      color = legend_title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_color_brewer(palette = "Dark2")
  
  if (!is.null(facet_by)) {
    facet_formula <- as.formula(paste("~", facet_by))
    p <- p + facet_wrap(facet_formula)
  }
  
  print(p)
  
  if (!is.null(filename)) {
    ggsave(
      filename = paste0("./plots/entropy/", filename, ".jpeg"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300
    )
  }
}




entropy_DNAm_BW_high_low_relatedness <- entropy_DNAm_BW %>%
  mutate(high_low_relatedness = ifelse(Relatedness > mean(Relatedness), 'high', 'low'))





###########Residual_Age_predicted_from_entropy

plot_boxplot_relatedness(
  data = entropy_DNAm_BW_high_low_relatedness,
  xvar = "high_low_relatedness",
  yvar = "Residual_Age_predicted_from_entropy",
  xlab = 'Relatedness',
  ylab = 'Residual Age',
  facet_by = 'Sex',
  title = "Residual Age by Relatedness",
  filename = "Residual_Age_predicted_from_entropy_by_Relatedness_BySex",
  legend_title = "Relatedness Level",
  show_legend = TRUE
)


plot_boxplot_relatedness(
  data = entropy_DNAm_BW_high_low_relatedness,
  xvar = "high_low_relatedness",
  yvar = "Residual_Age_predicted_from_entropy",
  xlab = 'Relatedness',
  ylab = 'Residual Age',
  facet_by = NULL,
  title = "Residual Age by Relatedness",
  filename = "Residual_Age_predicted_from_entropy_by_Relatedness",
  legend_title = "Relatedness Level",
  show_legend = TRUE
)









