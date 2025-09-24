############################################################
# Relatedness & Epigenetic Aging – COMPLETE ANALYSIS PIPELINE
# - Computes EAA (clock residuals)
# - Computes methylome entropy (Shannon & KDE)
# - Tests relatedness vs EAA (added-variable)
# - Tests entropy vs age / relatedness / EAA
# - Runs PCA on beta matrix; relates PCs to EAA & relatedness
# - Mediation (PCs as mediators of Relatedness -> EAA)
# - Robustness: outlier trimming, train/test split
# - Saves figures and tables in organized folders
############################################################

## 0) Libraries & setup ------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(broom)
  library(glue)
  library(ggpubr)
  library(entropy)
  library(wateRmelon)
  library(mediation)
})

theme_set(theme_bw(base_size = 12))

# Create output folders
dir.create("plots", showWarnings = FALSE)
dir.create("plots/fig1", showWarnings = FALSE)
dir.create("plots/fig2", showWarnings = FALSE)
dir.create("plots/pca", showWarnings = FALSE)
dir.create("plots/mediation", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

save_plot <- function(p, file, w = 7, h = 5, dpi = 300) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = dpi)
}

## Load data ---------------------------------------------
# Annotation (not used below but kept for completeness)
map  <- read.csv("./data/Peromyscus_maniculatus_bairdii.hu_pman_2.1.100.HorvathMammalMethylChip40.v1.csv")

# Beta matrix (sesame-normalized; rows = CpGs, columns = SIDs/arrays)
beta <- read.csv("./data/sesame_data_BW.csv")

# Sample-level metadata including DNAmAgeTailFinal, Age, Sex, Relatedness, IDs
DNAm_Age <- read.xlsx("./data/DNAm_Age.xlsx") %>%
  mutate(
    Sex = factor(Sex),
    Relatedness = as.numeric(Relatedness),
    Age = as.numeric(Age),
    ExternalSampleName = as.character(ExternalSampleName)
  )

# If you have a separate key with SID <-> ExternalSampleName mapping, load it:
key <- read.xlsx("./data/DNAm mice info -BW dataset.xlsx") %>%
  mutate(ExternalSampleName = as.character(ExternalSampleName),
         SID = as.character(SID))

## Build long/wide beta & join metadata ------------------
# Ensure unique rownames and capture CpG ids
rownames(beta) <- make.unique(rownames(beta))


beta_long <- as.data.frame(beta) %>%
  rownames_to_column("RowID") %>%         # Retain row names in a temporary column
  mutate(CGid = X) %>%                    # Assign `X` values to `CGid`
  dplyr:: select(-X) %>%                          # Remove the original `X` column
  pivot_longer(-c(RowID, CGid),           # Keep `RowID` and `CGid` as fixed columns
               names_to = "SID",
               values_to = "Beta") %>%
  dplyr:: select(CGid, SID, Beta)  




  beta_with_meta <- beta_long %>%
    inner_join(key, by = "SID") %>% distinct()


# Make the wide matrix (rows = samples; cols = CpGs)
beta_wide <- beta_with_meta %>%
  dplyr :: select(ExternalSampleName, CGid, Beta) %>%
  distinct() %>%
  pivot_wider(names_from = CGid, values_from = Beta) %>%
  arrange(ExternalSampleName)

sample_ids <- beta_wide$ExternalSampleName
beta_mat   <- as.matrix(beta_wide %>% dplyr :: select(-ExternalSampleName))
rownames(beta_mat) <- sample_ids

## Entropy per sample (two ways) -------------------------
# Shannon entropy (binned) per sample
entropy_shannon <- beta_with_meta %>%
  group_by(ExternalSampleName) %>%
  summarize(
    entropy_shannon = {
      v <- na.omit(Beta)
      bins  <- hist(v, breaks = seq(0, 1, length.out = 11), plot = FALSE)
      probs <- bins$counts / sum(bins$counts)
      entropy::entropy.empirical(probs, unit = "log2")
    },
    .groups = "drop"
  )

# Continuous KDE-based entropy per sample
compute_entropy_kde <- function(beta_values, n = 512) {
  beta_values <- na.omit(beta_values)
  if (length(beta_values) < 50) return(NA_real_)
  dens <- density(beta_values, from = 0, to = 1, n = n)
  p <- dens$y
  p[p == 0] <- 1e-12
  -sum(p * log2(p)) * (dens$x[2] - dens$x[1])
}
entropy_kde <- apply(beta_mat, 1, compute_entropy_kde)
entropy_kde_df <- tibble(ExternalSampleName = names(entropy_kde),
                         entropy_kde = as.numeric(entropy_kde))

# Merge entropy into DNAm_Age
DNAm_Age <- DNAm_Age %>%
  left_join(entropy_shannon, by = "ExternalSampleName") %>%
  left_join(entropy_kde_df, by = "ExternalSampleName")

## Compute EAA (clock residuals) -------------------------
m_clock <- lm(DNAmAgeTailFinal ~ Age + Sex, data = DNAm_Age)
DNAm_Age <- DNAm_Age %>%
  mutate(EAA = resid(m_clock))

# Save model summary
write.xlsx(broom::tidy(m_clock), "tables/model_clock_ols.xlsx", overwrite = TRUE)

## FIGURE 1 – Relatedness & EAA (adjusted) ---------------
# Added-variable: residualize both by Age+Sex
res_x <- resid(lm(Relatedness ~ Age + Sex, data = DNAm_Age))
res_y <- resid(lm(EAA ~ Age + Sex, data = DNAm_Age))  # equals EAA (sanity check)

df_add <- tibble(res_rel = res_x, res_eaa = res_y)
fit_add <- lm(res_eaa ~ res_rel, data = df_add)
b_add   <- coef(fit_add)[2]; r2_add <- summary(fit_add)$r.squared
p_add   <- summary(fit_add)$coefficients["res_rel", "Pr(>|t|)"]

p_added_variable_relatedness_EAA <- ggplot(df_add, aes(res_rel, res_eaa)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Added-variable plot: Relatedness vs EAA (adjusted for Age+Sex)",
       subtitle = glue("β = {round(b_add, 4)}, partial R² = {round(r2_add, 3)}, p = {formatC(p_add, format='e', digits=2)}"),
       x = "Relatedness (residualized by Age + Sex)",
       y = "EAA (residualized by Age + Sex)") +
  theme_bw(base_size = 12)

print(p_added_variable_relatedness_EAA)

save_plot(p_added_variable_relatedness_EAA, "plots/added_variable_relatedness_EAA.png")


# ===========================
# Leave-one-out (jackknife) – Relatedness → EAA

# ===========================

DNAm_Age$._idx <- seq_len(nrow(DNAm_Age))  # keep a simple index to label samples

jack_results <- purrr::map_dfr(DNAm_Age$._idx, function(i) {
  # leave one out
  df_loo <- df_add[-i, , drop = FALSE]
  fit_loo <- lm(res_eaa ~ res_rel, data = df_loo)
  tibble(
    left_out = DNAm_Age$ExternalSampleName[i],
    beta     = coef(fit_loo)[2],
    se       = coef(summary(fit_loo))["res_rel","Std. Error"],
    t        = coef(summary(fit_loo))["res_rel","t value"],
    p        = coef(summary(fit_loo))["res_rel","Pr(>|t|)"]
  )
})

# Save table for supplement
write.xlsx(jack_results, "tables/Table_S_LOO_relatedness_EAA.xlsx", overwrite = TRUE)

jack_extremes <- jack_results %>%
  summarize(min_beta = min(beta, na.rm = TRUE),
            max_beta = max(beta, na.rm = TRUE),
            median_beta = median(beta, na.rm = TRUE),
            prop_sig_p05 = mean(p < 0.05, na.rm = TRUE))
write.xlsx(jack_extremes, "tables/LOO_beta_summary.xlsx", overwrite = TRUE)


# Plot β across leave-one-out fits



beta_full <- coef(fit_add)[2]
p_loo <- ggplot(jack_results, aes(x = reorder(left_out, beta), y = beta)) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = beta_full, linetype = 2, color = "steelblue") +
  coord_flip() +
  labs(title = "Leave-one-out β(Relatedness) for EAA (adjusted for Age+Sex)",
       subtitle = "Dashed line = full-sample β",
       x = "Sample left out", y = "β (res_eaa ~ res_rel)") +
  theme_bw(base_size = 11)

print(p_loo)
ggsave("plots/fig1/fig1C_LOO_beta_relatedness_EAA.png", p_loo, width = 7.5, height = 10, dpi = 300)

# Compact distribution view (helpful for main figure inset or supplement)
p_loo_dist <- ggplot(jack_results, aes(beta)) +
  geom_density(fill = "grey85", color = "grey25") +
  geom_vline(xintercept = beta_full, color = "steelblue", size = 1) +
  labs(title = "Distribution of leave-one-out β(Relatedness) → EAA",
       subtitle = "Vertical line = full-sample β",
       x = "β", y = "Density") +
  theme_bw(base_size = 12)

print(p_loo_dist)

ggsave("plots/fig1/fig1C_LOO_beta_density.png", p_loo_dist, width = 6.5, height = 4.2, dpi = 300)


## FIGURE 2 – Entropy analyses ---------------------------
# Entropy vs Age
fit_ent_age <- lm(entropy_kde ~ Age, data = DNAm_Age)

p_ent_age <- ggplot(DNAm_Age, aes(Age, entropy_kde)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "Entropy (KDE) vs chronological age",
       subtitle = glue("β = {round(coef(fit_ent_age)[2], 4)}, p = {formatC(summary(fit_ent_age)$coefficients[2,4], format='e', digits=2)}"),
       x = "Age", y = "Entropy (KDE)") +
  theme_bw(base_size = 12)

print(p_ent_age)

save_plot(p_ent_age, "plots/fig2/fig2A_entropy_vs_age.png")

# Entropy vs Relatedness, adjusted Age+Sex
fit_ent_rel <- lm(entropy_kde ~ Relatedness + Age + Sex, data = DNAm_Age)
p_ent_rel <- ggplot(DNAm_Age, aes(Relatedness, entropy_kde)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray30") +
  labs(title = "Entropy (KDE) vs Relatedness (unadjusted view)",
       subtitle = glue("Adjusted test (Age+Sex): β(Relatedness) p = {signif(coef(summary(fit_ent_rel))['Relatedness','Pr(>|t|)'], 3)}"),
       x = "Relatedness", y = "Entropy (KDE)") +
  theme_bw(base_size = 12)

print(p_ent_rel)

save_plot(p_ent_rel, "plots/fig2/fig2B_entropy_vs_relatedness.png")

# EAA vs Entropy (partial: adjust Age+Sex+Relatedness)
res_EAA   <- resid(lm(EAA ~ Relatedness, data = DNAm_Age))
res_Ent   <- resid(lm(entropy_kde ~ Relatedness, data = DNAm_Age))
fit_eaa_ent <- lm(res_EAA ~ res_Ent)
p_eaa_ent <- ggplot(tibble(res_Ent, res_EAA), aes(res_Ent, res_EAA)) +
  geom_point(alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "EAA vs Entropy (adjusted for Relatedness)",
       subtitle = glue("β = {round(coef(fit_eaa_ent)[2],4)}, partial R² = {round(summary(fit_eaa_ent)$r.squared,3)}, p = {formatC(summary(fit_eaa_ent)$coefficients[2,4], format='e', digits=2)}"),
       x = "Entropy (residualized by Relatedness)",
       y = "EAA (residualized by Relatedness)") +
  theme_bw(base_size = 12)

print(p_eaa_ent)

save_plot(p_eaa_ent, "plots/fig2/fig2C_EAA_vs_entropy_partial.png")

# Residual entropy vs residual relatedness (Age+Sex adjusted)
res_rel_AS <- resid(lm(Relatedness ~ Age + Sex, data = DNAm_Age))
res_ent_AS <- resid(lm(entropy_kde ~ Age + Sex, data = DNAm_Age))

fit_resres <- lm(res_ent_AS ~ res_rel_AS)
p_resres <- ggplot(tibble(res_rel_AS, res_ent_AS), aes(res_rel_AS, res_ent_AS)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Residual Entropy vs Residual Relatedness (both adjusted for Age+Sex)",
       subtitle = glue("β = {round(coef(fit_resres)[2],4)}, R² = {round(summary(fit_resres)$r.squared,3)}, p = {signif(summary(fit_resres)$coefficients[2,4],3)}"),
       x = "Relatedness residual (Age+Sex)",
       y = "Entropy residual (Age+Sex)") +
  theme_bw(base_size = 12)

print(p_resres)

save_plot(p_resres, "plots/fig2/fig2D_residual_entropy_vs_residual_relatedness.png")


## PCA on beta matrix ------------------------------------
pca_res <- prcomp(beta_mat, center = TRUE, scale. = TRUE)
var_explained <- (pca_res$sdev)^2 / sum(pca_res$sdev^2)
cum_var <- cumsum(var_explained)

# Scree
# Build data frame of variance explained
df_scree <- tibble(
  PC      = seq_along(var_explained),
  PropVar = as.numeric(var_explained)
)

# How many PCs to show (cap at available length)
k_show <- min(10, nrow(df_scree))

# Scree plot
p_scree <- ggplot(dplyr::slice_head(df_scree, n = k_show), aes(PC, PropVar)) +
  geom_point() + geom_line() +
  labs(title = "Scree plot", y = "Proportion of variance explained", x = "Principal component") +
  theme_bw(base_size = 12)

print(p_scree)

# Cumulative variance plot
p_cum <- ggplot(dplyr::slice_head(df_scree, n = k_show),
                aes(PC, cumsum(PropVar))) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.80, linetype = 2, color = "red") +
  labs(title = "Cumulative variance explained", y = "Cumulative proportion", x = "Principal component") +
  theme_bw(base_size = 12)

print(p_cum)
# Save
ggsave("plots/pca/scree_plot.png", p_scree, width = 7, height = 5, dpi = 300)
ggsave("plots/pca/cumulative_variance.png", p_cum, width = 7, height = 5, dpi = 300)


# Merge PCs (first 5) into DNAm_Age
pc_df <- as.data.frame(pca_res$x[, 1:5]) %>%
  rownames_to_column("ExternalSampleName") %>%
  mutate(ExternalSampleName = as.character(ExternalSampleName))
DNAm_Age <- DNAm_Age %>% left_join(pc_df, by = "ExternalSampleName")

# PCA colored by relatedness & by sex
p_pca_rel <- ggplot(DNAm_Age, aes(PC1, PC2, color = Relatedness)) +
  geom_point(size = 2, alpha = 0.9) + scale_color_viridis_c() +
  labs(title = "PCA of methylation: colored by Relatedness") +
  theme_bw(base_size = 12)

print(p_pca_rel)
save_plot(p_pca_rel, "plots/pca/pca_PC1_PC2_by_relatedness.png")

p_pca_sex <- ggplot(DNAm_Age, aes(PC1, PC2, color = Sex)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(title = "PCA of methylation: colored by Sex") +
  theme_bw(base_size = 12)

print(p_pca_sex)


save_plot(p_pca_sex, "plots/pca/pca_PC1_PC2_by_sex.png")

# Residual Age (EAA) ~ PCs
DNAm_Age <- DNAm_Age %>%
  mutate(resid_age = resid(lm(DNAmAgeTailFinal ~ Age + Sex, data = .)))

fit_pc_resid <- lm(resid_age ~ PC1 + PC2 + PC3, data = DNAm_Age)
write.xlsx(broom::tidy(fit_pc_resid), "tables/model_residAge_vs_PCs.xlsx", overwrite = TRUE)

p_res_PC1 <- ggplot(DNAm_Age, aes(PC1, resid_age)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Residual DNAmAge ~ PC1") +
  theme_bw(base_size = 12)

print(p_res_PC1)
save_plot(p_res_PC1, "plots/pca/residAge_PC1.png")

p_res_PC2 <- ggplot(DNAm_Age, aes(PC2, resid_age)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Residual DNAmAge ~ PC2") +
  theme_bw(base_size = 12)

print(p_res_PC2)

save_plot(p_res_PC2, "plots/pca/residAge_PC2.png")

## Model comparison: do PCs improve EAA model ------------
m_base <- lm(resid_age ~ Relatedness + Age + Sex, data = DNAm_Age)
m_pc   <- lm(resid_age ~ Relatedness + Age + Sex + PC1 + PC2, data = DNAm_Age)
anova_res <- anova(m_base, m_pc)
write.xlsx(broom::tidy(m_base), "tables/model_base.xlsx", overwrite = TRUE)
write.xlsx(broom::tidy(m_pc),   "tables/model_plusPCs.xlsx", overwrite = TRUE)
write.xlsx(as.data.frame(anova_res), "tables/anova_base_vs_plusPCs.xlsx", overwrite = TRUE)

## Mediation: PCs mediating Relatedness -> EAA ----------
med_model1    <- lm(PC1 ~ Relatedness + Age + Sex, data = DNAm_Age)
outcome_mod1  <- lm(resid_age ~ PC1 + Relatedness + Age + Sex, data = DNAm_Age)
med_out1 <- mediate(med_model1, outcome_mod1, treat = "Relatedness", mediator = "PC1",
                    boot = TRUE, sims = 1000)
capture.output(summary(med_out1), file = "tables/mediation_PC1.txt")

med_model2    <- lm(PC2 ~ Relatedness + Age + Sex, data = DNAm_Age)
outcome_mod2  <- lm(resid_age ~ PC2 + Relatedness + Age + Sex, data = DNAm_Age)
med_out2 <- mediate(med_model2, outcome_mod2, treat = "Relatedness", mediator = "PC2",
                    boot = TRUE, sims = 1000)
capture.output(summary(med_out2), file = "tables/mediation_PC2.txt")

#---Tidy a mediate() object into one row with CIs -------

tidy_med <- function(mo, label = "PC") {
  # mediate objects typically expose *_avg values and *_sims vectors
  # Fallback to non-avg fields if needed.
  getv <- function(x, alt) if (!is.null(mo[[x]])) mo[[x]] else mo[[alt]]
  
  # Estimates
  acme  <- getv("z.avg", "z0")     # average ACME
  ade   <- getv("d.avg", "d0")     # average ADE
  te    <- mo$tau.coef             # total effect
  pm    <- getv("n.avg", "n0")     # proportion mediated
  
  # Sims for CIs (bootstrap)
  acme_s <- getv("z.avg.sims", "z0.sims")
  ade_s  <- getv("d.avg.sims", "d0.sims")
  te_s   <- if (!is.null(mo$tau.sims)) mo$tau.sims else NULL
  pm_s   <- getv("n.avg.sims", "n0.sims")
  
  ci <- function(v) if (is.null(v)) c(NA, NA) else stats::quantile(v, c(.025, .975), na.rm = TRUE)
  
  tibble(
    Mediator = label,
    term     = c("ACME (indirect)", "ADE (direct)", "Total effect", "Proportion mediated"),
    estimate = c(acme, ade, te, pm),
    ci_low   = c(ci(acme_s)[1], ci(ade_s)[1], ci(te_s)[1], ci(pm_s)[1]),
    ci_high  = c(ci(acme_s)[2], ci(ade_s)[2], ci(te_s)[2], ci(pm_s)[2])
  )
}

med_tab_PC1 <- tidy_med(med_out1, label = "PC1")
med_tab_PC2 <- tidy_med(med_out2, label = "PC2")

med_tab_both <- bind_rows(med_tab_PC1, med_tab_PC2)

# Save a tidy table for the supplement
write.xlsx(med_tab_both, "tables/Table_S_mediation_PC1_PC2.xlsx", overwrite = TRUE)

#---Plot ACME/ADE/TE side-by-side for PC1 vs PC2 --------

plot_med_bar <- function(df, terms = c("ACME (indirect)", "ADE (direct)", "Total effect")) {
  df %>%
    filter(term %in% terms) %>%
    ggplot(aes(x = term, y = estimate, color = Mediator)) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.15, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    labs(x = NULL, y = "Effect (units of EAA)",
         title = "Mediation components for Relatedness → EAA",
         subtitle = "Points = estimate; bars = 95% bootstrap CI") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

p_med_core <- plot_med_bar(med_tab_both)

print(p_med_core)

ggsave("plots/mediation/mediation_PC1_PC2_ACME_ADE_TE.png", p_med_core, width = 7.5, height = 5, dpi = 300)

# Proportion mediated
p_med_prop <- med_tab_both %>%
  filter(term == "Proportion mediated") %>%
  ggplot(aes(x = Mediator, y = estimate, color = Mediator)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion mediated",
       title = "Proportion of Relatedness → EAA mediated by PCs",
       subtitle = "Estimate with 95% bootstrap CI") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

print(p_med_prop)


ggsave("plots/mediation/mediation_PC1_PC2_proportion.png", p_med_prop, width = 5.5, height = 4.5, dpi = 300)

#--- Sensitivity plots for PC1 and PC2 -------------------

sens1 <- medsens(med_out1, rho.by = 0.1, sims = 100)
png("plots/mediation/mediation_sensitivity_PC1.png", width = 1000, height = 420, res = 140)
plot(sens1)
dev.off()

sens2 <- medsens(med_out2, rho.by = 0.1, sims = 100)
png("plots/mediation/mediation_sensitivity_PC2.png", width = 1000, height = 420, res = 140)
plot(sens2)
dev.off()

## Sex-stratified mediation (optional) -----------------
sex_levels <- levels(DNAm_Age$Sex)
all_lines <- list()
for (sx in sex_levels) {
  subdat <- DNAm_Age %>% filter(Sex == sx)
  med_m1 <- lm(PC1 ~ Relatedness + Age, data = subdat)
  out_m1 <- lm(resid_age ~ PC1 + Relatedness + Age, data = subdat)
  mo1    <- mediate(med_m1, out_m1, treat = "Relatedness", mediator = "PC1",
                    boot = TRUE, sims = 1000)
  
  med_m2 <- lm(PC2 ~ Relatedness + Age, data = subdat)
  out_m2 <- lm(resid_age ~ PC2 + Relatedness + Age, data = subdat)
  mo2    <- mediate(med_m2, out_m2, treat = "Relatedness", mediator = "PC2",
                    boot = TRUE, sims = 1000)
  
  capture.output(
    list(Sex = sx, PC1 = summary(mo1), PC2 = summary(mo2)),
    file = glue("tables/mediation_sex_{sx}.txt")
  )
}

## Outlier-trimmed mediation ---------------------------
q_low  <- quantile(DNAm_Age$resid_age, 0.01, na.rm = TRUE)
q_high <- quantile(DNAm_Age$resid_age, 0.99, na.rm = TRUE)
DNAm_trim <- DNAm_Age %>% filter(resid_age > q_low, resid_age < q_high)

med_model_trim   <- lm(PC1 ~ Relatedness + Age + Sex, data = DNAm_trim)
outcome_model_tm <- lm(resid_age ~ PC1 + Relatedness + Age + Sex, data = DNAm_trim)
med_out_trim <- mediate(med_model_trim, outcome_model_tm, treat = "Relatedness", mediator = "PC1",
                        boot = TRUE, sims = 1000)
capture.output(summary(med_out_trim), file = "tables/mediation_PC1_trimmed.txt")

## 12) Train/Test predictive check -------------------------
set.seed(123)
idx  <- sample(seq_len(nrow(DNAm_Age)), size = floor(0.7 * nrow(DNAm_Age)))
train <- DNAm_Age[idx, ]; test <- DNAm_Age[-idx, ]

m_train <- lm(resid_age ~ Relatedness + Age + Sex + PC1 + PC2, data = train)
test$pred <- predict(m_train, newdata = test)
cor_val   <- cor(test$pred, test$resid_age, use = "complete.obs")
writeLines(glue("Holdout correlation (pred vs true resid_age): {round(cor_val,3)}"),
           "tables/holdout_correlation.txt")



#--- Sex-stratified ACME plot -----------------


have_all <- all(c("med_F_PC1","med_F_PC2","med_M_PC1","med_M_PC2") %in% ls())

if (have_all) {
  med_sex_tab <- bind_rows(
    tidy_med(med_F_PC1, "PC1 (F)"),
    tidy_med(med_F_PC2, "PC2 (F)"),
    tidy_med(med_M_PC1, "PC1 (M)"),
    tidy_med(med_M_PC2, "PC2 (M)")
  )
  
  p_sex_acme <- med_sex_tab %>%
    filter(term == "ACME (indirect)") %>%
    ggplot(aes(x = Mediator, y = estimate, color = Mediator)) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    labs(x = NULL, y = "ACME (indirect effect)",
         title = "Sex-stratified ACME for Relatedness → EAA",
         subtitle = "Dots = estimate; bars = 95% CI") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave("plots/mediation/mediation_sex_stratified_ACME.png", p_sex_acme,
         width = 7.5, height = 4.5, dpi = 300)
}

#--- Show bootstrap distribution of ACME ------
# Helpful for supplement: density of indirect effects

plot_acme_density <- function(mo, label = "PC") {
  sims <- if (!is.null(mo$z.avg.sims)) mo$z.avg.sims else mo$z0.sims
  tibble(value = sims) %>%
    ggplot(aes(value)) +
    geom_density(fill = "grey85", color = "grey20") +
    geom_vline(xintercept = mean(sims, na.rm = TRUE), color = "steelblue", size = 1) +
    geom_vline(xintercept = quantile(sims, c(.025,.975), na.rm = TRUE),
               linetype = 2, color = "steelblue") +
    labs(title = glue("Bootstrap distribution of ACME ({label})"),
         x = "ACME (indirect effect)", y = "Density") +
    theme_bw(base_size = 12)
}

p_acme_PC1 <- plot_acme_density(med_out1, "PC1")
print(p_acme_PC1)

p_acme_PC2 <- plot_acme_density(med_out2, "PC2")

print(p_acme_PC2)



ggsave("plots/mediation/acme_density_PC1.png", p_acme_PC1, width = 6.5, height = 4.2, dpi = 300)
ggsave("plots/mediation/acme_density_PC2.png", p_acme_PC2, width = 6.5, height = 4.2, dpi = 300)



## Compact tables for manuscript -----------------------
# Table S1-like: Relatedness -> EAA (OLS, HC3, robust)
library(sandwich); library(lmtest); library(MASS)
ols <- lm(EAA ~ Relatedness, data = DNAm_Age)
hc3 <- coeftest(ols, vcov. = vcovHC(ols, type = "HC3"))

rlm_fit <- rlm(EAA ~ Relatedness, data = DNAm_Age)
tabS1 <- tibble(
  Method    = c("OLS", "OLS + HC3 SE", "Robust (rlm)"),
  estimate  = c(coef(ols)["Relatedness"],
                hc3["Relatedness","Estimate"],
                coef(rlm_fit)["Relatedness"]),
  std.error = c(coef(summary(ols))["Relatedness","Std. Error"],
                hc3["Relatedness","Std. Error"],
                coef(summary(rlm_fit))["Relatedness","Std. Error"]),
  statistic = c(coef(summary(ols))["Relatedness","t value"],
                hc3["Relatedness","t value"],
                coef(summary(rlm_fit))["Relatedness","t value"]),
  p.value   = c(coef(summary(ols))["Relatedness","Pr(>|t|)"],
                hc3["Relatedness","Pr(>|t|)"],
                NA_real_)
)
write.xlsx(tabS1, "tables/Table_S1_relatedness_EAA.xlsx", overwrite = TRUE)

# Save a model digest
write.xlsx(broom::glance(ols), "tables/ols_glance.xlsx", overwrite = TRUE)
