############################################################
# Relatedness & Epigenetic Aging – COMPLETE PIPELINE

# Includes:
#   - EAA residuals
#   - Relatedness vs EAA regression (robustness)
#   - Entropy vs Age/Relatedness/EAA
#   - PCA analysis of DNAm
#   - Mediation (PC1, PC2)
#   - Sex-stratified mediation

############################################################

## 0) Libraries & setup ------------------------------------
CORE_PKGS <- c(
  "tidyverse",   # dplyr, tidyr, tibble, purrr, ggplot2, readr, stringr, forcats
  "readxl",
  "openxlsx",
  "writexl",
  "broom",
  "glue",
  "ggpubr",
  "entropy",
  "mediation",
  "MASS",        # robust regression helpers
  "sandwich",    # HC3 robust SE
  "lmtest"       # coeftest, tests
)

ensure_packages <- function(pkgs, repos = getOption("repos")) {
  # Why: make this idempotent & CI-friendly; fail fast if install fails.
  if (is.null(repos) || identical(repos, c(CRAN = "@CRAN@"))) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, quiet = TRUE)
  }
  invisible(NULL)
}

load_quietly <- function(pkgs) {
  failed <- character(0)
  suppressPackageStartupMessages({
    for (p in pkgs) {
      ok <- require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      if (!ok) failed <- c(failed, p)
    }
  })
  if (length(failed)) {
    stop(sprintf("Failed to load packages: %s", paste(failed, collapse = ", ")), call. = FALSE)
  }
  invisible(pkgs)
}

# ---- Run once per session ----
ensure_packages(CORE_PKGS)
load_quietly(CORE_PKGS)



theme_set(theme_bw(base_size = 12))
set.seed(12345)  # reproducibility

# Create output folders
dir.create("plots", showWarnings = FALSE)
dir.create("plots/EAA_vs_Relatedness", showWarnings = FALSE)
dir.create("plots/entropy", showWarnings = FALSE)
dir.create("plots/pca", showWarnings = FALSE)
dir.create("plots/mediation", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

save_plot <- function(p, file, w = 7, h = 5, dpi = 600) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = dpi)
}

## 1) Load data ---------------------------------------------
beta <- read.csv("./data/sesame_data_BW.csv")



# Ensure CGids are unique
rownames(beta) <- make.unique(rownames(beta))

# Step 1: Convert `beta` to long format
beta_long <- as.data.frame(beta) %>%
  rownames_to_column("RowID") %>%         # Retain row names in a temporary column
  mutate(CGid = X) %>%                    # Assign `X` values to `CGid`
   dplyr :: select(-X) %>%                          # Remove the original `X` column
  pivot_longer(-c(RowID, CGid),           # Keep `RowID` and `CGid` as fixed columns
               names_to = "SID",
               values_to = "Beta") %>%
  dplyr :: select(CGid, SID, Beta)                 # Reorder columns for clarity


# Step 2: Merge with `key` to add metadata
beta_with_metadata <- beta_long %>%
  inner_join(key, by = c("SID" = "SID")) %>% distinct()




# Ensure Sex is a factor and Age is numeric
complete_data <- beta_with_metadata %>%
  mutate(
    Sex = as.factor(Sex),
    Age = as.numeric(Age)
  ) 



# Create beta matrix: rows = samples, columns = CpGs
beta_wide <- complete_data %>%
  dplyr :: select(ExternalSampleName, CGid, Beta) %>%
  distinct() %>%
  pivot_wider(names_from = CGid, values_from = Beta) %>%
  arrange(ExternalSampleName)  # Ensure rows are sorted by Sample ID




#Convert to matrix
beta_mat <- as.matrix(beta_wide %>% dplyr :: select(-ExternalSampleName))

# Store sample IDs
sample_ids <- beta_wide$ExternalSampleName

rownames(beta_mat) <- sample_ids


DNAm_Age <- read.xlsx("./data/DNAm_Age.xlsx") %>%
  mutate(
    Sex = factor(Sex),
    Relatedness = as.numeric(Relatedness),
    Age = as.numeric(Age),
    ExternalSampleName = as.character(ExternalSampleName)
  )

## 2) Compute EAA (epigenetic age acceleration) --------------
m_clock <- lm(DNAmAgeTailFinal ~ Age + Sex, data = DNAm_Age)
DNAm_Age <- DNAm_Age %>%
  mutate(EAA = resid(m_clock))

write.xlsx(broom::tidy(m_clock), "tables/model_clock_ols.xlsx", overwrite = TRUE)

## 3)Relatedness vs. EAA -------------------------
# Added-variable plot
res_x <- resid(lm(Relatedness ~ Age + Sex, data = DNAm_Age))
res_y <- DNAm_Age$EAA
df_add <- tibble(res_rel = res_x, res_eaa = res_y)

fit_add <- lm(res_eaa ~ res_rel, data = df_add)
summary(fit_add)

p_added <- ggplot(df_add, aes(res_rel, res_eaa)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(x = "Relatedness (residualized by Age+Sex)",
       y = "EAA (residualized by Age+Sex)",
       title = "Added-variable plot: Relatedness vs EAA") +
  theme_bw(base_size = 12)

save_plot(p_added, "plots/EAA_vs_Relatedness/added_variable.png")

# OLS, HC3 robust SEs, and Robust regression
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
write.xlsx(tabS1, "tables/Table_S_relatedness_EAA.xlsx", overwrite = TRUE)

# Leave-one-out analysis
DNAm_Age$._idx <- seq_len(nrow(DNAm_Age))
jack_results <- purrr::map_dfr(DNAm_Age$._idx, function(i) {
  df_loo <- df_add[-i, , drop = FALSE]
  fit_loo <- lm(res_eaa ~ res_rel, data = df_loo)
  tibble(
    left_out = DNAm_Age$ExternalSampleName[i],
    beta     = coef(fit_loo)[2],
    p        = coef(summary(fit_loo))["res_rel","Pr(>|t|)"]
  )
})
write.xlsx(jack_results, "tables/Table_S_LOO_relatedness_EAA.xlsx", overwrite = TRUE)

p_loo <- ggplot(jack_results, aes(x = reorder(left_out, beta), y = beta)) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = coef(fit_add)[2], linetype = 2, color = "steelblue") +
  coord_flip() +
  labs(title = "Leave-one-out β(Relatedness) for EAA",
       subtitle = "Dashed line = full-sample β",
       x = "Sample left out", y = "β") +
  theme_bw(base_size = 11)
save_plot(p_loo, "plots/EAA_vs_Relatedness/EAA_vs_Relatedness_LOO.png", w = 7.5, h = 10)

p_loo_dist <- ggplot(jack_results, aes(beta)) +
  geom_density(fill = "grey85", color = "grey25") +
  geom_vline(xintercept = coef(fit_add)[2], color = "steelblue", size = 1) +
  labs(title = "Distribution of leave-one-out β",
       subtitle = "Vertical line = full-sample β",
       x = "β", y = "Density") +
  theme_bw(base_size = 12)
save_plot(p_loo_dist, "plots/EAA_vs_Relatedness/LOO_density.png")

## 4) Entropy analysis ----------------------------

# Custom function to calculate entropy of beta values using KDE
compute_entropy_kde <- function(beta_values, n = 512) {
  beta_values <- na.omit(beta_values)
  if (length(beta_values) < 10) return(NA)  # Skip if too few CpGs
  
  # Estimate KDE of beta distribution
  dens <- density(beta_values, from = 0, to = 1, n = n)
  
  # Remove 0s to avoid log(0)
  p <- dens$y
  p[p == 0] <- 1e-10
  
  # Differential entropy: ∫ -p(x) log(p(x)) dx
  entropy <- -sum(p * log2(p)) * (dens$x[2] - dens$x[1])
  return(entropy)
}

# beta_mat: rows = samples, columns = CpGs (from previous step)
dim(beta_mat)  # Should be 96 x ~30,000+


# Compute entropy per sample
entropy_vec <- apply(beta_mat, 1, compute_entropy_kde)

# Store in your metadata
entropy_df <- data.frame(
  ExternalSampleName = sample_ids,
  entropy_kde = entropy_vec
) %>%
  mutate(ExternalSampleName = as.character(ExternalSampleName))



DNAm_Age <- DNAm_Age %>%
  left_join(entropy_df, by = "ExternalSampleName")



# Entropy vs Age
fit_ent_age <- lm(entropy_kde ~ Age, data = DNAm_Age)

p_ent_age <- ggplot(DNAm_Age, aes(Age, entropy_kde)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "Entropy (KDE) vs chronological age",
       subtitle = glue("β = {round(coef(fit_ent_age)[2],4)}, R² = {round(summary(fit_ent_age)$r.squared,3)}, p = {formatC(summary(fit_ent_age)$coefficients[2,4], format='e', digits=2)}"),
       x = "Age (years)", y = "Entropy (KDE)") +
  theme_bw(base_size = 12)

print(p_ent_age)

save_plot(p_ent_age, "plots/entropy/entropy_vs_age.png")



# EAA vs Entropy (partial: adjust Age+Sex+Relatedness)
res_EAA   <- resid(lm(EAA ~ Age+Sex+Relatedness, data = DNAm_Age))
res_Ent   <- resid(lm(entropy_kde ~ Age+Sex+Relatedness, data = DNAm_Age))
fit_eaa_ent <- lm(res_EAA ~ res_Ent)
p_eaa_ent <- ggplot(tibble(res_Ent, res_EAA), aes(res_Ent, res_EAA)) +
  geom_point(alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "EAA vs Entropy (both adjusted for Age+Sex+Relatedness)",
       subtitle = glue("β = {round(coef(fit_eaa_ent)[2],4)}, partial R² = {round(summary(fit_eaa_ent)$r.squared,3)}, p = {formatC(summary(fit_eaa_ent)$coefficients[2,4], format='e', digits=2)}"),
       x = "Entropy residual",
       y = "EAA residual") +
  theme_bw(base_size = 12)

print(p_eaa_ent)

save_plot(p_eaa_ent, "plots/entropy/EAA_vs_entropy_partial(adjust Age+Sex+Relatedness).png")

# Residual entropy vs residual relatedness (Age+Sex adjusted)
res_rel_AS <- resid(lm(Relatedness ~ Age + Sex, data = DNAm_Age))
res_ent_AS <- resid(lm(entropy_kde ~ Age + Sex, data = DNAm_Age))

fit_resres <- lm(res_ent_AS ~ res_rel_AS)
p_resres <- ggplot(tibble(res_rel_AS, res_ent_AS), aes(res_rel_AS, res_ent_AS)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Residual Entropy vs Residual Relatedness (both adjusted for Age+Sex)",
       subtitle = glue("β = {round(coef(fit_resres)[2],4)}, partial R² = {round(summary(fit_resres)$r.squared,3)}, p = {signif(summary(fit_resres)$coefficients[2,4],3)}"),
       x = "Relatedness residual",
       y = "Entropy residual") +
  theme_bw(base_size = 12)

print(p_resres)

save_plot(p_resres, "plots/entropy/residual_entropy_vs_residual_relatedness(Age+Sex adjusted).png")


## 5) PCA analysis ------------------------------

# Run PCA
pca_res <- prcomp(beta_mat, center = TRUE, scale. = TRUE)
pc_df <- as.data.frame( pca_res$x[, 1:5])  # Get PC1 to PC5
pc_df <- rownames_to_column(pc_df, var = "row_id")
pc_df <- pc_df %>%
  mutate(ExternalSampleName = sample_ids)


var_explained <- (pca_res$sdev)^2 / sum(pca_res$sdev^2)
df_scree <- tibble(PC = seq_along(var_explained), PropVar = var_explained)

p_scree <- ggplot(slice_head(df_scree, n = 10), aes(PC, PropVar)) +
  geom_point() + geom_line() +
  labs(title = "Scree plot") + theme_bw(base_size = 12)
save_plot(p_scree, "plots/pca/scree.png")

# Cumulative variance plot
p_cum <- ggplot(dplyr::slice_head(df_scree, n = k_show),
                aes(PC, cumsum(PropVar))) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.80, linetype = 2, color = "red") +
  labs(title = "Cumulative variance explained", y = "Cumulative proportion", x = "Principal component") +
  theme_bw(base_size = 12)

print(p_cum)
save_plot(p_cum, "plots/pca/cumulative_variance.png")


DNAm_Age <- left_join(DNAm_Age, pc_df, by = "ExternalSampleName")



# ---------- helpers ----------

save_plot <- function(p, path, width = 7, height = 5, dpi = 300) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)  # why: ensure dir exists
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("< 1e-4")
  sprintf("= %.3g", p)
}

is_discrete_like <- function(x) {
  is.factor(x) || is.character(x) || (is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) <= 10L)
}

pc_pvals <- function(df, predictor) {
  pred_quo  <- enquo(predictor)
  pred_name <- as_label(pred_quo)
  
  # select(): avoid tidy-eval inside tidyselect; be explicit
  cols <- c("PC1", "PC2", pred_name)
  d <- df %>%
    as.data.frame() %>%
    dplyr::select(all_of(cols)) %>%     # robust to column order & package conflicts
    filter(complete.cases(.))
  
  x <- d[[pred_name]]
  
  # per-PC p-values
  p_pc <- function(pc) {
    m <- lm(reformulate(termlabels = pred_name, response = pc), data = d)
    if (is.factor(x)) {
      an <- anova(m)
      as.numeric(an[1, "Pr(>F)"])      # why: single-term model → first row is predictor
    } else {
      tidy(m) %>% filter(term == pred_name) %>% pull(p.value) %>% as.numeric()
    }
  }
  
  p1 <- tryCatch(p_pc("PC1"), error = function(e) NA_real_)
  p2 <- tryCatch(p_pc("PC2"), error = function(e) NA_real_)
  
  # MANOVA (Pillai) — build formula via string; no tidy-eval in base
  p_mv <- NA_real_
  if (is.factor(x) && dplyr::n_distinct(x) > 1L) {
    f_mv <- as.formula(paste("cbind(PC1, PC2) ~", pred_name))
    mv   <- manova(f_mv, data = d)
    st   <- summary(mv, test = "Pillai")$stats
    p_mv <- suppressWarnings(as.numeric(st[1, "Pr(>F)"]))
  }
  
  list(p_pc1 = p1, p_pc2 = p2, p_manova = p_mv)
}

subtitle_for <- function(label, pvals) {
  paste0(
    label, ": ",
    "PC1: p", fmt_p(pvals$p_pc1), ", ",
    "PC2: p", fmt_p(pvals$p_pc2),
    ifelse(is.na(pvals$p_manova), "", paste0(" | MANOVA p", fmt_p(pvals$p_manova)))
  )
}

# ---------- input checks ----------
req_cols <- c("PC1", "PC2", "Relatedness", "Sex")
missing_cols <- setdiff(req_cols, names(DNAm_Age))


if (length(missing_cols)) {
  stop("Missing columns in DNAm_Age: ", paste(missing_cols, collapse = ", "))
}

DNAm_Age <- DNAm_Age %>%
  mutate(
    Relatedness = if (is_discrete_like(Relatedness)) factor(Relatedness) else Relatedness,
    Sex         = if (is_discrete_like(Sex))         factor(Sex)         else Sex
  )

# ---------- stats ----------
rel_p <- pc_pvals(DNAm_Age, Relatedness)
sex_p <- pc_pvals(DNAm_Age, Sex)

rel_sub <- subtitle_for("Relatedness", rel_p)
sex_sub <- subtitle_for("Sex", sex_p)

# ---------- plots ----------
p_pca_rel <- ggplot(DNAm_Age, aes(PC1, PC2, color = Relatedness)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(
    title = "PCA colored by Relatedness",
    subtitle = rel_sub,
    x = "PC1", y = "PC2", color = "Relatedness"
  ) +
  { if (is.factor(DNAm_Age$Relatedness)) scale_color_viridis_d() else scale_color_viridis_c() } +
  theme_bw(base_size = 12)

print(p_pca_rel)

save_plot(p_pca_rel, "plots/pca/pca_relatedness.png")

p_pca_sex <- ggplot(DNAm_Age, aes(PC1, PC2, color = Sex)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(
    title = "PCA of methylation: colored by Sex",
    subtitle = sex_sub,
    x = "PC1", y = "PC2", color = "Sex"
  ) +
  { if (is.factor(DNAm_Age$Sex)) scale_color_viridis_d() else scale_color_viridis_c() } +
  theme_bw(base_size = 12)

print(p_pca_sex)
save_plot(p_pca_sex, "plots/pca/pca_PC1_PC2_by_sex.png")

message("Relatedness — ", rel_sub)
message("Sex — ", sex_sub)


# Residual Age (EAA) ~ PCs
DNAm_Age <- DNAm_Age %>%
  mutate(resid_age = resid(lm(DNAmAgeTailFinal ~ Age + Sex, data = .)))


# PC1 vs residual DNAmAge
p_res_PC1 <- ggplot(DNAm_Age, aes(PC1, resid_age)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top", # position in plot
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  labs(title = "Residual DNAmAge ~ PC1",
       x = "PC1", y = "Residual DNAmAge") +
  theme_bw(base_size = 12)

print(p_res_PC1)


# PC2 vs residual DNAmAge
p_res_PC2 <- ggplot(DNAm_Age, aes(PC2, resid_age)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top",
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  labs(title = "Residual DNAmAge ~ PC2",
       x = "PC2", y = "Residual DNAmAge") +
  theme_bw(base_size = 12)

print(p_res_PC2)

# Save plots
ggsave("plots/pca/residAge_PC1_with_stats.png", p_res_PC1, width = 6.5, height = 4.5, dpi = 300)
ggsave("plots/pca/residAge_PC2_with_stats.png", p_res_PC2, width = 6.5, height = 4.5, dpi = 300)


## 6) Figure 5 – Mediation ----------------------------------
med_model1 <- lm(PC1 ~ Relatedness + Age + Sex, data = DNAm_Age)
outcome_model1 <- lm(EAA ~ PC1 + Relatedness + Age + Sex, data = DNAm_Age)
med_out1 <- mediate(med_model1, outcome_model1,
                    treat = "Relatedness", mediator = "PC1",
                    boot = TRUE, sims = 5000)
summary(med_out1)
capture.output(summary(med_out1), file = "tables/mediation_PC1.txt")

med_model2 <- lm(PC2 ~ Relatedness + Age + Sex, data = DNAm_Age)
outcome_model2 <- lm(EAA ~ PC2 + Relatedness + Age + Sex, data = DNAm_Age)
med_out2 <- mediate(med_model2, outcome_model2,
                    treat = "Relatedness", mediator = "PC2",
                    boot = TRUE, sims = 5000)
summary(med_out2)
capture.output(summary(med_out2), file = "tables/mediation_PC2.txt")

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

## 7) Sex-stratified mediation (Supplementary) ---------------
set.seed(123)  # reproducible bootstrap

# ---------- input checks ----------
req <- c("Sex","Relatedness","Age","EAA","PC1","PC2")
miss <- setdiff(req, names(DNAm_Age))
if (length(miss)) stop("Missing columns in DNAm_Age: ", paste(miss, collapse = ", "))
if (!is.factor(DNAm_Age$Sex)) DNAm_Age$Sex <- factor(DNAm_Age$Sex)
sex_levels <- levels(DNAm_Age$Sex)
med_levels <- c("PC1","PC2")

# ---------- helpers ----------
fmt_num <- function(x, d = 5) sprintf(paste0("%.", d, "f"), x)
fmt_p   <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, "<1e-4", sprintf("%.3g", p)))
fmt_p_vec <- function(p) paste0("p ", ifelse(is.na(p), "NA", ifelse(p < 1e-4, "<1e-4", sprintf("%.3g", p))))
stars   <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

# Why: avoid NSE capture in bootstrap by using a fixed mediator name.
run_mediation_cell <- function(dat, mediator, sims = 5000, min_n = 10L) {
  stopifnot(is.character(mediator), length(mediator) == 1L)
  sub <- dat %>%
    mutate(MED = .data[[mediator]]) %>%
    filter(complete.cases(Relatedness, Age, EAA, MED))
  n_cell <- nrow(sub)
  if (n_cell < min_n || isTRUE(var(sub$MED, na.rm = TRUE) == 0)) {
    return(tibble(Mediator = mediator, n = n_cell,
                  term = c("ACME","ADE","Total"),
                  estimate = NA_real_, lower = NA_real_, upper = NA_real_, p_value = NA_real_))
  }
  med_m <- lm(MED ~ Relatedness + Age, data = sub)
  out_m <- lm(EAA ~ MED + Relatedness + Age, data = sub)
  res <- tryCatch({
    mo <- mediation::mediate(model.m = med_m, model.y = out_m,
                             treat = "Relatedness", mediator = "MED",
                             boot = TRUE, sims = sims)
    sm <- summary(mo)
    tibble(
      Mediator = mediator,
      n        = n_cell,
      term     = c("ACME","ADE","Total"),
      estimate = c(sm$d0, sm$z0, sm$tau.coef),
      lower    = c(sm$d0.ci[1], sm$z0.ci[1], sm$tau.ci[1]),
      upper    = c(sm$d0.ci[2], sm$z0.ci[2], sm$tau.ci[2]),
      p_value  = c(sm$d0.p, sm$z0.p, sm$tau.p)
    )
  }, error = function(e) {
    warning(sprintf("Mediation failed for mediator=%s: %s", mediator, conditionMessage(e)))
    tibble(Mediator = mediator, n = n_cell,
           term = c("ACME","ADE","Total"),
           estimate = NA_real_, lower = NA_real_, upper = NA_real_, p_value = NA_real_)
  })
  res
}

# ---------- collect results: long format (Sex × Mediator × {ACME,ADE,Total}) ----------
res_long <- map_dfr(sex_levels, function(sx) {
  sub <- DNAm_Age %>% filter(Sex == sx)
  map_dfr(med_levels, ~run_mediation_cell(sub, .x)) %>%
    mutate(Sex = sx, .before = 1)
})

# ---------- enrich: q-values (by sex), stars, % mediated from SAME run ----------
res_long <- res_long %>%
  group_by(Sex) %>%                                       # sex-specific BH control
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    stars = case_when(is.na(p_value) ~ "",
                      p_value < .001 ~ "***",
                      p_value < .01  ~ "**",
                      p_value < .05  ~ "*",
                      TRUE ~ "")
  )

# --- build percent_mediated FROM THE SAME RUN (res_long must have ACME, Total) ---
prop_med <- res_long %>%
  dplyr::select(Sex, Mediator, term, estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
  dplyr::mutate(
    percent_mediated = dplyr::case_when(
      is.finite(ACME) & is.finite(Total) & Total != 0 ~ 100 * ACME / Total,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(Sex, Mediator, percent_mediated)

# --- table: join percent_mediated, then format/widen ---
table_supp <- res_long %>%
  dplyr::left_join(prop_med, by = c("Sex","Mediator")) %>%     # <-- ensures the column exists
  dplyr::mutate(
    estimate  = rnd(estimate, 6),
    lower     = rnd(lower, 6),
    upper     = rnd(upper, 6),
    p_value   = rnd(p_value, 4),
    q_value   = rnd(q_value, 4),
    percent_mediated = rnd(percent_mediated, 1)                # now found
  ) %>%
  dplyr::mutate(term = factor(term, levels = c("ACME","ADE","Total"))) %>%
  dplyr::arrange(Sex, Mediator, term) %>%
  dplyr::mutate(col_prefix = paste0(term, "_")) %>%
  tidyr::pivot_wider(
    id_cols = c(Sex, Mediator, n, percent_mediated),
    names_from = col_prefix,
    values_from = c(estimate, lower, upper, p_value, q_value, stars),
    names_vary = "slowest"
  ) %>%
  dplyr::rename(
    `ACME_est`  = `estimate_ACME_`, `ACME_CI_lo`= `lower_ACME_`,  `ACME_CI_hi`= `upper_ACME_`,
    `ACME_p`    = `p_value_ACME_`,  `ACME_q`    = `q_value_ACME_`, `ACME_sig`  = `stars_ACME_`,
    `ADE_est`   = `estimate_ADE_`,  `ADE_CI_lo` = `lower_ADE_`,   `ADE_CI_hi` = `upper_ADE_`,
    `ADE_p`     = `p_value_ADE_`,   `ADE_q`     = `q_value_ADE_`,  `ADE_sig`   = `stars_ADE_`,
    `Total_est` = `estimate_Total_`,`Total_CI_lo`=`lower_Total_`,  `Total_CI_hi`=`upper_Total_`,
    `Total_p`   = `p_value_Total_`, `Total_q`   = `q_value_Total_`,`Total_sig` = `stars_Total_`,
    `%_mediated`= `percent_mediated`
  ) %>%
  dplyr::select(
    Sex, Mediator, n,
    ACME_est, ACME_CI_lo, ACME_CI_hi, ACME_p, ACME_q, ACME_sig,
    ADE_est,  ADE_CI_lo,  ADE_CI_hi,  ADE_p,  ADE_q,  ADE_sig,
    Total_est, Total_CI_lo, Total_CI_hi, Total_p, Total_q, Total_sig,
    `%_mediated`
  )



writexl::write_xlsx(list(`Table S2 — Sex-stratified mediation` = table_supp),
                    path = "tables/Table_S_sex_stratified_mediation.xlsx")
write.csv(table_supp, "tables/Table_S_sex_stratified_mediation.csv", row.names = FALSE)


# ---------- plot with significance highlighting and labels ----------
# --- (re)build percent_mediated FROM THE SAME RUN (res_long must have ACME & Total) ---
prop_med <- res_long %>%
  dplyr::select(Sex, Mediator, term, estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
  dplyr::mutate(
    percent_mediated = dplyr::case_when(
      is.finite(ACME) & is.finite(Total) & Total != 0 ~ 100 * ACME / Total,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(Sex, Mediator, percent_mediated)

# --- build plotting df: JOIN percent_mediated BEFORE mutate() ---
res_plot <- res_long %>%
  dplyr::left_join(prop_med, by = c("Sex","Mediator")) %>%   # <-- key fix
  dplyr::mutate(
    sig     = !is.na(p_value) & p_value < 0.05,
    est_ci  = paste0(sprintf("%.4f", estimate), " [", sprintf("%.4f", lower), ", ", sprintf("%.4f", upper), "]"),
    label   = paste0(est_ci, "\n", fmt_p_vec(p_value), " | q ", fmt_p(q_value), " ", stars(p_value)),
    term    = factor(term, levels = c("ACME","ADE","Total")),
    panel_lab = glue::glue("{Sex} (n={n}) — {Mediator}\n% mediated ≈ {ifelse(is.na(percent_mediated), 'NA', paste0(round(percent_mediated,1),'%'))}")
  )

# --- plot ---
pal <- c(`TRUE` = "#b22222", `FALSE` = "grey40")
p_med_sig <- ggplot2::ggplot(res_plot, aes(x = term, y = estimate, color = as.character(sig))) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  ggplot2::geom_point(size = 2.6, alpha = 0.95) +
  ggplot2::geom_text(aes(label = label), vjust = -1.1, size = 3.1, lineheight = 0.95, na.rm = TRUE, show.legend = FALSE) +
  ggplot2::facet_wrap(~ panel_lab, ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = pal, breaks = c("TRUE","FALSE"), labels = c("Significant","Not significant")) +
  ggplot2::labs(
    title    = "Sex-stratified mediation: Relatedness → PC → EAA",
    subtitle = "ACME (indirect), ADE (direct), Total (from same mediation run) with 95% CI; labels show estimate [CI], p | q (BH) and stars",
    x = NULL, y = "Effect size", color = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(strip.text = element_text(face = "bold"),
                 panel.grid.minor = element_blank(),
                 legend.position = "bottom")

print(p_med_sig)

ggplot2::ggsave("plots/mediation/mediation_sex_facet_sig.png", p_med_sig, width = 10, height = 6.2, dpi = 600)


# ---------- auto-generate manuscript paragraph + caption ----------
pick <- function(df, sex, med, path) df %>% filter(Sex == sex, Mediator == med, term == path) %>% dplyr::slice(1)
sx_levels <- unique(res_plot$Sex)

f_PC1_ACME <- pick(res_plot, sx_levels[1], "PC1", "ACME")
f_PC1_ADE  <- pick(res_plot, sx_levels[1], "PC1", "ADE")
f_PC1_TOT  <- pick(res_plot, sx_levels[1], "PC1", "Total")
f_PC2_ACME <- pick(res_plot, sx_levels[1], "PC2", "ACME")

m_PC1_ACME <- pick(res_plot, sx_levels[2], "PC1", "ACME")
m_PC1_ADE  <- pick(res_plot, sx_levels[2], "PC1", "ADE")
m_PC1_TOT  <- pick(res_plot, sx_levels[2], "PC1", "Total")
m_PC2_ACME <- pick(res_plot, sx_levels[2], "PC2", "ACME")

para <- glue(
  "Sex-stratified mediation analyses revealed distinct patterns. In {sx_levels[1]} (n = {f_PC1_ACME$n}), ",
  "PC1 showed an indirect-only effect (ACME = {fmt_num(f_PC1_ACME$estimate, 5)}, ",
  "95% CI [{fmt_num(f_PC1_ACME$lower, 5)}, {fmt_num(f_PC1_ACME$upper, 5)}], p = {fmt_p(f_PC1_ACME$p_value)}), ",
  "while the direct (ADE = {fmt_num(f_PC1_ADE$estimate, 5)}, p = {fmt_p(f_PC1_ADE$p_value)}) ",
  "and total effects based on the same model (Total = {fmt_num(f_PC1_TOT$estimate, 5)}, p = {fmt_p(f_PC1_TOT$p_value)}) were not significant; ",
  "PC2 showed no significant mediation (ACME p = {fmt_p(f_PC2_ACME$p_value)}). ",
  "In {sx_levels[2]} (n = {m_PC1_ADE$n}), PC1 exhibited a small direct effect ",
  "(ADE = {fmt_num(m_PC1_ADE$estimate, 5)}, 95% CI [{fmt_num(m_PC1_ADE$lower, 5)}, {fmt_num(m_PC1_ADE$upper, 5)}], ",
  "p = {fmt_p(m_PC1_ADE$p_value)}), with negligible mediation (ACME p = {fmt_p(m_PC1_ACME$p_value)}); ",
  "by contrast, PC2 mediated the association (ACME = {fmt_num(m_PC2_ACME$estimate, 5)}, ",
  "95% CI [{fmt_num(m_PC2_ACME$lower, 5)}, {fmt_num(m_PC2_ACME$upper, 5)}], p = {fmt_p(m_PC2_ACME$p_value)}). ",
  "Benjamini–Hochberg q-values (by sex) appear in figure labels."
)

caption <- glue(
  "Sex-stratified mediation of relatedness → EAA via PC1/PC2. ",
  "Points = ACME (indirect), ADE (direct), Total (all from the same mediation run); bars = 95% CI. ",
  "Labels show estimate [CI], raw p | BH q (by sex) and significance (*** <0.001, ** <0.01, * <0.05)."
)

dir.create("text", showWarnings = FALSE)
writeLines(para,    "text/sex_stratified_mediation_paragraph.txt")
writeLines(caption, "text/sex_stratified_mediation_caption.txt")
cat("\n--- Manuscript paragraph ---\n", para, "\n\n")
cat("--- Figure caption ---\n", caption, "\n")



############################################################
# End of pipeline
############################################################
