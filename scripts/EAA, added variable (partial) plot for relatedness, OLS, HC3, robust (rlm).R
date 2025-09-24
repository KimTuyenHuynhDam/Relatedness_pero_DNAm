# --- Packages ---
# install.packages(c("broom","lmtest","sandwich","MASS","dplyr","writexl"))
library(broom)
library(lmtest)
library(sandwich)
library(MASS)
library(dplyr)
library(writexl)

# EAA–Relatedness analysis: end‑to‑end manuscript‑ready figures
# Author: <your name>
# Date: <auto>
# Description: Computes EAA, runs nested models, robust/sensitivity analyses,
#              and generates publication‑quality figures with in‑plot statistics.

# =============================
# 0) Packages & Theme
# =============================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, broom, ggplot2, viridis, ggpubr, ggrepel,
  car, lmtest, sandwich, MASS, rsq, effectsize,
  lme4, performance, modelsummary
)

theme_set(theme_minimal(base_size = 12))
update_geom_defaults("point", list(alpha = 0.85))

# Utility: save function (PDF and PNG at hi‑res)
save_fig <- function(p, filename, width = 6, height = 4.5, dpi = 600) {
  dir.create("figures", showWarnings = FALSE)
  ggsave(file.path("figures", paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path("figures", paste0(filename, ".png")), p, width = width, height = height, dpi = dpi)
}

# =============================
# 1) Data prep
# =============================
# Expect a data.frame named DNAm_Age with columns:
#   DNAmAgeTailFinal (numeric), Age (numeric), Sex (factor/char), Relatedness (numeric)
# Optional: Cohort (factor) for mixed models, SampleID (character) for labeling
DNAm_Age = read.xlsx("./data/DNAm_Age.xlsx")

DNAm_Age <- DNAm_Age %>%
  rename(DNAmAge = DNAmAgeTailFinal) %>%
  mutate(
    Sex = as.factor(Sex),
    Relatedness = as.numeric(Relatedness),
    Age = as.numeric(Age)
  ) %>%
  drop_na(DNAmAge, Age, Sex, Relatedness)

# Force SampleID to be ExternalSampleName for labeling
if (!"ExternalSampleName" %in% names(DNAm_Age)) {
  stop("DNAm_Age must contain 'ExternalSampleName'.")
}
DNAm_Age <- DNAm_Age %>% mutate(SampleID = as.character(ExternalSampleName))

# Use ExternalSampleName as SampleID
DNAm_Age$SampleID <- as.character(DNAm_Age$ExternalSampleName)

# =============================
# 2) Compute EAA (adjust for Age + Sex)
# =============================
m0 <- lm(DNAmAge ~ Age + Sex, data = DNAm_Age)          # base model used to compute EAA
DNAm_Age <- DNAm_Age %>% mutate(EAA = resid(m0))

# =============================
# 3) Test relatedness effect
# =============================
# 3a) Nested models on DNAmAge (does Relatedness add beyond Age + Sex?)
m1 <- lm(DNAmAge ~ Age + Sex + Relatedness, data = DNAm_Age)
nested <- anova(m0, m1)

# Partial R^2 uniquely explained by Relatedness
# (rsq.partial returns partial R^2 for each predictor in the full model)
partial_r2 <- tryCatch({ rsq::rsq.partial(m1, adj = FALSE) }, error = function(e) NULL)

# 3b) EAA ~ Relatedness (OLS and robust)
ml_eaa <- lm(EAA ~ Relatedness, data = DNAm_Age)
coefs_eaa <- broom::tidy(ml_eaa)
# Robust SEs
rob_eaa <- lmtest::coeftest(ml_eaa, vcov = sandwich::vcovHC(ml_eaa, type = "HC3"))
rob_eaa_df <- broom::tidy(rob_eaa) %>% rename(estimate = estimate, std.error = std.error, statistic = statistic, p.value = p.value)

# Robust regression (M‑estimator)
rlm_eaa <- MASS::rlm(EAA ~ Relatedness, data = DNAm_Age)
rlm_tidy <- broom::tidy(rlm_eaa)

# Nonlinearity check (quadratic)
ml_eaa_lin  <- lm(EAA ~ Relatedness, data = DNAm_Age)
ml_eaa_quad <- lm(EAA ~ poly(Relatedness, 2), data = DNAm_Age)
poly_cmp <- anova(ml_eaa_lin, ml_eaa_quad)

# Influence diagnostics
cooks_cut <- 4 / nrow(DNAm_Age)
inf_tbl <- tibble(
  SampleID   = DNAm_Age$SampleID,
  cooksd     = cooks.distance(ml_eaa),
  leverage   = hatvalues(ml_eaa),
  flagged    = cooksd > cooks_cut | leverage > 2*mean(leverage)
)

# =============================
# 4) Mixed model (optional if Cohort exists)
# =============================
mm_res <- NULL
if ("Cohort" %in% names(DNAm_Age)) {
  DNAm_Age <- DNAm_Age %>% mutate(Cohort = as.factor(Cohort))
  m0_lmm <- lmer(DNAm_Age ~ Age + Sex + (1|Cohort), data = DNAm_Age, REML = FALSE)
  m1_lmm <- lmer(DNAm_Age ~ Age + Sex + Relatedness + (1|Cohort), data = DNAm_Age, REML = FALSE)
  mm_cmp <- anova(m0_lmm, m1_lmm)
  mm_r2  <- performance::r2(m1_lmm)
  mm_res <- list(mm_cmp = mm_cmp, mm_r2 = mm_r2)
}

# =============================
# 5) Print key stats to console (for the record)
# =============================
cat("\n==== NESTED MODELS (DNAmAge) ====\n"); print(nested)
cat("\nPartial R^2 (full model)\n"); if (!is.null(partial_r2)) print(partial_r2)
cat("\n==== EAA ~ Relatedness (OLS) ====\n"); print(summary(ml_eaa))
cat("\n==== EAA ~ Relatedness (Robust SE) ====\n"); print(rob_eaa)
cat("\n==== EAA ~ Relatedness (rlm) ====\n"); print(summary(rlm_eaa))
cat("\n==== Linearity check (EAA) ====\n"); print(poly_cmp)
cat("\n==== Influence (top flagged) ====\n"); print(inf_tbl %>% filter(flagged) %>% arrange(desc(cooksd)) %>% head(10))
if (!is.null(mm_res)) { cat("\n==== Mixed model (optional) ====\n"); print(mm_res$mm_cmp); print(mm_res$mm_r2) }

# =============================
# 6) Figures (manuscript‑ready)
# =============================
# Helper: build annotation strings
fmt_p <- function(p) ifelse(p < 2.2e-16, "p < 2e-16", paste0("p = ", formatC(p, format = "e", digits = 2)))

# Extract stats
anova_p    <- nested$`Pr(>F)`[2]
r2_m0      <- summary(m0)$r.squared
r2_m1      <- summary(m1)$r.squared
beta_rel   <- broom::tidy(m1) %>% filter(term == "Relatedness") %>% pull(estimate)
p_rel_m1   <- broom::tidy(m1) %>% filter(term == "Relatedness") %>% pull(p.value)
r2_eaa     <- summary(ml_eaa)$r.squared
beta_eaa   <- coef(ml_eaa)["Relatedness"]
p_eaa      <- broom::tidy(ml_eaa) %>% filter(term == "Relatedness") %>% pull(p.value)
p_eaa_rob  <- rob_eaa_df %>% filter(term == "Relatedness") %>% pull(p.value)
beta_rlm   <- rlm_tidy %>% filter(term == "Relatedness") %>% pull(estimate)

# (Fig 1) DNAm Age vs Chronological Age (colored by Relatedness)
fig1 <- ggplot(DNAm_Age, aes(x = Age, y = DNAmAge, color = Relatedness)) +
  geom_point(size = 2.6) +
  geom_smooth(method = stats::lm, method.args = list(), se = TRUE, color = "black", formula = y ~ x) +
  scale_color_viridis(option = "C", end = 0.95) +
  labs(x = "Chronological Age", y = "DNAm Age",
       title = "DNAm Age vs Chronological Age",
       color = "Relatedness") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.6,
           label = paste0(
             "Model: DNAmAge ~ Age + Sex\n",
             "R^2 = ", sprintf("%.3f", r2_m0), "\n",
             "EAA = residuals of this model")
  )

print(fig1)
save_fig(fig1, "Fig1_DNAmAge_vs_Age")

# (Fig 2) EAA vs Relatedness with OLS and robust slope annotation
fig2 <- ggplot(DNAm_Age, aes(x = Relatedness, y = EAA)) +
  geom_point(size = 2.6) +
  geom_smooth(method = stats::lm, method.args = list(), se = TRUE, formula = y ~ x) +
  labs(x = "Parental Relatedness", y = "Epigenetic Age Acceleration (EAA)",
       title = "Higher Relatedness associates with lower EAA") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.6,
           label = paste0(
             "OLS: slope = ", sprintf("%.4f", beta_eaa), ", R^2 = ", sprintf("%.3f", r2_eaa), ", ", fmt_p(p_eaa), "\n",
             "Robust (HC3): ", fmt_p(p_eaa_rob), "; rlm slope = ", sprintf("%.4f", beta_rlm)
           ))


print(fig2)
save_fig(fig2, "Fig2_EAA_vs_Relatedness")

# (Fig 3) Added‑variable (partial) plot for Relatedness controlling Age + Sex
# Compute residuals
res_y <- resid(lm(DNAmAge ~ Age + Sex, data = DNAm_Age))
res_x <- resid(lm(Relatedness ~ Age + Sex, data = DNAm_Age))
partial_df <- tibble(res_y = res_y, res_x = res_x)
partial_fit <- lm(res_y ~ res_x)
partial_p   <- summary(partial_fit)$coefficients[2,4]
partial_b   <- coef(partial_fit)[2]
partial_r2  <- summary(partial_fit)$r.squared

fig3 <- ggplot(partial_df, aes(x = res_x, y = res_y)) +
  geom_point(size = 2.6) +
  geom_smooth(method = stats::lm, method.args = list(), se = TRUE) +
  labs(x = "Relatedness (residualized by Age + Sex)",
       y = "DNAm Age (residualized by Age + Sex)",
       title = "Added variable plot: unique effect of Relatedness") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.6,
           label = paste0(
             "slope = ", sprintf("%.4f", partial_b), ", R^2 = ", sprintf("%.3f", partial_r2), ", ", fmt_p(partial_p)
           ))


print(fig3)
save_fig(fig3, "Fig3_AddedVariable_Relatedness")

# (Fig 4) Influence (Cook's distance vs leverage) with labels for flagged points
infl_df <- inf_tbl %>% left_join(DNAm_Age %>% dplyr :: select(SampleID, Relatedness, EAA), by = "SampleID")
fig4 <- ggplot(infl_df, aes(x = leverage, y = cooksd, label = ifelse(flagged, SampleID, ""))) +
  geom_point() +
  geom_hline(yintercept = cooks_cut, linetype = "dashed") +
  geom_vline(xintercept = 2*mean(infl_df$leverage), linetype = "dashed") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 12) +
  labs(x = "Leverage (hat)", y = "Cook's distance",
       title = "Influence diagnostics for EAA ~ Relatedness") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.6,
           label = paste0("Dashed lines: Cook's cutoff = ", sprintf("%.3f", cooks_cut), "; leverage ≈ 2×mean"))

print(fig4)
save_fig(fig4, "Fig4_Influence_EAA")

# =============================
# 7) Model table (manuscript‑friendly)
# =============================
models <- list(
  "DNAmAge ~ Age + Sex" = m0,
  "DNAmAge ~ Age + Sex + Relatedness" = m1,
  "EAA ~ Relatedness (OLS)" = ml_eaa
)
msummary(models,
         statistic = c("std.error", "statistic", "p.value"),
         gof_omit = "IC|F|RMSE",
         stars = TRUE,
         output = file.path("figures", "Model_Table.html"))

cat("\nAll figures saved to ./figures and model table to ./figures/Model_Table.html\n")




# ============== TABLE S1: OLS vs HC3 vs Robust (rlm) ==============
ml_eaa  <- lm(EAA ~ Relatedness, data = DNAm_Age)
rob_eaa <- coeftest(ml_eaa, vcov = vcovHC(ml_eaa, type = "HC3"))
rlm_eaa <- rlm(EAA ~ Relatedness, data = DNAm_Age)

n_all  <- nobs(ml_eaa)
r2_all <- summary(ml_eaa)$r.squared

# tidy rows for the Relatedness coefficient
tbl_ols <- tidy(ml_eaa) |> filter(term == "Relatedness") |>
  transmute(Method = "OLS",
            Estimate = estimate, SE = std.error, t = statistic, p = p.value,
            N = n_all, R2 = r2_all)

tbl_hc3 <- tidy(rob_eaa) |> filter(term == "Relatedness") |>
  transmute(Method = "OLS + HC3 SE",
            Estimate = estimate, SE = std.error, t = statistic, p = p.value,
            N = n_all, R2 = r2_all)

# MASS::rlm doesn't return p-values; approximate from t with df = n - 2 (simple linear model)
rlm_row <- tidy(rlm_eaa) |> filter(term == "Relatedness")
t_rlm   <- rlm_row$statistic
p_rlm   <- 2 * pt(abs(t_rlm), df = n_all - 2, lower.tail = FALSE)

tbl_rlm <- rlm_row |>
  transmute(Method = "Robust (rlm)",
            Estimate = estimate, SE = std.error, t = statistic, p = p_rlm,
            N = n_all, R2 = r2_all)

table_S1 <- bind_rows(tbl_ols, tbl_hc3, tbl_rlm)

# ============== TABLE S2: Sensitivity (ALL vs EXCLUDING OUTLIERS) ==============
# Replace with your actual ExternalSampleName values:
outlier_ids <- c("38961","38962")

DNAm_all   <- DNAm_Age
DNAm_trim  <- DNAm_Age |> filter(!ExternalSampleName %in% outlier_ids)

make_rows <- function(df, label) {
  fit_ols  <- lm(EAA ~ Relatedness, data = df)
  fit_hc3  <- coeftest(fit_ols, vcov = vcovHC(fit_ols, type = "HC3"))
  fit_rlm  <- rlm(EAA ~ Relatedness, data = df)
  
  n  <- nobs(fit_ols)
  r2 <- summary(fit_ols)$r.squared
  
  r1 <- tidy(fit_ols) |> filter(term == "Relatedness") |>
    transmute(Dataset = label, Method = "OLS",
              Estimate = estimate, SE = std.error, t = statistic, p = p.value,
              N = n, R2 = r2)
  
  r2_ <- tidy(fit_hc3) |> filter(term == "Relatedness") |>
    transmute(Dataset = label, Method = "OLS + HC3 SE",
              Estimate = estimate, SE = std.error, t = statistic, p = p.value,
              N = n, R2 = r2)
  
  rlm_row <- tidy(fit_rlm) |> filter(term == "Relatedness")
  t_rlm   <- rlm_row$statistic
  p_rlm   <- 2 * pt(abs(t_rlm), df = n - 2, lower.tail = FALSE)
  
  r3 <- rlm_row |>
    transmute(Dataset = label, Method = "Robust (rlm)",
              Estimate = estimate, SE = std.error, t = statistic, p = p_rlm,
              N = n, R2 = r2)
  
  bind_rows(r1, r2_, r3)
}

table_S2 <- bind_rows(
  make_rows(DNAm_all,  "All samples"),
  make_rows(DNAm_trim, "Excluding 38961 & 38962")
)

# Optional: nicer rounding for presentation
round_df <- function(x, digits = 6) {
  nums <- sapply(x, is.numeric)
  x[nums] <- lapply(x[nums], function(col) signif(col, digits))
  x
}
table_S1 <- round_df(table_S1, 6)
table_S2 <- round_df(table_S2, 6)

# ============== WRITE TO EXCEL (two sheets) ==============
writexl::write_xlsx(
  list(
    "Table_S1_Robustness" = table_S1,
    "Table_S2_Sensitivity" = table_S2
  ),
  path = "figures/Tables_S1_S2.xlsx"
)

message("Saved: figures/Tables_S1_S2.xlsx")
