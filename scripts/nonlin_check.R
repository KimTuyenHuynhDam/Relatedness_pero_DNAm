
library(tidyverse)
library(splines) 
library(broom)
library(writexl)
library(openxlsx)


DNAm_Age = read.xlsx("./data/DNAm_Age.xlsx")

DNAm_Age <- DNAm_Age %>%
  mutate(
    Sex = as.factor(Sex),
    Relatedness = as.numeric(Relatedness),
    Age = as.numeric(Age)
  )


  m0 <- lm(DNAmAgeTailFinal ~ Age + Sex, data = DNAm_Age)
  DNAm_Age$EAA <- resid(m0)


m_lin <- lm(EAA ~ Relatedness, data = DNAm_Age)
m_ns  <- lm(EAA ~ ns(Relatedness, df=3), data = DNAm_Age)
nl_cmp <- anova(m_lin, m_ns)

# Save table
write_xlsx(list("Nonlinearity_ANOVA"=tidy(nl_cmp)), "tables/Nonlinearity_NS.xlsx")

# Plot spline shape
p <- ggplot(DNAm_Age, aes(Relatedness, EAA)) +
  geom_point() +
  geom_smooth(method=lm, formula = y ~ splines::ns(x, 3), se=TRUE) +
  labs(title="Spline fit (df=3)", subtitle=paste0("ANOVA p=", signif(nl_cmp$`Pr(>F)`[2],3))) +
  theme_minimal()
ggsave("figures/Fig_NS_shape.png", p, width=6, height=5, dpi=600)
