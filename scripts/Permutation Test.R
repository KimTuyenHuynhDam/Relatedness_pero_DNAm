
library(tidyverse)
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
t_obs <- summary(m_lin)$coef["Relatedness","t value"]

set.seed(123)
B <- 5000
t_perm <- replicate(B, {
  perm <- DNAm_Age
  perm$EAA <- sample(perm$EAA)
  summary(lm(EAA ~ Relatedness, data = perm))$coef["Relatedness","t value"]
})
p_perm <- mean(abs(t_perm) >= abs(t_obs))

write.xlsx(list("Permutation"=tibble(t_obs=t_obs, p_perm=p_perm)), "tables/Permutation.xlsx")
