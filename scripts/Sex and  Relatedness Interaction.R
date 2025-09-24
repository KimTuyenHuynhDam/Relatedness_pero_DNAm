
library(tidyverse); library(broom); library(openxlsx)


DNAm_Age = read.xlsx("./data/DNAm_Age.xlsx")

DNAm_Age <- DNAm_Age %>%
  mutate(
    Sex = as.factor(Sex),
    Relatedness = as.numeric(Relatedness),
    Age = as.numeric(Age)
  )


m0 <- lm(DNAmAgeTailFinal ~ Age + Sex, data = DNAm_Age)

DNAm_Age$EAA <- resid(m0)


m_int <- lm(EAA ~ Age + Sex*Relatedness, data = DNAm_Age)
tab <- tidy(m_int)

write.xlsx(list("Interaction_Sex_Relatedness"=tab), "tables/Interaction.xlsx")


