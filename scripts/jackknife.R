
library(tidyverse)
library(broom)
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

jack <- sapply(1:nrow(DNAm_Age), function(i) {
  coef(lm(EAA ~ Relatedness, data = DNAm_Age[-i, ]))["Relatedness"]
})
jack_df <- tibble(Sample=DNAm_Age$ExternalSampleName, beta=jack)

write.xlsx(list("Jackknife_Beta"=jack_df), "tables/Jackknife_Beta.xlsx")

p <- ggplot(jack_df, aes(reorder(Sample,beta), beta)) +
  geom_point() + coord_flip() +
  geom_hline(yintercept=coef(lm(EAA ~ Relatedness, data = DNAm_Age))["Relatedness"], lty=2) +
  labs(title="Leave-one-out β(Relatedness)", x="Sample left out", y="β") +
  theme_minimal()
ggsave("figures/Fig_Jackknife.png", p, width=6, height=6, dpi=600)
