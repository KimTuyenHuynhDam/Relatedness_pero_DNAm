
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



DNAm_Age <- DNAm_Age %>%
  mutate(RelGroup = cut(Relatedness, quantile(Relatedness, c(0,.33,.66,1)), include.lowest=TRUE))

by_ter <- DNAm_Age %>%
  group_by(RelGroup) %>%
  do(tidy(lm(EAA ~ Relatedness, data=.))) %>%
  filter(term=="Relatedness")

write.xlsx(list("Tertile_Effects"=by_ter), "tables/Tertile_Effects.xlsx")

p <- ggplot(by_ter, aes(RelGroup, estimate)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=estimate-2*std.error, ymax=estimate+2*std.error), width=.1) +
  geom_hline(yintercept=0, lty=2) +
  labs(title="β(Relatedness) by tertile", x="Relatedness tertile", y="β") +
  theme_minimal()

print(p)
ggsave("figures/Fig_Tertiles.png", p, width=6, height=5, dpi=600)
