


library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv('build/results.csv')
df <- df %>% mutate(fitness=factor(fitness))

print(head(df))
ggplot(df, aes(model, omega, color=fitness)) + 
  theme_bw() +
  geom_point(position = position_jitter(w = 0.25), size=3) +
  scale_y_log10() +
  scale_x_discrete("model", limits=c("noselection","empiricalvalues","purifyingchem","homoresidue"),
                   labels = c("noselection" = "Neutral",
                                        "empiricalvalues" = "EmpiricalValue",
                                       "purifyingchem" = "ChemicalAffinity","homoresidue" = "Empirical\nHomoresidue"))
sessionInfo()

df <- df %>% filter(model=='homoresidue')
df <- df %>% filter(generation==20000)
ggplot(df, aes(model, omega, color=fitness)) + 
  theme_bw() +
  geom_point(position = position_jitter(w = 0.25), size=3) +
  scale_color_discrete(breaks=c("0.01", "0.0105"),
                      labels=c("neutral", "purifying"))+
  scale_x_discrete("model", 
                   labels = c("noselection" = "Neutral",
                              "empiricalvalues" = "EmpiricalValue",
                              "purifyingchem" = "ChemicalAffinity","homoresidue" = "Empirical\nHomoresidue"))
