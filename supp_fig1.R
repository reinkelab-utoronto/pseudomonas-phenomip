library(plyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggdendro)
library(reshape2)
library(grid)
library(cluster) 
library(gplots)
library(tidyr)
library(stringr)
library(pracma)
library(rstatix)
library(readr)

percent_males <- read.csv("pheno8_percent_males.csv")
df <- pivot_longer(percent_males, cols = c("G5", "G7", "G9"), names_to = "Generation")
head(df, n = 20)
df$Condition <- factor(df$Condition, levels = c("MO1", "MO2", "PA", "PV", "M395", "M114", "M117"))

supp_fig1 <- ggplot(df, aes(x = Generation, y = value, group = interaction(Pool, Condition))) +
  geom_point(aes(color = Condition)) +
  geom_line(aes(color = Condition)) +
  facet_grid(. ~ Pool) +
  ylab("% male") +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

# View the result
supp_fig1

# Export — single-panel figure, so height is scaled to ~1/3 of the
# original 220mm used for the 3-row main Figure 1 (width unchanged)
ggsave("supp_fig1.pdf", 
       plot = supp_fig1, 
       width = 180, 
       height = 220/3, 
       units = "mm", 
       device = cairo_pdf)
