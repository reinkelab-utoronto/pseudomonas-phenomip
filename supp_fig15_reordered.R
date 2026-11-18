library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

library(plyr)
library(data.table)
library(gplots)
library(ggdendro)
library(reshape2)
library(grid)
library(cluster) 
library(tidyr)
library(stringr)
library(pracma)
library(rstatix)
library(writexl)
library(readr)
library(cowplot)
library(pheatmap)
library(corrplot)

final_freezing_plan <- read.csv("final_freezing_plan.csv")
pool_one <- c(1, 4, 9, 11, 16)
pool_two <- c(2,3,6,8, 15)
pool_three <- c(10,12,14)
pool_four <- c(5,7,13)
control <- c("CB4856", "JU1400", "VC20019")

pool <- c()
for (i in 1:nrow(final_freezing_plan)) {
  subpool <- final_freezing_plan[i, "Subpool"]
  
  if (subpool %in% pool_one) {
    pool <- append(pool, "pool1")
  } else if (subpool %in% pool_two) {
    pool <- append(pool, "pool2")
  } else if (subpool %in% pool_three) {
    pool <- append(pool, "pool3")
  } else if (subpool %in% pool_four) {
    pool <- append(pool, "pool4")
  } else if (identical(subpool, "Control")) {
    pool <- append(pool, "control")
  } else {
    pool <- append(pool, "removed")
  }
  
}
final_freezing_plan$pool <- pool

bac <- read.csv("bac_all3data.csv")
mic <- read.csv("mic_all3data.csv")

pathogens <- rbind(bac, mic)

pathogens <- pathogens %>%
  filter(! strain == "VC20019") %>%
  filter(! strain == "CB4856") %>%
  filter(! strain == "JU1400")

AUC <- pathogens %>%
  filter(method == "AUC")
slope <- pathogens %>%
  filter(method == "slope")
PC1 <- pathogens %>%
  filter(method == "PC1")

big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")

# ==============================================================================
# 1. SHARED PREP
# ==============================================================================

strain_order <- c("ECA722", "ECA730", "ECA744", "ECA347", "ECA742", "ECA733", "XZ2019", 
                  "ECA703", "ECA760", "ECA807", "CX11264", "CX11315", "QG556", "QX1793", 
                  "JU2526", "NIC258", "NIC265", "JU782", "NIC262", "NIC251", "WN2001", 
                  "ECA369", "ECA928", "QX1212", "ECA923", "ECA738", "XZ1515", "ECA705", 
                  "ECA743", "ECA592", "ECA594", "ECA930", "ED3005", "MY18", "MY795", 
                  "JU2829", "JU2016", "MY518", "JU2017", "JU2513", "MY2713", "JU1896", 
                  "JU2841", "NIC272", "NIC276", "QG2841", "JU323", "QG2813", "QG2811", 
                  "QG2832", "QG2837", "QG2836", "QG2825", "BRC20263", "JU2811", "QG2846", 
                  "JU3169", "NIC269", "JU2522", "JU2001", "JU2619", "JU2838", "JU3144", 
                  "XZ1735", "QW947", "NIC528", "QG2855", "QG2873", "NIC1119", "ED3077", 
                  "ED3040", "CB4852", "QG2843", "XZ2020", "JU2825", "JU2257", "JU394", 
                  "NIC1049", "JU2575", "WN2066", "JU2853", "JU1249", "JU367", "WN2064", 
                  "JU2907", "JU2572", "JU1792", "MY2147", "NIC511", "JU2586", "NIC515", 
                  "NIC266", "NIC268", "EG4349", "JU1666", "MY772", "JU1530", "NIC166", 
                  "NIC236", "ED3012", "JU847", "JU1652", "NIC195", "JU1212", "JU1246", 
                  "EG4347", "PX179", "JU1200", "CB4854", "CX11292", "NIC231", "LKC34", 
                  "MY1", "JU310", "JU3128", "PB306", "CB4857", "EG4946", "ED3011", 
                  "XZ1756", "ED3046", "JU2464", "QG536", "RC301", "GXW1", "JU3225", 
                  "MY2585", "NIC267", "ED3048", "JU393", "NIC501", "JU2106", "WN2050", 
                  "JU2610", "JU3291", "JU2250", "JU792", "JU1934", "JU3132", "JU1491", 
                  "NIC207", "CB4932", "NIC242", "NIC259", "JU2576", "WN2063", "JU2592", 
                  "JU311", "NIC529", "NIC526", "NIC527")

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "pool4" = "violet")
condition_levels <- c("PA", "PV", "M395", "M114", "M117", "LUA", "LVS", "M1", "M2", "M5", "MA", "PP")

# ==============================================================================
# 2. DATA PREP
# ==============================================================================

combined_heat_df <- bind_rows(
  AUC %>% mutate(Metric = "AUC"),
  slope %>% mutate(Metric = "regression")
) %>%
  filter(!(condition %in% c("MO-mic", "MO-bac"))) %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    strain = factor(strain, levels = rev(strain_order)),
    pool = factor(pool, levels = c("pool1", "pool2", "pool3", "pool4"))
  )

big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")
my_strains <- unique(as.character(combined_heat_df$strain))
tips_to_drop <- setdiff(big_tree$tip.label, my_strains)
trimmed_tree <- drop.tip(big_tree, tips_to_drop)

strain_pool_lookup <- combined_heat_df %>%
  distinct(strain, pool) %>%
  mutate(label = as.character(strain)) %>%
  dplyr::select(label, pool)

# ==============================================================================
# 3. BUILD PLOTS
# ==============================================================================

tree <- ggtree(trimmed_tree, branch.length = "none") %<+% strain_pool_lookup +
  geom_tiplab(aes(color = pool), size = 5 / .pt, align = TRUE, linesize = 0.2) +
  scale_color_manual(
    values = pool_colors, name = "Pool",
    guide = guide_legend(override.aes = list(label = "\u25CF", size = 8 / .pt))
  ) +
  # REDUCED upper expansion to tightly constrain the bounding box
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + 
  coord_cartesian(clip = "off") + 
  theme(
    # INCREASED negative right margin to -25 to overlap with heatmap panel
    plot.margin = margin(0, -25, 0, 0),
    legend.position = "top",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm")
  )

faceted_heatmaps <- ggplot(combined_heat_df, aes(x = condition, y = strain, fill = value)) +
  geom_tile() +  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    limits = c(-10, 10), oob = squish, name = "Average modified Z-score"
  ) +
  facet_wrap(~ Metric, ncol = 2) +
  labs(x = "Condition", y = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length.y = unit(0, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
    axis.title.x = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top", 
    legend.text = element_text(size = 6), 
    legend.title = element_text(size = 8),
    legend.key.height = unit(3, "mm"), 
    legend.key.width = unit(10, "mm"),
    # INCREASED negative left margin to -25 to meet the tree
    plot.margin = margin(0, 0, -15, -25)
  )

supp_fig15 <- tree + faceted_heatmaps +
  plot_layout(widths = c(1, 2), guides = "collect") & 
  theme(legend.position = "top")

# ==============================================================================
# 4. EXPORT
# ==============================================================================

ggsave("supp_fig15.pdf", 
       plot = supp_fig15, 
       width = 180,   
       height = 220,  
       units = "mm", 
       device = cairo_pdf)


tree_s15_df <- data.frame(newick = write.tree(trimmed_tree, file = ""))


tip_order <- function(phy) ggtree(phy)$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

identical(tip_order(trimmed_tree), levels(droplevels(combined_heat_df$strain)))

bac_mic_heatmap <- combined_heat_df
