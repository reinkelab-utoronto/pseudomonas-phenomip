library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "control" = "pink")

# ==============================================================================
# 1. SHARED PREP
# ==============================================================================
big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")
my_strains <- unique(total_correlation_df_auc_no_control$strain)

tips_to_drop <- setdiff(big_tree$tip.label, my_strains)
trimmed_tree <- drop.tip(big_tree, tips_to_drop)

strain_order <- c("ECA722", "ECA744", "ECA347", "ECA742", "ECA733", "XZ2019", "ECA703", 
                  "ECA807", "NIC251", "JU782", "NIC265", "WN2001", "ECA928", "QX1212", 
                  "ECA923", "ECA738", "XZ1515", "ECA705", "ECA743", "ECA592", "ECA594", 
                  "ECA930", "ED3005", "MY18", "MY795", "JU2829", "JU2016", "MY518", 
                  "JU2017", "JU2513", "MY2713", "JU2841", "NIC272", "NIC276", "JU323", 
                  "JU3169", "NIC269", "JU2522", "JU2001", "JU2619", "JU2838", "BRC20263", 
                  "JU2811", "QG2846", "QG2825", "QG2836", "JU3144", "XZ1735", "QW947", 
                  "NIC528", "ED3077", "NIC1119", "ED3040", "CB4852", "QG2843", "JU2825", 
                  "JU2257", "JU394", "NIC1049", "JU2575", "WN2066", "JU2853", "JU1249", 
                  "JU367", "WN2064", "JU2907", "JU2572", "JU1792", "MY2147", "NIC511", 
                  "JU2586", "NIC515", "NIC266", "EG4349", "JU1666", "MY772", "JU1530", 
                  "NIC166", "NIC236", "JU847", "NIC195", "JU1212", "JU1246", "JU1200", 
                  "PX179", "CB4854", "NIC231", "MY1", "JU310", "JU3128", "ED3046", 
                  "JU2464", "ED3011", "PB306", "RC301", "GXW1", "JU3225", "MY2585", 
                  "NIC267", "ED3048", "JU393", "NIC501", "JU2106", "WN2050", "JU2610", 
                  "JU2250", "JU792", "JU1934", "JU3132", "CB4932", "NIC242", "NIC207", 
                  "JU2576", "WN2063", "JU2592", "JU311", "NIC526", "NIC527")

# ==============================================================================
# 2. DATA PREP
# ==============================================================================
# Filter: keep strains where BOTH AUC and regression have at least 2/5 conditions
# with |average| > 5
combined_heat_df <- bind_rows(
  total_correlation_df_auc_no_control %>% filter(condition != "MO") %>% mutate(Metric = "AUC"),
  total_correlation_df_regression_no_control %>% filter(condition != "MO") %>% mutate(Metric = "regression")
) %>%
  mutate(
    condition = factor(condition, levels = c("PA", "PV", "M395", "M114", "M117")),
    strain = factor(strain, levels = rev(strain_order))
  )

# Strains remaining after filtering
filtered_strains <- unique(as.character(combined_heat_df$strain))

# Drop tips not in the filtered strain set
tips_to_drop2 <- setdiff(trimmed_tree$tip.label, filtered_strains)
trimmed_tree2 <- drop.tip(trimmed_tree, tips_to_drop2)

# --- Strain -> pool lookup for coloring tree tips ---
# NOTE: combined_heat_df must contain a `pool` column carried over from
# total_correlation_df_auc_no_control / total_correlation_df_regression_no_control.
# If it doesn't, build this lookup from strains_analysed.csv instead:
#   strain_pool_lookup <- read.csv("strains_analysed.csv") %>%
#     distinct(strain, pool) %>% rename(label = strain)
strain_pool_lookup <- combined_heat_df %>%
  distinct(strain, pool) %>%
  mutate(label = as.character(strain)) %>%
  dplyr::select(label, pool)

# Tree (with tip labels, stripped of right padding), tips colored by pool
# Tip label size matched to the 6pt axis-text size used elsewhere (Fig 1 convention:
# size given in points via /.pt, since geom_tiplab's `size` arg is in mm like other
# ggplot2 text geoms)
tree <- ggtree(trimmed_tree2, branch.length = "none") %<+% strain_pool_lookup +
  geom_tiplab(aes(color = pool), size = 6 / .pt, align = TRUE, linesize = 0.2) +
  scale_color_manual(
    values = pool_colors, name = "Pool",
    # geom_tiplab is a text geom, so its default legend key shows the
    # letter glyph + dashed alignment line; override it to draw a plain
    # solid dot instead, matching the look of a normal color legend
    guide = guide_legend(override.aes = list(label = "\u25CF", size = 8 / .pt))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(
    plot.margin = margin(0, -5, 0, 0),
    legend.position = "top",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm")
  )

# ==============================================================================
# 3. BUILD PLOTS
# ==============================================================================

# Panel: Tree + Heatmaps
# Formatting matched to main Figure 1 (theme_bw, base_size = 8, 8pt axis titles,
# 6pt axis text, 8pt bold strip text, 8pt legend title, 6pt legend text)
faceted_heatmaps <- ggplot(combined_heat_df, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-10, 10), oob = squish, name = "Z-score") +
  facet_wrap(~ Metric, ncol = 2) +
  labs(x = "Condition", y = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length.y = unit(0, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.title.x = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top", legend.text = element_text(size = 6), legend.title = element_text(size = 8),
    legend.key.height = unit(3, "mm"), legend.key.width = unit(10, "mm"),
    plot.margin = margin(0, 0, -15, -10)
  )

# Collect legends (Pool + Z-score) into one row instead of two separate ones
supp_6 <- tree + faceted_heatmaps +
  plot_layout(widths = c(1.5, 4), guides = "collect") &
  theme(legend.position = "top")
supp_6

# ==============================================================================
# 4. EXPORT
# ==============================================================================
# Full-page figure (long strain list needs the extra vertical room to keep
# tip labels legible at the matched 6pt size)
ggsave("supp_fig6.pdf",
       plot = supp_6,
       width = 180,
       height = 220,
       units = "mm",
       device = cairo_pdf)


worm_species_newick <- write.tree(trimmed_tree2, file = "")
worm_species_tree_df <- data.frame(newick = worm_species_newick)


combined_heat_df_allstrains <- combined_heat_df
normality_S6 <- combined_heat_df_allstrains %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarize(
    n         = sum(!is.na(average)),
    shapiro_W = unname(shapiro.test(average)$statistic),
    shapiro_p = shapiro.test(average)$p.value,
    .groups = "drop"
  )

normality_S6