library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(ggdendro)

# ==============================================================================
# 0. CONDITION LABELS (full names; species in italics, strain names upright)
# ==============================================================================
fig2_condition_labels <- c(
  MO   = "'Uninfected control'",
  PA   = "italic('P. aeruginosa')~'PA14'",
  PV   = "italic('P. vranovensis')~'BIGb446'",
  M395 = "italic('P. canadensis')~'MYb395'",
  M114 = "italic('P. poae')~'MYb114'",
  M117 = "italic('P. poae')~'MYb117'"
)
# Turns condition codes into plotmath labels, in whatever order ggplot asks for
label_conditions <- function(x) parse(text = fig2_condition_labels[as.character(x)])

# ==============================================================================
# 1. SHARED PREP
# ==============================================================================
# NOTE: big_tree / trimmed_tree / strain_order below are no longer used to
# order Panel A (that's now driven by hierarchical clustering in Section 2).
# Left here commented out in case you still need the phylogeny or this vector
# elsewhere. Delete if not needed.

# big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")
# my_strains <- unique(total_correlation_df_auc_no_control$strain)
#
# tips_to_drop <- setdiff(big_tree$tip.label, my_strains)
# trimmed_tree <- drop.tip(big_tree, tips_to_drop)
#
# strain_order <- c("ECA722", "ECA744", "ECA347", "ECA742", "ECA733", "XZ2019", "ECA703",
#                   "ECA807", "NIC251", "JU782", "NIC265", "WN2001", "ECA928", "QX1212",
#                   "ECA923", "ECA738", "XZ1515", "ECA705", "ECA743", "ECA592", "ECA594",
#                   "ECA930", "ED3005", "MY18", "MY795", "JU2829", "JU2016", "MY518",
#                   "JU2017", "JU2513", "MY2713", "JU2841", "NIC272", "NIC276", "JU323",
#                   "JU3169", "NIC269", "JU2522", "JU2001", "JU2619", "JU2838", "BRC20263",
#                   "JU2811", "QG2846", "QG2825", "QG2836", "JU3144", "XZ1735", "QW947",
#                   "NIC528", "ED3077", "NIC1119", "ED3040", "CB4852", "QG2843", "JU2825",
#                   "JU2257", "JU394", "NIC1049", "JU2575", "WN2066", "JU2853", "JU1249",
#                   "JU367", "WN2064", "JU2907", "JU2572", "JU1792", "MY2147", "NIC511",
#                   "JU2586", "NIC515", "NIC266", "EG4349", "JU1666", "MY772", "JU1530",
#                   "NIC166", "NIC236", "JU847", "NIC195", "JU1212", "JU1246", "JU1200",
#                   "PX179", "CB4854", "NIC231", "MY1", "JU310", "JU3128", "ED3046",
#                   "JU2464", "ED3011", "PB306", "RC301", "GXW1", "JU3225", "MY2585",
#                   "NIC267", "ED3048", "JU393", "NIC501", "JU2106", "WN2050", "JU2610",
#                   "JU2250", "JU792", "JU1934", "JU3132", "CB4932", "NIC242", "NIC207",
#                   "JU2576", "WN2063", "JU2592", "JU311", "NIC526", "NIC527")

# ==============================================================================
# 2. DATA PREP
# ==============================================================================
# Filter: keep strains where BOTH AUC and regression have at least 2/5 conditions
# with |average| > 5
combined_heat_df <- bind_rows(
  total_correlation_df_auc_no_control %>% filter(condition != "MO") %>% mutate(Metric = "AUC"),
  total_correlation_df_regression_no_control %>% filter(condition != "MO") %>% mutate(Metric = "regression")
) %>%
  group_by(strain, Metric) %>%
  mutate(passes = sum(average > 5 | average < -5) >= 2) %>%
  group_by(strain) %>%
  filter(all(passes)) %>%
  ungroup() %>%
  dplyr::select(-passes)

heat_2a_df <- combined_heat_df   # unique name for Data S7; other scripts reuse "combined_heat_df"

# ------------------------------------------------------------------------------
# Hierarchical clustering (replaces phylogeny-based ordering)
# ------------------------------------------------------------------------------
# Build a strain x condition matrix, averaged across Metric, purely to derive
# ONE consistent row/column order shared by both the AUC and regression facets.
clust_matrix <- combined_heat_df %>%
  group_by(strain, condition) %>%
  summarise(average = mean(average, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = average) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Cluster strains (rows -> y axis) and conditions (columns -> x axis)
strain_hclust    <- hclust(dist(clust_matrix), method = "average")
condition_hclust <- hclust(dist(t(clust_matrix)), method = "average")

# Convert the strain dendrogram into a ggtree-compatible tree object
strain_dendrogram <- ape::as.phylo(strain_hclust)

# Extract the tip order exactly as ggtree will draw it, so the heatmap's
# y-axis lines up with the dendrogram (don't just trust strain_hclust$order -
# ggtree's plotting order can differ slightly in y-position bookkeeping)
tree_tip_order <- ggtree(strain_dendrogram)$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

condition_order_clustered <- colnames(clust_matrix)[condition_hclust$order]

# Condition dendrogram data, duplicated once per Metric so its facet panel
# widths line up with the AUC / regression facets in faceted_heatmaps
condition_dendro_data <- dendro_data(as.dendrogram(condition_hclust), type = "rectangle")

condition_dendro_segments <- bind_rows(
  condition_dendro_data$segments %>% mutate(Metric = "AUC"),
  condition_dendro_data$segments %>% mutate(Metric = "regression")
) %>%
  mutate(Metric = factor(Metric, levels = c("AUC", "regression")))

condition_dendrogram <- ggplot(condition_dendro_segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) +
  facet_wrap(~ Metric, ncol = 2) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_reverse(expand = expansion(mult = c(0.05, 0))) +
  theme_void(base_size = 8) +
  theme(
    strip.text = element_blank(),
    plot.margin = margin(-5, 0, 0, -10)
  )

# Apply clustered orders as factor levels
combined_heat_df <- combined_heat_df %>%
  mutate(
    strain    = factor(strain, levels = tree_tip_order),
    condition = factor(condition, levels = condition_order_clustered)
  )

# Strains remaining after filtering
filtered_strains <- unique(as.character(combined_heat_df$strain))

# Tree (dendrogram from hierarchical clustering, with tip labels)
tree <- ggtree(strain_dendrogram, branch.length = "none") +
  geom_tiplab(size = 1.5, align = TRUE, linesize = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme(plot.margin = margin(0, -5, 0, 0))

newick_2a <- write.tree(strain_dendrogram, file = "")
tree_2a_df <- data.frame(newick = newick_2a)

# ==============================================================================
# 3. BUILD PLOTS
# ==============================================================================

# Panel A: Tree + Heatmaps
faceted_heatmaps <- ggplot(combined_heat_df, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-10, 10), oob = squish, name = "Z-score") +
  scale_x_discrete(labels = label_conditions) +    # full pathogen names
  facet_wrap(~ Metric, ncol = 2) +
  labs(x = NULL, y = NULL) +
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

# Spacer to sit below the tree, matching the tree's column width
tree_corner_spacer <- plot_spacer()

plot_2a <- (tree + faceted_heatmaps) / (tree_corner_spacer + condition_dendrogram) +
  plot_layout(widths = c(0.9, 4), heights = c(6, 0.7))
plot_2a

normality_2a <- combined_heat_df %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarize(
    n         = sum(!is.na(average)),
    shapiro_W = unname(shapiro.test(average)$statistic),
    shapiro_p = shapiro.test(average)$p.value,
    .groups = "drop"
  )

normality_2a

# Panel B: Averaged Growth Curves
four_strain_averaged <- abundance_df_withg4 %>%
  filter(strain %in% c("JU2610", "ED3005", "JU2522", "WN2050", "JU2250")) %>%
  group_by(strain, condition, generation) %>%
  summarise(avg_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(strain = factor(strain, levels = c("JU2610", "ED3005", "JU2522", "WN2050", "JU2250")))

plot_2b_faceted <- ggplot(four_strain_averaged, aes(x = generation, y = avg_abundance, color = condition, group = condition)) +
  geom_point() + geom_line() +
  scale_colour_discrete(breaks = c("MO", "PA", "PV", "M395", "M114", "M117"),   # legend order
                        labels = label_conditions, name = "Condition") +   # full names; MO = uninfected control
  facet_wrap(~strain, ncol = 3, scales = "free_y") +
  labs(x = "Generation", y = "Average strain abundance (log₂ fold change vs. generation 4)") +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8), axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8), legend.text = element_text(size = 6),
    plot.margin = margin(-15, 0, 0, 0)
  ) +
  guides(colour = guide_legend(ncol = 1))

# ==============================================================================
# 4. FINAL ASSEMBLY (Aligned & Font-Synchronized)
# ==============================================================================

# 1. Wrap the plots to prevent facet-tagging issues
p1 <- wrap_elements(plot_2a)
p2 <- wrap_elements(plot_2b_faceted)

# 2. Assemble with matching font size and plain (non-bold) weight
final_figure_2 <- (p1 / p2) +
  plot_layout(heights = c(2.4, 1.6), guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = "topleft"
  )
# Export (Consider setting height to 220mm to match Figure 1's canvas ratio if they still look off)
ggsave("Figure_2_Final.pdf", 
       plot = final_figure_2, 
       width = 180, 
       height = 220, 
       units = "mm", 
       device = cairo_pdf)