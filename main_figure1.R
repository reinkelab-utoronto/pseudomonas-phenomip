library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)   # needed for Panel B (species tree), moved in from supp_fig6.R
library(ape)      # needed for read.tree()

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "control" = "pink")

# Shapes and full names for each pathogen, shared by panels D and E
# (species names in italics; strain names upright)
cond_shapes <- c(PA = 7, PV = 8, M395 = 3, M114 = 17, M117 = 16)
cond_labels <- parse(text = c(
  PA   = "italic('P. aeruginosa')~'PA14'",
  PV   = "italic('P. vranovensis')~'BIGb446'",
  M395 = "italic('P. canadensis')~'MYb395'",
  M114 = "italic('P. poae')~'MYb114'",
  M117 = "italic('P. poae')~'MYb117'"
))
scale_condition_shape <- scale_shape_manual(values = cond_shapes, breaks = names(cond_shapes),
                                            labels = cond_labels, name = "Condition")

# Figure 1a - World Map
map <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "white",
               color = "grey60",
               linewidth = 0.2) +
  geom_point(data = joint_coordinates_divergence,
             aes(x = longitude, y = latitude, color = pool),
             size = 1.5,
             alpha = 0.8) +
  scale_color_manual(values = pool_colors, name = "Pool") + 
  coord_cartesian(ylim = c(-50, 60)) +
  theme_void(base_size = 8) +                 
  theme(
    legend.position = "bottom",               
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

coords <- joint_coordinates_divergence %>%
  dplyr::select(strain, pool, longitude, latitude, origin)

#addWorksheet(supp_data, "Table S2 (Figure 1A)")
#writeData(supp_data, "Table S2 (Figure 1A)", coords)


# ==============================================================================
# Figure 1b - Species Tree (moved in from supp_fig6.R, was Supp Fig 6B)
# ==============================================================================
# Read the tree from your file path
bac_species_tree <- read.tree("SpeciesTree_rooted copy.tree")

# 1. Replace all underscores with spaces
bac_species_tree$tip.label <- gsub("_", " ", bac_species_tree$tip.label)

# 2. Change "Psuedomonas" to "P." to save space
bac_species_tree$tip.label <- gsub("Psuedomonas", "P.", bac_species_tree$tip.label)

# 3. Format the labels into plotmath syntax to italicize Genus/abbr. and species
bac_species_tree$tip.label <- gsub("^(\\S+)\\s+(\\S+)\\s+(.*)$", "italic('\\1 \\2')~'\\3'", bac_species_tree$tip.label)

# Plot the tree using ggtree
panel_1b_species_tree <- ggtree(bac_species_tree) +
  geom_tiplab(size = 6 / .pt, align = FALSE, parse = TRUE) +  
  geom_nodelab(size = 6 / .pt, vjust = -0.5, hjust = -0.2) + 
  theme_tree2(base_size = 8) +                  
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5),
    plot.tag = element_text(size = 10, face = "plain")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.8)))

panel_1b_species_tree_wrapped <- wrap_elements(panel_1b_species_tree) +
  theme(plot.tag = element_text(size = 10, face = "plain"))

#addWorksheet(supp_data, "Table S2 (Figure 1A)")
#writeData(supp_data, "Table S2 (Figure 1A)", coords)

newick <- write.tree(bac_species_tree, file = "")

tree_df <- data.frame(newick = newick)

# Figure 1d - Correlation 
total_correlation_df$method[total_correlation_df$method == "regression z-score"] <- "regression"
total_correlation_df$method[total_correlation_df$method == "AUC z-score"] <- "AUC"

r2_data_overall <- total_correlation_df %>%
  group_by(method) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

total_correlation_p_overall <- ggplot(
  total_correlation_df,
  aes(x = z_R1, y = z_R2, color = pool, shape = condition)
) +
  scale_color_manual(values = pool_colors, name = "Pool") + 
  scale_condition_shape +                          # full pathogen names, shapes shared with panel E
  guides(shape = "none") +                         # shown once, centred below D and E (see assembly)
  geom_point() +
  facet_grid(~method) +
  geom_smooth(
    data = total_correlation_df,
    aes(x = z_R1, y = z_R2),
    method = "lm",
    se = FALSE,
    color = "red",
    fullrange = TRUE,
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)",
    shape = "Condition"
  ) +
  geom_text(
    data = r2_data_overall,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.5,
    size = 8 / .pt,
    color = "black"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",  
    legend.box = "vertical",
    axis.title = element_text(size = 8),      
    axis.text = element_text(size = 6),       
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"), 
    plot.subtitle = element_text(size = 8),   
    legend.title = element_text(size = 8),    
    legend.text = element_text(size = 6)
  )

normality_df_1d <- total_correlation_df %>%
  group_by(method) %>%
  mutate(residual = resid(lm(z_R2 ~ z_R1, na.action = na.exclude))) %>%
  ungroup() %>%
  pivot_longer(c(z_R1, z_R2, residual),
               names_to = "variable", values_to = "value") %>%
  group_by(method, variable) %>%
  summarize(
    n           = sum(!is.na(value)),
    shapiro_W   = unname(shapiro.test(value)$statistic),
    shapiro_p   = shapiro.test(value)$p.value,
    .groups = "drop"
  )

normality_df_1d


# Figure 1e
# 1. Source the scripts into isolated environments
env_auc <- new.env()
env_reg <- new.env()

source("minipheno/correlation_allstatsAUC.R", local = env_auc)
source("minipheno/correlation_allstatsregression.R", local = env_reg)

# 2. Extract Data
df_auc <- env_auc$plot_df %>% mutate(Metric = "AUC")
df_reg <- env_reg$plot_df %>% mutate(Metric = "regression")
combined_plot_df <- bind_rows(df_auc, df_reg)
combined_df <- combined_plot_df
colnames(combined_df) <- c("strain", "condition", "phenoMIP z-score", "mini-phenoMIP z-score", "Metric")

# 3. Extract Labels
label_df <- data.frame(
  Metric = c("AUC", "regression"),
  label = c(env_auc$global_label_avg, env_reg$global_label_avg)
)

# 4. Build Figure 1e
fig_1e <- ggplot(combined_plot_df, aes(x = pheno6_z, y = mini_z)) +
  geom_point(aes(colour = strain, shape = condition), size = 2, alpha = 0.8) +
  scale_condition_shape +                          # same shapes as panel D (legend hidden below)
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", fullrange = TRUE, linewidth = 0.5) +
  geom_text(data = label_df, 
            aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, 
            hjust = 1.1,
            vjust = -0.15,
            size = 8 / .pt) +
  facet_wrap(~ Metric) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(x = "PhenoMIP mean z-score", y = "Mini-phenoMIP mean z-score", colour = "Strain") +  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),          
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  guides(colour = guide_legend(nrow = 4, byrow = TRUE), shape = "none")

get_stats <- function(env, metric) {
  d <- env$plot_df
  data.frame(
    Metric       = metric,
    n            = sum(complete.cases(d$pheno6_z, d$mini_z)),
    pearson_r    = env$r_global_avg,
    pearson_p    = cor.test(d$pheno6_z, d$mini_z, method = "pearson")$p.value,
    r2           = env$r2_global_avg,
    spearman_rho = env$rho_global_avg,
    spearman_p   = cor.test(d$pheno6_z, d$mini_z, method = "spearman", exact = FALSE)$p.value
  )
}

r2_data_1e <- dplyr::bind_rows(
  get_stats(env_auc, "AUC"),
  get_stats(env_reg, "regression")
)

# Per-condition breakdown (label column dropped: its line breaks look messy in Excel)
stats_s5_by_condition <- dplyr::bind_rows(
  env_auc$stats_df_avg %>% dplyr::mutate(Metric = "AUC"),
  env_reg$stats_df_avg %>% dplyr::mutate(Metric = "regression")
) %>%
  dplyr::select(Metric, condition, R, R2, rho)


normality_df_1e <- combined_df %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(residual = resid(lm(`mini-phenoMIP z-score` ~ `phenoMIP z-score`,
                                    na.action = na.exclude))) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(c(`phenoMIP z-score`, `mini-phenoMIP z-score`, residual),
                      names_to = "variable", values_to = "value") %>%
  dplyr::group_by(Metric, variable) %>%
  dplyr::summarize(
    n         = sum(!is.na(value)),
    shapiro_W = unname(shapiro.test(value)$statistic),
    shapiro_p = shapiro.test(value)$p.value,
    .groups = "drop"
  )

normality_df_1e

# --- Shared condition legend for panels D and E ---
# Drawn from a small helper plot and placed centred below both panels
condition_legend_plot <- ggplot(data.frame(condition = names(cond_shapes), x = 1, y = 1),
                                aes(x, y, shape = condition)) +
  geom_point() +
  scale_condition_shape +
  guides(shape = guide_legend(nrow = 1)) +         # one row; use nrow = 2 if it is too wide
  theme_void(base_size = 8) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

legend_grob <- ggplotGrob(condition_legend_plot)
keep <- grepl("guide-box", legend_grob$layout$name) &
  !vapply(legend_grob$grobs, inherits, logical(1), "zeroGrob")
condition_legend <- wrap_elements(full = legend_grob$grobs[[which(keep)[1]]]) +
  theme(plot.margin = margin(t = -8, r = 0, b = 0, l = 0))   # negative top margin = closer to D and E
# --- Patchwork Assembly ---

# 1. Create a true empty plot so patchwork recognizes it and assigns the "C" tag
biorender_placeholder <- ggplot() + theme_void(base_size = 8)

# 2. Layout for 3 rows:
# Row 1: A (Map) and B (Species Tree) split the columns
# Row 2: C (BioRender Spacer) spans both columns
# Row 3: D (Correlation 1) and E (Correlation 2) split the columns
# Row 4: F, the shared condition legend, centred across both columns
layout_design <- "
  AB
  CC
  DE
  FF
"

# 3. Assembly
final_figure <- map + panel_1b_species_tree_wrapped + biorender_placeholder + total_correlation_p_overall + fig_1e +
  condition_legend +
  plot_layout(design = layout_design, heights = c(1, 1, 1, 0.04)) +   # last row: legend only
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", ""))) &   # "" = no letter for the legend
  theme(
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = "topleft"
  )
# View the result
final_figure

# Export with a larger height to accommodate the 3 full rows
ggsave("Figure_1_Final.pdf", 
       plot = final_figure, 
       width = 180, 
       height = 220, 
       units = "mm", 
       device = cairo_pdf)