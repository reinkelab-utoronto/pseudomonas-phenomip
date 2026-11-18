library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(scales)
library(dendextend)

# ==============================================================================
# 1. SHARED SETTINGS
# ==============================================================================
# (Assuming total_auc_df, total_regression_df, and final_freezing_plan are loaded in your environment)

# Defined first: used in the data prep below
controls <- c("JU1400", "CB4856", "VC20019")

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange",
                 "control" = "pink", "did_not_freeze" = "grey50")

# Explicit condition order from the left, enforced across both panels
condition_levels <- c("PA", "PV", "M395", "M114", "M117")

# ==============================================================================
# 2. LOAD & PREP DATA
# ==============================================================================
auc_avg_zscore_df <- total_auc_df %>%
  group_by(pool, strain, condition) %>%
  summarise(avg_z_score = mean(z_score, na.rm = TRUE), .groups = "drop") %>%
  filter(!strain %in% controls, condition != "MO")

regression_avg_zscore_df <- total_regression_df %>%
  group_by(pool, strain, condition) %>%
  summarise(avg_z_score = mean(z_score, na.rm = TRUE), .groups = "drop") %>%
  filter(!strain %in% controls, condition != "MO")

# Dendrogram-derived strain vectors (heatmap row order)
# Heatmap row orders are read from the dendrograms below (section 5)


# ==============================================================================
# 3. CLUSTERING HELPER (shared by the figure and the tree export)
# ==============================================================================
cluster_trees <- function(avg_zscore_df) {
  wide_matrix <- avg_zscore_df %>%
    dplyr::select(strain, condition, avg_z_score) %>%
    pivot_wider(names_from = condition, values_from = avg_z_score) %>%
    column_to_rownames("strain") %>%
    as.matrix()
  wide_matrix <- wide_matrix[complete.cases(wide_matrix), ]
  
  # Strains (y-axis): Euclidean distance, complete linkage
  hc <- hclust(dist(wide_matrix), method = "complete")
  
  # Bacteria/conditions (x-axis), rotated to condition_levels as drawn
  col_hc   <- hclust(dist(t(wide_matrix)), method = "complete")
  col_dend <- dendextend::rotate(as.dendrogram(col_hc), order = condition_levels)
  
  list(wide_matrix    = wide_matrix,
       strain_tree    = as.phylo(hc),
       condition_tree = as.phylo(as.hclust(col_dend)))
}

# Tip order as ggtree draws it (bottom to top), to compare with the heatmap rows
tip_order <- function(phy) {
  ggtree(phy)$data %>% filter(isTip) %>% arrange(y) %>% pull(label)
}

# ==============================================================================
# 4. HELPER: Build One Tree + Heatmap Panel
# ==============================================================================
build_panel <- function(avg_zscore_df, strain_names, show_legend = TRUE) {
  
  # Clustering (same function used for the tree export)
  cl          <- cluster_trees(avg_zscore_df)
  wide_matrix <- cl$wide_matrix
  tree        <- cl$strain_tree
  col_tree    <- cl$condition_tree
  
  # Re-shape for ggplot
  tidy_heatmap_df <- as.data.frame(wide_matrix) %>%
    rownames_to_column("strain") %>%
    pivot_longer(-strain, names_to = "condition", values_to = "average") %>%
    mutate(
      strain    = factor(strain, levels = strain_names),
      condition = factor(condition, levels = condition_levels)
    )
  
  # Assign pool per strain from final_freezing_plan
  all_strains_pool <- vapply(as.character(tidy_heatmap_df$strain), function(s) {
    if (s %in% final_freezing_plan$Strain) {
      as.character(final_freezing_plan$pool[final_freezing_plan$Strain == s][1])
    } else {
      "did_not_freeze"
    }
  }, character(1))
  
  tidy_heatmap_df$pool <- factor(all_strains_pool,
                                 levels = c("pool1", "pool2", "pool3", "did_not_freeze"))
  
  # Align y-axis label colors with factor levels
  y_axis_colors <- pool_colors[
    as.character(tidy_heatmap_df$pool)[match(levels(tidy_heatmap_df$strain), tidy_heatmap_df$strain)]
  ]
  
  # ==============================
  # PLOT ADJUSTMENTS HERE
  # ==============================
  # Row Tree Plot (Left) - expand=c(0,0) and extreme negative right margin
  tree_plot <- ggtree(tree, branch.length = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.5)) +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(0, -80, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank()
    )
  
  # Column Tree Plot (Bottom, flipped & reversed to face the heatmap)
  col_tree_plot <- ggtree(col_tree, branch.length = "none") +
    coord_flip() +
    scale_y_reverse(expand = c(0, 0.5)) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Heatmap Plot (Right) - Removed invisible y-axis ticks and lengths, bumped negative left margin
  heatmap_plot <- ggplot(tidy_heatmap_df, aes(x = condition, y = strain, fill = average)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-10, 10),
      oob = squish,
      name = "Average Z-score"
    ) +
    scale_x_discrete(expand = c(0, 0.5)) +
    scale_y_discrete(guide = guide_axis(n.dodge = 1), expand = c(0, 0.5)) +
    labs(x = "", y = "") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 4, color = y_axis_colors),
      axis.ticks.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt"),
      axis.title  = element_text(size = 8),
      legend.position = "top",
      legend.text  = element_text(size = 6),
      legend.title = element_text(size = 8),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(10, "mm"),
      plot.title    = element_text(size = 8, face = "bold"),
      plot.margin = margin(0, 0, 0, -30),
      plot.background = element_blank()
    )
  
  # Patchwork logic to combine the 3 structural elements
  if (show_legend) {
    design <- "
      #D
      AC
      #B
    "
    panel <- tree_plot + col_tree_plot + heatmap_plot + guide_area() +
      plot_layout(design = design, widths = c(0.5, 2), heights = c(0.2, 2, 0.1), guides = "collect") &
      theme(plot.margin = margin(0, 0, 0, 0))
  } else {
    design <- "
      AC
      #B
    "
    heatmap_plot <- heatmap_plot + theme(legend.position = "none")
    panel <- tree_plot + col_tree_plot + heatmap_plot +
      plot_layout(design = design, widths = c(0.5, 2), heights = c(2, 0.1)) &
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  wrap_elements(panel)
}

# ==============================================================================
# 5. TREES + ALIGNMENT CHECK
# ==============================================================================
trees_A <- cluster_trees(auc_avg_zscore_df)          # Panel A (AUC)
trees_B <- cluster_trees(regression_avg_zscore_df)   # Panel B (regression)

# Heatmap rows follow the dendrogram tips, so rows and tree always line up
strain_names_auc <- tip_order(trees_A$strain_tree)
strain_names_regression <- tip_order(trees_B$strain_tree)

# ==============================================================================
# 6. ASSEMBLE PLOTS
# ==============================================================================
# Panel A - AUC (Legend generated at the top)
panel_A <- build_panel(auc_avg_zscore_df, strain_names_auc, show_legend = TRUE) &
  theme(plot.margin = margin(t = 2, r = 0, b = -6, l = 0))

# Panel B - Regression (No duplicate legend)
panel_B <- build_panel(regression_avg_zscore_df, strain_names_regression, show_legend = FALSE) &
  theme(plot.margin = margin(t = -6, r = 0, b = 0, l = 0))

# Combine Panel A on top of Panel B
supp_fig8 <- panel_A / panel_B +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold"),
      plot.tag.position = "topleft"
    )
  ) &
  theme(
    plot.margin = margin(0, 2, 0, 2)
  )

supp_fig8

# ==============================================================================
# 7. EXPORT FIGURE
# ==============================================================================
ggsave("supp_fig8.pdf",
       plot = supp_fig8,
       width = 180,
       height = 220,
       units = "mm",
       device = cairo_pdf)

# ==============================================================================
# 8. TREES FOR SUPPLEMENTARY DATA (Newick; clustering dendrograms,
#    Euclidean distance, complete linkage)
# ==============================================================================
tree_s8a_strain_df    <- data.frame(newick = write.tree(trees_A$strain_tree,    file = ""))  # Panel A, y-axis
tree_s8a_condition_df <- data.frame(newick = write.tree(trees_A$condition_tree, file = ""))  # Panel A, x-axis
tree_s8b_strain_df    <- data.frame(newick = write.tree(trees_B$strain_tree,    file = ""))  # Panel B, y-axis
tree_s8b_condition_df <- data.frame(newick = write.tree(trees_B$condition_tree, file = ""))  # Panel B, x-axis