library(ape)
library(ggtree)
big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")

#total_correlation_df_regression_no_control <- read.csv("bac_total_correlation_df_regression_no_control.csv")
my_strains <- unique(total_correlation_df_regression_no_control$strain)

total_correlation_df_regression_no_control$condition <- factor(
  total_correlation_df_regression_no_control$condition,
  levels = c("MO", "PA", "PV", "M395", "M114", "M117")
)
# Find tips that are NOT in your list of strains
tips_to_drop <- setdiff(big_tree$tip.label, my_strains)

# Prune the tree to keep only the desired strains
trimmed_tree <- drop.tip(big_tree, tips_to_drop)

# Save the trimmed tree to a new file
#write.tree(trimmed_tree, file = "trimmed_tree.newick")

# Optional: Visualize the trimmed tree
tree <- ggtree(trimmed_tree, branch.length = "none") 

tree2 <- ggtree(trimmed_tree, branch.length = "none") +
  geom_tiplab(size = 2)
ggsave("regression_trimmed_tree_no_controls.pdf", tree2, units = 'mm', width = 297, height = 210)

# Row order = tree tips from top to bottom, read from the pruned tree
strain_order <- rev(dplyr::pull(dplyr::arrange(dplyr::filter(ggtree(trimmed_tree)$data, isTip), y), label))

total_correlation_df_regression_no_control$strain <- factor(
  total_correlation_df_regression_no_control$strain, 
  levels = rev(strain_order)  # Reverse order to match top-down view
)

total_correlation_df_regression_no_control <- total_correlation_df_regression_no_control %>%
  mutate(pool = factor(pool, levels = c("pool1", "pool2", "pool3")))

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "brown")
y_axis_colors <- pool_colors[total_correlation_df_regression_no_control$pool]

heatmap_plot <- ggplot(total_correlation_df_regression_no_control, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    limits = c(-10, 10),  
    oob = squish,  
    name = "Average modified Z-score"
  ) +
  labs(
    x = "Condition", 
    y = "", 
    title = "Heatmap of average z-scores", 
    subtitle = "Calculated using regression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size =15), 
    axis.title.x = element_text(size = 18), 
    legend.position = "top",
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 18),  
    plot.title = element_text(size = 22, face = "bold"),  
    plot.subtitle = element_text(size = 18)    
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1)) + 
  theme(
    axis.text.y = element_text(color = y_axis_colors[match(levels(total_correlation_df_regression_no_control$strain), total_correlation_df_regression_no_control$strain)])
  )

tree <- ggtree(trimmed_tree, branch.length = "none") 
heatmap_plot <- heatmap_plot + theme(plot.margin = margin(0, 0, 0, -22))


final_plot_regression <- plot_grid(tree, heatmap_plot, align = "h", rel_widths = c(0.5, 2))



ggsave(paste0("phenomip6_regression_bac_zscore_heatmap_dendro.pdf"), 
       plot =final_plot_regression, units = 'mm', width = 300, height = 350)



heatmap_plot_no_strainnames <- ggplot(
  filter(total_correlation_df_regression_no_control, condition != "MO"), 
  aes(x = condition, y = strain, fill = average)
) +
  geom_tile() +  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    limits = c(-10, 10),  
    oob = squish,  
    name = "Average modified Z-score"
  ) +
  labs(
    x = "Condition", 
    y = "", 
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
    axis.title.x = element_text(size = 18),
    legend.position = "top",
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 18),  
    plot.title = element_text(size = 22, face = "bold"),  
    plot.subtitle = element_text(size = 18)  
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1))

tree <- ggtree(trimmed_tree, branch.length = "none") 
heatmap_plot_nostrainnames <- heatmap_plot_no_strainnames + theme(plot.margin = margin(0, 0, 0, -22))


final_plot_regression2 <- plot_grid(tree, heatmap_plot_nostrainnames, align = "h", rel_widths = c(0.5, 2))



ggsave(paste0("phenomip6_regression_bac_zscore_heatmap_dendro_nonames.pdf"), 
       plot =final_plot_regression2, units = 'mm', width = 300, height = 350)




