controls <- c("JU1400", "CB4856", "VC20019")

auc_avg_zscore_df <- total_auc_df %>%
  group_by(pool, strain, condition) %>%
  summarise(avg_z_score = mean(z_score, na.rm = TRUE), .groups = "drop") %>%
  filter(! strain %in% controls)

library(tibble)
wide_matrix <- auc_avg_zscore_df %>%
  select(strain, condition, avg_z_score) %>%
  pivot_wider(names_from = condition, values_from = avg_z_score) %>%
  column_to_rownames("strain")

wide_matrix <- as.matrix(wide_matrix)


tidy_heatmap_df <- as.data.frame(wide_matrix) %>%
  rownames_to_column("strain") %>%
  pivot_longer(-strain, names_to = "condition", values_to = "average")


# Compute distance matrix and hierarchical clustering
distance_matrix <- dist(wide_matrix)  # Compute distance between strains
hc <- hclust(distance_matrix, method = "complete")  # Use same method as pheatmap()

# Convert to dendrogram object
dend <- as.dendrogram(hc)

# Extract strain order from the dendrogram
strain_order <- order.dendrogram(dend)  
#strain_names <- rownames(wide_matrix)[strain_order]
tree <- as.phylo(hc)  # Convert hclust to phylo for ggtree
# Heatmap row order = dendrogram tips, bottom to top as ggtree draws them
strain_names <- dplyr::pull(dplyr::arrange(dplyr::filter(ggtree(tree)$data, isTip), y), label)

tree_plot <- ggtree(tree, branch.length = "none") +
  geom_tiplab(size = 2)+
  theme(plot.margin = margin(0, 0, 0, -22))
ggsave("auc_hc_tree_plot.pdf", tree_plot, width = 40, height = 21, 
       units = "cm", device='pdf')


# Ensure the strain column is a factor with levels ordered according to the dendrogram
tidy_heatmap_df <- as.data.frame(wide_matrix) %>%
  rownames_to_column("strain") %>%
  pivot_longer(-strain, names_to = "condition", values_to = "average") %>%
  mutate(strain = factor(strain, levels = strain_names))  # Ensure strain is ordered correctly

# Hierarchical clustering for conditions (columns)
col_distance_matrix <- dist(t(wide_matrix))  # Transpose to cluster columns
col_hc <- hclust(col_distance_matrix, method = "complete")
col_dend <- as.dendrogram(col_hc)
col_tree <- as.phylo(col_hc)

# Generate condition dendrogram plot
col_tree_plot2 <- ggtree(col_tree, branch.length = "none") +
  geom_tiplab(hjust = 1) +
  geom_nodelab(aes(label = node), hjust = -0.2, vjust = -0.5) +  # Add node labels
  scale_y_reverse() +  # Flip the tree vertically
  theme(plot.margin = margin(10, 0, 0, 0))+
  coord_flip()

ggsave(paste0("auc_hclust_condition_dendrogram_pre-rotate.pdf"), 
       plot =col_tree_plot2, units = 'mm', width = 300, height = 350)

col_tree_plot <- ggtree(col_tree, branch.length = "none") +
  scale_y_reverse() +  # Flip the tree vertically
  theme(plot.margin = margin(10, 0, 0, 0))+
  coord_flip()

col_tree_plot3 <- col_tree_plot +
  geom_tiplab()

rotated_tree <- col_tree_plot %>%
  ggtree::rotate(node = 8) %>%
  ggtree::rotate(node = 9)

rotated_tree3 <- col_tree_plot3 %>%
  ggtree::rotate(node = 8) %>%
  ggtree::rotate(node = 9)

ggsave(paste0("auc_hclust_condition_dendrogram_after-rotate.pdf"), 
       plot =rotated_tree3, units = 'mm', width = 300, height = 350)



# Ensure the condition column is a factor with the correct levels
tidy_heatmap_df$condition <- factor(
  tidy_heatmap_df$condition,
  levels = c("MO", "PA", "PV", "M395", "M114", "M117")
)

all_strains_pool <- c()
for (i in 1:nrow(tidy_heatmap_df)) {
  strain <- tidy_heatmap_df$strain[i]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
tidy_heatmap_df$pool <- all_strains_pool


tidy_heatmap_df <- tidy_heatmap_df %>%
  mutate(pool = factor(pool, levels = c("pool1", "pool2", "pool3")))

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "brown")
y_axis_colors <- pool_colors[tidy_heatmap_df$pool]

# Create heatmap
heatmap_plot <- ggplot(tidy_heatmap_df, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-10, 10),
    oob = scales::squish,
    name = "Average modified Z-score"
  ) +
  labs(
    x = "",
    y = "",
    title = "Heatmap of Average Z-scores",
    subtitle = "Using Hierarchical Clustering of Z-scores calculated using AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1)) + 
  theme(
    axis.text.y = element_text(color = y_axis_colors[match(levels(tidy_heatmap_df$strain), tidy_heatmap_df$strain)])
  )

# Print the plot
print(heatmap_plot)



heatmap_plot <- heatmap_plot + theme(plot.margin = margin(0, 0, 0, -18))

tree_plot2 <- ggtree(tree, branch.length = "none") +
  theme(plot.margin = margin(0, 0, 0, 10))

top_row <- plot_grid(tree_plot2, heatmap_plot, align = "h", rel_widths = c(0.5, 2))

rotated_tree <- rotated_tree +
  theme(plot.margin = margin(0, 0, 0, 0))

rotated_tree2 <- rotated_tree +
  theme(plot.margin = margin(0, 0, 0, 250))


final_plot_auc <- plot_grid(top_row, rotated_tree2, ncol = 1, rel_heights = c(2, 0.1))

#final_plot_auc <- plot_grid(top_row, rotated_tree, ncol = 1, rel_heights = c(2, 0.5), align = "v")



#final_plot_auc <- plot_grid(tree_plot2, heatmap_plot, align = "h", rel_widths = c(0.5, 2))

#tree_plot2 + heatmap_plot + rotated_tree


ggsave(paste0("phenomip6_auc_bac_zscore_heatmap_dendro_hclust.pdf"), 
       plot =final_plot_auc, units = 'mm', width = 300, height = 350)

# Ensure column order matches the desired levels
condition_order <- c("MO", "M114", "M117", "M395", "PA", "PV")
wide_matrix <- wide_matrix[, condition_order]  # Reorder columns

# Define the color breaks to cap at -10 and 10
breaks <- seq(-10, 10, length.out = 51)  # Creates 50 color bins from -10 to 10

# Generate pheatmap with ordered conditions
auc_pheatmap <- pheatmap(
  wide_matrix,
  cluster_rows = TRUE,    
  cluster_cols = FALSE,   
  scale = "none",        
  color = colorRampPalette(c("blue", "white", "red"))(50), 
  breaks = breaks,  
  main = "Hierarchical Clustering of Strains by Condition"
)

ggsave(paste0("phenomip6_auc_pheatmap.pdf"), 
       plot =auc_pheatmap, units = 'mm', width = 300, height = 350)
ggdendrogram(hc, rotate = TRUE)



