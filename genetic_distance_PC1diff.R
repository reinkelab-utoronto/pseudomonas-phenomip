big_tree <- read.tree("WI.20231213.hard-filter.isotype.min4.tree")

PC1_df_novc <- non_divergent_PC1_df_avg_zscore_df %>%
  filter(! strain == "VC20019")

my_strains <- unique(PC1_df_novc$strain)

# Find tips that are NOT in your list of strains
tips_to_drop <- setdiff(big_tree$tip.label, my_strains)

# Prune the tree to keep only the desired strains
trimmed_tree <- drop.tip(big_tree, tips_to_drop)

# Compute the pairwise patristic distance matrix
distance_matrix <- cophenetic.phylo(trimmed_tree)

# Reorder the distance matrix by row and column names alphabetically
distance_matrix_sorted <- distance_matrix[order(rownames(distance_matrix)), order(colnames(distance_matrix))]

# View part of the matrix
head(distance_matrix_sorted)

for (i in 2:nrow(conditions)) {
  each_condition <- conditions$types[i]
  condition_data <- PC1_df_novc[PC1_df_novc$condition == each_condition, ]
  head(condition_data)
  
  condition_data <- condition_data %>%
    arrange(strain)
  
  # Get the vector of averages and strain names
  averages <- condition_data$average
  strain_names <- condition_data$strain
  
  # Compute the absolute differences
  diff_matrix <- outer(averages, averages, FUN = function(x, y) abs(x - y))
  
  # Set row and column names
  rownames(diff_matrix) <- strain_names
  colnames(diff_matrix) <- strain_names
  
  # View the matrix
  diff_matrix
  
  identical(rownames(diff_matrix), sort(strain_names))
  identical(colnames(diff_matrix), sort(strain_names))
  
  identical(colnames(diff_matrix), colnames(distance_matrix_sorted))
  identical(rownames(diff_matrix), rownames(distance_matrix_sorted))
  
  
  
  distance_long <- melt(as.matrix(distance_matrix_sorted), varnames = c("Strain1", "Strain2"), value.name = "Distance")
  distance_upper <- distance_long[as.vector(upper.tri(distance_matrix_sorted)), ]
  
  diff_matrix_long <- melt(as.matrix(diff_matrix), varnames = c("Strain1", "Strain2"), value.name = "Zscore_difference")
  diff_upper <- diff_matrix_long[as.vector(upper.tri(diff_matrix)), ]
  
  
  nrow(distance_upper)
  nrow(diff_upper)
  
  if (!identical(distance_upper$Strain1, diff_upper$Strain1) || 
      !identical(distance_upper$Strain2, diff_upper$Strain2)) {
    print("Error!")
  } else {
    print("Both matrices same strain names yay!")
  }
  
  
  head(distance_upper)
  head(diff_upper)
  
  combined_df <- cbind(distance_upper, Zscore_difference = diff_upper$Zscore_difference)
  
  head(combined_df)
  
  
  r2 <- summary(lm(Zscore_difference ~ Distance, data = combined_df))$r.squared
  r2_label <- paste0("R² = ", round(r2, 3))
  
  # Step 1: Label the rows based on strains
  combined_df$highlight <- case_when(
    combined_df$Strain1 == "CB4856" | combined_df$Strain2 == "CB4856" ~ "CB4856",
    combined_df$Strain1 == "QX1212" | combined_df$Strain2 == "QX1212" ~ "QX1212",
    TRUE ~ "Other"
  )
  
  # Step 2: Plot with color mapping
  distance_diffp <- ggplot(combined_df, aes(x = Distance, y = Zscore_difference)) +
    geom_point(aes(color = highlight), alpha = 0.7, size = 2) +
    scale_color_manual(
      values = c("CB4856" = "forestgreen", "QX1212" = "orange", "Other" = "grey60")
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      title = "Correlation between Genetic Distance and PC1 Difference",
      x = "Patristic Distance",
      y = "PC1 score Difference",
      subtitle = each_condition,
      color = "Highlighted Strains"
    ) +
    annotate(
      "text",
      x = max(combined_df$Distance),
      y = max(combined_df$Zscore_difference),
      label = r2_label,
      hjust = 1.1,
      vjust = 1.2,
      size = 5,
      color = "blue"
    ) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 15),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.position = "bottom"
    )
  
  ggsave(paste0("pheno6_", each_condition, "_PC1_distance_diff_correlation.pdf"), plot =distance_diffp, units = 'mm', width = 297, height = 210)
  
  binned_df <- combined_df %>%
    mutate(distance_bin = cut(Distance,
                              breaks = seq(0, max(Distance) + 0.01, by = 0.05),
                              include.lowest = TRUE)) %>%
    group_by(distance_bin) %>%
    summarise(avg_zscore = median(Zscore_difference, na.rm = TRUE),
              count = n())
  
  distance_diff_hisp <- ggplot(binned_df, aes(x = distance_bin, y = avg_zscore)) +
    geom_col(fill = "steelblue") +
    labs(x = "Patristic Distance (binned)",
         y = "Median Z-score Difference",
         title = "Median Z-score Difference by Genetic Distance Bin",
         subtitle = each_condition) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("pheno6_", each_condition, "pc1_diff_median_histogram.pdf"), plot =distance_diff_hisp, units = 'mm', width = 297, height = 210)
  
}


