PC_condition <- c()
G7_loading <- c()
G9_loading <- c()
total_PC1_df <- c()
#no average controls
#average controls
mean_zscore_df <- merged_pc %>%
  group_by(strain, condition, generation, pool) %>%
  summarise(mean_z_score = mean(z_score_value, na.rm = TRUE), .groups = "drop")

merged_pc_G7G9 <- mean_zscore_df %>%
  pivot_wider(
    names_from = generation,  # Use generation values (7, 9) as column names
    values_from = mean_z_score,  # Fill new columns with mean_z_score values
    names_prefix = "mean_zscore_G"  # Prefix column names with "mean_z_scoreG"
  )

variance_byPC1 <- c()
for (i in 2:length(condition_list)) {
  each_condition <- condition_list[i]
  #PCA_condition <- append(PCA_condition, each_condition)
  print(each_condition)
  
  condition_data <- merged_pc_G7G9 %>%
    filter(condition == each_condition)  
  
  pca_result <- prcomp(condition_data[, 4:5], center = TRUE, scale. = TRUE)
  condition_data$PC1 <- pca_result$x[, 1]
  
  total_PC1_df <- rbind(total_PC1_df, condition_data)
  # Standard deviations from PCA result
  std_devs <- pca_result$sdev
  
  loadings <- pca_result$rotation[, 1] # Loadings for PC1
  
  # Variance explained by PC1
  var_explained_pc1 <- (std_devs[1]^2) / sum(std_devs^2)
  var_explained_pc1_percent <- var_explained_pc1 * 100 # Print as percentage
  print(var_explained_pc1_percent)
  
  
  PC_condition <- append(PC_condition, each_condition)
  G7_loading <- append(G7_loading, loadings[1][[1]])
  G9_loading <- append(G9_loading, loadings[2][[1]])
  variance_byPC1 <- c(variance_byPC1, var_explained_pc1_percent)
}
PCA_stats <- data.frame(condition = PC_condition,
                        G7_loading,
                        G9_loading,
                        variance_explained_by_PC1 = variance_byPC1)


varianceby_PC1p <- ggplot(PCA_stats, aes(x = condition, y = variance_explained_by_PC1, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  geom_text(aes(label = round(variance_explained_by_PC1, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 6) +  # Increase annotation text size
  labs(
    y = "Variance explained by PC1", 
    title = paste("Variance explained by PC1 by condition (PhenoMIP 6.0)", pathogen_name)
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20),  # Base text size
    axis.title = element_text(size = 22),  # Axis labels
    axis.text = element_text(size = 18),   # Tick labels
    legend.text = element_text(size = 18), # Legend labels (if applicable)
    legend.title = element_text(size = 20), # Legend title (if applicable)
    plot.title = element_text(size = 24, face = "bold") # Title
  )
ggsave(paste0("variance explained by PC1 by condition (PhenoMIP 6.0) ", pathogen_name, ".pdf"), plot =varianceby_PC1p, units = 'mm', width = 297, height = 210)

  
