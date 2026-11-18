cluster <- c()
for (i in 1:nrow(merged_pc_G7G9)) {
  each_condition <- merged_pc_G7G9[i, 2] 
  
  if (each_condition == "PA" | each_condition == "PV") {
    cluster <- append(cluster, "non-MYBs")
  } else {
    cluster <- append(cluster, "MYBs")
  }
}
merged_pc_G7G9$cluster <- cluster

merged_pc_G7G9_mean <- merged_pc_G7G9 %>%
  group_by(strain, cluster) %>%
  summarise(
    mean_zscore_G7 = mean(mean_zscore_G7, na.rm = TRUE),
    mean_zscore_G9 = mean(mean_zscore_G9, na.rm = TRUE),
    .groups = "drop"
  )

PC1_bycluster_df <- data.frame()
PC_condition <- c()
G7_loading <- c()
G9_loading <- c()
cluster_name <- c()
variance_byPC1 <- c()
for (i in 1:length(unique(cluster))) {
  
  each_cluster <- unique(cluster)[i]
  
  cluster_data <- merged_pc_G7G9_mean %>%
    filter(cluster == each_cluster)
  
  pca_result <- prcomp(cluster_data[, 3:4], center = TRUE, scale. = TRUE)
  cluster_data$PC1 <- pca_result$x[, 1]
  
  PC1_bycluster_df <- rbind(PC1_bycluster_df, cluster_data)
  # Standard deviations from PCA result
  std_devs <- pca_result$sdev
  
  loadings <- pca_result$rotation[, 1] # Loadings for PC1
  
  # Variance explained by PC1
  var_explained_pc1 <- (std_devs[1]^2) / sum(std_devs^2)
  var_explained_pc1_percent <- var_explained_pc1 * 100 # Print as percentage
  print(var_explained_pc1_percent)
  
  cluster_name <- append(cluster_name, each_cluster)
  G7_loading <- append(G7_loading, loadings[1][[1]])
  G9_loading <- append(G9_loading, loadings[2][[1]])
  variance_byPC1 <- c(variance_byPC1, var_explained_pc1_percent)
  
}
PCA_stats <- data.frame(cluster_name,
                        G7_loading,
                        G9_loading,
                        variance_explained_by_PC1 = variance_byPC1)

varianceby_PC1p <- ggplot(PCA_stats, aes(x = cluster_name, y = variance_explained_by_PC1, fill = cluster_name)) +
  geom_bar(stat = "identity", position = "dodge") +  
  geom_text(aes(label = round(variance_explained_by_PC1, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 6) +  # Increase annotation text size
  labs(
    y = "Variance explained by PC1", 
    title = paste("Variance explained by PC1 by cluster (PhenoMIP 6.0)", pathogen_name)
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
ggsave(paste0("variance explained by PC1 by cluster mean (PhenoMIP 6.0) ", pathogen_name, ".pdf"), plot =varianceby_PC1p, units = 'mm', width = 297, height = 210)


MYBs <- PC1_bycluster_df %>%
  filter(cluster == "MYBs")

MYBs_summary <- MYBs %>%
  select(strain, PC1)
colnames(MYBs_summary) <- c("strain", "MYBs_mean_PC1")


nonMYBs <- PC1_bycluster_df %>%
  filter(cluster == "non-MYBs")
nonMYBs_summary <- nonMYBs %>%
  select(strain, PC1)
colnames(nonMYBs_summary) <- c("strain", "non-MYBs_mean_PC1")


dir.create("non_divergent_nemascan_input/MYBs_mean", recursive = TRUE, showWarnings = FALSE)
dir.create("non_divergent_nemascan_input/non-MYBs_mean", recursive = TRUE, showWarnings = FALSE)
# NemaScan needs CaeNDR strain names (to_caendr_strain_name() is defined in pheno6_bac_v2.R)
MYBs_summary <- to_caendr_strain_name(MYBs_summary)
nonMYBs_summary <- to_caendr_strain_name(nonMYBs_summary)
write_tsv(MYBs_summary, "non_divergent_nemascan_input/MYBs_mean/Pheno6_MYBs_mean_PC1_nemascan_input.tsv")
write_tsv(nonMYBs_summary, "non_divergent_nemascan_input/non-MYBs_mean/Pheno6_non-MYBs_mean_PC1_nemascan_input.tsv")





