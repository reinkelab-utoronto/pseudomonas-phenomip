library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# NOTE: This script assumes `merged_pc`, `h2_df`, and `pathogen_name` 
# are already loaded in your environment prior to running.

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "control" = "pink")

# ==============================================================================
# PANEL A: overall_PCA
# ==============================================================================
merged_pc_noMO <- merged_pc %>%
  filter(!condition == "MO")

R1_normalised_df <- merged_pc_noMO %>%
  filter(replicate_group == "R1") %>%
  group_by(strain, condition, generation) %>%
  summarise(mean_z_score = mean(z_score_value), .groups = "drop") %>%
  mutate(replicate_group = "R1")

R2_normalised_df <- merged_pc_noMO %>%
  filter(replicate_group == "R2") %>%
  group_by(strain, condition, generation) %>%
  summarise(mean_z_score = mean(z_score_value), .groups = "drop") %>%
  mutate(replicate_group = "R2")

mean_normalised_df <- rbind(R1_normalised_df, R2_normalised_df)

wide_df <- mean_normalised_df %>%
  unite("group", condition, generation, replicate_group, sep = "_") %>%
  pivot_wider(names_from = group, values_from = mean_z_score, values_fill = list(z_score_value = 0))

#pca_data <- wide_df %>% select(-strain)
pca_result_overall <- prcomp(pca_data, scale. = TRUE)

variance_explained <- pca_result_overall$sdev^2 / sum(pca_result_overall$sdev^2)

metadata <- data.frame(
  group = colnames(pca_data),
  condition = gsub("^(.*?)_\\d+_R\\d+", "\\1", colnames(pca_data)),
  generation = as.numeric(gsub("^.*?_(\\d+)_R\\d+", "\\1", colnames(pca_data))),
  replicate = gsub("^.*?_\\d+_(R\\d+)", "\\1", colnames(pca_data))
)

pca_loadings_metadata <- cbind(as.data.frame(pca_result_overall$rotation), metadata)

overall_PCA <- ggplot(pca_loadings_metadata, aes(x = PC1, y = PC2, shape = replicate, colour = condition)) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ generation, labeller = as_labeller(function(x) paste("generation", x))) +
  labs(
    x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"),
    shape = "Replicate",
    color = "Condition"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold")
  )

# ==============================================================================
# DATA PREP
# ==============================================================================
#mean_zscore_df <- merged_pc %>%
#  group_by(strain, condition, generation, pool) %>%
#  summarise(mean_z_score = mean(z_score_value, na.rm = TRUE), .groups = "drop")
#
#merged_pc_G7G9 <- mean_zscore_df %>%
#  pivot_wider(
#    names_from = generation,
#    values_from = mean_z_score,
#    names_prefix = "mean_zscore_G"
#  )
#
## ==============================================================================
## varianceby_PC1p (by condition)
## ==============================================================================
#PC_condition <- c()
#variance_byPC1_cond <- c()
#
#for (each_condition in unique(merged_pc_G7G9$condition)) {
#  condition_data <- merged_pc_G7G9 %>% filter(condition == each_condition)  
#  
#  if(nrow(condition_data) > 0) {
#    pca_result_cond <- prcomp(condition_data[, c("mean_zscore_G7", "mean_zscore_G9")], center = TRUE, scale. = TRUE)
#    
#    std_devs <- pca_result_cond$sdev
#    var_explained_pc1_percent <- ((std_devs[1]^2) / sum(std_devs^2)) * 100
#    
#    PC_condition <- append(PC_condition, each_condition)
#    variance_byPC1_cond <- c(variance_byPC1_cond, var_explained_pc1_percent)
#  }
#}
#
#PCA_stats_cond <- data.frame(
#  condition = PC_condition,
#  variance_explained_by_PC1 = variance_byPC1_cond
#)
#
#varianceby_PC1p_condition <- ggplot(PCA_stats_cond, aes(x = condition, y = variance_explained_by_PC1)) +
#  geom_bar(stat = "identity", position = "dodge") +  
#  geom_text(aes(label = round(variance_explained_by_PC1, 2)), 
#            position = position_dodge(width = 0.9), 
#            vjust = -0.5, size = 5 / .pt) +  
#  labs(
#    y = "Variance explained by PC1"
#  ) +
#  theme_bw(base_size = 8) +
#  theme(
#    legend.position = "none",
#    axis.title = element_text(size = 8),
#    axis.text = element_text(size = 6),
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    plot.title = element_text(size = 8, face = "bold")
#  )
#
## ==============================================================================
##  varianceby_PC1p (by cluster)
## ==============================================================================
#merged_pc_G7G9$cluster <- ifelse(merged_pc_G7G9$condition %in% c("PA", "PV"), "non-MYBs", "MYBs")
#
#merged_pc_G7G9_median_clust <- merged_pc_G7G9 %>%
#  group_by(strain, cluster) %>%
#  summarise(
#    median_zscore_G7 = median(mean_zscore_G7, na.rm = TRUE),
#    median_zscore_G9 = median(mean_zscore_G9, na.rm = TRUE),
#    .groups = "drop"
#  )
#
#cluster_name_list <- c()
#variance_byPC1_clust <- c()
#
#for (each_cluster in unique(merged_pc_G7G9_median_clust$cluster)) {
#  cluster_data <- merged_pc_G7G9_median_clust %>% filter(cluster == each_cluster)
#  
#  if(nrow(cluster_data) > 0) {
#    pca_result_clust <- prcomp(cluster_data[, c("median_zscore_G7", "median_zscore_G9")], center = TRUE, scale. = TRUE)
#    
#    std_devs <- pca_result_clust$sdev
#    var_explained_pc1_percent <- ((std_devs[1]^2) / sum(std_devs^2)) * 100
#    
#    cluster_name_list <- append(cluster_name_list, each_cluster)
#    variance_byPC1_clust <- c(variance_byPC1_clust, var_explained_pc1_percent)
#  }
#}
#
#PCA_stats_clust <- data.frame(
#  cluster_name = cluster_name_list,
#  variance_explained_by_PC1 = variance_byPC1_clust
#)
#
#varianceby_PC1p_cluster <- ggplot(PCA_stats_clust, aes(x = cluster_name, y = variance_explained_by_PC1)) +
#  geom_bar(stat = "identity", position = "dodge") +  
#  geom_text(aes(label = round(variance_explained_by_PC1, 2)), 
#            position = position_dodge(width = 0.9), 
#            vjust = -0.5, size = 5 / .pt) +  
#  labs(
#    y = "Variance explained by PC1",
#    x = "Cluster"
#  ) +
#  theme_bw(base_size = 8) +
#  theme(
#    legend.position = "none",
#    axis.title = element_text(size = 8),
#    axis.text = element_text(size = 6),
#    plot.title = element_text(size = 8, face = "bold")
#  )

# ==============================================================================
# PANEL B: Broad-sense Heritability (formerly Panel E; retagged to B since the
# species tree panel that used to sit here has moved to main_figure1.R)
# ==============================================================================
h2_df_supp <- h2_df %>%
  filter(! condition == "MO") %>%
  mutate(broad_sense_heritability = broad_sense_heritability * 100)

h2_df_supp$method[h2_df_supp$method == "AUC z-score"] <- "AUC"
h2_df_supp$method[h2_df_supp$method == "regression slope z-score"] <- "regression"

total_h2_plot <- ggplot(h2_df_supp, aes(x = condition, y = broad_sense_heritability, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(broad_sense_heritability, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 5 / .pt) +  
  facet_grid(~method) +
  scale_fill_manual(values = pool_colors, name = "Pool") + 
  labs(
    y = "Broad-sense Heritability",
    x = "Condition"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    plot.title = element_text(size = 8, face = "bold"),  
    axis.title = element_text(size = 8),              
    axis.text = element_text(size = 6),                 
    legend.text = element_text(size = 6),               
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold")
  )

# ==============================================================================
# ASSEMBLE - Patchwork (2 full-width rows: A on top, B on bottom)
# ==============================================================================
# Panels B (varianceby_PC1p_condition) and C (varianceby_PC1p_cluster) removed
# from the layout below. Their generating code above is left untouched in case
# you want them back later. The species tree panel (formerly tagged B here)
# has moved to main_figure1.R.

# New layout: A (overall_PCA) spans the full top row,
# B (total_h2_plot) spans the full bottom row
layout_design <- "
  A
  B
"

supp_fig7 <- overall_PCA + total_h2_plot +
  # varianceby_PC1p_condition +   # <- commented out
  # varianceby_PC1p_cluster +     # <- commented out
  plot_layout(design = layout_design) +
  plot_annotation(
    tag_levels = 'A', 
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )

# Export with settings matching Main Figure 1
ggsave("supp_fig7.pdf", 
       plot = supp_fig7, 
       width = 180, 
       height = 220*2/3, 
       units = "mm", 
       device = cairo_pdf)

