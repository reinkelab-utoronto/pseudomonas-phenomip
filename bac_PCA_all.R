#overall PCA

merged_pc <- merged_pc %>%
  filter(! condition == "MO")

R1_normalised_df <- merged_pc %>%
  filter(replicate_group == "R1") %>%
  group_by(strain, condition, generation) %>%
  summarise(mean_z_score = mean(z_score_value), .groups = "drop")

R2_normalised_df <- merged_pc %>%
  filter(replicate_group == "R2") %>%
  group_by(strain, condition, generation) %>%
  summarise(mean_z_score = mean(z_score_value), .groups = "drop")

R1_normalised_df <- R1_normalised_df %>% mutate(replicate_group = "R1")
R2_normalised_df <- R2_normalised_df %>% mutate(replicate_group = "R2")
mean_normalised_df <- rbind(R1_normalised_df, R2_normalised_df)

wide_df <- mean_normalised_df %>%
  unite("group", condition, generation, replicate_group, sep = "_") %>%
  pivot_wider(names_from = group, values_from = mean_z_score, values_fill = list(z_score_value = 0))

pca_data <- wide_df %>% dplyr::select(-strain)

#pca_data_imputed <- pca_data %>%
#  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) #Every NA in a column is replaced with that column’s mean.

pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract the PCA scores (coordinates of the samples in the PCA space)
pca_scores <- as.data.frame(pca_result$x)

# Extract the variance explained by each principal component
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)


# Extract metadata from the column names of pca_data
metadata <- data.frame(
  group = colnames(pca_data),
  condition = gsub("^(.*?)_\\d+_R\\d+", "\\1", colnames(pca_data)),  # Extract condition
  generation = as.numeric(gsub("^.*?_(\\d+)_R\\d+", "\\1", colnames(pca_data))),  # Extract generation
  replicate = gsub("^.*?_\\d+_(R\\d+)", "\\1", colnames(pca_data))  # Extract replicate
)


# Extract PCA loadings (contributions of each variable to the principal components)
pca_loadings <- as.data.frame(pca_result$rotation)

# Add metadata to the PCA loadings
pca_loadings_metadata <- cbind(pca_loadings, metadata)

# Plot PCA results for the loadings (samples)
overall_PCA <- ggplot(pca_loadings_metadata, aes(x = PC1, y = PC2, shape = replicate,
                                                 colour = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  facet_wrap(~ generation) +  # Facet by generation
  labs(
    x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"),
    title = "PCA Loadings (Samples)",
    shape = "Replicate",
    color = "Condition"
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 20),  # Increases overall text size
    axis.title = element_text(size = 22),  # Axis labels
    axis.text = element_text(size = 14),   # Tick labels
    legend.text = element_text(size = 18), # Legend labels
    legend.title = element_text(size = 20), # Legend title
    strip.text = element_text(size = 20),  # Facet labels
    plot.title = element_text(size = 24, face = "bold") # Title
  )
ggsave(paste0("bac_overall_PCA.pdf"), plot =overall_PCA, units = 'mm', width = 297, height = 210)


cor_matrix <- cor(pca_data, use = "complete.obs")
# Open a PDF device
pdf("bac_overall_PCA_corr_matrix.pdf", width = 11.69, height = 8.27)

# Generate the heatmap
heatmap(cor_matrix)

# Close the PDF device
dev.off()


############################PCA by all

#average controls
mean_zscore_df <- merged_pc %>%
  group_by(strain, condition, generation) %>%
  summarise(mean_z_score = mean(z_score_value, na.rm = TRUE), .groups = "drop")

merged_pc_G7G9 <- mean_zscore_df %>%
  pivot_wider(
    names_from = generation,  # Use generation values (7, 9) as column names
    values_from = mean_z_score,  # Fill new columns with mean_z_score values
    names_prefix = "mean_zscore_G"  # Prefix column names with "mean_z_scoreG"
  )

merged_pc_G7G9_median_byall <- merged_pc_G7G9 %>%
  group_by(strain) %>%
  summarise(
    median_zscore_G7 = median(mean_zscore_G7, na.rm = TRUE),
    median_zscore_G9 = median(mean_zscore_G9, na.rm = TRUE),
    .groups = "drop"
  )

PC1_byall_df <- data.frame()
PC_condition <- c()
G7_loading <- c()
G9_loading <- c()
variance_byPC1 <- c()

pca_result <- prcomp(merged_pc_G7G9_median_byall[, 2:3], center = TRUE, scale. = TRUE)
merged_pc_G7G9_median_byall$PC1 <- pca_result$x[, 1]

PC1_byall_df <- rbind(PC1_byall_df, merged_pc_G7G9_median_byall)
# Standard deviations from PCA result
std_devs <- pca_result$sdev

loadings <- pca_result$rotation[, 1] # Loadings for PC1

# Variance explained by PC1
var_explained_pc1 <- (std_devs[1]^2) / sum(std_devs^2)
var_explained_pc1_percent <- var_explained_pc1 * 100 # Print as percentage
print(var_explained_pc1_percent)

G7_loading <- append(G7_loading, loadings[1][[1]])
G9_loading <- append(G9_loading, loadings[2][[1]])
variance_byPC1 <- c(variance_byPC1, var_explained_pc1_percent)

PCA_stats <- data.frame(G7_loading,
                        G9_loading,
                        variance_explained_byall_PC1 = variance_byPC1)
write.csv(merged_pc_G7G9, "merged_pc_G7G9_noMO.csv")

pca_byall_nemascan_input <- merged_pc_G7G9_median_byall %>%
  dplyr::select(strain, PC1)
colnames(pca_byall_nemascan_input) <- c("strain", "PC1_byall")

#varianceby_PC1p <- ggplot(PCA_stats, aes(x = condition, y = variance_explained_by_PC1)) +
#  geom_bar(stat = "identity", position = "dodge") +  
#  geom_text(aes(label = round(variance_explained_by_PC1, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +  # Add text annotations
#  labs(y = "Variance explained by PC1", title = paste("Variance explained by PC1 (PhenoMIP 6.0)", pathogen_name))
#ggsave(paste0("variance explained by PC1 (PhenoMIP 6.0) ", pathogen_name, ".pdf"), plot =varianceby_PC1p, units = 'mm', width = 297, height = 210)


#write_tsv(pca_byall_nemascan_input, "non_divergent_nemascan_input/PCA_byall/Pheno6_PC1_byall_nemascan_input.tsv")




