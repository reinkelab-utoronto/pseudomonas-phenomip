library(plyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggdendro)
library(reshape2)
library(grid)
library(cluster) 
library(gplots)
library(tidyr)
library(stringr)
library(pracma)
library(rstatix)
library(writexl)
library(readr)
library(scales)
library(cowplot)
library(pheatmap)
library(purrr)
library(broom)

pathogen_name <- "bac"
pool1 <- read.csv("upstream/PhenoMIP6.0_bac_P1/by.rep/growth.curve.data/GC.abundance.all.data.txt", 
                  sep = "\t")
pool2 <- read.csv("upstream/PhenoMIP6.0_bac_P2/by.rep/growth.curve.data/GC.abundance.all.data.txt", 
                  sep = "\t")
pool3 <- read.csv("upstream/PhenoMIP6.0_bac_P3/by.rep/growth.curve.data/GC.abundance.all.data.txt", 
                  sep = "\t")
pool_names <- c("pool1", "pool2", "pool3")

# Define the columns to remove
cols_to_remove <- c("PV.R08.E05", "M395.R14.E05", "MO.R03.E09", "MO.R02.E09", 
                    "MO.R03.E05", "MO.R02.E05", "MO.R03.E07", "MO.R02.E07",
                    "M395.R13.E05", "M114.R09.E05", "M117.R12.E05", "M117.R11.E05", 
                    "M114.R10.E05", "PV.R07.E05")

# Remove the specified columns
pool1 <- pool1 %>% 
  select(-all_of(cols_to_remove))

# Strains are named by their CaeNDR isotype throughout (e.g. "PB306", obtained as CaeNDR
# strain ECA259). NemaScan needs CaeNDR strain names, so in its input files any name that is
# an isotype but not itself a CaeNDR strain is replaced by that isotype's reference strain.
caendr_strain_data <- read.csv("20231213_c_elegans_strain_data.csv")
isotype_reference <- caendr_strain_data %>%
  filter(as.logical(isotype_ref_strain)) %>%
  dplyr::select(isotype, reference_strain = strain)
to_caendr_strain_name <- function(df) {
  df %>%
    left_join(isotype_reference, by = c("strain" = "isotype")) %>%
    mutate(strain = if_else(!strain %in% caendr_strain_data$strain & !is.na(reference_strain),
                            reference_strain, strain)) %>%
    dplyr::select(-reference_strain)
}
  


pools <- list(pool1, pool2, pool3)

final_freezing_plan <- read.csv("final_freezing_plan.csv")


pool_one <- c(1, 4, 9, 11, 16)
pool_two <- c(2,3,6,8, 15)
pool_three <- c(10,12,14)
pool_four <- c(5,7,13)
control <- c("CB4856", "JU1400", "VC20019")

pool <- c()
for (i in 1:nrow(final_freezing_plan)) {
  subpool <- final_freezing_plan[i, "Subpool"]
  
  if (subpool %in% pool_one) {
    pool <- append(pool, "pool1")
  } else if (subpool %in% pool_two) {
    pool <- append(pool, "pool2")
  } else if (subpool %in% pool_three) {
    pool <- append(pool, "pool3")
  } else if (subpool %in% pool_four) {
    pool <- append(pool, "pool4")
  } else if (identical(subpool, "Control")) {
    pool <- append(pool, "control")
  } else {
    pool <- append(pool, "removed")
  }
  
}
final_freezing_plan$pool <- pool

all_strains_pool <- c()
for (i in 1:nrow(pool1)) {
  strain <- pool1[i, 2]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
pool1_strains <- data.frame(pool1$strain.name, all_strains_pool)
pool1$pool <- all_strains_pool
write.csv(pool1, "GC.abundance.all.pool1.csv")

all_strains_pool <- c()
for (i in 1:nrow(pool2)) {
  strain <- pool2[i, 2]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
pool2_strains <- data.frame(pool2$strain.name, all_strains_pool)
pool2$pool <- all_strains_pool
write.csv(pool2, "GC.abundance.all.pool2.csv")

all_strains_pool <- c()
for (i in 1:nrow(pool3)) {
  strain <- pool3[i, 2]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
pool3_strains <- data.frame(pool3$strain.name, all_strains_pool)
pool3$pool <- all_strains_pool
write.csv(pool3, "GC.abundance.all.pool3.csv")


check_pool_strains <- list(pool1_strains, pool2_strains, pool3_strains)

divergent_strains <- c("CB4856",
                       "ECA347",
                       "ECA592",
                       "ECA594",
                       "ECA703",
                       "ECA705",
                       "ECA722",
                       "ECA733",
                       "ECA738",
                       "ECA742",
                       "ECA743",
                       "ECA744",
                       "ECA807",
                       "ECA923",
                       "ECA928",
                       "ECA930",
                       "JU3225",
                       "XZ2019")

good_strains <- read.csv("good_strains.csv")

source("./making_zscores.R")

# Compute R^2 values for auc vs abundance sum
auc_vs_abundancesum <- total_auc_df %>%
  group_by(pool) %>%
  summarize(r2 = summary(lm(auc_value ~ abundance_sum))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

auc_vs_abundancesum_p <- ggplot(total_auc_df, aes(x = auc_value, y = abundance_sum, color = condition)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) + # Regression line with full range
  facet_wrap(~ pool) +
  geom_text(data = auc_vs_abundancesum, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +
  labs(
    title = "Correlation between AUC and Sum of Abundance",
    x = "AUC Value",
    y = "Abundance Sum",
    color = "Condition"
  ) +
    theme(
      axis.title = element_text(size = 20),   # Increase axis title font size
      axis.text = element_text(size = 20),    # Increase axis text font size
      strip.text = element_text(size = 20),   # Increase facet label font size
      plot.title = element_text(size = 20, face = "bold"),  
      legend.position = "bottom", # Increase plot title size
      plot.subtitle = element_text(size = 20) # Increase subtitle size
    )
ggsave("auc_vs_abundancesum.pdf", plot =auc_vs_abundancesum_p, units = 'mm', width = 297, height = 210)

merged_pc <- merged_pc %>%
  filter(!strain %in% divergent_strains) %>%
  filter(!strain == "VC20019") #because VC20019 not in cendr

auc_median_sd_summary <- auc_all_conditions_sds %>%
  group_by(strain, pool) %>%
  summarize(median_sd = first(median_sd), .groups = 'drop')

#auc_sd <- ggplot(auc_all_conditions_sds, aes(x = strain, y = sd_of_auc_values)) +
#  geom_violin(alpha = 0.5) +  # Violin plot with transparency
#  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = condition)) +  # Overlay points with jitter
#  geom_point(data = auc_median_sd_summary, aes(y = median_sd), color = "red", size = 2.5) +  # Add median_sd as red points
#  labs(x = "Strain", y = "SD of AUC of each condition", title = "Variation in AUC Values by Strain") +
#  facet_grid(~pool, scale = "free") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#        legend.position = "bottom")  # Rotate x-axis labels for readability
#
#ggsave(paste0("auc_sd", ".pdf"), plot =auc_sd, units = 'mm', width = 297, height = 210)

#regression_median_sd_summary <- regression_all_conditions_sds %>%
#  group_by(strain, pool) %>%
#  summarize(median_sd = first(median_sd), .groups = 'drop')
#
#regression_sd <- ggplot(regression_all_conditions_sds, aes(x = strain, y = sd_of_slope_values)) +
#  geom_violin(alpha = 0.5) +  # Violin plot with transparency
#  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = condition)) +  # Overlay points with jitter
#  geom_point(data = regression_median_sd_summary, aes(y = median_sd), color = "red", size = 2.5) +  # Add median_sd as red points
#  labs(x = "Strain", y = "SD of regression slope of each condition", title = "Variation in Slope Values by Strain") +
#  facet_grid(~pool, scale ="free")+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#        legend.position = "bottom")  # Rotate x-axis labels for readability
#ggsave(paste0("regression_sd", ".pdf"), plot =regression_sd, units = 'mm', width = 297, height = 210)


#PCA
source("./bac_PCA_all.R") #PC1_byall_df is condition level PCA, not by pool
source("./bac_PCA_bycondition.R") #total_PC1_df is by condition
source("./bac_PCA_bycluster.R") #median between conditions
source("./bac_PCA_bycluster_mean.R") #mean between conditions



# Compute R^2 values for each condition and method (move to supp fig3)
#r2_data <- total_correlation_df %>%
#  group_by(method, condition) %>%
#  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
#  mutate(label = paste0("R² = ", round(r2, 3)))
#
#total_correlation_p <- ggplot(total_correlation_df, aes(x = z_R1, y = z_R2, color= pool)) +
#  geom_point() +
#  geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) + # Regression line with full range
#  facet_grid(condition ~ method) +
#  labs(
#    x = "Z-score (Replicate 1)",
#    y = "Z-score (Replicate 2)",
#    title = paste("Comparison of Z-scores calculated from AUC between Replicates (PhenoMIP 6.0)")) +
#  # Add R² labels
#  geom_text(data = r2_data, aes(x = Inf, y = -Inf, label = label),
#            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +
#  theme(
#    legend.position = "bottom",
#    theme(
#      axis.title = element_text(size = 20),   # Increase axis title font size
#      axis.text = element_text(size = 20),    # Increase axis text font size
#      strip.text = element_text(size = 20),   # Increase facet label font size
#      plot.title = element_text(size = 20, face = "bold"),  # Increase plot title size
#      plot.subtitle = element_text(size = 20) # Increase subtitle size
#    ))
#ggsave(paste0("phenomip6_", pathogen_name, "_","auc_regression_correlation.pdf"), plot =total_correlation_p, units = 'mm', width = 297, height = 210)


#broad sense heritability
total_h2_plot <- ggplot(h2_df, aes(x = condition, y = broad_sense_heritability, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(broad_sense_heritability, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 2.5) +  # Add text annotations
  facet_grid(~method) +
  labs(
    title = paste("Broad-sense Heritability by Condition for PhenoMIP 6.0", pathogen_name, sep = " "),
    y = "Broad-sense Heritability",
    x = "Condition"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 20, face = "bold"),  # Increase title text size
    axis.title.x = element_text(size = 20),              # Increase x-axis title text size
    axis.title.y = element_text(size = 20),              # Increase y-axis title text size
    axis.text = element_text(size = 20),                 # Increase axis tick label text size
    legend.text = element_text(size = 20),               # Increase legend label text size
    legend.title = element_text(size = 20)               # Increase legend title text size
  )
ggsave(paste0("phenomip6_", pathogen_name, "_","broadsense.pdf"), plot =total_h2_plot, units = 'mm', width = 297, height = 210)

total_correlation_df$average <- (total_correlation_df$z_R1 + total_correlation_df$z_R2) / 2

total_correlation_df_auc <- total_correlation_df %>%
  filter(method == "AUC z-score")
total_correlation_df_regression <- total_correlation_df %>%
  filter(method == "regression z-score")


#to combine with other pathogen data
auc_avg_zscore_df <- total_auc_df %>%
  group_by(pool, strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop") 

regression_avg_zscore_df <- total_regression_df %>%
  group_by(pool, strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop") 

auc_avg_zscore_df2 <- total_auc_df %>%
  group_by(pool, strain, condition) %>%
  summarise(value = mean(z_score, na.rm = TRUE), .groups = "drop") 
auc_avg_zscore_df2$method <- rep("AUC", nrow(auc_avg_zscore_df2))

regression_avg_zscore_df2 <- total_regression_df %>%
  group_by(pool, strain, condition) %>%
  summarise(value = mean(z_score, na.rm = TRUE), .groups = "drop") 
regression_avg_zscore_df2$method <- rep("slope", nrow(regression_avg_zscore_df2))

total_PC1_dfx <- data.frame(pool= total_PC1_df$pool, 
                            strain= total_PC1_df$strain, 
                            condition= total_PC1_df$condition, 
                            value= total_PC1_df$PC1, 
                            method= rep("PC1", nrow(total_PC1_df)))
all_3data <- rbind(rbind(auc_avg_zscore_df2, regression_avg_zscore_df2), total_PC1_dfx)
all_3data$condition[all_3data$condition == "MO"] <- "MO-bac"
write.csv(all_3data, "bac_all3data.csv", row.names = F)

# Create a list of data frames
data_list <- list(AUC = total_auc_df, average_AUC= total_correlation_df_auc, Regression = total_regression_df, 
                  average_regression = total_correlation_df_regression,
                  PC1_bypool = total_PC1_df)

# Write to an Excel file
write_xlsx(data_list, paste0("Pheno6_", pathogen_name, "_output.xlsx"))

# phenoMIP data for the 9 strains also tested in the mini-phenoMIP, read by the
# Figure 1E and Supp. Fig. S5 correlation scripts. These two files used to be made by hand.
#  - M117 is dropped (not tested in the mini-phenoMIP)
#  - M395 replicates R13/R14 are renamed R11/R12 to match the mini-phenoMIP numbering
#  - the control strain VC20019 is in all three pools: each replicate's values are averaged
#    across the pools, so VC20019 has one value per replicate like the other strains
minipheno_strains <- c("ED3005", "JU1530", "JU2250", "JU2522", "JU2610",
                       "JU2907", "JU782", "VC20019", "WN2050")
make_9strain_table <- function(df) {
  df <- df %>%
    filter(strain %in% minipheno_strains, condition != "M117") %>%
    mutate(replicate = case_when(
      condition == "M395" & replicate == "R13" ~ "R11",
      condition == "M395" & replicate == "R14" ~ "R12",
      TRUE ~ as.character(replicate)))
  vc20019 <- df %>%
    filter(strain == "VC20019") %>%
    group_by(strain, condition, replicate) %>%
    summarise(across(where(is.numeric), mean),
              pool = "pool1-3 mean", method = first(method), .groups = "drop")
  bind_rows(filter(df, strain != "VC20019"), vc20019) %>%
    select(all_of(names(df)))
}
write.csv(make_9strain_table(total_auc_df), "pheno6_9strains_auc.csv", row.names = FALSE)
write.csv(make_9strain_table(total_regression_df), "pheno6_9strains_regression.csv", row.names = FALSE)

total_PC1_df3 <- total_PC1_df
all_strains_pool <- c()
for (i in 1:nrow(total_PC1_df)) {
  strain <- total_PC1_df$strain[i]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
total_PC1_df$pool <- all_strains_pool

#ranked z-scores and PC1
for (i in 2:length(condition)) {
  condition_name <- condition[i]
  print(condition[i])
  auc_each_condition_zscore <- total_correlation_df_auc %>%
    filter(condition == condition_name)
  
  auc_long <- auc_each_condition_zscore %>%
    pivot_longer(cols = c(z_R1, z_R2, average), names_to = "metric", values_to = "value")
  
  auc_p <- ggplot(auc_long, aes(x = reorder(strain, -value), y = value, color = metric)) +
    geom_point(size = 2) +  # Plot points
    geom_hline(yintercept = c(10, -10), linetype = "dashed", color = "red") +  # Add horizontal lines
    facet_wrap(~pool, scales = "free_x") +  # Facet by pool
    labs(
      x = "Strain",
      y = "Z-Score",
      color = "Metric",
      title = "AUC z-scores by Pool"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10), 
          legend.position = "bottom")  # Rotate x-axis labels for better readability
  
  
  auc_plot_name <- paste0("phenomip6.0_nemascan_input/", condition_name, "/Pheno6_", condition_name, "_auc_zscore_ordered.pdf")
  dir.create(dirname(auc_plot_name), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(auc_p, filename = auc_plot_name, width = 45, height = 21, 
         units = "cm", device='pdf')  
  
  regression_each_condition_zscore <- total_correlation_df_regression %>%
    filter(condition == condition_name)
  
  regression_long <- regression_each_condition_zscore %>%
    pivot_longer(cols = c(z_R1, z_R2, average), names_to = "metric", values_to = "value")
  
  regression_p <- ggplot(regression_long, aes(x = reorder(strain, -value), y = value, color = metric)) +
    geom_point(size = 2) +  # Plot points
    facet_wrap(~pool, scales = "free_x") +  # Facet by pool
    labs(
      x = "Strain",
      y = "Z-Score",
      color = "Metric",
      title = "regression z-scores by Pool"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1,  size = 10), 
          legend.position = "bottom")  # Rotate x-axis labels for better readability
  
  regression_plot_name <- paste0("phenomip6.0_nemascan_input/", condition_name, "/Pheno6_", condition_name, "_regression_zscore_ordered.pdf")
  dir.create(dirname(regression_plot_name), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(regression_p, filename = regression_plot_name, width = 45, height = 21, 
         units = "cm", device='pdf')
  
  PC1_each_condition <- total_PC1_df %>%
    filter(condition == condition_name)
  
  PC1_each_condition_long <- PC1_each_condition %>%
    pivot_longer(cols = c(mean_zscore_G7, mean_zscore_G9, PC1), names_to = "metric", values_to = "value")
  
  pc1_p <- ggplot(PC1_each_condition_long, aes(x = reorder(strain, -value), y = value, color = metric)) +
    geom_point(size = 2) +  # Plot points
    labs(
      x = "Strain",
      y = "PC1",
      color = "Metric",
      title = "PC1"
    ) +
    facet_grid(~pool, scales = "free_x") +  # Facet by pool
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 9), 
          legend.position = "bottom")  # Rotate x-axis labels for better readability
  
  pc1_plot_name <- paste0("phenomip6.0_ranked_plots/", condition_name, "/Pheno6_", condition_name, "_PC1_ordered.pdf")
  dir.create(dirname(pc1_plot_name), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(pc1_p, filename = pc1_plot_name, width = 40, height = 21, 
         units = "cm", device='pdf')  
  
}

#control strains correlation
#write.csv(total_correlation_df, "total_correlation_df.csv", row.names = F)
#write.csv(total_auc_df, "total_auc_df.csv", row.names = F)
#write.csv(total_regression_df, "total_regression_df.csv", row.names = F)

source("./pheno6_controls_correlation.R")




 #auc_heatmap <- ggplot(total_correlation_df_auc, aes(x = condition, y = strain, fill = average)) +
#  geom_tile() +  
#  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name = "Average modified Z-score") +
#  labs(x = "Condition", y = "Strain", title = "Heatmap of Average Z-Scores", subtitle = "Calculated using AUC") +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    legend.position = "bottom"
#  )
#
#final_plot <- tree_plot | auc_heatmap
#final_plot


# For creating the heatmaps
total_auc_df_no_control <- total_auc_df %>%
  filter(!strain %in% c("CB4856", "VC20019", "JU1400"))
#write.csv(total_correlation_df_auc_no_control, "bac_total_correlation_df_auc_no_control.csv", row.names = F)

total_correlation_df_auc_no_control <- total_auc_df_no_control %>%
  group_by(pool, strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop")

source("./pheno6_bac_auc_heatmap_dendro.R")



auc_heatmap <- ggplot(total_correlation_df_auc_no_control, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    limits = c(-10, 10),  # Set explicit scale limits
    oob = squish,  # Cap values beyond ±10
    name = "Average modified Z-score"
  ) +
  facet_wrap(~pool, scales = "free") +  # Facet by the 'pool' column
  labs(
    x = "Condition", 
    y = "Strain", 
    title = "Heatmap of average z-scores", 
    subtitle = "Calculated using AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    legend.position = "bottom"  # Move legend to the bottom
  )

total_regression_df_no_control <- total_regression_df %>%
  filter(!strain %in% c("CB4856", "VC20019", "JU1400"))
#write.csv(total_correlation_df_regression_no_control, "bac_total_correlation_df_regression_no_control.csv", row.names = F)

total_correlation_df_regression_no_control <- total_regression_df_no_control %>%
  group_by(pool, strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop")
#write.csv(total_correlation_df_regression_no_control, "bac_total_correlation_df_regression_no_control.csv", row.names = F)

source("./pheno6_bac_regression_heatmap_dendro.R")

regression_heatmap <- ggplot(total_correlation_df_regression, aes(x = condition, y = strain, fill = average)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "Average modified Z-score") +  # Adjust the color gradient and legend title
  facet_wrap(~pool, scale = "free") +  # Facet by the 'pool' column
  labs(x = "Condition", y = "Strain", title = "Heatmap of average z-scores", subtitle = "Calculated using regression slope") +  # Label the axes
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    legend.position = "bottom"  # Move legend to the bottom
  )



total_PC1_df_no_control <- total_PC1_df %>%
  filter(!strain %in% c("CB4856", "VC20019", "JU1400"))
#write.csv(total_PC1_df_no_control, "bac_total_PC1_df_no_control.csv", row.names = F)

source("./pheno6_bac_PC1_heatmap_dendro.R")

pc1_heatmap <- ggplot(total_PC1_df, aes(x = condition, y = strain, fill = PC1)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "PC1") +  # Adjust the color gradient and legend title
  facet_wrap(~pool, scale = "free") +  # Facet by the 'pool' column
  labs(x = "Condition", y = "Strain", title = "Heatmap of PC1 values") +  # Label the axes
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    legend.position = "bottom"  # Move legend to the bottom
  )

total_PC1_df_control <- total_PC1_df3 %>%
  filter(strain %in% c("CB4856", "VC20019", "JU1400")) %>%
  select(pool, strain, condition, avg_z_score = PC1)
total_PC1_df_control$method <- rep("PC1", nrow(total_PC1_df_control))

filtered_auc_avg_zscore_df$method <- rep("AUC", nrow(filtered_auc_avg_zscore_df))
filtered_regression_avg_zscore_df$method <- rep("regression slope", nrow(filtered_regression_avg_zscore_df))
filtered_df <- rbind(filtered_auc_avg_zscore_df, filtered_regression_avg_zscore_df)
#filtered_df <- rbind(filtered_df, total_PC1_df_control)


colnames(filtered_df) <- c("pool", "strain", "condition", "metric", "method")

filtered_df$condition <- factor(
  filtered_df$condition,
  levels = c("MO", "M114", "M117", "M395", "PA", "PV")
)

filtered_df$method <- factor(
  filtered_df$method,
  levels = c("AUC", "regression slope")
)





#hclust

#AUC
source("./auc_hclust.R") #removed controls in auc_avg_zscore_df

#regression
source("./regression_hclust.R") #removed controls in regression_avg_zscore_df


auc_avg_zscore_df <- total_auc_df %>%
  group_by(strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop") 
  
regression_avg_zscore_df <- total_regression_df %>%
  group_by(strain, condition) %>%
  summarise(average = mean(z_score, na.rm = TRUE), .groups = "drop") 

PC1_df_avg_zscore_df <- total_PC1_df %>%
  group_by(strain, condition) %>%
  summarise(average = mean(PC1, na.rm = TRUE), .groups = "drop") 


#nemascan input
#for (i in 1:length(condition)) {
#  condition_name <- condition[i]
#  
#  auc_avg_zscore_df <- auc_avg_zscore_df[!auc_avg_zscore_df$strain == "VC20019",]
#  
#  auc_each_condition_zscore <- auc_avg_zscore_df %>%
#    filter(condition == condition_name) %>%
#    select(strain, average)
#  
#  colnames(auc_each_condition_zscore) <- c("strain", paste0(condition_name, "_avg_auc_zscore"))
#  
#  regression_avg_zscore_df <- regression_avg_zscore_df[!regression_avg_zscore_df$strain == "VC20019",]
#  regression_each_condition_zscore <- regression_avg_zscore_df %>%
#    filter(condition == condition_name) %>%
#    select(strain, average)
#  
#  colnames(regression_each_condition_zscore) <- c("strain", paste0(condition_name, "_avg_slope_zscore"))
#  
#  
#  PC1_df_avg_zscore_df <- PC1_df_avg_zscore_df[!PC1_df_avg_zscore_df$strain == "VC20019",]
#
#  PC1_input <- PC1_df_avg_zscore_df %>%
#    filter(condition == condition_name) %>%
#    select(strain, average)
#  
#  colnames(PC1_input) <- c("strain", paste0(condition_name, "_avg_PC1"))
#  
#  
#  output_folder <- paste0("phenomip6.0_nemascan_input/", condition_name)
#  dir.create(output_folder, recursive = TRUE)
#  
#  auc_file_name <- paste0(output_folder,"/Pheno6_", condition_name, "_all_pools_", "avg_auc_zscore_nemascan_input.tsv")
#  regression_file_name <- paste0(output_folder, "/Pheno6_", condition_name, "_all_pools_", "avg_regression_zscore_nemascan_input.tsv")
#  PC1_file_name <- paste0(output_folder, "/Pheno6_", condition_name, "_all_pools_", "PC1_nemascan_input.tsv")
#  
#  write_tsv(auc_each_condition_zscore, auc_file_name)
#  write_tsv(regression_each_condition_zscore, regression_file_name)
#  write_tsv(PC1_input, PC1_file_name)
#  
#}




non_divergent_auc_avg_zscore_df <- auc_avg_zscore_df %>%
  filter(!strain %in% divergent_strains) %>%
  filter(! strain == "VC20019")

non_divergent_auc_avg_zscore_df_noMO <- auc_avg_zscore_df %>%
  filter(!strain %in% divergent_strains) %>%
  filter(! strain == "VC20019") %>%
  filter(! condition == "MO")

non_divergent_regression_avg_zscore_df <- regression_avg_zscore_df %>%
  filter(!strain %in% divergent_strains) %>%
  filter(! strain == "VC20019")

non_divergent_regression_avg_zscore_df_noMO <- regression_avg_zscore_df %>%
  filter(!strain %in% divergent_strains) %>%
  filter(! strain == "VC20019") %>%
  filter(! condition == "MO")

non_divergent_PC1_df_avg_zscore_df <- PC1_df_avg_zscore_df %>%
  filter(!strain %in% divergent_strains)

#non divergent nemascan input
for (i in 2:length(condition)) {
  condition_name <- condition[i]
  
  auc_each_condition_zscore <- non_divergent_auc_avg_zscore_df %>%
    filter(condition == condition_name) %>%
    select(strain, average)
  
  colnames(auc_each_condition_zscore) <- c("strain", paste0(condition_name, "_avg_auc_zscore"))
  
  regression_each_condition_zscore <- non_divergent_regression_avg_zscore_df %>%
    filter(condition == condition_name) %>%
    select(strain, average)
  
  colnames(regression_each_condition_zscore) <- c("strain", paste0(condition_name, "_avg_slope_zscore"))
  

  PC1_input <- non_divergent_PC1_df_avg_zscore_df %>%
    filter(condition == condition_name) %>%
    select(strain, average)
  
  colnames(PC1_input) <- c("strain", paste0(condition_name, "_avg_PC1"))
  
  
  output_folder <- paste0("non_divergent_nemascan_input/", condition_name)
  dir.create(output_folder, recursive = TRUE)
  
  auc_file_name <- paste0(output_folder,"/Pheno6_", condition_name, "_all_pools_", "avg_auc_zscore_nemascan_input.tsv")
  regression_file_name <- paste0(output_folder, "/Pheno6_", condition_name, "_all_pools_", "avg_regression_zscore_nemascan_input.tsv")
  PC1_file_name <- paste0(output_folder, "/Pheno6_", condition_name, "_all_pools_", "PC1_nemascan_input.tsv")
  
  # NemaScan needs CaeNDR strain names (see to_caendr_strain_name() at the top)
  write_tsv(to_caendr_strain_name(auc_each_condition_zscore), auc_file_name)
  write_tsv(to_caendr_strain_name(regression_each_condition_zscore), regression_file_name)
  write_tsv(to_caendr_strain_name(PC1_input), PC1_file_name)
  
}


#heritability input
source("./auc_heritability_input.R")
source("./slope_heritability_input.R")
#no input for PC1 because lack of "technical replicates"- PCA will group strains so each strain only has 1 PC1 value

#worldmap
source("./worldmapv2.R")


#correlation between genetic distance and z-scores
source("./genetic_distance_zscorediff_auc.R")
source("./genetic_distance_zscorediff_regression.R")
source("./genetic_distance_PC1diff.R")



#main figures
# The mini-phenoMIP analysis writes minipheno_auc.csv and minipheno_regression.csv,
# which Figure 1E reads, so it has to run before Figure 1 (supp_fig4 runs it again later)
source("minipheno/mini_pheno_analysis.R", local = new.env())
source("main_figure1.R")
source("main_figure2.R")
source("main_figure3.R") #M114 BSA full stats
source("main_figure4.R") #PA14 BSA half stats
source("main_figure5.R") #JU2250-M395 BSA stats, chr X RILs
source("main_figure6.R") #strain level FISH for M114, PA14 and M395

#supp figures
source("supp_fig1.R") #% male
source("supp_fig2.R") #correlation of MIPs
source("supp_fig3.R") #correlation of z-scores
source("supp_fig4_reordered.R") #mini pheno correlation between replicates
source("supp_fig5_reordered.R") #mini pheno correlation between pheno6 and minipheno, by condition
source("supp_fig6_reodered.R") #full trees
source("supp_fig7_reordered.R") #PCA, H^2
source("supp_fig8_reordered.R") #hclust
source("supp_fig9_reorderd.R") #location
#supp 10-14 in the worldmap.R

source("supp_fig15_reordered.R") #bac+mic full tree
source("supp_fig16_reordered.R") #bac+mic hclust
source("supp_fig17_reordered.R") #spearman matrix
setwd("gwas")
source("supp_fig18.R") #supp 18 is nemascan results in gwas folder
setwd("..")

#supp 19 is BSA biorender
source("supp_fig20_reordered.R") #BSA stats low for M114
source("supp_fig21_reordered.R") #M114 BSA N2 allele f and log10 odds ratio
source("supp_fig22.R") #wild isolates tested for M114 candidate gene
source("supp_fig23.R") #BSA stats low of PA14
source("supp_fig24.R") #PA14 BSA N2 allele f and log10 odds ratio
source("supp_fig25.R") #M114 and PA14 RILs all chr 
source("supp_fig26_reordered.R") #specificity
source("supp_fig27_reordered.R") #BSA stats for EG4946, chr X RILs
source("supp_fig28_reordered.R") #M395 BSA N2 allele f and log10 odds ratio
source("supp_fig29_reordered.R") #M395 JU2250 and EG4946 RILs all chr
#supp 30 is FISH stain bacteria
source("supp_fig31.R") #FISH body size
#supp 32 is PA14 FISH
#supp 33 is M114 FISH
#supp 34 is M395 FISH

source("supp_data.R")
