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

ab <- read.csv("minipheno/minipheno_graphing/by.rep/growth.curve.data/GC.abundance.all.data.txt",
               sep = "\t")
strain_list <- ab$strain.name
conditions <- data.frame(types= c("MO", "PA", "PV", "M114", "M395"),
                         replicate = c(4, 2, 2, 2, 2))
generations <- c("E02", "E04", "E06")


strain <- c()
sample <- c()
abundance <- c()
for (k in 1:length(generations)) {
  generation_name <- generations[k]
  each_generation <- ab %>%
    dplyr::select(2, contains(generation_name))
  each_generation$strain.name <- ab$strain.name
  
  for (j in 1:length(ab$strain.name)) {
    each_strain <- ab$strain.name[j]
    each_strain_gen <- each_generation %>%
      filter(strain.name == each_strain) 
    
    each_strain_gen_abundance <- each_strain_gen %>%
      dplyr::select(2:ncol(each_strain_gen))
    
    
    strain <- append(strain, rep(each_strain, length(as.numeric(each_strain_gen_abundance))))
    sample <- append(sample, (colnames(each_strain_gen)[2:length(colnames(each_strain_gen))]))
    abundance <- append(abundance, as.numeric(each_strain_gen_abundance))
    
  }
}
final_df <- data.frame(strain, name = sample, abundance)


final_df$condition <- sapply(final_df$name, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
final_df$generation <- sapply(final_df$name, function(x) strsplit(x, ".", fixed = TRUE)[[1]][3])
final_df$replicate <- sapply(final_df$name, function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])

generation_name <- c()
for (i in 1:nrow(final_df)) {
  gen <- final_df[i, 5]
  if (identical(gen, "E04")) {
    generation_name <- append(generation_name, 4)
  } else if (identical(gen, "E06")) {
    generation_name <- append(generation_name, 6)
  } else if (identical(gen, "E02")) {
    generation_name <- append(generation_name, 2)
  } else if (identical(gen, "E00")) {
    generation_name <- append(generation_name, 0)
  }
}
final_df$generation_name <- generation_name
#total_abundance_df <- rbind(total_abundance_df, final_df)


rep_list <- unique(final_df$replicate)
condition_list <- unique(final_df$condition)


#AUC
abundance_df <- data.frame()
for (i in 1:length(ab$strain.name)) {
  each_strain <- ab$strain.name[i]
  each_strain_df <- final_df %>%
    filter(strain == each_strain)
  
  for (j in 1:length(rep_list)) {
    each_rep <- rep_list[j]
    each_strain_rep <- each_strain_df %>%
      filter(replicate == each_rep) 
    
    each_strain_rep$condition <- sapply(each_strain_rep$name, 
                                        function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
    
    df <- data.frame(strain = c(rep(each_strain, (nrow(each_strain_rep)+1))),
                     replicate = c(rep(unique(each_strain_rep$replicate), (nrow(each_strain_rep)+1))),
                     condition = c(rep(unique(each_strain_rep$condition), (nrow(each_strain_rep)+1))),
                     generation = c(0, each_strain_rep$generation_name), 
                     abundance = c(0, each_strain_rep$abundance))
    
    
    abundance_df <- rbind(abundance_df, df)
  }
}


strain_auc <- c()
condition_auc <- c()
replicate_auc <- c()
sumy <- c()
auc <- c()
for (i in 1:length(ab$strain.name)) {
  each_strain <- ab$strain.name[i]
  each_strain_df <- final_df %>%
    filter(strain == each_strain)
  
  for (j in 1:length(rep_list)) {
    each_rep <- rep_list[j]
    each_strain_rep <- each_strain_df %>%
      filter(replicate == each_rep) 
    
    df <- data.frame(generation = c(0, each_strain_rep$generation_name), 
                     abundance = c(0, each_strain_rep$abundance))
    
    auc_value <- trapz(y =df$abundance, 
                       x=df$generation)  
    
    strain_auc <- append(strain_auc, unique(each_strain_rep$strain))
    condition_auc <- append(condition_auc, unique(each_strain_rep$condition))
    replicate_auc <- append(replicate_auc, unique(each_strain_rep$replicate))
    auc <- append(auc, auc_value)
    sumy <- append(sumy, sum(df$abundance))
    
  }  
}
auc_df <- data.frame(strain = strain_auc,
                     condition = condition_auc,
                     replicate = replicate_auc,
                     auc_value = auc,
                     abundance_sum = sumy)

custom_color_palette <- c(
  "M114" = "#F8766D",
  "M117" = "#B79F00",
  "M395" = "#00BA38",
  "MO"   = "#00BFC4",
  "PA"   = "#619CFF",
  "PV"   = "#F564E3"
)


output_folder <- "minipheno_growth_curves"
dir.create(output_folder, showWarnings = FALSE)
# Define a custom set of shapes for your replicates
custom_shapes <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

for (i in 1:length(ab$strain.name)) {
  each_strain <- ab$strain.name[i]
  each_strain_plot <- abundance_df %>%
    filter(strain == each_strain)
  
  growth_curve <- ggplot(each_strain_plot, aes(x = generation, y = abundance, group = replicate, 
                                               color = condition, shape = replicate)) +
    geom_point() +
    geom_line()+
    scale_shape_manual(values = custom_shapes) +
    scale_color_manual(values = custom_color_palette) + 
    labs(title = paste(each_strain, "growth curve", sep = " "))
  
  plot_name <- paste(output_folder, "/", each_strain, "_growth_curve.pdf", sep = "")
  
  ggsave(plot_name, 
         plot =growth_curve, units = 'mm', 
         width = 297, height = 210)
  
}










z_score <- c()
median_sd <- c()
auc_all_conditions_sds <- data.frame()
for (i in 1:length(strain_list)) {
  each_strain_auc <- auc_df %>%
    filter(strain == strain_list[i])
  
  std_dev_per_condition <- aggregate(auc_value ~ condition, each_strain_auc, sd)
  colnames(std_dev_per_condition) <- c("condition", "sd_of_auc_values")
  std_dev_per_condition$strain <- rep(strain_list[i], nrow(std_dev_per_condition))

  median_sd_all_conditions <- median(std_dev_per_condition$sd_of_auc_values)
  
  std_dev_per_condition$median_sd <- rep(median_sd_all_conditions, nrow(std_dev_per_condition))
  auc_all_conditions_sds <- rbind(auc_all_conditions_sds, std_dev_per_condition)
  
  # Extract data for the control group (MO)
  control_auc <- each_strain_auc$auc_value[each_strain_auc$condition == "MO"]
  
  # Calculate the mean and standard deviation of the control group
  control_median <- median(control_auc)
  #control_sd <- sd(control_auc)
  
  # Calculate z-scores for each value in the data frame
  z_score_value <- (each_strain_auc$auc_value - control_median) / median_sd_all_conditions
  median_sd_value <- rep(median_sd_all_conditions, length(z_score_value))
  z_score <- append(z_score, z_score_value)
  median_sd <- append(median_sd, median_sd_value)
}
auc_df$z_score <- z_score
auc_df$median_auc_sd <- median_sd
auc_df$method <- rep("AUC z-score", nrow(auc_df))

auc_df <- auc_df %>%
  filter(is.finite(z_score))


#AUC replicates
auc_df_rep <- auc_df %>%
  mutate(replicate_group = case_when(
    replicate %in% c("R01", "R02", "R05", "R07", "R09", "R11") ~ "R1",
    replicate %in% c("R03", "R04", "R06", "R08", "R10", "R12") ~ "R2",
    TRUE ~ NA_character_  # If there are other replicates not specified
  ))


auc_df_R1 <- auc_df_rep %>% filter(replicate_group == "R1") %>% dplyr::select(strain, condition, auc_value, z_score)
auc_df_R2 <- auc_df_rep %>% filter(replicate_group == "R2") %>% dplyr::select(strain, condition, auc_value, z_score)

if (identical(auc_df_R1$strain, auc_df_R2$strain)) {
  print("strain okay!")
} else {
  print("strain not okay!")
}

if (identical(auc_df_R1$condition, auc_df_R2$condition)) {
  print("condition okay!")
} else {
  print("condition not okay!")
}

colnames(auc_df_R1) <- c("strain", "condition", "measurement_R1", "z_R1")
colnames(auc_df_R2) <- c("strain.2", "condition.2", "measurement_R2", "z_R2")

auc_merged_data <- cbind(auc_df_R1, auc_df_R2)

if (identical(auc_merged_data$strain, auc_merged_data$strain.2)) {
  print("strain okay after merge!")
} else {
  print("strain not okay after merge!")
}

if (identical(auc_merged_data$condition, auc_merged_data$condition.2)) {
  print("condition okay after merge!")
} else {
  print("condition not okay after merge!")
}

total_correlation_df <- data.frame()
auc_merged_data$method <- rep("AUC z-score", nrow(auc_merged_data))
total_correlation_df <- rbind(total_correlation_df, auc_merged_data)

# Compute R^2 values for each condition
auc_r2_data <- auc_merged_data %>%
  group_by(condition) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))


# Create a correlation plot faceted by condition
auc_correlation_p <- ggplot(auc_merged_data, aes(x = z_R1, y = z_R2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue", fullrange = TRUE, size = 0.2) + # Regression line with full range
  facet_wrap(~ condition) +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)",
    title = "Comparison of Z-scores calculated from AUC between Replicates (PhenoMIP 6.0)") +
  # Add R² labels
  geom_text(data = auc_r2_data, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +
  theme(
    axis.title = element_text(size = 16),   # Increase axis title font size
    axis.text = element_text(size = 14),    # Increase axis text font size
    strip.text = element_text(size = 15),   # Increase facet label font size
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size
    plot.subtitle = element_text(size = 16) # Increase subtitle size
  )
ggsave("mini_pheno_auc_correlation.pdf", plot =auc_correlation_p, units = 'mm', width = 297, height = 210)



#AUC violin
mo_data <- auc_df[auc_df$condition == "MO", ]


auc_violin_p <- ggplot(mo_data, aes(x = strain, y = z_score)) +
  geom_violin(fill = "lightblue") +
  geom_jitter(data = auc_df, aes(x = strain, y = z_score, 
                                 color = condition), 
              width = 0.15, size = 2) +
  labs(x = "Strain", y = "Z-Score",
       title = "mini-PhenoMIP Violin Plot") +
  scale_color_manual(values = custom_color_palette)+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title=element_text(size=14, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggsave("mini_pheno_auc_violin.pdf", plot =auc_violin_p, units = 'mm', width = 297, height = 210)



#regression
regression_strain <- c()
regression_rep <- c()
regression_condition <- c()
regression_slope <- c()
regression_intercept <- c()
for (i in 1:length(ab$strain.name)) {
  each_strain <- ab$strain.name[i]
  each_strain_df <- final_df %>%
    filter(strain == each_strain)
  
  for (j in 1:length(rep_list)) {
    each_rep <- rep_list[j]
    each_strain_rep <- each_strain_df %>%
      filter(replicate == each_rep) 
    condition_name <- unique(each_strain_rep$condition)
    
    df <- data.frame(generation = c(0, each_strain_rep$generation_name), 
                     abundance = c(0, each_strain_rep$abundance))
    
    linear_model <- lm(df$abundance ~ df$generation)
    intercept <- coef(summary(linear_model))[1, 1]
    slope <- coef(summary(linear_model))[2, 1]
    
    regression_strain <- append(regression_strain, each_strain)
    regression_rep <- append(regression_rep, each_rep)
    regression_condition <- append(regression_condition, condition_name)
    regression_slope <- append(regression_slope, slope)
    regression_intercept <- append(regression_intercept, intercept)
    
    
    #regression_value <- trapz(y = df$abundance, 
    #                   x=df$generation)  
    #
    #strain_regression <- append(strain_regression, unique(each_strain_rep$strain))
    #condition_regression <- append(condition_regression, unique(each_strain_rep$condition))
    #replicate_regression <- append(replicate_regression, unique(each_strain_rep$replicate))
    #regression <- append(regression, regression_value)
  }  
}
regression_df <- data.frame(strain = regression_strain,
                            replicate= regression_rep,
                            condition = regression_condition,
                            slope = regression_slope,
                            intercept = regression_intercept)

z_score <- c()
median_sd <- c()
for (i in 1:length(strain_list)) {
  each_strain_regression <- regression_df %>%
    filter(strain == strain_list[i])
  
  std_dev_per_condition <- aggregate(slope ~ condition, each_strain_regression, sd)
  colnames(std_dev_per_condition) <- c("condition", "sd_of_slope_values")
  std_dev_per_condition$strain <- rep(strain_list[i], nrow(std_dev_per_condition))

  median_sd_all_conditions <- median(std_dev_per_condition$sd_of_slope_values)
  
  std_dev_per_condition$median_sd <- rep(median_sd_all_conditions, nrow(std_dev_per_condition))

  # Extract data for the control group (MO)
  control_regression <- each_strain_regression$slope[each_strain_regression$condition == "MO"]
  
  # Calculate the mean and standard deviation of the control group
  control_median <- median(control_regression)
  #control_sd <- sd(control_regression)
  
  # Calculate z-scores for each value in the data frame
  z_score_value <- (each_strain_regression$slope - control_median) / median_sd_all_conditions
  median_sd_value <- rep(median_sd_all_conditions, length(z_score_value))
  z_score <- append(z_score, z_score_value)
  median_sd <- append(median_sd, median_sd_value)
}
regression_df$z_score <- z_score
regression_df$median_auc_sd <- median_sd
regression_df$method <- rep("regression slope z-score", nrow(regression_df))


regression_df <- regression_df %>%
  filter(is.finite(z_score))

#regression replicates
regression_df_rep <- regression_df %>%
  mutate(replicate_group = case_when(
    replicate %in% c("R01", "R02", "R05", "R07", "R09", "R11") ~ "R1",
    replicate %in% c("R03", "R04", "R06", "R08", "R10", "R12") ~ "R2",
    TRUE ~ NA_character_  # If there are other replicates not specified
  ))


regression_df_R1 <- regression_df_rep %>% filter(replicate_group == "R1") %>% dplyr::select(strain, condition, slope, z_score)
regression_df_R2 <- regression_df_rep %>% filter(replicate_group == "R2") %>% dplyr::select(strain, condition, slope, z_score)

if (identical(regression_df_R1$strain, regression_df_R2$strain)) {
  print("strain okay!")
} else {
  print("strain not okay!")
}

if (identical(regression_df_R1$condition, regression_df_R2$condition)) {
  print("condition okay!")
} else {
  print("condition not okay!")
}

colnames(regression_df_R1) <- c("strain", "condition", "measurement_R1", "z_R1")
colnames(regression_df_R2) <- c("strain.2", "condition.2", "measurement_R2", "z_R2")

regression_merged_data <- cbind(regression_df_R1, regression_df_R2)

if (identical(regression_merged_data$strain, regression_merged_data$strain.2)) {
  print("strain okay after merge!")
} else {
  print("strain not okay after merge!")
}

if (identical(regression_merged_data$condition, regression_merged_data$condition.2)) {
  print("condition okay after merge!")
} else {
  print("condition not okay after merge!")
}


regression_merged_data$method <- rep("regression z-score", nrow(regression_merged_data))

# Compute R^2 values for each condition
regression_r2_data <- regression_merged_data %>%
  group_by(condition) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

#weird_PA <- regression_merged_data[regression_merged_data$condition == "PA",] %>%
#  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
#  mutate(label = paste0("R² = ", round(r2, 3)))
#
#weird_PV <- regression_merged_data[regression_merged_data$condition == "PV",] %>%
#  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
#  mutate(label = paste0("R² = ", round(r2, 3)))




# Create a correlation plot faceted by condition
regression_correlation_p <- ggplot(regression_merged_data, aes(x = z_R1, y = z_R2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue", fullrange = TRUE, size = 0.5) + # Regression line with full range
  facet_wrap(~ condition) +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)",
    title = "Comparison of Z-scores calculated from regression slope between Replicates") +
  # Add R² labels
  geom_text(data = regression_r2_data, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +
  theme(
    axis.title = element_text(size = 16),   # Increase axis title font size
    axis.text = element_text(size = 14),    # Increase axis text font size
    strip.text = element_text(size = 15),   # Increase facet label font size
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size
    plot.subtitle = element_text(size = 16) # Increase subtitle size
  )
ggsave("mini_pheno_regression_correlation.pdf", plot =regression_correlation_p, units = 'mm', width = 297, height = 210)


mo_data_regression <- regression_df[regression_df$condition == "MO", ]

regression_violin_p <- ggplot(mo_data_regression, aes(x = strain, y = z_score)) +
  geom_violin(fill = "lightblue") +
  geom_jitter(data = regression_df, aes(x = strain, y = z_score, 
                                        color = condition), 
              width = 0.15, size = 2) +
  labs(x = "Strain", y = "Z-Score") +
  scale_color_manual(values = custom_color_palette)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("mini_pheno_regression_vioin.pdf", plot =regression_violin_p, units = 'mm', width = 297, height = 210)




for (i in 1:length(strain_list)) {
  each_strain <- strain_list[i]
  each_strain_abundance <- abundance_df %>%
    filter(strain == each_strain)
  each_strain_regression <- regression_df %>%
    filter(strain == each_strain)
  
  
  # Create a dataframe for fitted lines
  fitted_lines <- each_strain_regression %>%
    group_by(strain, replicate, condition) %>%
    do({
      data.frame(
        generation = seq(min(each_strain_abundance$generation), max(each_strain_abundance$generation), by = 1),
        abundance = .$slope * seq(min(each_strain_abundance$generation), max(each_strain_abundance$generation), by = 1) + .$intercept
      )
    }) %>%
    ungroup()
  
  
  # Plot the data
  regression_growth_curve <- ggplot(each_strain_abundance, aes(x = generation, y = abundance, group = replicate, color = condition, shape= replicate)) +
    geom_point() +  # Scatter plot for abundance points
    geom_line(data = fitted_lines, aes(x = generation, y = abundance, group = replicate)) +  # Overlay regression lines
    scale_shape_manual(values = custom_shapes) +
    scale_color_manual(values = custom_color_palette) + 
    labs(
      x = "Generation",
      y = "Abundance",
      title = paste("Growth Curve of", each_strain, "with Linear Regression Lines"),
      color = "Condition"
    ) +
    theme(
#      axis.title = element_text(size = 14),
#      axis.text = element_text(size = 12),
#      legend.title = element_text(size = 12),
#      legend.text = element_text(size = 10),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  dir.create("mini_pheno_regression_growth_curves", showWarnings = FALSE)
  plot_name <- paste("mini_pheno_regression_growth_curves/", each_strain, "_regression_growth_curve.pdf", sep = "")
  
  ggsave(plot_name, 
         plot =regression_growth_curve, units = 'mm', 
         width = 297, height = 210)
  
}

library(writexl)

write_xlsx(
  list(
    "AUC" = auc_df,
    "average_AUC" = auc_merged_data,
    "Regression" = regression_df,
    "average_regression" = regression_merged_data
  ),
  path = "miniPheno_output.xlsx"
)

write.csv(auc_df, "minipheno_auc.csv", row.names = F)
write.csv(regression_df, "minipheno_regression.csv", row.names = F)
