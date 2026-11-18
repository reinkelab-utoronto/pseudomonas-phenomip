#total_correlation_df <- read.csv("total_correlation_df.csv")
#total_auc_df <-  read.csv("total_auc_df.csv")
#total_regression_df<- read.csv("total_regression_df.csv")

controls <- c("JU1400", "CB4856", "VC20019")

total_correlation_df$condition <- factor(
  total_correlation_df$condition,
  levels = c("MO", "M114", "M117", "M395", "PA", "PV")
)
filtered_df <- total_correlation_df[total_correlation_df$strain %in% controls, ]



#AUC
auc_avg_zscore_df <- total_auc_df %>%
  group_by(pool, strain, condition) %>%
  summarise(avg_z_score = mean(z_score, na.rm = TRUE), .groups = "drop")

filtered_auc_avg_zscore_df <- auc_avg_zscore_df[auc_avg_zscore_df$strain %in% controls, ]

wide_data <- filtered_auc_avg_zscore_df %>%
  pivot_wider(names_from = pool, values_from = avg_z_score)

p1p2_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool2 ~ pool1))$r.squared  # Fit linear model and extract R²
  )

p1p2_wide_data_with_r2 <- wide_data %>%
  left_join(p1p2_r_squared_values, by = "strain")


correlation_plot_pool1_vs_pool2 <- ggplot(p1p2_wide_data_with_r2, aes(x = pool1, y = pool2)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 1 Avg Z-Score", 
    y = "Pool 2 Avg Z-Score", 
    title = "Correlation Plot: Pool 1 vs Pool 2",
    subtitle = "Average of z-scores calculated using AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool1_vs_pool2, filename = "controls_auc_zscore_correlation_plot_pool1_vs_pool2.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')

p1p3_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool3 ~ pool1))$r.squared
  )

p1p3_wide_data_with_r2 <- wide_data %>%
  left_join(p1p3_r_squared_values, by = "strain")

correlation_plot_pool1_vs_pool3 <- ggplot(p1p3_wide_data_with_r2, aes(x = pool1, y = pool3)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 1 Avg Z-Score", 
    y = "Pool 3 Avg Z-Score", 
    title = "Correlation Plot: Pool 1 vs Pool 3",
    subtitle = "Average of z-scores calculated using AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool1_vs_pool3, filename = "controls_auc_zscore_correlation_plot_pool1_vs_pool3.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')


p2p3_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool2 ~ pool3))$r.squared
  )

p2p3_wide_data_with_r2 <- wide_data %>%
  left_join(p2p3_r_squared_values, by = "strain")

correlation_plot_pool2_vs_pool3 <- ggplot(p2p3_wide_data_with_r2, aes(x = pool2, y = pool3)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 2 Avg Z-Score", 
    y = "Pool 3 Avg Z-Score", 
    title = "Correlation Plot: Pool 2 vs Pool 3",
    subtitle = "Average of z-scores calculated using AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool2_vs_pool3, filename = "controls_auc_zscore_correlation_plot_pool2_vs_pool3.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')

#regression
controls <- c("JU1400", "CB4856", "VC20019")

regression_avg_zscore_df <- total_regression_df %>%
  group_by(pool, strain, condition) %>%
  summarise(avg_z_score = mean(z_score, na.rm = TRUE), .groups = "drop")

filtered_regression_avg_zscore_df <- regression_avg_zscore_df[regression_avg_zscore_df$strain %in% controls, ]

wide_data <- filtered_regression_avg_zscore_df %>%
  pivot_wider(names_from = pool, values_from = avg_z_score)

p1p2_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool2 ~ pool1))$r.squared  # Fit linear model and extract R²
  )

p1p2_wide_data_with_r2 <- wide_data %>%
  left_join(p1p2_r_squared_values, by = "strain")


correlation_plot_pool1_vs_pool2 <- ggplot(p1p2_wide_data_with_r2, aes(x = pool1, y = pool2)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 1 Avg Z-Score", 
    y = "Pool 2 Avg Z-Score", 
    title = "Correlation Plot: Pool 1 vs Pool 2",
    subtitle = "Average of z-scores calculated using regression slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool1_vs_pool2, filename = "controls_regression_zscore_correlation_plot_pool1_vs_pool2.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')

p1p3_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool3 ~ pool1))$r.squared
  )

p1p3_wide_data_with_r2 <- wide_data %>%
  left_join(p1p3_r_squared_values, by = "strain")

correlation_plot_pool1_vs_pool3 <- ggplot(p1p3_wide_data_with_r2, aes(x = pool1, y = pool3)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 1 Avg Z-Score", 
    y = "Pool 3 Avg Z-Score", 
    title = "Correlation Plot: Pool 1 vs Pool 3",
    subtitle = "Average of z-scores calculated using regression slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool1_vs_pool3, filename = "controls_regression_zscore_correlation_plot_pool1_vs_pool3.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')


p2p3_r_squared_values <- wide_data %>%
  group_by(strain) %>%
  summarise(
    r_squared = summary(lm(pool2 ~ pool3))$r.squared
  )

p2p3_wide_data_with_r2 <- wide_data %>%
  left_join(p2p3_r_squared_values, by = "strain")

correlation_plot_pool2_vs_pool3 <- ggplot(p2p3_wide_data_with_r2, aes(x = pool2, y = pool3)) +
  geom_point(aes(color = condition), size = 3) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +     # Add regression line
  geom_text(
    aes(label = paste("R² =", round(r_squared, 2))),  # Add R² annotation
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 8, color = "blue"
  ) +
  facet_wrap(~ strain, scales = "free") +                     # Facet by strain
  labs(
    x = "Pool 2 Avg Z-Score", 
    y = "Pool 3 Avg Z-Score", 
    title = "Correlation Plot: Pool 2 vs Pool 3",
    subtitle = "Average of z-scores calculated using regression slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(correlation_plot_pool2_vs_pool3, filename = "controls_regression_zscore_correlation_plot_pool2_vs_pool3.pdf", width = 45, height = 21, 
       units = "cm", device='pdf')

####others

control_auc_r2_data_all <- filtered_df %>%
  group_by(strain) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))


control_z_correlation_all <- ggplot(filtered_df, aes(x = z_R1, y = z_R2, color = pool, shape = condition)) +
  geom_point(size = 2) +          
  facet_grid(~strain) +   
  geom_smooth(aes(group = strain), method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) +  # One regression line per facet
  geom_text(data = control_auc_r2_data_all, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +  
  labs(
    x = "Z-score (Replicate 1)",  
    y = "Z-score (Replicate 2)",  
    color = "Condition",          
    title = "Comparison of Z-scores (Replicate 1 and Replicate 2)",
    subtitle = "Control strains: JU1400, CB4856, VC20019"
  ) +
  theme(
    legend.position = "bottom",   
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12),    
    strip.text = element_text(size = 12),   
    plot.title = element_text(size = 16, face = "bold"),  
    plot.subtitle = element_text(size = 14) 
  )
ggsave(paste0("phenomip6_controls_auc_zscore_correlation_all.pdf"), plot =control_z_correlation_all, units = 'mm', width = 297, height = 210)


control_auc_r2_data <- filtered_df %>%
  group_by(strain, pool) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))


control_z_correlation <- ggplot(filtered_df, aes(x = z_R1, y = z_R2, color = pool, shape = condition)) +
  geom_point(size = 2) +          
  facet_grid(pool~strain) +   
  geom_smooth(aes(group = strain), method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) +  # One regression line per facet
  geom_text(data = control_auc_r2_data, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +  
  labs(
    x = "Z-score (Replicate 1)",  
    y = "Z-score (Replicate 2)",  
    color = "Condition",          
    title = "Comparison of Z-scores (Replicate 1 and Replicate 2)",
    subtitle = "Control strains: JU1400, CB4856, VC20019"
  ) +
  theme(
    legend.position = "bottom",   
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12),    
    strip.text = element_text(size = 12),   
    plot.title = element_text(size = 16, face = "bold"),  
    plot.subtitle = element_text(size = 14) 
  )
ggsave(paste0("phenomip6_controls_auc_zscore_correlation.pdf"), plot =control_z_correlation, units = 'mm', width = 297, height = 210)


filtered_df <- total_correlation_df_regression[total_correlation_df_regression$strain %in% controls, ]

control_regression_r2_data_all <- filtered_df %>%
  group_by(strain) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))


control_z_correlation_all_regression <- ggplot(filtered_df, aes(x = z_R1, y = z_R2, color = pool, shape = condition)) +
  geom_point(size = 2) +          
  facet_grid(~strain) +   
  geom_smooth(aes(group = strain), method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) +  # One regression line per facet
  geom_text(data = control_regression_r2_data_all, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +  
  labs(
    x = "Z-score (Replicate 1)",  
    y = "Z-score (Replicate 2)",  
    color = "Condition",          
    title = "Comparison of Z-scores (Replicate 1 and Replicate 2)",
    subtitle = "Control strains: JU1400, CB4856, VC20019"
  ) +
  theme(
    legend.position = "bottom",   
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12),    
    strip.text = element_text(size = 12),   
    plot.title = element_text(size = 16, face = "bold"),  
    plot.subtitle = element_text(size = 14) 
  )
ggsave(paste0("phenomip6_controls_regression_zscore_correlation_all.pdf"), plot =control_z_correlation_all_regression, units = 'mm', width = 297, height = 210)


control_regression_r2_data <- filtered_df %>%
  group_by(strain, pool) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))


control_z_correlation_regression <- ggplot(filtered_df, aes(x = z_R1, y = z_R2, color = pool, shape = condition)) +
  geom_point(size = 2) +          
  facet_grid(pool~strain) +   
  geom_smooth(aes(group = strain), method = "lm", se = FALSE, color = "black", fullrange = TRUE, size = 0.5) +  # One regression line per facet
  geom_text(data = control_regression_r2_data, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 5, color = "red") +  
  labs(
    x = "Z-score (Replicate 1)",  
    y = "Z-score (Replicate 2)",  
    color = "Condition",          
    title = "Comparison of Z-scores (Replicate 1 and Replicate 2)",
    subtitle = "Control strains: JU1400, CB4856, VC20019"
  ) +
  theme(
    legend.position = "bottom",   
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12),    
    strip.text = element_text(size = 12),   
    plot.title = element_text(size = 16, face = "bold"),  
    plot.subtitle = element_text(size = 14) 
  )
ggsave(paste0("phenomip6_controls_regression_zscore_correlation.pdf"), plot =control_z_correlation_regression, units = 'mm', width = 297, height = 210)

