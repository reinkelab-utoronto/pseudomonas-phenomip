library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Source the analysis script into an isolated environment
env_corr <- new.env()
source("minipheno/mini_pheno_analysis.R", local = env_corr)

# 2. Extract and filter data frames to remove the "MO" condition
auc_data_filtered <- env_corr$auc_merged_data %>% filter(condition != "MO")
reg_data_filtered <- env_corr$regression_merged_data %>% filter(condition != "MO")

# 3. Recalculate R² labels for the filtered data
auc_r2_filtered <- auc_data_filtered %>%
  group_by(condition) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

reg_r2_filtered <- reg_data_filtered %>%
  group_by(condition) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

# 4. Build auc_corr_plot (Panel A) matching main figure formatting
auc_corr_plot <- ggplot(auc_data_filtered, aes(x = z_R1, y = z_R2)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = TRUE, linewidth = 0.5) +
  facet_wrap(~ condition) +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)"  ) +
  geom_text(
    data = auc_r2_filtered,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.5,
    size = 8 / .pt,
    color = "black"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold")
  )

# 5. Build reg_corr_plot (Panel B) matching main figure formatting
reg_corr_plot <- ggplot(reg_data_filtered, aes(x = z_R1, y = z_R2)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = TRUE, linewidth = 0.5) +
  facet_wrap(~ condition) +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)"  ) +
  geom_text(
    data = reg_r2_filtered,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.5,
    size = 8 / .pt,
    color = "black"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold")
  )

# 6. Combine side-by-side using Patchwork with automatic A/B tags
combined_correlations <- auc_corr_plot + reg_corr_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )

# View the result
combined_correlations

# 7. Export the final combined correlation figure
ggsave("supp_fig4.pdf", 
       plot = combined_correlations, 
       width = 180, 
       height = 220/2, 
       units = "mm", 
       device = cairo_pdf)



# ==============================================================================
# Stats for Supp Fig 4 (mini-phenoMIP replicate correlation, per condition, no MO)
# ==============================================================================
supp4_data <- dplyr::bind_rows(
  auc_data_filtered %>% dplyr::mutate(Metric = "AUC"),        # Panel A
  reg_data_filtered %>% dplyr::mutate(Metric = "regression")  # Panel B
)

# Stats table: one row per facet
stats_supp4 <- supp4_data %>%
  dplyr::group_by(Metric, condition) %>%
  dplyr::summarize(
    n            = sum(complete.cases(z_R1, z_R2)),
    r2           = summary(lm(z_R2 ~ z_R1))$r.squared,   # same R² as the plot labels
    pearson_r    = cor(z_R1, z_R2, use = "complete.obs"),
    pearson_p    = cor.test(z_R1, z_R2)$p.value,
    spearman_rho = cor(z_R1, z_R2, method = "spearman", use = "complete.obs"),
    spearman_p   = cor.test(z_R1, z_R2, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

# Normality: z_R1, z_R2 and regression residuals, per facet
normality_supp4 <- supp4_data %>%
  dplyr::group_by(Metric, condition) %>%
  dplyr::mutate(residual = resid(lm(z_R2 ~ z_R1, na.action = na.exclude))) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(c(z_R1, z_R2, residual),
                      names_to = "variable", values_to = "value") %>%
  dplyr::group_by(Metric, condition, variable) %>%
  dplyr::summarize(
    n         = sum(!is.na(value)),
    shapiro_W = unname(safe_shapiro(value)["W"]),
    shapiro_p = unname(safe_shapiro(value)["p"]),
    .groups = "drop"
  )

stats_supp4
