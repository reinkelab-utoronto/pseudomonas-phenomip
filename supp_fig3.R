library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

pool_colors <- c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "control" = "pink")

# --- Panel A ---
total_correlation_df_forplotting <- total_correlation_df
total_correlation_df_forplotting <- total_correlation_df_forplotting %>%
  filter(!condition == "MO")

r2_data <- total_correlation_df_forplotting %>%
  group_by(method, condition) %>%
  summarize(r2 = summary(lm(z_R2 ~ z_R1))$r.squared) %>%
  mutate(label = paste0("R² = ", round(r2, 3)))

total_correlation_p <- ggplot(total_correlation_df_forplotting, aes(x = z_R1, y = z_R2, color = pool)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE, linewidth = 0.5) +
  facet_grid(condition ~ method) +
  scale_color_manual(values = pool_colors, name = "Pool") +
  labs(
    x = "Z-score (Replicate 1)",
    y = "Z-score (Replicate 2)"
  ) +
  geom_text(data = r2_data, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 8 / .pt, color = "red") +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

# --- Panel B ---
filtered_df_noMO <- filtered_df %>%
  filter(!condition == "MO")

filtered_df_noMO$method <- as.character(filtered_df_noMO$method)
filtered_df_noMO$method[filtered_df_noMO$method == "regression slope"] <- "regression"
controls_heatmap <- ggplot(filtered_df_noMO, aes(x = condition, y = pool, fill = metric)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    oob = squish,
    name = "Metric"
  ) +
  facet_grid(strain ~ method, scale = "free") +
  labs(
    x = "Condition",
    y = ""
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )

# --- Patchwork assembly ---
# Panel A takes 2/3 of vertical space, panel B takes 1/3
final_supp_fig3 <- total_correlation_p / controls_heatmap +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )

# View the result
final_supp_fig3

# Export
ggsave("supp_fig3.pdf",
       plot = final_supp_fig3,
       width = 180,
       height = 220,
       units = "mm",
       device = cairo_pdf)


# ==============================================================================
# Stats for Supp Fig 3A (replicate correlation per method x condition, no MO)
# ==============================================================================

# Stats table: one row per facet
stats_supp3a <- total_correlation_df_forplotting %>%
  dplyr::group_by(method, condition) %>%
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
normality_supp3a <- total_correlation_df_forplotting %>%
  dplyr::group_by(method, condition) %>%
  dplyr::mutate(residual = resid(lm(z_R2 ~ z_R1, na.action = na.exclude))) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(c(z_R1, z_R2, residual),
                      names_to = "variable", values_to = "value") %>%
  dplyr::group_by(method, condition, variable) %>%
  dplyr::summarize(
    n         = sum(!is.na(value)),
    shapiro_W = unname(safe_shapiro(value)["W"]),
    shapiro_p = unname(safe_shapiro(value)["p"]),
    .groups = "drop"
  )

stats_supp3a
normality_supp3a