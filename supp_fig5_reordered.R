library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Source both scripts into isolated environments
env_auc <- new.env()
env_reg <- new.env()

source("minipheno/correlation_allstatsAUC.R", local = env_auc)
source("minipheno/correlation_allstatsregression.R", local = env_reg)

# 2. Extract plot data and stats data from both environments
df_auc_avg <- env_auc$plot_df
stats_auc_avg <- env_auc$stats_df_avg

df_reg_avg <- env_reg$plot_df
stats_reg_avg <- env_reg$stats_df_avg

# 3. Build Panel A (AUC averaged correlation)
p_auc_supp5 <- ggplot(df_auc_avg, aes(x = pheno6_z, y = mini_z)) +
  geom_point(aes(colour = strain), size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = TRUE, linewidth = 0.5) +
  facet_wrap(~ condition) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    x = "PhenoMIP mean z-score", 
    y = "Mini-phenoMIP mean z-score",
    colour = "Strain"
  ) +
  geom_text(
    data = stats_auc_avg,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.1,
    size = 8 / .pt,
    color = "black"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  # Force the title to sit perfectly centered ABOVE the legend items
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))

# 4. Build Panel B (Regression averaged correlation)
p_reg_supp5 <- ggplot(df_reg_avg, aes(x = pheno6_z, y = mini_z)) +
  geom_point(aes(colour = strain), size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = TRUE, linewidth = 0.5) +
  facet_wrap(~ condition) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    x = "PhenoMIP mean z-score", 
    y = "Mini-phenoMIP mean z-score",
    colour = "Strain"
  ) +
  geom_text(
    data = stats_reg_avg,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.1,
    size = 8 / .pt,
    color = "black"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  # Force the title to sit perfectly centered ABOVE the legend items
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))

# 5. Combine side-by-side using Patchwork 
supp_fig_5_combined <- p_auc_supp5 + p_reg_supp5 +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = 'A') & # The '&' ensures the theme below applies to the collected patchwork wrapper
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = "topleft"
  )

# View the result
supp_fig_5_combined

# 6. Export Supplementary Figure 5
ggsave("supp_fig5.pdf", 
       plot = supp_fig_5_combined, 
       width = 180, 
       height = 220/2, 
       units = "mm", 
       device = cairo_pdf)


# ==============================================================================
# Stats for Supp Fig 5 (phenoMIP vs mini-phenoMIP mean z-scores, per condition)
# ==============================================================================
supp5_data <- dplyr::bind_rows(
  df_auc_avg %>% dplyr::mutate(Metric = "AUC"),        # Panel A
  df_reg_avg %>% dplyr::mutate(Metric = "regression")  # Panel B
)

# n and p-values per facet (computed from the same plot_df the figure uses)
supp5_pvals <- supp5_data %>%
  dplyr::group_by(Metric, condition) %>%
  dplyr::summarize(
    n          = sum(complete.cases(pheno6_z, mini_z)),
    pearson_p  = cor.test(pheno6_z, mini_z)$p.value,
    spearman_p = suppressWarnings(cor.test(pheno6_z, mini_z, method = "spearman")$p.value),
    .groups = "drop"
  )

# Stats table: R, R², rho from the sourced scripts (= plot labels) + n and p-values
stats_supp5 <- dplyr::bind_rows(
  stats_auc_avg %>% dplyr::mutate(Metric = "AUC"),
  stats_reg_avg %>% dplyr::mutate(Metric = "regression")
) %>%
  dplyr::select(Metric, condition, pearson_r = R, r2 = R2, spearman_rho = rho) %>%
  dplyr::left_join(supp5_pvals, by = c("Metric", "condition")) %>%
  dplyr::relocate(n, .after = condition)

# Normality: phenoMIP, mini-phenoMIP and regression residuals, per facet
normality_supp5 <- supp5_data %>%
  dplyr::group_by(Metric, condition) %>%
  dplyr::mutate(residual = resid(lm(mini_z ~ pheno6_z, na.action = na.exclude))) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(c(pheno6_z, mini_z, residual),
                      names_to = "variable", values_to = "value") %>%
  dplyr::group_by(Metric, condition, variable) %>%
  dplyr::summarize(
    n         = sum(!is.na(value)),
    shapiro_W = unname(safe_shapiro(value)["W"]),
    shapiro_p = unname(safe_shapiro(value)["p"]),
    .groups = "drop"
  )

stats_supp5
