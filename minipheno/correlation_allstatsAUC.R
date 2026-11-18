library(plyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggdendro)
library(reshape2)
library(grid)
library(cluster)
library(tidyr)
library(stringr)
library(pracma)
library(rstatix)
library(ggpubr)
library(broom)

# ==============================================================================
# 0. Shared Plot Styling (matches my_plot_theme used in the growth-curve /
#    regression script, so all figures across both scripts look consistent)
# ==============================================================================
# Canvas: full A4 landscape (297 x 210 mm)
my_plot_theme <- theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box       = "horizontal",
    legend.direction = "horizontal",
    plot.title      = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 20),
    axis.text       = element_text(size = 18),
    strip.text      = element_text(size = 20, face = "bold"),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.2, "cm")
  )

# Legend guide spec (added separately to each plot, since guides() can't be
# combined with a theme object via `+`)
my_legend_guides <- guides(colour = guide_legend(nrow = 2, byrow = TRUE),
                           shape  = guide_legend(nrow = 2, byrow = TRUE))

fig_width  <- 297
fig_height <- 210
fig_units  <- "mm"

pheno6_auc <- read.csv("pheno6_9strains_auc.csv")
minipheno_auc <- read.csv("minipheno_auc.csv")

head(pheno6_auc)
head(minipheno_auc)

# ============================================================
# DATA PREPARATION
# ============================================================

correlation_df <- pheno6_auc %>%
  dplyr::select(strain, condition, replicate, pheno6_z = z_score) %>%
  inner_join(
    minipheno_auc %>%
      dplyr::select(strain, condition, replicate, minipheno_z = z_score),
    by = c("strain", "condition", "replicate")
  ) %>%
  mutate(replicate_group = case_when(
    replicate %in% c("R01", "R02", "R05", "R07", "R09", "R11") ~ "R1",
    replicate %in% c("R03", "R04", "R06", "R08", "R10", "R12") ~ "R2",
    TRUE ~ NA_character_
  ))

# Present in pheno6 but missing in minipheno
anti_join(pheno6_auc, minipheno_auc, by = c("strain", "condition", "replicate"))

# Present in minipheno but missing in pheno6
anti_join(minipheno_auc, pheno6_auc, by = c("strain", "condition", "replicate"))

correlation_df[correlation_df$strain == 'VC20019',]

correlation_df <- correlation_df %>%
  filter(condition != "MO")

# Averaged z-scores
pheno6_avg_z <- pheno6_auc %>%
  group_by(strain, condition) %>%
  summarise(mean_z_score = mean(z_score, na.rm = TRUE), .groups = "drop")

minipheno_avg_z <- minipheno_auc %>%
  group_by(strain, condition) %>%
  summarise(mean_z_score = mean(z_score, na.rm = TRUE), .groups = "drop")

head(pheno6_avg_z)
head(minipheno_avg_z)
identical(pheno6_avg_z$strain, minipheno_avg_z$strain)
identical(pheno6_avg_z$condition, minipheno_avg_z$condition)
identical(pheno6_avg_z$mean_z_score, minipheno_avg_z$mean_z_score)

plot_df <- pheno6_avg_z %>%
  rename(pheno6_z = mean_z_score) %>%
  inner_join(
    minipheno_avg_z %>% rename(mini_z = mean_z_score),
    by = c("strain", "condition")
  ) %>%
  filter(condition != "MO")

# ============================================================
# SUMMARY STATS: R, R^2, and Spearman rho per condition
# (replicate-level)
# ============================================================

stats_df <- correlation_df %>%
  group_by(condition) %>%
  summarise(
    R   = cor(pheno6_z, minipheno_z, use = "complete.obs", method = "pearson"),
    R2  = R^2,
    rho = cor(pheno6_z, minipheno_z, use = "complete.obs", method = "spearman"),
    # Combined multi-line label, e.g.:
    # R = 0.850
    # R² = 0.723
    # rho = 0.812
    label = sprintf("R = %.3f\nR\u00b2 = %.3f\nrho = %.3f", R, R2, rho),
    .groups = "drop"
  )

# Global Pearson R, R^2, and Spearman rho (replicate-level)
r_global   <- cor(correlation_df$pheno6_z, correlation_df$minipheno_z,
                  use = "complete.obs", method = "pearson")
r2_global  <- r_global^2
rho_global <- cor(correlation_df$pheno6_z, correlation_df$minipheno_z,
                  use = "complete.obs", method = "spearman")

global_label_reps <- sprintf("R = %.3f\nR\u00b2 = %.3f\nrho = %.3f",
                             r_global, r2_global, rho_global)

# ============================================================
# SUMMARY STATS: R, R^2, and Spearman rho per condition
# (averaged)
# ============================================================

stats_df_avg <- plot_df %>%
  group_by(condition) %>%
  summarise(
    R   = cor(pheno6_z, mini_z, use = "complete.obs", method = "pearson"),
    R2  = R^2,
    rho = cor(pheno6_z, mini_z, use = "complete.obs", method = "spearman"),
    label = sprintf("R = %.3f\nR\u00b2 = %.3f\nrho = %.3f", R, R2, rho),
    .groups = "drop"
  )

# Global Pearson R, R^2, and Spearman rho (averaged)
r_global_avg   <- cor(plot_df$pheno6_z, plot_df$mini_z,
                      use = "complete.obs", method = "pearson")
r2_global_avg  <- r_global_avg^2
rho_global_avg <- cor(plot_df$pheno6_z, plot_df$mini_z,
                      use = "complete.obs", method = "spearman")

global_label_avg <- sprintf("R = %.3f\nR\u00b2 = %.3f\nrho = %.3f",
                            r_global_avg, r2_global_avg, rho_global_avg)

# ============================================================
# PLOTS: Faceted by condition (replicate-level)
# ============================================================

p_faceted_reps <- ggplot(correlation_df,
                         aes(x = pheno6_z,
                             y = minipheno_z)) +
  geom_point(aes(colour = strain, shape = replicate_group), size = 3.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", fullrange = TRUE) +
  geom_text(data = stats_df,
            aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 5.5) +
  facet_wrap(~condition) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_shape_manual(values = c("R1" = 16, "R2" = 17)) +
  labs(x = "Pheno6 z-score", y = "MiniPheno z-score") +
  my_plot_theme +
  my_legend_guides

ggsave("plot_AllStats_faceted_replicatesauc.pdf", plot = p_faceted_reps,
       units = fig_units, width = fig_width, height = fig_height)

# ============================================================
# PLOTS: Faceted by condition (averaged)
# ============================================================

p_faceted_avg <- ggplot(plot_df, aes(x = pheno6_z, y = mini_z)) +
  geom_point(aes(colour = strain), size = 3.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", fullrange = TRUE) +
  geom_text(data = stats_df_avg,
            aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 5) +
  facet_wrap(~condition) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(x = "Pheno6 mean z-score", y = "Minipheno mean z-score") +
  my_plot_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

ggsave("plot_AllStats_faceted_averagedauc.pdf", plot = p_faceted_avg,
       units = fig_units, width = fig_width, height = fig_height)

# ============================================================
# PLOTS: Global (replicate-level)
# ============================================================

p_global_reps <- ggplot(correlation_df,
                        aes(x = pheno6_z,
                            y = minipheno_z)) +
  geom_point(aes(colour = strain, shape = replicate_group), size = 3.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", fullrange = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = global_label_reps,
           hjust = -0.1, vjust = 1.1, size = 6.5) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_shape_manual(values = c("R1" = 16, "R2" = 17)) +
  labs(x = "Pheno6 z-score", y = "MiniPheno z-score") +
  my_plot_theme +
  my_legend_guides

ggsave("plot_AllStats_global_replicatesauc.pdf", plot = p_global_reps,
       units = fig_units, width = fig_width, height = fig_height)

# ============================================================
# PLOTS: Global (averaged)
# ============================================================

p_global_avg <- ggplot(plot_df, aes(x = pheno6_z, y = mini_z)) +
  geom_point(aes(colour = strain), size = 3.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", fullrange = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = global_label_avg,
           hjust = -0.1, vjust = 1.1, size = 6.5) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(x = "Pheno6 mean z-score", y = "Minipheno mean z-score") +
  my_plot_theme

ggsave("plot_AllStats_global_averagedauc.pdf", plot = p_global_avg,
       units = fig_units, width = fig_width, height = fig_height)