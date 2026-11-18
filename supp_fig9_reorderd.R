library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)    # geom_bracket()
library(ggrepel)   # geom_text_repel()
library(FSA)       # dunnTest()

# Needs total_location_df from the location/divergence loop script

# Shapiro-Wilk that returns NA instead of erroring (n < 3, n > 5000, or all values identical)
safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

# Define consistent colors
pool_colors <- c(
  "pool1" = "purple",
  "pool2" = "green",
  "pool3" = "orange",
  "control" = "pink",
  "did_not_freeze" = "grey60"
)

# ============================================================
# 1. Clean Data (keep AUC and regression, drop MO)
# ============================================================
total_location_df_clean <- total_location_df %>%
  filter(! condition == "MO") %>%
  filter(method %in% c("AUC", "regression")) %>%
  mutate(
    condition = factor(condition, levels = c("M114", "M117", "M395", "PA", "PV")),
    method    = factor(method, levels = c("AUC", "regression"))
  ) %>%
  droplevels()

# ============================================================
# 2. Stats (these tables also drive the brackets on the plot)
# ============================================================

# 2a. Descriptives per location (what the boxplots show)
desc_s9 <- total_location_df_clean %>%
  dplyr::group_by(method, condition, location) %>%
  dplyr::summarise(
    n      = sum(!is.na(average)),
    mean   = mean(average, na.rm = TRUE),
    sd     = sd(average, na.rm = TRUE),
    median = median(average, na.rm = TRUE),
    q1     = quantile(average, 0.25, na.rm = TRUE),
    q3     = quantile(average, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# 2b. Kruskal-Wallis per facet
kw_s9 <- total_location_df_clean %>%
  dplyr::group_by(method, condition) %>%
  dplyr::summarise(
    n       = sum(!is.na(average)),
    kw_chi2 = unname(kruskal.test(average, factor(location))$statistic),
    kw_df   = unname(kruskal.test(average, factor(location))$parameter),
    kw_p    = kruskal.test(average, factor(location))$p.value,
    .groups = "drop"
  )

# 2c. Dunn's test (two-sided, Benjamini-Hochberg), only where Kruskal-Wallis p < 0.05
sig_facets <- kw_s9 %>%
  dplyr::filter(kw_p < 0.05) %>%
  dplyr::select(method, condition)

if (nrow(sig_facets) > 0) {
  dunn_s9 <- total_location_df_clean %>%
    dplyr::semi_join(sig_facets, by = c("method", "condition")) %>%
    dplyr::group_by(method, condition) %>%
    dplyr::group_modify(~ {
      dt <- FSA::dunnTest(average ~ factor(location), data = .x, method = "bh")$res
      tibble::tibble(comparison = dt$Comparison,
                     Z          = dt$Z,
                     p_unadj    = dt$P.unadj,
                     p_adjusted = dt$P.adj)
    }) %>%
    dplyr::ungroup()
} else {
  dunn_s9 <- tibble::tibble(method = factor(), condition = factor(),
                            comparison = character(), Z = numeric(),
                            p_unadj = numeric(), p_adjusted = numeric())
}

dunn_s9 <- dunn_s9 %>%
  dplyr::mutate(signif = dplyr::case_when(
    p_adjusted < 0.001 ~ "***",
    p_adjusted < 0.01  ~ "**",
    p_adjusted < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# 2d. Normality per location group (context for using a non-parametric test)
normality_s9 <- total_location_df_clean %>%
  dplyr::group_by(method, condition, location) %>%
  dplyr::summarise(
    n         = sum(!is.na(average)),
    shapiro_W = unname(safe_shapiro(average)["W"]),
    shapiro_p = unname(safe_shapiro(average)["p"]),
    .groups = "drop"
  )

kw_s9
print(dunn_s9, n = Inf)

# ============================================================
# 3. Significance Brackets (every significant Dunn comparison, any facet)
# ============================================================
y_max     <- max(total_location_df_clean$average, na.rm = TRUE)
step_size <- 0.1 * y_max

sig_comparisons <- dunn_s9 %>%
  dplyr::filter(p_adjusted < 0.05) %>%
  separate(comparison, into = c("group1", "group2"), sep = " - ") %>%
  dplyr::group_by(method, condition) %>%
  dplyr::arrange(desc(p_adjusted), .by_group = TRUE) %>%
  dplyr::mutate(y.position = y_max * 1.1 + row_number() * step_size) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p.adj.signif = signif)

# geom_bracket() errors on an empty data frame, and ggplot ignores NULL layers
bracket_layer <- if (nrow(sig_comparisons) > 0) {
  geom_bracket(
    data = sig_comparisons,
    aes(xmin = group1, xmax = group2, label = p.adj.signif, y.position = y.position),
    tip.length = 0.02,
    vjust = -0.2,
    inherit.aes = FALSE
  )
} else {
  NULL
}

# ============================================================
# 4. Identify Outliers for Labels
# ============================================================
labelled_points <- total_location_df_clean %>%
  filter(average > 25 | average < -25)

# ============================================================
# 5. Create the Plot
# ============================================================
location_p <- ggplot(total_location_df_clean, aes(x = location, y = average)) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    aes(color = pool),
    width = 0.15,
    height = 0,
    alpha = 0.5,
    size = 1.5
  ) +
  # Strain labels
  geom_text_repel(
    data = labelled_points,
    aes(label = strain),
    size = 6 / .pt,
    max.overlaps = 20
  ) +
  # Significance brackets (none drawn if no comparison is significant)
  bracket_layer +
  facet_grid(method ~ condition) +
  scale_color_manual(
    values = pool_colors,
    name = "Pool"
  ) +
  labs(
    x = "Location",
    y = "Value (Z-score)",
    color = "Pool"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  ) +
  coord_cartesian(clip = "off")

# View result
location_p

# ============================================================
# 6. Export
# ============================================================
ggsave(
  "supp_fig9.pdf",
  plot = location_p,
  width = 180,
  height = 220 * (2/3),
  units = "mm",
  device = cairo_pdf
)