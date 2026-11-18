#world map
library("ggplot2")
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(purrr)
library(dunn.test)
library(sf)
library(ggpubr)
library(ggrepel)  # For nice labels

# Shapiro-Wilk that returns NA instead of erroring (n < 3, n > 5000, or all values identical)
safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

all_locations <- read.csv("20231213_c_elegans_strain_data.csv")
all_divergence <- read.csv("wild_isolate_variants.csv")
world_map <- map_data("world")

#sf::sf_use_s2(FALSE)  # Temporarily disable s2 geometry engine
# Get world map with both country and continent information
world <- ne_countries(scale = "large", returnclass = "sf")

# Convert your data to spatial points (keeping your NA handling)
all_locations_sf <- st_as_sf(all_locations,
                             coords = c("longitude", "latitude"),
                             crs = 4326,
                             na.fail = FALSE)

# Perform spatial join to get BOTH country and continent information
sf::sf_use_s2(FALSE)
result <- st_join(all_locations_sf, world[c("name", "continent", "iso_a3")])

# Add new columns to original data frame
all_locations$country <- result$name        # Country name
all_locations$country_code <- result$iso_a3 # 3-letter country code
all_locations$continent <- result$continent # Continent

# Optional: Reorder columns to group the new geographic info together
col_order <- c(setdiff(names(all_locations), c("country", "country_code", "continent")),
               "country", "country_code", "continent")
all_locations <- all_locations[, col_order]

all_locations$location <- all_locations$continent
all_locations$location[all_locations$country == "New Zealand"] <- "New Zealand"
all_locations$location[all_locations$country == "Australia"] <- "Australia"

# Check if any cell in a row contains "Hawaii" (case-insensitive)
hawaii_rows <- apply(all_locations, 1, function(row) {
  any(grepl("hawaii", row, ignore.case = TRUE))
})

# Update the location column to "Hawaii" where the condition is TRUE
all_locations$location[hawaii_rows] <- "Hawaii"

# Our strains are named by CaeNDR isotype, so take each isotype's coordinates from its
# reference strain. Strains with no CaeNDR coordinates (VC20019, CB4852, PB306) get NA
# and are dropped below by filter(!is.na(location)).
isotype_locations <- all_locations %>%
  filter(as.logical(isotype_ref_strain)) %>%
  dplyr::select(strain = isotype, latitude, longitude, location)

all_strains_pool <- c()
for (i in 1:nrow(auc_avg_zscore_df)) {
  strain <- as.data.frame(auc_avg_zscore_df)[i, 1]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
auc_avg_zscore_df$pool <- all_strains_pool

all_strains_pool <- c()
for (i in 1:nrow(regression_avg_zscore_df)) {
  strain <- as.data.frame(regression_avg_zscore_df)[i, 1]
  if (strain %in% final_freezing_plan$Strain) {
    pool_assigned <- final_freezing_plan$pool[final_freezing_plan$Strain == strain]
  } else {
    pool_assigned <- "did_not_freeze"
  }
  all_strains_pool <- append(all_strains_pool, pool_assigned)
}
regression_avg_zscore_df$pool <- all_strains_pool

# Start with the right columns, so filtering works even if no Dunn's test runs
empty_dunn <- data.frame(comparison = character(), Z = numeric(), p = numeric(),
                         p_adjusted = numeric(), pathogen = character(),
                         method = character())
dunn_df_sig_total <- empty_dunn   # significant comparisons only (for brackets)
dunn_df_total     <- empty_dunn   # ALL comparisons (for the supplementary table)
total_location_df <- data.frame()
kw_total <- data.frame()

# Collected across the loop for the supplementary tables (divergence figures)
div_stats_total <- data.frame()
div_norm_total  <- data.frame()
div_data_total  <- data.frame()

# Initialize the figure number counter starting at 10
fig_num <- 10

for (i in 2:nrow(conditions)) {
  each_condition <- conditions[i, 1]
  print(each_condition)
  
  auc_avg_zscore_df_condition <- auc_avg_zscore_df %>%
    filter(condition == each_condition)
  
  regression_avg_zscore_df_condition <- regression_avg_zscore_df %>%
    filter(condition == each_condition)
  
  print("Strains in our set but not in CENDR geography locations")
  print(setdiff(auc_avg_zscore_df_condition$strain, isotype_locations$strain))
  print(setdiff(regression_avg_zscore_df_condition$strain, isotype_locations$strain))
  
  print("Strains in our set but not in divergence paper:")
  print(setdiff(auc_avg_zscore_df_condition$strain, all_divergence$Isotype.reference.strain))
  print(setdiff(regression_avg_zscore_df_condition$strain, all_divergence$Isotype.reference.strain))
  
  auc_avg_zscore_with_location <- auc_avg_zscore_df_condition %>%
    left_join(isotype_locations,
              by = "strain") %>%
    filter(!is.na(average) & !is.na(location))
  auc_avg_zscore_with_location$method <- rep("AUC", nrow(auc_avg_zscore_with_location))
  
  regression_avg_zscore_with_location <- regression_avg_zscore_df_condition %>%
    left_join(isotype_locations,
              by = "strain") %>%
    filter(!is.na(average) & !is.na(location))
  regression_avg_zscore_with_location$method <- rep("regression", nrow(regression_avg_zscore_with_location))
  
  total_location_df <- rbind(total_location_df, auc_avg_zscore_with_location)
  total_location_df <- rbind(total_location_df, regression_avg_zscore_with_location)
  
  auc_kw_result <- kruskal.test(average ~ location, data = auc_avg_zscore_with_location)
  
  if (auc_kw_result$p.value < 0.05) {
    message("Kruskal-Wallis p < 0.05: Running Dunn's test...")
    # altp = TRUE: two-sided p-values, so the p < 0.05 filter below is correct
    dunn_result <- dunn.test(auc_avg_zscore_with_location$average, auc_avg_zscore_with_location$location,
                             method = "bh", altp = TRUE)
    print(dunn_result)
    
    dunn_summary <- data.frame(
      comparison = dunn_result$comparisons,
      Z = dunn_result$Z,
      p = dunn_result$altP,
      p_adjusted = dunn_result$altP.adjusted
    )
    
    dunn_summary$pathogen <- rep(each_condition, nrow(dunn_summary))
    dunn_summary$method <- rep("AUC", nrow(dunn_summary))
    
    dunn_df_total <- rbind(dunn_df_total, dunn_summary)   # all comparisons
    
    # View only significant comparisons
    dunn_summary_sig <- dunn_summary %>%
      filter(p_adjusted < 0.05) %>%
      arrange(p_adjusted)
    
    dunn_df_sig_total <- rbind(dunn_df_sig_total, dunn_summary_sig)
    
  } else {
    message("AUC Kruskal-Wallis p >= 0.05: No significant differences between groups. Dunn's test skipped.")
  }
  
  regression_kw_result <- kruskal.test(average ~ location, data = regression_avg_zscore_with_location)
  
  if (regression_kw_result$p.value < 0.05) {
    message("Kruskal-Wallis p < 0.05: Running Dunn's test...")
    # altp = TRUE: two-sided p-values, so the p < 0.05 filter below is correct
    dunn_result <- dunn.test(regression_avg_zscore_with_location$average, regression_avg_zscore_with_location$location,
                             method = "bh", altp = TRUE)
    print(dunn_result)
    
    dunn_summary <- data.frame(
      comparison = dunn_result$comparisons,
      Z = dunn_result$Z,
      p = dunn_result$altP,
      p_adjusted = dunn_result$altP.adjusted
    )
    
    dunn_summary$pathogen <- rep(each_condition, nrow(dunn_summary))
    dunn_summary$method <- rep("regression", nrow(dunn_summary))
    
    dunn_df_total <- rbind(dunn_df_total, dunn_summary)   # all comparisons
    
    # View only significant comparisons
    dunn_summary_sig <- dunn_summary %>%
      filter(p_adjusted < 0.05) %>%
      arrange(p_adjusted)
    
    dunn_df_sig_total <- rbind(dunn_df_sig_total, dunn_summary_sig)
    
  } else {
    message("Regression Kruskal-Wallis p >= 0.05: No significant differences between groups. Dunn's test skipped.")
  }
  
  kw_log <- data.frame(
    condition = each_condition,
    AUC_p = auc_kw_result$p.value,
    regression_p = regression_kw_result$p.value
  )
  
  kw_total <- rbind(kw_total, kw_log)
  
  AUC_merged_data <- auc_avg_zscore_with_location %>%
    left_join(all_divergence, by = c("strain" = "Isotype.reference.strain"))
  regression_merged_data <- regression_avg_zscore_with_location %>%
    left_join(all_divergence, by = c("strain" = "Isotype.reference.strain"))
  #CB4952 not in worldmap so not here
  
  print("Are 2 data frames identical?")
  print(identical(nrow(AUC_merged_data), nrow(regression_merged_data)))
  
  AUC_merged_data$method <- rep("AUC", nrow(auc_avg_zscore_with_location))
  regression_merged_data$method <- rep("regression", nrow(regression_avg_zscore_with_location))
  
  joint_coordinates_divergence <- rbind(AUC_merged_data, regression_merged_data)
  
  joint_coordinates_divergence$method <- factor(
    joint_coordinates_divergence$method,
    levels = c("AUC", "regression")
  )
  
  # Plot with strain locations
  map_p <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(map_id = region),
             fill = "white", color = "grey", size = 0.3) +
    geom_point(data = joint_coordinates_divergence,
               aes(x = longitude, y = latitude, fill = average, color = average),
               shape = 21, size = 1., stroke = 0.4) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-10, 10),
      oob = scales::squish
    ) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-10, 10),
      oob = scales::squish
    ) +
    coord_cartesian(ylim = c(-50, 60)) +
    facet_wrap(~method) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(title = "Geographic distributions of C. elegans wild isolates",
         subtitle = each_condition,
         fill = "Z-Score") +
    guides(fill = guide_colorbar(title = "Value (z-score)"), color = "none")
  
  ggsave(paste0("phenomip6_", each_condition, "_", "worldmap.pdf"), plot = map_p, units = 'mm', width = 500, height = 210)
  
  #correlation of z-scores with divergence
  long_data <- joint_coordinates_divergence %>%
    pivot_longer(
      cols = c(Hyper.divergent.genome..kb.,
               Fraction.hyper.divergent,
               Genome.wide.variants,
               Divergent.region.variants,
               Fraction.divergent.region.variants),
      names_to = "divergence_metric",
      values_to = "divergence_value"
    )
  
  long_data$divergence_metric <- factor(
    long_data$divergence_metric,
    levels = c(
      "Genome.wide.variants",
      "Hyper.divergent.genome..kb.",
      "Fraction.hyper.divergent",
      "Divergent.region.variants",
      "Fraction.divergent.region.variants"
    ),
    labels = c(
      "Genome-wide\nVariants",
      "Hyper-divergent\nGenome (kb)",
      "Fraction\nHyper-divergent",
      "Divergent Region\nVariants",
      "Fraction Divergent\nRegion Variants"
    )
  )
  
  long_data$method <- factor(
    long_data$method,
    levels = c("AUC", "regression")
  )
  
  r2_labels <- long_data %>%
    group_by(method, divergence_metric) %>%
    summarise(
      model = list(lm(average ~ divergence_value, data = cur_data())),
      .groups = "drop"
    ) %>%
    mutate(r2 = map_dbl(model, ~summary(.x)$r.squared)) %>%
    mutate(label = paste0("R² = ", round(r2, 2)))
  
  # ----------------------------------------------------------------------------
  # Supplementary tables for this condition (tagged with this figure's number)
  # ----------------------------------------------------------------------------
  div_data <- long_data %>%
    dplyr::filter(!is.na(divergence_value), !is.na(average))
  
  # Stats table: one row per facet (method x divergence metric)
  div_stats <- div_data %>%
    dplyr::group_by(method, divergence_metric) %>%
    dplyr::summarize(
      n            = dplyr::n(),
      r2           = summary(lm(average ~ divergence_value))$r.squared,  # = plot label
      pearson_r    = cor(divergence_value, average),
      pearson_p    = cor.test(divergence_value, average)$p.value,
      spearman_rho = cor(divergence_value, average, method = "spearman"),
      spearman_p   = suppressWarnings(
        cor.test(divergence_value, average, method = "spearman")$p.value),
      .groups = "drop"
    ) %>%
    dplyr::mutate(supp_fig = fig_num, condition = each_condition, .before = 1)
  
  # Normality of the regression residuals, per facet
  div_norm <- div_data %>%
    dplyr::group_by(method, divergence_metric) %>%
    dplyr::summarize(
      n         = dplyr::n(),
      shapiro_W = unname(safe_shapiro(resid(lm(average ~ divergence_value)))["W"]),
      shapiro_p = unname(safe_shapiro(resid(lm(average ~ divergence_value)))["p"]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(supp_fig = fig_num, condition = each_condition, .before = 1)
  
  # Plotted data for the data workbook
  div_plot_data <- div_data %>%
    dplyr::select(strain, pool, location, method, average,
                  divergence_metric, divergence_value) %>%
    dplyr::mutate(supp_fig = fig_num, condition = each_condition, .before = 1)
  
  div_stats_total <- rbind(div_stats_total, div_stats)
  div_norm_total  <- rbind(div_norm_total,  div_norm)
  div_data_total  <- rbind(div_data_total,  div_plot_data)
  
  divergence_correlation_p <- ggplot(long_data, aes(x = divergence_value, y = average, color = pool)) +
    geom_point(alpha = 0.8, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
    facet_grid(method ~ divergence_metric, scales = "free_x") +
    scale_color_manual(values = c("pool1" = "purple", "pool2" = "green", "pool3" = "orange", "control" = "pink")) +
    geom_text(
      data = r2_labels,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.5,
      inherit.aes = FALSE,
      size = 8 / .pt,
      color = "black"
    ) +
    labs(
      title = "Correlation of Divergence Metrics with Average Z-score",
      subtitle = each_condition,
      x = "Divergence Metric Value",
      y = "Value (z-score)",
      color = "Pool"
    ) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      panel.grid.minor = element_blank()
    )
  
  # Updated saving line using the figure number counter
  ggsave(paste0("supp_fig", fig_num, ".pdf"), plot = divergence_correlation_p, units = 'mm',
         width = 180, height = 220/2)
  
  # Increment the counter for the next condition (keep this AFTER the tables above)
  fig_num <- fig_num + 1
}

# ==============================================================================
# One data frame per divergence figure: stats_suppNN, normality_suppNN, data_suppNN
# ==============================================================================
# Line breaks in the facet labels look messy in Excel, so replace them with spaces
div_stats_total$divergence_metric <- gsub("\n", " ", div_stats_total$divergence_metric)
div_norm_total$divergence_metric  <- gsub("\n", " ", div_norm_total$divergence_metric)
div_data_total$divergence_metric  <- gsub("\n", " ", div_data_total$divergence_metric)

for (f in unique(div_stats_total$supp_fig)) {
  assign(paste0("stats_supp", f),     dplyr::filter(div_stats_total, supp_fig == f))
  assign(paste0("normality_supp", f), dplyr::filter(div_norm_total,  supp_fig == f))
  assign(paste0("data_supp", f),      dplyr::filter(div_data_total,  supp_fig == f))
}

# Which condition went into which figure
print(dplyr::distinct(div_stats_total, supp_fig, condition))

# Kruskal-Wallis results for every condition (Supp Fig 9)
kw_total

#The Kruskal-Wallis test indicated a significant difference in average z-scores among locations (p = 0.04).
#However, no individual pairwise comparisons reached statistical significance after adjustment for multiple testing,
#suggesting that differences may be subtle or distributed across multiple locations rather than driven by a few specific pairs.

total_location_df <- na.omit(total_location_df)

total_location_df <- total_location_df %>%
  filter(! condition == "MO")

dunn_df_sig_total

# 1. Prepare the significance data
sig_comparisons_regression <- dunn_df_sig_total %>%
  filter(pathogen == "PV", method == "regression") %>%
  separate(comparison, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    p.adj = p_adjusted,
    p.adj.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    condition = "PV",
    method = "regression"
  )

sig_comparisons <- sig_comparisons_regression

# 2. Calculate optimal y-positions
y_max <- max(total_location_df$average, na.rm = TRUE)
step_size <- 0.1 * y_max
sig_comparisons <- sig_comparisons %>%
  arrange(method, desc(p.adj)) %>%
  group_by(method) %>%  # separate brackets per facet row
  mutate(y.position = y_max * 1.1 + (row_number() * step_size)) %>%
  ungroup()

labelled_points <- total_location_df %>%
  filter(average > 25 | average < -25)  # Adjust threshold as needed

total_location_df$condition <- factor(
  total_location_df$condition,
  levels = c("M114", "M117", "M395", "PA", "PV")
)

total_location_df$method <- factor(
  total_location_df$method,
  levels = c("AUC", "regression")
)

sig_comparisons$method <- factor(
  sig_comparisons$method,
  levels = c("AUC", "regression")
)

# Only draw brackets if at least one comparison is significant;
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

# 3. Create the plot
location_p <- ggplot(total_location_df, aes(x = location, y = average)) +
  geom_boxplot(
    alpha = 0.7, width = 0.6, outlier.shape = NA, color = "black"
  ) +
  geom_jitter(
    aes(color = pool),
    width = 0.15, height = 0,
    alpha = 0.5, size = 1.5
  ) +
  geom_text_repel(
    data = labelled_points,
    aes(label = strain),
    size = 3,
    max.overlaps = 20
  ) +
  # Significance brackets (NULL when nothing is significant)
  bracket_layer +
  facet_grid(method ~ condition) +
  scale_color_manual(values = c(
    "pool1" = "purple",
    "pool2" = "green",
    "pool3" = "brown",
    "control" = "pink",
    "did_not_freeze" = "grey60"
  )) +
  labs(x = "Location", y = "Value (Z-score)", color = "Pool",
       title = "Geographic Variation in wild isolates phenotypes",
       subtitle = "Statistical comparison across wild isolate locations using Kruskal-Wallis and Dunn's post-hoc tests") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")

ggsave("location_comparison.pdf", plot = location_p, units = 'mm', width = 297, height = 210)

#kruskal-wallis is like ANOVA but use rank sum, compares rank sums between every pair of groups
#because >3 groups, non normal data, unequal sample sizes
#dunn test to see which differ

#one world map
world_map <- map_data("world")

map <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "white",
               color = "grey70",
               linewidth = 0.2) +
  geom_point(data = joint_coordinates_divergence,
             aes(x = longitude, y = latitude, color = pool),
             size = 1.5,
             alpha = 0.8) +
  coord_cartesian(ylim = c(-50, 60)) +
  theme_void() +
  theme(legend.position = "bottom")

ggsave("world_map.pdf", plot = map, units = 'mm', width = 297, height = 130)