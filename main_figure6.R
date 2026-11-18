library(tidyverse)
library(patchwork)
library(DescTools)

# ==============================================================================
# 0. Shared Plot Styling (same as Figure 5)
# ==============================================================================
my_plot_theme <- theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 8, face = "plain", hjust = 0),
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 6),
    strip.text      = element_text(size = 8, face = "bold"),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    plot.margin     = margin(5, 5, 5, 5)
  )

Y_LABEL <- "Average bacterial burden (% fluorescence)"

# ==============================================================================
# 1. Helpers
# ==============================================================================

# Read one CSV (one column per strain) into long format, keeping column order
read_fish <- function(file) {
  wide <- readr::read_csv(file, show_col_types = FALSE)   # handles the BOM
  wide %>%
    pivot_longer(everything(), names_to = "strain", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(strain = factor(strain, levels = names(wide)))
}

stars_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Dunnett's test of every strain vs N2
dunnett_vs_N2 <- function(d) {
  test <- DunnettTest(value ~ strain, data = d, control = "N2")[[1]]
  tibble(strain = sub("-N2$", "", rownames(test)),
         p_adj  = test[, "pval"]) %>%
    mutate(label = stars_from_p(p_adj))
}

# x-axis labels: italicise the gene name in pmk-1(km25)
strain_labels <- function(x) {
  lab <- ifelse(x == "pmk-1(km25)",
                'italic("pmk-1")*"(km25)"',
                paste0('"', x, '"'))
  parse(text = lab)
}

# Number of strains significantly different from N2 (adjusted p < 0.05)
n_significant <- function(d) sum(dunnett_vs_N2(d)$p_adj < 0.05)

# Highest point or SD bar in one dataset
panel_top <- function(d) {
  summ <- d %>%
    group_by(strain) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  max(c(d$value, summ$mean + replace_na(summ$sd, 0)))
}

# One panel: replicate points, mean +/- SD, Dunnett brackets vs N2
# top = where the brackets start; pass the same value to every panel
# to give them a shared y-axis
fish_panel <- function(d, title, top = panel_top(d), n_brackets = NULL) {
  sig <- dunnett_vs_N2(d)
  print(sig)   # all p-values in the console for the figure legend
  
  # Only draw brackets for significant comparisons
  sig <- sig %>% filter(p_adj < 0.05)
  
  summ <- d %>%
    group_by(strain) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  
  # Stack brackets above the top, nearest strain lowest
  if (is.null(n_brackets)) n_brackets <- nrow(sig)
  step <- top * 0.10
  sig  <- sig %>%
    mutate(strain = factor(strain, levels = levels(d$strain))) %>%
    arrange(strain) %>%
    mutate(y = top + step * row_number())
  tip <- step * 0.2   # bracket tip length, same in every panel
  
  ggplot(d, aes(x = strain, y = value)) +
    geom_jitter(width = 0.15, height = 0, size = 0.8, alpha = 0.6, color = "gray50") +
    geom_errorbar(data = summ, inherit.aes = FALSE,
                  aes(x = strain, ymin = pmax(mean - sd, 0), ymax = mean + sd),
                  width = 0, color = "gray20", linewidth = 0.35) +
    geom_point(data = summ, inherit.aes = FALSE,
               aes(x = strain, y = mean), size = 1.5, color = "black") +
    # Brackets: horizontal bar from N2 to each strain, with short end tips
    geom_segment(data = sig, inherit.aes = FALSE,
                 aes(x = 1, xend = as.numeric(strain), y = y, yend = y),
                 linewidth = 0.3) +
    geom_segment(data = sig, inherit.aes = FALSE,
                 aes(x = 1, xend = 1, y = y, yend = y - tip),
                 linewidth = 0.3) +
    geom_segment(data = sig, inherit.aes = FALSE,
                 aes(x = as.numeric(strain), xend = as.numeric(strain),
                     y = y, yend = y - tip),
                 linewidth = 0.3) +
    geom_text(data = sig, inherit.aes = FALSE,
              aes(x = (1 + as.numeric(strain)) / 2, y = y, label = label),
              vjust = 0.2, size = 2.3) +
    scale_x_discrete(labels = strain_labels) +
    scale_y_continuous(limits = c(0, top + step * (n_brackets + 1.5)),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = title, x = NULL, y = Y_LABEL) +
    my_plot_theme +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# 2. Load Data and Build Panels
# ==============================================================================
PA14   <- read_fish("PA14_FISH_area.csv")   # <-- EDIT path if needed
MYb114 <- read_fish("M114_FISH_area.csv")   # <-- EDIT path if needed
MYb395 <- read_fish("M395_FISH_area.csv")   # <-- EDIT path if needed

# Shared y-axis for the MYb114 and MYb395 panels (B and C); PA14 in panel A
# keeps its own scale. Room is left for the larger number of significant brackets.
shared_top   <- max(panel_top(MYb114), panel_top(MYb395))
max_brackets <- max(n_significant(MYb114), n_significant(MYb395))

# Panel titles: full names, species in italics
p6a <- fish_panel(PA14,   expression(italic("P. aeruginosa")~"PA14"))
p6b <- fish_panel(MYb114, expression(italic("P. poae")~"MYb114"),
                  top = shared_top, n_brackets = max_brackets)
p6c <- fish_panel(MYb395, expression(italic("P. canadensis")~"MYb395"),
                  top = shared_top, n_brackets = max_brackets)

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================
final_figure_6 <- p6a + p6b + p6c +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = "topleft"
  )

final_figure_6   # view in RStudio

ggsave("Figure_6_Final.pdf",
       plot   = final_figure_6,
       width  = 180,
       height = 220/3,
       units  = "mm",
       device = cairo_pdf)

# ==============================================================================
# 4. Stats tables for supplementary data
# ==============================================================================
panels_6 <- list(PA14 = PA14, MYb114 = MYb114, MYb395 = MYb395)

# Shapiro-Wilk that returns NA instead of erroring (n < 3, or all values identical)
safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

# Stats table: descriptives (match the plotted mean +/- SD) + Dunnett vs N2
fish_stats <- function(d, panel) {
  desc <- d %>%
    group_by(strain) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") %>%
    mutate(strain = as.character(strain))
  
  dun <- DunnettTest(value ~ strain, data = d, control = "N2")[[1]]
  dun_df <- tibble(
    strain     = sub("-N2$", "", rownames(dun)),
    diff_vs_N2 = dun[, "diff"],
    lwr_ci     = dun[, "lwr.ci"],
    upr_ci     = dun[, "upr.ci"],
    p_adj      = dun[, "pval"]
  ) %>%
    mutate(signif = stars_from_p(p_adj))
  
  desc %>%
    left_join(dun_df, by = "strain") %>%   # N2 row gets NA comparison columns
    mutate(panel = panel, .before = 1)
}

# Normality: Shapiro-Wilk per strain + on the pooled ANOVA residuals
fish_normality <- function(d, panel) {
  per_strain <- d %>%
    group_by(strain) %>%
    summarise(n         = n(),
              shapiro_W = unname(safe_shapiro(value)["W"]),
              shapiro_p = unname(safe_shapiro(value)["p"]),
              .groups = "drop") %>%
    mutate(strain = as.character(strain))
  
  res <- safe_shapiro(residuals(aov(value ~ strain, data = d)))
  pooled <- tibble(strain = "All strains (ANOVA residuals)", n = nrow(d),
                   shapiro_W = unname(res["W"]), shapiro_p = unname(res["p"]))
  
  bind_rows(per_strain, pooled) %>% mutate(panel = panel, .before = 1)
}

stats_fig6     <- bind_rows(imap(panels_6, fish_stats))
normality_fig6 <- bind_rows(imap(panels_6, fish_normality))

stats_fig6
normality_fig6

# ==============================================================================
# 5. Raw data for supplementary data (same layout as fig5de$raw)
# ==============================================================================
panels_6_raw <- list("A: PA14 FISH"   = PA14,
                     "B: MYb114 FISH" = MYb114,
                     "C: MYb395 FISH" = MYb395)

fig6_raw <- bind_rows(imap(panels_6_raw, ~ .x %>%
                             arrange(strain) %>%
                             mutate(strain = as.character(strain), panel = .y, .before = 1)))

fig6_raw