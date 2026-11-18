library(tidyverse)
library(patchwork)
library(DescTools)

# ==============================================================================
# 0. Shared Plot Styling (same as Figure 6)
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


# ==============================================================================
# 1. Helpers
# ==============================================================================

# Read one CSV (one column per strain) into long format.
# x-axis order: N2 (the control), then the two parents (the cross's own parent
# first), then the RILs in the order they appear in the file
PARENTS <- c("JU2610", "JU2522")

read_pheno <- function(file, parent) {
  wide   <- readr::read_csv(file, show_col_types = FALSE)   # handles the BOM
  others <- setdiff(names(wide), c("N2", PARENTS))
  order  <- c("N2", parent, setdiff(PARENTS, parent), others)
  order  <- order[order %in% names(wide)]
  wide %>%
    pivot_longer(everything(), names_to = "strain", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(strain = factor(strain, levels = order))
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

# x-axis labels: RILs shown as "RIL JU4", "RIL J21", as in Figures 3 and 4
strain_labels <- function(x) ifelse(x %in% c("N2", PARENTS), x, paste("RIL", x))

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
#   y_breaks = tick positions on the y-axis (default: chosen automatically)
fish_panel <- function(d, y_label, title = NULL, top = panel_top(d), n_brackets = NULL,
                       y_breaks = waiver()) {
  sig <- dunnett_vs_N2(d)
  cat("\n==", title, "-", y_label, "\n"); print(sig)   # all p-values for the legend
  
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
                       breaks = y_breaks,
                       expand = expansion(mult = c(0, 0))) +
    labs(title = title, x = NULL, y = y_label) +
    my_plot_theme +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# 2. Load Data and Build Panels
# ==============================================================================
JU2610_embryo <- read_pheno("JU2610_RILs_specificity_embryo.csv",  "JU2610")   # <-- EDIT path if needed
JU2610_gravid <- read_pheno("JU2610_RILs_specificity_gravid.csv",  "JU2610")
JU2522_embryo <- read_pheno("JU2522_RILs_specificity_embryos.csv", "JU2522")
JU2522_gravid <- read_pheno("JU2522_RILs_specificity_gravid.csv",  "JU2522")

GRAVID_BREAKS <- seq(0, 100, 25)   # ticks stop at 100 (space above holds brackets)

s26a <- fish_panel(JU2610_embryo, "Average normalised embryo", title = "JU2610 RILs on PA14")
s26b <- fish_panel(JU2610_gravid, "% gravid",                  title = "JU2610 RILs on PA14",
                   y_breaks = GRAVID_BREAKS)
s26c <- fish_panel(JU2522_embryo, "Average normalised embryo", title = "JU2522 RILs on MYb114")
s26d <- fish_panel(JU2522_gravid, "% gravid",                  title = "JU2522 RILs on MYb114",
                   y_breaks = GRAVID_BREAKS)

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================
supp_fig26 <- s26a + s26b + s26c + s26d +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(size = 10, face = "plain"))
  ) &
  theme(plot.tag = element_text(size = 10, face = "plain"))

supp_fig26   # view in RStudio

ggsave("supp_fig26.pdf",
       plot   = supp_fig26,
       width  = 180,
       height = 220/2,
       units  = "mm",
       device = cairo_pdf)

# ==============================================================================
# 4. Stats tables + raw data for supplementary data
# ==============================================================================
# Panel labels in the same style as Figures 3E, 4E, 5D-E, Supp Figs 22 and 26
# (cross + bacterium in brackets, since A/C and B/D share a measure)
panels_s26 <- list(
  "A: Average normalised embryo (JU2610 RILs on PA14)"   = JU2610_embryo,
  "B: % gravid (JU2610 RILs on PA14)"                    = JU2610_gravid,
  "C: Average normalised embryo (JU2522 RILs on MYb114)" = JU2522_embryo,
  "D: % gravid (JU2522 RILs on MYb114)"                  = JU2522_gravid
)

# Shapiro-Wilk that returns NA instead of erroring (n < 3, or all values identical)
safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

# Raw data: every replicate, as plotted (panel, strain, value)
pheno_raw <- function(d, panel) {
  d %>%
    arrange(strain) %>%                                  # plot order (N2, parents, RILs)
    transmute(panel  = panel,
              strain = strain_labels(as.character(strain)),   # as on the x-axis
              value)
}

# Stats table: descriptives (match the plotted mean +/- SD) + Dunnett vs N2
pheno_stats <- function(d, panel) {
  desc <- d %>%
    group_by(strain) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") %>%
    mutate(strain = as.character(strain))
  
  dun <- DunnettTest(value ~ strain, data = d, control = "N2")[[1]]
  dun_df <- tibble(
    strain       = sub("-N2$", "", rownames(dun)),
    diff_vs_N2   = dun[, "diff"],
    lwr_ci_vs_N2 = dun[, "lwr.ci"],
    upr_ci_vs_N2 = dun[, "upr.ci"],
    p_adj_vs_N2  = dun[, "pval"]
  ) %>%
    mutate(signif_vs_N2 = stars_from_p(p_adj_vs_N2))
  
  desc %>%
    left_join(dun_df, by = "strain") %>%          # N2 row gets NA comparison columns
    mutate(strain = strain_labels(strain),
           panel  = panel, .before = 1)
}

# Normality: Shapiro-Wilk per strain + on the pooled ANOVA residuals
pheno_normality <- function(d, panel) {
  per_strain <- d %>%
    group_by(strain) %>%
    summarise(n         = n(),
              shapiro_W = unname(safe_shapiro(value)["W"]),
              shapiro_p = unname(safe_shapiro(value)["p"]),
              .groups = "drop") %>%
    mutate(strain = strain_labels(as.character(strain)))
  
  res <- safe_shapiro(residuals(aov(value ~ strain, data = d)))
  pooled <- tibble(strain = "All strains (ANOVA residuals)", n = nrow(d),
                   shapiro_W = unname(res["W"]), shapiro_p = unname(res["p"]))
  
  bind_rows(per_strain, pooled) %>% mutate(panel = panel, .before = 1)
}

# Levene's test (median-centred) for equal variances across strains
pheno_levene <- function(d, panel) {
  lev <- DescTools::LeveneTest(value ~ strain, data = d)
  tibble(panel    = panel,
         F_value  = lev[1, "F value"],
         df1      = lev[1, "Df"],
         df2      = lev[2, "Df"],
         levene_p = lev[1, "Pr(>F)"])
}

raw_s26       <- bind_rows(imap(panels_s26, pheno_raw))
stats_s26     <- bind_rows(imap(panels_s26, pheno_stats))
normality_s26 <- bind_rows(imap(panels_s26, pheno_normality))
levene_s26    <- bind_rows(imap(panels_s26, pheno_levene))

raw_s26
stats_s26
normality_s26
levene_s26

# ==============================================================================
# Tables for the supplementary workbooks
# ==============================================================================
#   Data:  raw_s26
#   Stats: stats_s26, normality_s26, levene_s26