library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# 1. Define a shared, compact theme matching your previous figures
supp_theme <- theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0) # Pull legend up slightly
  )

# 2. Define a function to process a single TSV and return combined G8+G6 plots
process_combined_bsa <- function(strain_name, file_path) {
  
  # A. Read the TSV
  df <- read_tsv(file_path, col_names = c(
    "chrom", "pos", 
    "G6_PA1", "G6_PA2", "G6_uni1", "G6_uni2", 
    "G8_PA1", "G8_PA2", "G8_uni1", "G8_uni2"
  ), show_col_types = FALSE)
  
  # B. Separate the comma-delimited counts and remove MtDNA
  df_sep <- df %>%
    separate(G6_PA1,  c("g6_m1_r", "g6_m1_a"), sep=",", convert=TRUE) %>%
    separate(G6_PA2,  c("g6_m2_r", "g6_m2_a"), sep=",", convert=TRUE) %>%
    separate(G6_uni1, c("g6_u1_r", "g6_u1_a"), sep=",", convert=TRUE) %>%
    separate(G6_uni2, c("g6_u2_r", "g6_u2_a"), sep=",", convert=TRUE) %>%
    separate(G8_PA1,  c("g8_m1_r", "g8_m1_a"), sep=",", convert=TRUE) %>%
    separate(G8_PA2,  c("g8_m2_r", "g8_m2_a"), sep=",", convert=TRUE) %>%
    separate(G8_uni1, c("g8_u1_r", "g8_u1_a"), sep=",", convert=TRUE) %>%
    separate(G8_uni2, c("g8_u2_r", "g8_u2_a"), sep=",", convert=TRUE) %>%
    filter(chrom != "MtDNA")
  
  # C. Apply depth filter (>= 10 for ALL samples) and calculate AF
  df_filt <- df_sep %>%
    filter(
      (g6_m1_r + g6_m1_a) >= 10, (g6_m2_r + g6_m2_a) >= 10,
      (g6_u1_r + g6_u1_a) >= 10, (g6_u2_r + g6_u2_a) >= 10,
      (g8_m1_r + g8_m1_a) >= 10, (g8_m2_r + g8_m2_a) >= 10,
      (g8_u1_r + g8_u1_a) >= 10, (g8_u2_r + g8_u2_a) >= 10
    ) %>%
    mutate(
      g6_m1_AF = g6_m1_r / (g6_m1_r + g6_m1_a),
      g6_m2_AF = g6_m2_r / (g6_m2_r + g6_m2_a),
      g6_u1_AF = g6_u1_r / (g6_u1_r + g6_u1_a),
      g6_u2_AF = g6_u2_r / (g6_u2_r + g6_u2_a),
      g8_m1_AF = g8_m1_r / (g8_m1_r + g8_m1_a),
      g8_m2_AF = g8_m2_r / (g8_m2_r + g8_m2_a),
      g8_u1_AF = g8_u1_r / (g8_u1_r + g8_u1_a),
      g8_u2_AF = g8_u2_r / (g8_u2_r + g8_u2_a),
      pos_Mb = pos / 1e6,
      chrom = factor(chrom, levels=c("I","II","III","IV","V","X"))
    )
  
  # D. Format Allele Frequency data for plotting
  df_af <- df_filt %>%
    dplyr::select(chrom, pos_Mb, ends_with("_AF")) %>%
    pivot_longer(cols = ends_with("_AF"), names_to = "sample", values_to = "AF")
  
  # E. Calculate Odds Ratios
  clamp <- function(p) { ifelse(p == 1, 1 - 1e-6, ifelse(p == 0, 1e-6, p)) }
  
  df_or <- df_filt %>%
    mutate(
      o_g6_m1 = clamp(g6_m1_AF) / (1 - clamp(g6_m1_AF)),
      o_g6_m2 = clamp(g6_m2_AF) / (1 - clamp(g6_m2_AF)),
      o_g6_u1 = clamp(g6_u1_AF) / (1 - clamp(g6_u1_AF)),
      o_g6_u2 = clamp(g6_u2_AF) / (1 - clamp(g6_u2_AF)),
      o_g8_m1 = clamp(g8_m1_AF) / (1 - clamp(g8_m1_AF)),
      o_g8_m2 = clamp(g8_m2_AF) / (1 - clamp(g8_m2_AF)),
      o_g8_u1 = clamp(g8_u1_AF) / (1 - clamp(g8_u1_AF)),
      o_g8_u2 = clamp(g8_u2_AF) / (1 - clamp(g8_u2_AF)),
      
      log10_OR_g6_r1 = log10(o_g6_m1 / o_g6_u1),
      log10_OR_g6_r2 = log10(o_g6_m2 / o_g6_u2),
      log10_OR_g8_r1 = log10(o_g8_m1 / o_g8_u1),
      log10_OR_g8_r2 = log10(o_g8_m2 / o_g8_u2)
    ) %>%
    dplyr::select(chrom, pos_Mb, starts_with("log10_OR_")) %>%
    pivot_longer(cols = starts_with("log10_OR_"), names_to = "group", values_to = "log10_OR")
  
  # F. Plot Allele Frequency (8 lines mapped to Supp Fig 18 colors)
  p_af <- ggplot(df_af, aes(x = pos_Mb, y = AF, color = sample)) +
    geom_smooth(se = FALSE, linewidth = 0.8) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(
      values = c(
        "g8_m1_AF" = "darkgreen",  "g8_m2_AF" = "green3",    # 50% Pathogen mapping to high
        "g6_m1_AF" = "orange",     "g6_m2_AF" = "gold",      # 10% Pathogen mapping to low
        "g8_u1_AF" = "darkblue",   "g8_u2_AF" = "blue",      # 50% Ctrl
        "g6_u1_AF" = "purple4",    "g6_u2_AF" = "purple"     # 10% Ctrl
      ),
      labels = c(
        "g8_m1_AF" = "50% Pathogen (Rep 1)",
        "g8_m2_AF" = "50% Pathogen (Rep 2)",
        "g6_m1_AF" = "10% Pathogen (Rep 1)",
        "g6_m2_AF" = "10% Pathogen (Rep 2)",
        "g8_u1_AF" = "Uninfected 50% Ctrl (Rep 1)",
        "g8_u2_AF" = "Uninfected 50% Ctrl (Rep 2)",
        "g6_u1_AF" = "Uninfected 10% Ctrl (Rep 1)",
        "g6_u2_AF" = "Uninfected 10% Ctrl (Rep 2)"
      )
    ) +
    labs(
      title = paste0(strain_name, ": LOESS Allele Frequency"), 
      x = "Genomic Position (Mb)", 
      y = "N2 Allele Frequency", 
      color = "Sample Group"
    ) +
    supp_theme +
    guides(color = guide_legend(ncol = 2)) # Forces exactly 2 columns to prevent width cutoff
  
  # G. Plot Odds Ratio (4 lines mapped to Supp Fig 18 colors)
  p_or <- ggplot(df_or, aes(x = pos_Mb, y = log10_OR, color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_smooth(se = FALSE, linewidth = 0.8) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(
      values = c(
        "log10_OR_g8_r1" = "darkgreen", "log10_OR_g8_r2" = "green3",
        "log10_OR_g6_r1" = "orange",    "log10_OR_g6_r2" = "gold"
      ),
      labels = c(
        "log10_OR_g8_r1" = "50% vs Ctrl (Rep 1)",
        "log10_OR_g8_r2" = "50% vs Ctrl (Rep 2)",
        "log10_OR_g6_r1" = "10% vs Ctrl (Rep 1)",
        "log10_OR_g6_r2" = "10% vs Ctrl (Rep 2)"
      )
    ) +
    labs(
      title = paste0(strain_name, ": LOESS Log10 Odds Ratio"), 
      x = "Genomic Position (Mb)", 
      y = "log10 Odds Ratio (Infected/Control)", 
      color = "Sample Group"
    ) +
    supp_theme +
    guides(color = guide_legend(ncol = 2)) # Forces exactly 2 columns to prevent width cutoff
  
  return(list(af = p_af, or = p_or))
}

# 3. Process the data for both strains directly from TSVs
plots_ju2522 <- process_combined_bsa("JU2522", "BSA_results/JU2522/JU2522_allele_depths.tsv")
plots_wn2050 <- process_combined_bsa("WN2050", "BSA_results/WN2050/WN2050_allele_depths.tsv")

# 4. Combine into a 2x2 grid with automatic A, B, C, D tagging and individual legends
supp_fig_24 <- (plots_ju2522$af | plots_ju2522$or) /
  (plots_wn2050$af | plots_wn2050$or) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )


# 5. Export Supplementary Figure 19
ggsave("supp_fig24.pdf", 
       plot = supp_fig_24, 
       width = 180, 
       height = 220, 
       units = "mm", 
       device = cairo_pdf)

bsa24_cols <- c("chrom", "pos",
                "G6_10pct_pathogen_rep1_REF_ALT", "G6_10pct_pathogen_rep2_REF_ALT",
                "G6_uninfected_rep1_REF_ALT",     "G6_uninfected_rep2_REF_ALT",
                "G8_50pct_pathogen_rep1_REF_ALT", "G8_50pct_pathogen_rep2_REF_ALT",
                "G8_uninfected_rep1_REF_ALT",     "G8_uninfected_rep2_REF_ALT")

raw_s24_ju2522 <- readr::read_tsv("BSA_results/JU2522/JU2522_allele_depths.tsv", col_names = bsa24_cols, show_col_types = FALSE)
raw_s24_wn2050 <- readr::read_tsv("BSA_results/WN2050/WN2050_allele_depths.tsv", col_names = bsa24_cols, show_col_types = FALSE)

# Must be below 1,048,576 (Excel's row limit)
nrow(raw_s24_ju2522); nrow(raw_s24_wn2050)