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

# 2. Define a function to process a single TSV and return combined AF + OR plots
#    (plus the per-SNP tables behind them, for supplementary data)
process_shared_ctrl_bsa <- function(strain_name, file_path) {
  
  # A. Read the TSV
  df <- read_tsv(file_path, col_names = c(
    "chrom", "pos",
    "full1", "full2",
    "low1", "low2",
    "uni1", "uni2"
  ), show_col_types = FALSE)
  
  # B0. Keep only biallelic sites: every sample must have exactly "REF,ALT" counts.
  #     Multiallelic sites (e.g. "12,8,3") would otherwise lose the extra
  #     allele's reads and give a biased N2 allele frequency.
  n_before <- nrow(df)
  df <- df %>%
    filter(if_all(full1:uni2, ~ grepl("^[0-9]+,[0-9]+$", .x)))
  message(sprintf("[%s] Removed %d non-biallelic sites (%.1f%% of %d)",
                  strain_name, n_before - nrow(df),
                  100 * (n_before - nrow(df)) / n_before, n_before))
  
  # B. Separate the comma-delimited counts and remove MtDNA
  df_sep <- df %>%
    separate(full1, c("f1_r", "f1_a"), sep=",", convert=TRUE) %>%
    separate(full2, c("f2_r", "f2_a"), sep=",", convert=TRUE) %>%
    separate(low1,  c("l1_r", "l1_a"), sep=",", convert=TRUE) %>%
    separate(low2,  c("l2_r", "l2_a"), sep=",", convert=TRUE) %>%
    separate(uni1,  c("u1_r", "u1_a"), sep=",", convert=TRUE) %>%
    separate(uni2,  c("u2_r", "u2_a"), sep=",", convert=TRUE) %>%
    filter(chrom != "MtDNA")
  
  # C. Apply depth filter (>= 10 for ALL 6 samples) and calculate AF
  df_filt <- df_sep %>%
    filter(
      (f1_r + f1_a) >= 10, (f2_r + f2_a) >= 10,
      (l1_r + l1_a) >= 10, (l2_r + l2_a) >= 10,
      (u1_r + u1_a) >= 10, (u2_r + u2_a) >= 10
    ) %>%
    mutate(
      f1_AF = f1_r / (f1_r + f1_a),
      f2_AF = f2_r / (f2_r + f2_a),
      l1_AF = l1_r / (l1_r + l1_a),
      l2_AF = l2_r / (l2_r + l2_a),
      u1_AF = u1_r / (u1_r + u1_a),
      u2_AF = u2_r / (u2_r + u2_a),
      pos_Mb = pos / 1e6,
      chrom = factor(chrom, levels=c("I","II","III","IV","V","X"))
    )
  
  # D. Format Allele Frequency data for plotting
  df_af <- df_filt %>%
    dplyr::select(chrom, pos_Mb, ends_with("_AF")) %>%
    pivot_longer(cols = ends_with("_AF"), names_to = "sample", values_to = "AF")
  
  # E. Calculate Odds Ratios against the shared controls
  clamp <- function(p) { ifelse(p == 1, 1 - 1e-6, ifelse(p == 0, 1e-6, p)) }
  
  df_or_wide <- df_filt %>%
    mutate(
      o_f1 = clamp(f1_AF) / (1 - clamp(f1_AF)),
      o_f2 = clamp(f2_AF) / (1 - clamp(f2_AF)),
      o_l1 = clamp(l1_AF) / (1 - clamp(l1_AF)),
      o_l2 = clamp(l2_AF) / (1 - clamp(l2_AF)),
      o_u1 = clamp(u1_AF) / (1 - clamp(u1_AF)),
      o_u2 = clamp(u2_AF) / (1 - clamp(u2_AF)),
      
      log10_OR_full_r1 = log10(o_f1 / o_u1),
      log10_OR_full_r2 = log10(o_f2 / o_u2),
      log10_OR_low_r1  = log10(o_l1 / o_u1),
      log10_OR_low_r2  = log10(o_l2 / o_u2)
    )
  
  df_or <- df_or_wide %>%
    dplyr::select(chrom, pos_Mb, starts_with("log10_OR_")) %>%
    pivot_longer(cols = starts_with("log10_OR_"), names_to = "group", values_to = "log10_OR")
  
  # F. Plot Allele Frequency (6 lines mapping full, low, and uni)
  p_af <- ggplot(df_af, aes(x = pos_Mb, y = AF, color = sample)) +
    geom_smooth(se = FALSE, linewidth = 0.8) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(
      values = c(
        "f1_AF" = "darkgreen",  "f2_AF" = "green3",
        "l1_AF" = "orange",     "l2_AF" = "gold",
        "u1_AF" = "dodgerblue", "u2_AF" = "blue"
      ),
      labels = c(
        "f1_AF" = "100% Pathogen (Rep 1)",
        "f2_AF" = "100% Pathogen (Rep 2)",
        "l1_AF" = "33% Pathogen (Rep 1)",
        "l2_AF" = "33% Pathogen (Rep 2)",
        "u1_AF" = "Uninfected (Rep 1)",
        "u2_AF" = "Uninfected (Rep 2)"
      )
    ) +
    labs(
      title = paste0(strain_name, ": LOESS Allele Frequency"),
      x = "Genomic Position (Mb)",
      y = "N2 Allele Frequency",
      color = "Sample Group"
    ) +
    supp_theme +
    guides(color = guide_legend(ncol = 2))
  
  # G. Plot Odds Ratio (4 lines mapping full and low against shared uni)
  p_or <- ggplot(df_or, aes(x = pos_Mb, y = log10_OR, color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_smooth(se = FALSE, linewidth = 0.8) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(
      values = c(
        "log10_OR_full_r1" = "darkgreen", "log10_OR_full_r2" = "green3",
        "log10_OR_low_r1"  = "orange",    "log10_OR_low_r2"  = "gold"
      ),
      labels = c(
        "log10_OR_full_r1" = "100% vs Ctrl (Rep 1)",
        "log10_OR_full_r2" = "100% vs Ctrl (Rep 2)",
        "log10_OR_low_r1"  = "33% vs Ctrl (Rep 1)",
        "log10_OR_low_r2"  = "33% vs Ctrl (Rep 2)"
      )
    ) +
    labs(
      title = paste0(strain_name, ": LOESS Log10 Odds Ratio"),
      x = "Genomic Position (Mb)",
      y = "log10 Odds Ratio (Infected/Control)",
      color = "Sample Group"
    ) +
    supp_theme +
    guides(color = guide_legend(ncol = 2))
  
  # H. Supplementary data (one row per SNP that passed the depth filter)
  # AF panel: read counts + N2 (REF) allele frequency for all 6 samples
  af_out <- df_filt %>%
    transmute(
      chrom = as.character(chrom), pos,
      path100_rep1_N2_reads  = f1_r, path100_rep1_other_reads  = f1_a, path100_rep1_N2_AF  = f1_AF,
      path100_rep2_N2_reads  = f2_r, path100_rep2_other_reads  = f2_a, path100_rep2_N2_AF  = f2_AF,
      path33_rep1_N2_reads   = l1_r, path33_rep1_other_reads   = l1_a, path33_rep1_N2_AF   = l1_AF,
      path33_rep2_N2_reads   = l2_r, path33_rep2_other_reads   = l2_a, path33_rep2_N2_AF   = l2_AF,
      uninfected_rep1_N2_reads = u1_r, uninfected_rep1_other_reads = u1_a, uninfected_rep1_N2_AF = u1_AF,
      uninfected_rep2_N2_reads = u2_r, uninfected_rep2_other_reads = u2_a, uninfected_rep2_N2_AF = u2_AF
    )
  
  # OR panel: log10 odds ratio of each infected sample vs its matching uninfected replicate
  or_out <- df_or_wide %>%
    transmute(
      chrom = as.character(chrom), pos,
      log10_OR_path100_vs_uninfected_rep1 = log10_OR_full_r1,
      log10_OR_path100_vs_uninfected_rep2 = log10_OR_full_r2,
      log10_OR_path33_vs_uninfected_rep1  = log10_OR_low_r1,
      log10_OR_path33_vs_uninfected_rep2  = log10_OR_low_r2
    )
  
  return(list(af = p_af, or = p_or, af_data = af_out, or_data = or_out))
}

# 3. Process the data for both strains directly from TSVs
plots_ju2610 <- process_shared_ctrl_bsa("JU2610", "BSA_results/JU2610/JU2610_allele_depths.tsv")
plots_ed3005 <- process_shared_ctrl_bsa("ED3005", "BSA_results/ED3005/ED3005_allele_depths.tsv")

# 4. Combine into a 2x2 grid with automatic A, B, C, D tagging and individual legends
supp_fig_21 <- (plots_ju2610$af | plots_ju2610$or) /
  (plots_ed3005$af | plots_ed3005$or) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  ) &
  theme(legend.position = "bottom")

# 5. Export Supplementary Figure 21
ggsave("supp_fig21.pdf",
       plot = supp_fig_21,
       width = 180,
       height = 220,
       units = "mm",
       device = cairo_pdf)


bsa_cols <- c("chrom", "pos",
              "path100_rep1_REF_ALT", "path100_rep2_REF_ALT",
              "path33_rep1_REF_ALT",  "path33_rep2_REF_ALT",
              "uninfected_rep1_REF_ALT", "uninfected_rep2_REF_ALT")

raw_s21_ju2610 <- readr::read_tsv("BSA_results/JU2610/JU2610_allele_depths.tsv", col_names = bsa_cols, show_col_types = FALSE)
raw_s21_ed3005 <- readr::read_tsv("BSA_results/ED3005/ED3005_allele_depths.tsv", col_names = bsa_cols, show_col_types = FALSE)

# Must be below 1,048,576 (Excel's row limit, including the header)
nrow(raw_s21_ju2610); nrow(raw_s21_ed3005)
