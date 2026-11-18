library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# ============================================================
# 1. Shared theme
# ============================================================

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
    legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)
  )


# ============================================================
# 2. Function to process one BSA dataset
# ============================================================

process_bsa <- function(strain_name, file_path) {
  
  # ----------------------------------------------------------
  # A. Read TSV
  # ----------------------------------------------------------
  
  # A. Read TSV
  df <- read_tsv(
    file_path,
    col_names = c("chrom", "pos", "m1", "m2", "uni1", "uni2"),
    show_col_types = FALSE
  )
  
  
  # ----------------------------------------------------------
  # B. Split allele-depth fields
  #
  # AD can occasionally contain more than two values if a site
  # has multiple ALT alleles. We retain only REF and first ALT.
  # ----------------------------------------------------------
  
  split_ad <- function(x) {
    
    pieces <- strsplit(as.character(x), ",", fixed = TRUE)
    
    ref <- sapply(pieces, function(z) {
      if (length(z) >= 1) z[1] else NA
    })
    
    alt <- sapply(pieces, function(z) {
      if (length(z) >= 2) z[2] else NA
    })
    
    data.frame(
      ref = as.numeric(ref),
      alt = as.numeric(alt)
    )
  }
  
  
  m1_ad <- split_ad(df$m1)
  m2_ad <- split_ad(df$m2)
  u1_ad <- split_ad(df$uni1)
  u2_ad <- split_ad(df$uni2)
  
  
  df_sep <- df %>%
    mutate(
      m1_r = m1_ad$ref,
      m1_a = m1_ad$alt,
      m2_r = m2_ad$ref,
      m2_a = m2_ad$alt,
      u1_r = u1_ad$ref,
      u1_a = u1_ad$alt,
      u2_r = u2_ad$ref,
      u2_a = u2_ad$alt
    ) %>%
    filter(chrom != "MtDNA")
  
  
  # ----------------------------------------------------------
  # C. Depth filter
  #
  # Require >=10 reads in ALL FOUR samples
  # ----------------------------------------------------------
  
  df_filt <- df_sep %>%
    filter(
      (m1_r + m1_a) >= 10,
      (m2_r + m2_a) >= 10,
      (u1_r + u1_a) >= 10,
      (u2_r + u2_a) >= 10
    ) %>%
    mutate(
      m1_AF = m1_r / (m1_r + m1_a),
      m2_AF = m2_r / (m2_r + m2_a),
      u1_AF = u1_r / (u1_r + u1_a),
      u2_AF = u2_r / (u2_r + u2_a),
      
      pos_Mb = pos / 1e6,
      
      chrom = factor(
        chrom,
        levels = c("I", "II", "III", "IV", "V", "X")
      )
    )
  
  
  # ----------------------------------------------------------
  # D. Format allele-frequency data
  # ----------------------------------------------------------
  
  df_af <- df_filt %>%
    dplyr::select(
      chrom,
      pos_Mb,
      m1_AF,
      m2_AF,
      u1_AF,
      u2_AF
    ) %>%
    pivot_longer(
      cols = ends_with("_AF"),
      names_to = "sample",
      values_to = "AF"
    )
  
  
  # ----------------------------------------------------------
  # E. Calculate log10 odds ratios
  #
  # Rep 1: m1 vs uni1
  # Rep 2: m2 vs uni2
  # ----------------------------------------------------------
  
  clamp <- function(p) {
    ifelse(
      p == 1,
      1 - 1e-6,
      ifelse(p == 0, 1e-6, p)
    )
  }
  
  
  df_or <- df_filt %>%
    mutate(
      
      odds_m1 = clamp(m1_AF) / (1 - clamp(m1_AF)),
      odds_m2 = clamp(m2_AF) / (1 - clamp(m2_AF)),
      
      odds_u1 = clamp(u1_AF) / (1 - clamp(u1_AF)),
      odds_u2 = clamp(u2_AF) / (1 - clamp(u2_AF)),
      
      log10_OR_rep1 = log10(odds_m1 / odds_u1),
      log10_OR_rep2 = log10(odds_m2 / odds_u2)
      
    ) %>%
    dplyr::select(
      chrom,
      pos_Mb,
      log10_OR_rep1,
      log10_OR_rep2
    ) %>%
    pivot_longer(
      cols = starts_with("log10_OR"),
      names_to = "group",
      values_to = "log10_OR"
    )
  
  
  # ==========================================================
  # F. Allele-frequency plot
  # ==========================================================
  
  p_af <- ggplot(
    df_af,
    aes(
      x = pos_Mb,
      y = AF,
      color = sample
    )
  ) +
    
    geom_smooth(
      se = FALSE,
      linewidth = 0.8
    ) +
    
    facet_wrap(
      ~chrom,
      scales = "free_x"
    ) +
    
    scale_color_manual(
      values = c(
        "m1_AF" = "darkgreen",
        "m2_AF" = "green3",
        "u1_AF" = "dodgerblue",
        "u2_AF" = "blue"
      ),
      
      labels = c(
        "m1_AF" = "100% Pathogen (Rep 1)",
        "m2_AF" = "100% Pathogen (Rep 2)",
        "u1_AF" = "Uninfected (Rep 1)",
        "u2_AF" = "Uninfected (Rep 2)"
      )
    ) +
    
    labs(
      title = paste0(
        strain_name,
        ": LOESS Allele Frequency"
      ),
      x = "Genomic Position (Mb)",
      y = "N2 Allele Frequency",
      color = "Sample Group"
    ) +
    
    supp_theme +
    
    guides(
      color = guide_legend(ncol = 2)
    )
  
  
  # ==========================================================
  # G. Odds-ratio plot
  # ==========================================================
  
  p_or <- ggplot(
    df_or,
    aes(
      x = pos_Mb,
      y = log10_OR,
      color = group
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey"
    ) +
    
    geom_smooth(
      se = FALSE,
      linewidth = 0.8
    ) +
    
    facet_wrap(
      ~chrom,
      scales = "free_x"
    ) +
    
    scale_color_manual(
      values = c(
        "log10_OR_rep1" = "darkgreen",
        "log10_OR_rep2" = "green3"
      ),
      
      labels = c(
        "log10_OR_rep1" = "100% Pathogen vs Control (Rep 1)",
        "log10_OR_rep2" = "100% Pathogen vs Control (Rep 2)"
      )
    ) +
    
    labs(
      title = paste0(
        strain_name,
        ": LOESS Log10 Odds Ratio"
      ),
      x = "Genomic Position (Mb)",
      y = "log10 Odds Ratio (Infected/Control)",
      color = "Sample Group"
    ) +
    
    supp_theme +
    
    guides(
      color = guide_legend(ncol = 1)
    )
  
  
  return(
    list(
      af = p_af,
      or = p_or
    )
  )
}


# ============================================================
# 3. Process JU2250 and EG4946
# ============================================================

plots_ju2250 <- process_bsa(
  "JU2250",
  "BSA_results/JU2250/JU2250_allele_depths.tsv"
)

plots_eg4946 <- process_bsa(
  "EG4946",
  "BSA_results/EG4946/EG4946_allele_depths.tsv"
)


# ============================================================
# 4. Combine into 2 x 2 supplementary figure
# ============================================================

supp_fig_28 <-
  (plots_ju2250$af | plots_ju2250$or) /
  (plots_eg4946$af | plots_eg4946$or) +
  
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 10,
        face = "plain"
      ),
      plot.tag.position = "topleft"
    )
  ) &
  
  theme(
    legend.position = "bottom"
  )


# ============================================================
# 5. Export Supplementary Figure 24
# ============================================================

ggsave(
  "supp_fig28.pdf",
  plot = supp_fig_28,
  width = 180,
  height = 220,
  units = "mm",
  device = cairo_pdf
)

read_raw_s28 <- function(file) {
  readr::read_tsv(file, show_col_types = FALSE,
                  col_names = c("chrom", "pos",
                                "path100_rep1_REF_ALT", "path100_rep2_REF_ALT",
                                "uninfected_rep1_REF_ALT", "uninfected_rep2_REF_ALT"))
}

raw_s28_ju2250 <- read_raw_s28("BSA_results/JU2250/JU2250_allele_depths.tsv")
raw_s28_eg4946 <- read_raw_s28("BSA_results/EG4946/EG4946_allele_depths.tsv")

# Must be below 1,048,576 (Excel's row limit)
nrow(raw_s28_ju2250); nrow(raw_s28_eg4946)