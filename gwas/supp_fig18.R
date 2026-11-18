# ==============================================================================
# Combined 2-Trait Figure: Manhattan + Violin (PxG) panels
# Layout: 2 columns x 2 rows, tagged A-D
#   A = M114 avg_auc_zscore Manhattan     B = M114 avg_auc_zscore Violin
#   C = M395 avg_slope_zscore Manhattan   D = M395 avg_slope_zscore Violin
# Only traits with a genome-wide significant QTL are plotted.
# Phenotypes are read directly from each NemaScan report (.html).
# ==============================================================================
# Install once if needed:
# install.packages(c("dplyr", "tidyr", "ggplot2", "jsonlite", "ggbeeswarm",
#                     "ggrepel", "patchwork", "ggtext"))
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
library(ggbeeswarm)
library(ggrepel)
library(patchwork)
library(ggtext)

chr_levels <- c("I", "II", "III", "IV", "V", "X")   # MtDNA excluded throughout

base_txt_size <- 8

house_theme <- theme(
  axis.title      = element_text(size = 8),
  axis.text       = element_text(size = 6),
  strip.text      = element_text(size = 8, face = "bold"),
  plot.title      = ggtext::element_markdown(size = 8, face = "bold"),
  plot.subtitle   = element_text(size = 8),
  legend.title    = element_text(size = 8),
  legend.text     = element_text(size = 6),
  legend.key.size = unit(3, "mm"),
  legend.position = "bottom"
)

## =============================================================================
## ---- USER CONFIG: one block per trait ---------------------------------------
## =============================================================================

traits <- list(
  
  M114 = list(
    trait_name        = "M114_avg_auc_zscore",
    display_title     = "*P. poae* MYb114",
    method            = "loco",
    gwa_file          = "M114_avg_auc_zscore_loco.gwa",
    qtl_file          = "M114_avg_auc_zscore_loco_qtl.tsv",
    n_tests_file      = "M114_independent_tests.txt",
    genotype_matrix   = "M114_genotype_matrix.tsv",
    heritability_file = "M114_trait_heritability.tsv",
    report_html       = "NemaScan_Report_M114_avg_auc_zscore.html",
    which_qtl         = 1,
    n_label_each_end  = 5,
    violin_y_label    = "z-score"
  ),
  
  M395 = list(
    trait_name        = "M395_avg_slope_zscore",
    display_title     = "*P. canadensis* MYb395",
    method            = "inbred",
    gwa_file          = "M395_avg_slope_zscore_inbred.gwa",
    qtl_file          = "M395_avg_slope_zscore_inbred_qtl.tsv",
    n_tests_file      = "M395_independent_tests.txt",
    genotype_matrix   = "M395_genotype_matrix.tsv",
    heritability_file = "M395_trait_heritability.tsv",
    report_html       = "NemaScan_Report_M395_avg_slope_zscore.html",
    which_qtl         = 1,
    n_label_each_end  = 5,
    violin_y_label    = "z-score"
  )
)

out_file   <- "../supp_fig18.pdf"
out_width  <- 180   # mm
out_height <- 220*(2/3)   # mm

## =============================================================================
## ---- Shared functions --------------------------------------------------------
## =============================================================================

extract_all_data_json <- function(html_path) {
  html_text <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  
  start_marker <- "const all_data="
  start_idx <- regexpr(start_marker, html_text, fixed = TRUE)
  start_pos <- start_idx + attr(start_idx, "match.length")
  
  chars <- strsplit(html_text, "")[[1]]
  n <- length(chars)
  
  depth <- 0
  in_string <- FALSE
  escape <- FALSE
  end_pos <- NA
  
  i <- start_pos
  while (i <= n) {
    c <- chars[i]
    if (in_string) {
      if (escape) {
        escape <- FALSE
      } else if (c == "\\") {
        escape <- TRUE
      } else if (c == "\"") {
        in_string <- FALSE
      }
    } else {
      if (c == "\"") {
        in_string <- TRUE
      } else if (c == "{") {
        depth <- depth + 1
      } else if (c == "}") {
        depth <- depth - 1
        if (depth == 0) {
          end_pos <- i
          break
        }
      }
    }
    i <- i + 1
  }
  
  json_txt <- substr(html_text, start_pos, end_pos)
  fromJSON(json_txt, simplifyVector = FALSE)
}

get_phenotype_df <- function(cfg) {
  # Always read the phenotypes from the HTML report, so they always match the
  # report's GWA results (no cached copy that could go out of date).
  message("[", cfg$trait_name, "] Extracting phenotypes from report: ", cfg$report_html)
  all_data <- extract_all_data_json(cfg$report_html)
  pheno_df <- data.frame(
    strain    = unlist(all_data$broad$strain),
    phenotype = unlist(all_data$broad$phenotype)
  )
  
  pheno_df
}

standardize_gwa_cols <- function(df) {
  # .gwa files vary in header casing/naming depending on inbred vs loco
  # (fastGWA-style vs GCTA mlma-style output), e.g.:
  #   inbred: CHR   SNP  POS  A1  A2  N     AF1   BETA  SE  P
  #   loco:   Chr   SNP  bp   A1  A2  Freq  b     se    p
  # Match case-insensitively against the ORIGINAL header (snapshot taken
  # once, before any renaming) so lookups don't get confused mid-loop.
  orig_cols <- colnames(df)
  upper_cols <- toupper(orig_cols)
  
  rename_if_missing <- function(df, canonical, candidates) {
    if (!(canonical %in% colnames(df))) {
      idx <- which(upper_cols %in% toupper(candidates))
      if (length(idx) > 0) colnames(df)[idx[1]] <- canonical
    }
    df
  }
  
  df <- rename_if_missing(df, "CHR", c("CHR", "CHROM"))
  df <- rename_if_missing(df, "POS", c("POS", "BP"))
  df <- rename_if_missing(df, "P",   c("P", "PVAL", "P_WALD"))
  df
}

build_manhattan_plot <- function(cfg) {
  gwa_raw <- read.delim(cfg$gwa_file, stringsAsFactors = FALSE)
  gwa_raw <- standardize_gwa_cols(gwa_raw)
  n_tests <- as.numeric(readLines(cfg$n_tests_file))
  qtl <- read.delim(cfg$qtl_file, stringsAsFactors = FALSE)
  h2_df <- read.delim(cfg$heritability_file, header = FALSE, stringsAsFactors = FALSE) %>%
    distinct()
  h2_value <- as.numeric(h2_df[1, 2])
  
  gwa <- gwa_raw %>% filter(CHR != "MtDNA")
  qtl <- qtl %>% filter(CHROM != "MtDNA")
  
  gwa$CHR <- factor(gwa$CHR, levels = chr_levels)
  gwa <- gwa %>% filter(!is.na(CHR), !is.na(P)) %>% arrange(CHR, POS)
  
  chr_offsets <- gwa %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POS)) %>%
    mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))
  
  gwa <- gwa %>%
    left_join(chr_offsets, by = "CHR") %>%
    mutate(BPcum = POS + offset, logP = -log10(P))
  
  axis_df <- gwa %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  sig_threshold <- -log10(0.05 / n_tests)
  
  qtl$CHR <- factor(qtl$CHROM, levels = chr_levels)
  qtl <- qtl %>%
    left_join(chr_offsets, by = "CHR") %>%
    mutate(peak_BPcum = peakPOS + offset,
           interval_bp = endPOS - startPOS,
           size_label = case_when(
             interval_bp < 1e6 ~ paste0(round(interval_bp / 1e3, 1), " kb"),
             TRUE              ~ paste0(round(interval_bp / 1e6, 2), " Mb")
           ),
           h2 = h2_value,
           annotation = paste0(size_label, "\nh\u00b2 = ", round(h2, 3))) %>%
    left_join(gwa %>% dplyr::select(CHR, POS, logP), by = c("CHR" = "CHR", "peakPOS" = "POS"))
  
  p <- ggplot(gwa, aes(x = BPcum, y = logP, color = CHR)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = rep(c("#3B4992", "#8491B4"), length.out = length(chr_levels))) +
    geom_hline(yintercept = sig_threshold, color = "red", linetype = "dashed", linewidth = 0.5) +
    geom_point(data = gwa %>% filter(logP >= sig_threshold), color = "#EE0000", size = 1.4) +
    geom_vline(data = qtl, aes(xintercept = peak_BPcum), color = "black", linetype = "dotted", linewidth = 0.4) +
    geom_text(data = qtl, aes(x = peak_BPcum, y = logP, label = annotation),
              inherit.aes = FALSE, vjust = -0.3, hjust = 0.5,
              size = 8 / .pt, color = "black", lineheight = 0.9) +
    scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$center) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))),
         title = paste0(if (!is.null(cfg$display_title)) cfg$display_title else cfg$trait_name,
                        " (", cfg$method, ")")) +
    theme_bw(base_size = base_txt_size) +
    house_theme +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  message(sprintf("[%s] sig threshold -log10P = %.3f | independent tests = %.0f | SNPs above threshold = %d",
                  cfg$trait_name, sig_threshold, n_tests, sum(gwa$logP >= sig_threshold)))
  
  # Supplementary-data versions (plotting helper columns removed)
  gwa_out <- gwa %>%
    dplyr::select(-any_of(c("chr_len", "offset", "BPcum"))) %>%
    mutate(CHR = as.character(CHR),
           above_threshold = logP >= sig_threshold)
  
  qtl_out <- qtl %>%
    dplyr::select(-any_of(c("CHR", "chr_len", "offset", "peak_BPcum", "size_label", "annotation"))) %>%
    rename(peak_logP = logP)
  
  summary_out <- data.frame(
    n_markers             = nrow(gwa),
    n_independent_tests   = n_tests,
    threshold_alpha       = 0.05 / n_tests,
    threshold_neglog10p   = sig_threshold,
    n_markers_above       = sum(gwa$logP >= sig_threshold),
    n_qtl                 = nrow(qtl),
    h2                    = h2_value
  )
  
  list(plot = p, qtl = qtl, gwa_out = gwa_out, qtl_out = qtl_out, summary_out = summary_out)
}

build_violin_plot <- function(cfg, qtl) {
  pheno_df <- get_phenotype_df(cfg)
  peak_marker <- qtl %>% slice(cfg$which_qtl)
  
  geno_row <- read.delim(cfg$genotype_matrix, stringsAsFactors = FALSE) %>%
    filter(CHROM == as.character(peak_marker$CHR), POS == peak_marker$peakPOS)
  ref_base <- geno_row$REF[1]   # reference nucleotide at the peak marker
  alt_base <- geno_row$ALT[1]   # alternative nucleotide at the peak marker
  
  geno <- geno_row %>%
    dplyr::select(-CHROM, -POS, -REF, -ALT) %>%
    pivot_longer(everything(), names_to = "strain", values_to = "geno_call") %>%
    mutate(allele = case_when(
      geno_call == -1 ~ "REF",
      geno_call == 1  ~ "ALT",
      TRUE ~ NA_character_
    ))
  
  pxg <- pheno_df %>%
    inner_join(geno, by = "strain") %>%
    filter(!is.na(allele)) %>%
    mutate(allele = factor(allele, levels = c("REF", "ALT")))
  
  top_strains <- pxg %>% arrange(desc(phenotype)) %>% slice_head(n = cfg$n_label_each_end) %>% pull(strain)
  bottom_strains <- pxg %>% arrange(phenotype) %>% slice_head(n = cfg$n_label_each_end) %>% pull(strain)
  
  pxg <- pxg %>%
    mutate(
      strain_category = case_when(
        strain %in% top_strains    ~ "Most Positive",
        strain %in% bottom_strains ~ "Most Negative",
        TRUE ~ "Other"
      ),
      label = ifelse(strain %in% c(top_strains, bottom_strains), strain, "")
    )
  
  message(sprintf("[%s] Peak marker %s:%s | strains plotted: %d (REF: %d, ALT: %d)",
                  cfg$trait_name, peak_marker$CHR, peak_marker$peakPOS,
                  nrow(pxg), sum(pxg$allele == "REF"), sum(pxg$allele == "ALT")))
  
  bee_pos <- position_beeswarm(cex = 2)
  violin_title <- paste0("Genotype at ", peak_marker$CHR, ":", format(peak_marker$peakPOS, big.mark = ","))
  
  p <- ggplot(pxg, aes(x = allele, y = phenotype)) +
    geom_violin(aes(fill = allele), alpha = 0.5, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = c(REF = "#A79F92", ALT = "mediumpurple4"), guide = "none") +
    geom_point(aes(color = strain_category, size = strain_category),
               position = bee_pos, shape = 16, alpha = 0.8) +
    scale_color_manual(values = c("Most Positive" = "blue", "Most Negative" = "red", "Other" = "grey40")) +
    scale_size_manual(values = c("Most Positive" = 1.6, "Most Negative" = 1.6, "Other" = 0.9), guide = "none") +
    geom_text_repel(aes(label = label), position = bee_pos,
                    size = 6 / .pt, color = "black", max.overlaps = 100, segment.size = 0.3) +
    labs(x = NULL, y = cfg$violin_y_label, title = violin_title, color = "Strain") +
    theme_bw(base_size = base_txt_size) +
    house_theme
  
  # Supplementary-data version: EVERY phenotyped strain at the peak marker,
  # REF and ALT (labelled or not), plus any strain without a REF/ALT call
  # (not plotted). Sorted by allele, then highest to lowest phenotype.
  ranks <- pxg %>%
    mutate(phenotype_rank = rank(-phenotype, ties.method = "min")) %>%   # 1 = highest, among plotted strains
    dplyr::select(strain, phenotype_rank, strain_category)
  
  pxg_out <- pheno_df %>%
    left_join(geno, by = "strain") %>%
    left_join(ranks, by = "strain") %>%
    mutate(
      chrom         = as.character(peak_marker$CHR),
      peak_pos      = peak_marker$peakPOS,
      peak_marker   = paste0(chrom, ":", peak_pos),
      ref_base      = ref_base,
      alt_base      = alt_base,
      allele        = ifelse(is.na(allele), "No call (not plotted)", allele),
      genotype_base = case_when(allele == "REF" ~ ref_base,
                                allele == "ALT" ~ alt_base,
                                TRUE ~ NA_character_),
      plotted            = allele %in% c("REF", "ALT"),
      labelled_in_figure = strain %in% c(top_strains, bottom_strains),
      label_group = case_when(
        strain_category == "Most Positive" ~ paste0("Top ", cfg$n_label_each_end, " (highest phenotype, blue)"),
        strain_category == "Most Negative" ~ paste0("Bottom ", cfg$n_label_each_end, " (lowest phenotype, red)"),
        plotted ~ "Not labelled (grey)",
        TRUE ~ NA_character_
      ),
      phenotype_measure = paste0(cfg$trait_name, " (", cfg$violin_y_label, ")")
    ) %>%
    arrange(factor(allele, levels = c("REF", "ALT", "No call (not plotted)")), desc(phenotype)) %>%
    transmute(strain, peak_marker, chrom, peak_pos, ref_base, alt_base,
              allele, genotype_base, geno_call, phenotype_measure, phenotype,
              phenotype_rank, plotted, labelled_in_figure, label_group)
  
  list(plot = p, pxg_out = pxg_out)
}

## =============================================================================
## ---- Build all four panels ----------------------------------------------------
## =============================================================================

m114_man <- build_manhattan_plot(traits$M114)
m114_vio <- build_violin_plot(traits$M114, m114_man$qtl)

m395_man <- build_manhattan_plot(traits$M395)
m395_vio <- build_violin_plot(traits$M395, m395_man$qtl)

## =============================================================================
## ---- Assemble 2-column x 2-row figure, tagged A-D ----------------------------
## =============================================================================

layout_design <- "
AB
CD
"

final_figure <- m114_man$plot + m114_vio$plot +
  m395_man$plot + m395_vio$plot +
  plot_layout(design = layout_design) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )

print(final_figure)

ggsave(out_file, final_figure, width = out_width, height = out_height, units = "mm",
       device = if (grepl("\\.pdf$", out_file)) cairo_pdf else NULL)

cat(sprintf("\nSaved combined 4-panel figure to: %s\n", out_file))

## =============================================================================
## ---- Supplementary data tables -----------------------------------------------
## =============================================================================

# GWAS summary: threshold, markers above it, number of QTL, heritability (all traits)
gwas_summary_s18 <- bind_rows(
  M114 = m114_man$summary_out,
  M395 = m395_man$summary_out,
  .id = "trait"
) %>%
  mutate(trait_name = c(traits$M114$trait_name, traits$M395$trait_name),
         method     = c(traits$M114$method,     traits$M395$method),
         .after = trait)

# QTL intervals (all traits)
qtl_s18 <- bind_rows(
  M114 = m114_man$qtl_out,
  M395 = m395_man$qtl_out,
  .id = "trait"
)

# Phenotype x genotype at the plotted peak marker (panels B, D)
# Long version: one row per strain (kept for reference / per-strain values)
pxg_long_s18 <- bind_rows(
  M114 = m114_vio$pxg_out,
  M395 = m395_vio$pxg_out,
  .id = "trait"
)

# Compact version: one row per trait/peak marker, strains listed by genotype
# (alphabetical, comma-separated)
pxg_s18 <- pxg_long_s18 %>%
  group_by(trait, peak_marker, ref_base, alt_base) %>%
  summarise(
    n_REF           = sum(allele == "REF"),
    REF_strains     = paste(sort(strain[allele == "REF"]), collapse = ", "),
    n_ALT           = sum(allele == "ALT"),
    ALT_strains     = paste(sort(strain[allele == "ALT"]), collapse = ", "),
    n_no_call       = sum(!plotted),
    no_call_strains = paste(sort(strain[!plotted]), collapse = ", "),
    labelled_top_strains    = paste(strain[grepl("^Top", label_group)],    collapse = ", "),   # blue labels
    labelled_bottom_strains = paste(strain[grepl("^Bottom", label_group)], collapse = ", "),   # red labels
    .groups = "drop"
  ) %>%
  arrange(factor(trait, levels = c("M114", "M395")))

# Per-allele summary of what each violin shows (plotted strains only)
pxg_summary_s18 <- pxg_long_s18 %>%
  filter(plotted) %>%
  group_by(trait, peak_marker, allele, genotype_base) %>%
  summarise(
    n_strains     = n(),
    mean          = mean(phenotype),
    sd            = sd(phenotype),
    median        = median(phenotype),
    q1            = quantile(phenotype, 0.25),
    q3            = quantile(phenotype, 0.75),
    min           = min(phenotype),
    max           = max(phenotype),
    n_labelled    = sum(labelled_in_figure),
    .groups = "drop"
  ) %>%
  arrange(trait, factor(allele, levels = c("REF", "ALT")))

# Full GWA results, one data frame per trait (panels A, C)
gwa_s18_m114 <- m114_man$gwa_out
gwa_s18_m395 <- m395_man$gwa_out

# Excel sheets hold at most 1,048,576 rows
sapply(list(M114 = gwa_s18_m114, M395 = gwa_s18_m395), nrow)
