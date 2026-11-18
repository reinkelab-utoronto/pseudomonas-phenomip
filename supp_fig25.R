library(tidyverse)
library(patchwork)

# ==============================================================================
# 0. Shared Plot Styling
# ==============================================================================
my_plot_theme <- theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title      = element_blank(),
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 6),
    strip.text      = element_text(size = 8, face = "bold"),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    plot.margin     = margin(5, 5, 5, 5)
  )

DONOR       <- "JYIS8"      # --A strain in vk hmm (gt = 1)
DONOR_LABEL <- "jyIs8"   # how the donor is labelled on the plot (as in Figs 3 and 4)

# Samples shown in panel B (names as they appear in the JU2522 genotype file)
KEEP_JU2522 <- c("JYIS8", "JU2522", "JU2522_RIL_J21", "JU2522_RIL_J24")

# ==============================================================================
# 1. Load Data
# ==============================================================================
# These TSVs come from the newer pipeline, so sample names are already clean
# (e.g. "JU2610_RIL_JU13", "JU2522_RIL_J11", "JU2610", "JYIS8")
read_geno <- function(file, parent) {
  readr::read_tsv(file, show_col_types = FALSE) %>%
    filter(chrom != "MtDNA") %>%
    mutate(gt = ifelse(gt == 1, DONOR, parent))
}

df_JU2610 <- read_geno("JU2610_RILs_joint.filtered.tsv", "JU2610")   # <-- EDIT path if needed
df_JU2522 <- read_geno("JU2522_RILs_joint.filtered.tsv", "JU2522") %>%   # <-- EDIT path if needed
  filter(sample %in% KEEP_JU2522)

# Stop early if any selected sample is missing from the file (e.g. a typo)
missing_JU2522 <- setdiff(KEEP_JU2522, unique(df_JU2522$sample))
if (length(missing_JU2522) > 0) {
  stop("Not found in the JU2522 genotype file: ", paste(missing_JU2522, collapse = ", "))
}

# ==============================================================================
# 2. Helper function: clean sample names/order + build genome-wide genotype plot
# ==============================================================================
## NOTE: "JYIS8" is kept as-is in all underlying data (gt, sample matching).
## Only the DISPLAYED text (y-axis labels + legend) is relabelled.

display_label <- function(x) ifelse(x == DONOR, DONOR_LABEL, x)

build_genotype_plot <- function(df, parent, show_x_axis_title = TRUE) {
  
  # "JU2610_RIL_JU13" -> "RIL JU13"
  df$sample <- gsub("_", " ", df$sample)
  df$sample <- sub(paste0("^", parent, " RIL "), "RIL ", df$sample)
  
  # Row order (top to bottom): donor, parent, then RILs by number
  # (JU4, JU13, ... / J5, J11, ...), set by name rather than file order
  rils   <- unique(df$sample[grepl("^RIL ", df$sample)])
  ril_id <- sub("^RIL ", "", rils)
  rils   <- rils[order(gsub("[0-9]", "", ril_id), as.numeric(gsub("[^0-9]", "", ril_id)))]
  
  line_order <- rev(c(DONOR, parent, rils))   # reversed: first strain at the top
  df$sample  <- factor(df$sample, levels = line_order)
  df$index   <- as.numeric(df$sample)
  
  strain_idx <- setNames(display_label(line_order), seq_along(line_order))
  
  # chromosomes in natural order
  df$chrom <- factor(df$chrom, levels = c("I", "II", "III", "IV", "V", "X"))
  
  ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
    scale_fill_manual(values = setNames(c("#0080FF", "#FF8000"), c(parent, DONOR)),
                      limits = c(parent, DONOR),
                      labels = display_label,
                      name   = "Genotype") +
    facet_grid(. ~ chrom, scales = "free", space = "free") +
    scale_x_continuous(labels = \(x) x / 1e6, expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq_along(strain_idx) + 0.5,
      labels = strain_idx,
      expand = c(0, 0)
    ) +
    labs(x = if (show_x_axis_title) "chromosomal position (Mb)" else NULL,
         y = NULL) +
    my_plot_theme +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 6),
          axis.title.x = if (show_x_axis_title) element_text(size = 8) else element_blank(),
          axis.text.x = if (show_x_axis_title) element_text(size = 6) else element_blank(),
          axis.ticks.x = if (show_x_axis_title) element_line() else element_blank())
}

# ==============================================================================
# 3. Build Panels A (JU2610) and B (JU2522)
# ==============================================================================
panel_a <- build_genotype_plot(df_JU2610, parent = "JU2610",
                               show_x_axis_title = FALSE)   # top panel: bottom one supplies the x-axis

panel_b <- build_genotype_plot(df_JU2522, parent = "JU2522",
                               show_x_axis_title = TRUE)

# ==============================================================================
# 4. Assemble: A on top, B on bottom
# ==============================================================================
# Panel heights proportional to the number of strains (now 6 vs 4), so every
# strain row is the same thickness in both panels
n_strains <- c(n_distinct(df_JU2610$sample), n_distinct(df_JU2522$sample))

supp_fig25 <- panel_a / panel_b +
  plot_layout(heights = n_strains) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

supp_fig25

# With 10 strain rows instead of 14, a shorter page may look better
# (e.g. height = 160)
ggsave("supp_fig25.pdf",
       plot = supp_fig25,
       width = 180,
       height = 220*2/3,
       units = "mm",
       device = cairo_pdf)
