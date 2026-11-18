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

# ==============================================================================
# 1. Load Data
# ==============================================================================
df_JU <- readr::read_tsv("JU2250_RILs_joint.filtered.tsv") %>%
  filter(chrom != "MtDNA") %>%
  mutate(gt = ifelse(gt == 1, "JYIS8", "JU2250"),
         sample = sub("_S\\d+_L001\\.sorted\\.bam", "", sample)) %>%
  filter(!sample == "JU2250_CAENDR")

df_EG <- readr::read_tsv("EG4946_RILs_joint.filtered.tsv") %>%
  filter(chrom != "MtDNA") %>%
  mutate(gt = ifelse(gt == 1, "JYIS8", "EG4946"),
         sample = sub("_S\\d+_L001\\.sorted\\.bam", "", sample)) %>%
  filter(!sample == "EG4946_CAENDR")

s29_breakpoints <- rbind(df_JU, df_EG)
# ==============================================================================
# 2. Helper function: clean sample names/order + build genome-wide genotype plot
# ==============================================================================
## NOTE: "JYIS8" is kept as-is in all underlying data (gt, sample matching).
## Only the DISPLAYED text (y-axis labels + legend) is relabelled to "jyIs8".

build_genotype_plot <- function(df, wild_name, strain_prefix, fill_values, show_x_axis_title = TRUE) {
  
  df$sample <- gsub("NIL", "RIL", df$sample)
  df$sample <- gsub("_", " ", df$sample)
  df$sample <- gsub(paste0(strain_prefix, " RIL"), "RIL", df$sample)
  
  sample_names <- unique(df$sample)
  df$sample <- factor(df$sample, levels = c(sample_names[12], sample_names[1], sample_names[2:11]))
  df$sample <- factor(df$sample, levels = rev(levels(df$sample)))
  df$index  <- as.numeric(df$sample)
  
  df <- df %>%
    arrange(sample) %>%
    mutate(index = as.numeric(factor(sample, levels = levels(sample))))
  
  # display-only relabel: "JYIS8" -> "jyIs8" for the y-axis strain labels
  strain_idx <- setNames(
    gsub("JYIS8", "jyIs8", levels(df$sample)),
    1:length(levels(df$sample))
  )
  
  # chromosomes in natural order (not filtered to X)
  df$chrom <- factor(df$chrom, levels = c("I", "II", "III", "IV", "V", "X"))
  
  p <- ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
    scale_fill_manual(values = fill_values, name = "Genotype",
                      labels = function(x) gsub("JYIS8", "jyIs8", x)) +
    facet_grid(. ~ chrom, scales = "free", space = "free") +
    scale_x_continuous(labels = \(x) x / 1e6, expand = c(0, 0)) +
    scale_y_continuous(
      breaks = 1:length(strain_idx) + 0.5,
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
  
  p
}

# ==============================================================================
# 3. Build Panels A (JU2250) and B (EG4946)
# ==============================================================================

panel_a <- build_genotype_plot(
  df_JU, wild_name = "JU2250", strain_prefix = "JU",
  fill_values = c("JU2250" = "#0080FF", "JYIS8" = "#FF8000"),
  show_x_axis_title = FALSE   # top panel: hide x-axis, bottom panel supplies it
)

panel_b <- build_genotype_plot(
  df_EG, wild_name = "EG4946", strain_prefix = "EG",
  fill_values = c("EG4946" = "#0080FF", "JYIS8" = "#FF8000"),
  show_x_axis_title = TRUE
)

# ==============================================================================
# 4. Assemble: A on top, B on bottom
# ==============================================================================

supp_fig29 <- panel_a / panel_b +
  plot_layout(heights = c(1, 1)) +
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

supp_fig29

ggsave("supp_fig29.pdf",
       plot = supp_fig29,
       width = 180,
       height = 220,
       units = "mm",
       device = cairo_pdf)