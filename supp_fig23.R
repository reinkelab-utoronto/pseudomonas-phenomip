library(tidyverse)
library(QTLseqr)
library(patchwork)

# ==============================================================================
# 0. Shared Plot Styling
# ==============================================================================
my_plot_theme <- theme_bw(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title      = element_blank(), # Removes default QTLseqr titles
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 6),
    strip.text      = element_text(size = 8, face = "bold"),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    plot.margin     = margin(5, 5, 5, 5)
  )

# ==============================================================================
# 1. Load Pre-Calculated Data (G8: Uninfected vs 50% Pathogen)
# ==============================================================================
# Loading the saved CSVs from your pipeline
ju2522_low_data <- read.csv("JU2522_G6_uni_vs_low_results.csv")
wn2050_low_data <- read.csv("WN2050_G6_uni_vs_low_results.csv")

# ==============================================================================
# 2. Build the Plots
# ==============================================================================

# JU2522 (G8 half vs uni)
p4a <- plotQTLStats(SNPset = ju2522_low_data, var = "Gprime", plotThreshold = TRUE, q = 0.01) + 
  my_plot_theme +
  ggtitle("JU2522") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

p4b <- plotQTLStats(SNPset = ju2522_low_data, var = "deltaSNP", plotIntervals = TRUE) + 
  my_plot_theme + 
  ggtitle("JU2522") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0)) 

# WN2050 (G8 half vs uni)
p4c <- plotQTLStats(SNPset = wn2050_low_data, var = "Gprime", plotThreshold = TRUE, q = 0.01) + 
  my_plot_theme +
  ggtitle("WN2050") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

p4d <- plotQTLStats(SNPset = wn2050_low_data, var = "deltaSNP", plotIntervals = TRUE) + 
  my_plot_theme + 
  ggtitle("WN2050") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

# Placeholder for E (The empty bottom row)
placeholder_e <- ggplot() + theme_void(base_size = 8)

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================
# Layout design: A 3x2 grid where 'EE' makes the final plot span the bottom row.
layout_design <- "
  AB
  CD
"
# Assemble directly without wrap_elements()
supp_fig23 <- p4a + p4b + p4c + p4d + placeholder_e +
  plot_layout(design = layout_design) +
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
# Export
ggsave("supp_fig23.pdf", 
       plot = supp_fig23, 
       width = 180, 
       height = 220*2/3, 
       units = "mm", 
       device = cairo_pdf)


# JU2522 (G8, uninfected vs half pathogen)
ju2522_low_qtl <- getQTLTable(
  SNPset   = ju2522_low_data,
  method   = "Gprime",
  alpha    = 0.01,
  export   = TRUE,
  fileName = "JU2522_G8_QTL_low_regions.csv"
)

# WN2050 (G8, uninfected vs half pathogen)
wn2050_low_qtl <- getQTLTable(
  SNPset   = wn2050_low_data,
  method   = "Gprime",
  alpha    = 0.01,
  export   = TRUE,
  fileName = "WN2050_G8_QTL_low_regions.csv"
)

