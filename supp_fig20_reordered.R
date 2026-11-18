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
# 1. Load Pre-Calculated Data
# ==============================================================================
ju2610_low_data <- read.csv("JU2610_uni_vs_low_results.csv")
ed3005_low_data <- read.csv("ED3005_uni_vs_low_results.csv")

# ==============================================================================
# 2. Build the Plots
# ==============================================================================

# JU2610
p3a <- plotQTLStats(SNPset = ju2610_low_data, var = "Gprime", plotThreshold = TRUE, q = 0.01) + 
  my_plot_theme +
  ggtitle("JU2610") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

p3b <- plotQTLStats(SNPset = ju2610_low_data, var = "deltaSNP", plotIntervals = TRUE) + 
  my_plot_theme +
  ggtitle("JU2610") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

# ED3005
p3c <- plotQTLStats(SNPset = ed3005_low_data, var = "Gprime", plotThreshold = TRUE, q = 0.01) + 
  my_plot_theme +
  ggtitle("ED3005") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

p3d <- plotQTLStats(SNPset = ed3005_low_data, var = "deltaSNP", plotIntervals = TRUE) + 
  my_plot_theme +
  ggtitle("ED3005") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

# Placeholder for E
placeholder_e <- ggplot() + theme_void(base_size = 8)

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================
# Layout design: 'EE' makes the final plot span the bottom row. 
layout_design <- "
  AB
  CD
"
# Assemble directly without wrap_elements() so the plots fill their designated space naturally
final_figure_3 <- p3a + p3b + p3c + p3d + placeholder_e +
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
ggsave("supp_fig20.pdf", 
       plot = final_figure_3, 
       width = 180, 
       height = 220*2/3, 
       units = "mm", 
       device = cairo_pdf)




library(QTLseqr)

# JU2610
ju2610_low_qtl <- getQTLTable(
  SNPset   = ju2610_low_data,
  method   = "Gprime",
  alpha    = 0.01,        # same q you used in plotQTLStats
  export   = TRUE,
  fileName = "ju2610_low_qtl_low_regions.csv"
)

# ED3005
ed3005_low_qtl <- getQTLTable(
  SNPset   = ed3005_low_data,
  method   = "Gprime",
  alpha    = 0.01,
  export   = TRUE,
  fileName = "ed3005_low_qtl_low_regions.csv"
)

