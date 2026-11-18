library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(psych)
library(corrplot)
library(grid)

# Needs `slope` and `AUC` from the Supp Fig 16 script

custom_order <- c("PP", "M2", "LVS", "M1", "MA", "LUA", "M5", "PV", "PA", "M117", "M114", "M395")

# ==============================================================================
# 0. HELPERS
# ==============================================================================
# psych::corr.test() returns a p matrix with UNADJUSTED p-values below the
# diagonal and ADJUSTED p-values above it. Reordering rows/columns moves values
# between the two triangles, so build symmetric matrices before reordering.
sym_p <- function(ct) {
  p_raw <- ct$p
  p_raw[upper.tri(p_raw)] <- t(p_raw)[upper.tri(p_raw)]   # unadjusted, both triangles
  p_adj <- p_raw
  p_adj[lower.tri(p_adj)] <- p.adjust(p_raw[lower.tri(p_raw)], method = "BH")
  p_adj[upper.tri(p_adj)] <- t(p_adj)[upper.tri(p_adj)]   # adjusted, both triangles
  diag(p_raw) <- 0
  diag(p_adj) <- 0
  
  # Check against corr.test's own adjusted values (original order, upper triangle)
  if (!isTRUE(all.equal(p_adj[upper.tri(p_adj)], ct$p[upper.tri(ct$p)]))) {
    warning("Rebuilt adjusted p-values do not match corr.test(); check before using.")
  }
  list(raw = p_raw, adj = p_adj)
}

# Long table: one row per condition pair, in the figure's order
cor_pairs <- function(ct, p, metric) {
  r     <- ct$r[custom_order, custom_order]
  p_raw <- p$raw[custom_order, custom_order]
  p_adj <- p$adj[custom_order, custom_order]
  n     <- if (length(ct$n) == 1) {
    matrix(ct$n, nrow(ct$r), ncol(ct$r), dimnames = dimnames(ct$r))
  } else {
    ct$n
  }
  n <- n[custom_order, custom_order]
  
  idx <- which(upper.tri(r), arr.ind = TRUE)
  tibble::tibble(
    Metric       = metric,
    condition_1  = rownames(r)[idx[, 1]],
    condition_2  = colnames(r)[idx[, 2]],
    n            = n[idx],
    spearman_rho = r[idx],
    p_unadj      = p_raw[idx],
    p_adj_BH     = p_adj[idx]
  ) %>%
    dplyr::mutate(signif = dplyr::case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    ))
}

# ==============================================================================
# 1. DATA PREPARATION: Slope Data
# ==============================================================================
slope_wide <- slope %>%
  dplyr::select(strain, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value)

cor_test_results_slope <- corr.test(
  slope_wide %>% dplyr::select(-strain),
  use = "pairwise",
  method = "spearman",
  adjust = "BH"
)

p_slope <- sym_p(cor_test_results_slope)

cor_matrix_slope <- cor_test_results_slope$r[custom_order, custom_order]
p_matrix_slope   <- p_slope$adj[custom_order, custom_order]   # BH-adjusted, symmetric

# ==============================================================================
# 2. DATA PREPARATION: AUC Data
# ==============================================================================
AUC_wide <- AUC %>%
  dplyr::select(strain, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value)

cor_test_results_AUC <- corr.test(
  AUC_wide %>% dplyr::select(-strain),
  use = "pairwise",
  method = "spearman",
  adjust = "BH"
)

p_AUC <- sym_p(cor_test_results_AUC)

cor_matrix_AUC <- cor_test_results_AUC$r[custom_order, custom_order]
p_matrix_AUC   <- p_AUC$adj[custom_order, custom_order]       # BH-adjusted, symmetric

# Both should be TRUE (symmetric, so no triangle mixing)
isSymmetric(unname(p_matrix_slope))
isSymmetric(unname(p_matrix_AUC))

# ==============================================================================
# 3. CREATE PANELS (Wrapping Base R corrplot for Patchwork)
# ==============================================================================
panel_A <- wrap_elements(panel = ~ {
  par(mar = c(4, 4, 3, 2))
  corrplot(
    cor_matrix_slope,
    method = "color",
    type = "upper",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    p.mat = p_matrix_slope,
    sig.level = 0.05,
    insig = "blank",
    diag = FALSE,
    number.cex = 0.5,
    col = colorRampPalette(c("blue", "white", "red"))(200),
    addgrid.col = "gray20",
    tl.cex = 0.6,
    cl.cex = 0.6
  )
})

panel_B <- wrap_elements(panel = ~ {
  par(mar = c(4, 4, 3, 2))
  corrplot(
    cor_matrix_AUC,
    method = "color",
    type = "upper",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    p.mat = p_matrix_AUC,
    sig.level = 0.05,
    insig = "blank",
    diag = FALSE,
    number.cex = 0.5,
    col = colorRampPalette(c("blue", "white", "red"))(200),
    addgrid.col = "gray20",
    tl.cex = 0.6,
    cl.cex = 0.6
  )
})

# ==============================================================================
# 4. ASSEMBLE PANELS AND EXPORT TO PDF
# ==============================================================================
combined_corr_fig <- panel_A + panel_B +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 10, face = "plain"),
      plot.tag.position = "topleft"
    )
  )

ggsave(
  filename = "supp_fig17.pdf",
  plot = combined_corr_fig,
  width = 180,
  height = 220/2,
  units = "mm",
  device = cairo_pdf
)

# ==============================================================================
# 5. STATS TABLE FOR SUPPLEMENTARY DATA
# ==============================================================================
# One row per condition pair; cells shown on the figure are signif != "ns".
# No normality test needed: Spearman's rho is rank-based.
stats_s17 <- dplyr::bind_rows(
  cor_pairs(cor_test_results_slope, p_slope, "regression"),   # Panel A
  cor_pairs(cor_test_results_AUC,   p_AUC,   "AUC")            # Panel B
)

stats_s17