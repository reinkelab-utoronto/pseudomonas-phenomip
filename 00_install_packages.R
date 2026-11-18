# Installs the R packages used by this repository.
# Run once before running the pipeline:  source("00_install_packages.R")

cran_pkgs <- c("ape", "broom", "cluster", "corrplot", "cowplot", "data.table", "dendextend",
               "DescTools", "dplyr", "dunn.test", "e1071", "fitdistrplus", "FSA", "ggbeeswarm",
               "ggdendro", "ggplot2", "ggpubr", "ggrepel", "ggtext", "gplots", "jsonlite", "lme4",
               "nortest", "openxlsx", "patchwork", "pheatmap", "plyr", "pracma", "psych", "purrr",
               "readr", "reshape2", "rnaturalearth", "rnaturalearthdata", "rstatix", "scales", "sf", "slider",
               "stringr", "tibble", "tidyr", "tidyverse", "writexl",
               "BiocManager", "remotes")

missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

# Bioconductor
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree", update = FALSE)

# GitHub
if (!requireNamespace("QTLseqr", quietly = TRUE)) remotes::install_github("bmansfeld/QTLseqr")

# High-resolution map data used by worldmapv2.R (rnaturalearth scale = "large")
if (!requireNamespace("rnaturalearthhires", quietly = TRUE))
  install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

# Report anything that still failed to install
all_pkgs <- c(setdiff(cran_pkgs, c("BiocManager", "remotes")), "ggtree", "QTLseqr", "rnaturalearthhires")
failed <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed) > 0) {
  stop("These packages could not be installed: ", paste(failed, collapse = ", "))
}
message("All packages installed.")
