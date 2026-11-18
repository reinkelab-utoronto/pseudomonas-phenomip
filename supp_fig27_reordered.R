library(tidyverse)
library(QTLseqr)
library(patchwork)
library(DescTools)

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
# 1. Load Data
# ==============================================================================
EG4946_data <- read.csv("EG4946_uni_vs_m_results.csv")

EG4946_embryo <- read.csv("EG4946_RILs_avg_embryos.csv", check.names = FALSE)   # <-- EDIT path if needed
EG4946_gravid <- read.csv("EG4946_RILs_gravid.csv",      check.names = FALSE)   # <-- EDIT path if needed
geno_raw      <- readr::read_tsv("EG4946_RILs_joint.filtered.tsv", show_col_types = FALSE)

# ---- Settings for panel C ----
PARENT <- "EG4946"   # --B strain in vk hmm (gt = 2)
DONOR  <- "JYIS8"    # --A strain in vk hmm (gt = 1)
DONOR_LABEL <- "jyIs8"   # how the donor is labelled on the plot (as in Figures 3 and 4)

# Phenotype files may call the donor "N2" (same strain as JYIS8 in sequencing)
PHENO_ALIASES <- c(N2 = DONOR)

# Genotyped samples to leave out of the plot (raw names in the TSV)
GENO_EXCLUDE <- c("EG4946_CAENDR")

# Chromosomes shown (in this order) in the genotype panel
GENO_CHROM <- c("X")

# Panel labels used in the supplementary tables (same style as Figures 3E, 4E, 5D-E)
panel_labels_s27c <- c(Embryo = "C: Average normalised embryo",
                       Gravid = "C: % gravid")

# ==============================================================================
# 2. Build the Plots
# ==============================================================================

# --- Panels A/B: QTL plots ---
s27a <- plotQTLStats(SNPset = EG4946_data, var = "Gprime", plotThreshold = TRUE, q = 0.01) +
  my_plot_theme +
  ggtitle("EG4946") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

s27b <- plotQTLStats(SNPset = EG4946_data, var = "deltaSNP", plotIntervals = TRUE) +
  my_plot_theme +
  ggtitle("EG4946") +
  theme(plot.title = element_text(size = 8, face = "plain", hjust = 0))

# --- Panel C: Embryo + Gravid + Genotype ---

## One naming scheme for both the genotype and phenotype files:
## "EG_NIL_A1_S9_L001.sorted.bam", "EG.RIL.A1", "EG RIL A1" -> "RIL A1";
## "EG4946_S19_L001.sorted.bam" -> "EG4946"; "jyIs8" / "JYIS8" -> "JYIS8"
clean_name <- function(x) {
  x <- sub("_S[0-9]+_L[0-9]+\\.sorted\\.bam$", "", x)   # old BAM-style names
  x <- gsub("NIL", "RIL", x)
  x <- trimws(gsub("[._-]+", " ", x))
  x <- sub("^\\S+ RIL ", "RIL ", x)                     # "EG RIL A1" -> "RIL A1"
  ifelse(toupper(x) == toupper(DONOR), DONOR, x)
}

## Genotype data
df_geno <- geno_raw %>%
  filter(chrom != "MtDNA", !sample %in% GENO_EXCLUDE) %>%
  mutate(gt = ifelse(gt == 1, DONOR, PARENT),
         sample = clean_name(sample))

## Row order (top to bottom): donor, parent, then RILs by letter then number (A1 ... B5)
rils   <- unique(df_geno$sample[grepl("^RIL ", df_geno$sample)])
ril_id <- sub("^RIL ", "", rils)
rils   <- rils[order(gsub("[0-9]", "", ril_id), as.numeric(gsub("[^0-9]", "", ril_id)))]

line_order <- rev(c(DONOR, PARENT, rils))   # reversed: first strain at the top

df_geno$sample <- factor(df_geno$sample, levels = line_order)
df_geno$index  <- as.numeric(df_geno$sample)

## Display-only relabel: "JYIS8" -> "N2"
display_label <- function(x) ifelse(x == DONOR, DONOR_LABEL, x)

## Reshape embryo + gravid to long format (individual replicates)
to_long <- function(df, measure_name) {
  df %>%
    pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(measure = measure_name)
}

combined_long <- bind_rows(to_long(EG4946_embryo, "Embryo"),
                           to_long(EG4946_gravid, "Gravid")) %>%
  mutate(name = clean_name(name),
         # "N2" -> "JYIS8"
         name = ifelse(name %in% names(PHENO_ALIASES), unname(PHENO_ALIASES[name]), name),
         # bare RIL IDs in the phenotype files (e.g. "A1") -> "RIL A1"
         name = ifelse(paste("RIL", name) %in% line_order, paste("RIL", name), name))

## Sanity check: every phenotype line must match a genotyped line
unmatched <- setdiff(unique(combined_long$name), line_order)
if (length(unmatched) > 0) {
  stop("These phenotype columns don't match any genotyped strain: ",
       paste(unmatched, collapse = ", "),
       "\nGenotyped strains: ", paste(rev(line_order), collapse = ", "),
       "\nAdd them to PHENO_ALIASES if they're another name for a genotyped strain.")
}
combined_long <- combined_long %>% mutate(name = factor(name, levels = line_order))

## Summarise mean/SD per line per measure (for SD bars + mean point)
summary_df <- combined_long %>%
  group_by(measure, name) %>%
  summarise(mean = mean(value), sd = sd(value), n = n(), .groups = "drop")

## --- Dunnett's test: every line vs the N2 control, run separately per measure ---
## NOTE: keep "JYIS8" (not "N2") in all logic below; it must match the
## strain name used in the data for relevel()/joins to work.
dunnett_result <- function(df, measure_name, control_strain) {
  d <- df %>% filter(measure == measure_name) %>% droplevels()
  d$name <- relevel(d$name, ref = control_strain)
  test <- DunnettTest(value ~ name, data = d, control = control_strain)[[1]]
  tibble(
    name    = sub(paste0("-", control_strain, "$"), "", rownames(test)),
    measure = measure_name,
    p_adj   = test[, "pval"]
  )
}

# "" (not "ns") so non-significant lines show no stars on the figure
stars_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

sig_donor <- bind_rows(
  dunnett_result(combined_long, "Embryo", DONOR),
  dunnett_result(combined_long, "Gravid", DONOR)
) %>% mutate(stars_donor = stars_from_p(p_adj)) %>% dplyr::select(measure, name, stars_donor)

summary_df <- summary_df %>%
  mutate(name = as.character(name)) %>%
  left_join(sig_donor, by = c("measure", "name")) %>%
  mutate(stars_donor = replace_na(stars_donor, ""),
         name = factor(name, levels = line_order))

## Phenotype panel (embryo or gravid); stars = significance vs N2
star_labels <- setNames(display_label(c(PARENT, DONOR)), c(PARENT, DONOR))   # genotype legend labels

pheno_panel <- function(measure_name, x_label, y_labels_right = FALSE) {
  pts  <- combined_long %>% filter(measure == measure_name)
  summ <- summary_df    %>% filter(measure == measure_name)
  
  p <- ggplot() +
    geom_jitter(data = pts, aes(x = value, y = name),
                height = 0.2, width = 0, size = 0.6, alpha = 0.5, color = "gray50") +
    geom_errorbarh(data = summ, aes(xmin = mean - sd, xmax = mean + sd, y = name),
                   height = 0, color = "gray20", linewidth = 0.35) +
    geom_point(data = summ, aes(x = mean, y = name), size = 1.5, color = "black") +
    geom_text(data = summ,
              aes(x = mean + sd, y = name, label = stars_donor),
              hjust = -0.2, vjust = 0.75, size = 2.3) +
    # limits pin the row order; add = 0.5 lines rows up with the genotype panel
    scale_y_discrete(limits = line_order, expand = expansion(add = 0.5),
                     position = if (y_labels_right) "right" else "left",
                     labels = display_label) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = x_label, y = NULL) +
    my_plot_theme +
    theme(panel.grid.minor = element_blank())
  
  if (!y_labels_right) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

embryo_p <- pheno_panel("Embryo", "Average normalised embryo")
gravid_p <- pheno_panel("Gravid", "% gravid", y_labels_right = TRUE)

## Genotype panels: one small plot per chromosome (no facet strips), each
## with its chromosome in the x-axis title. Without strips, every panel in
## row C has the same top and bottom, so the x-axes stay attached and the
## "C" tag sits beside the panels.
geno_fill <- scale_fill_manual(
  values = setNames(c("#0080FF", "#FF8000"), c(PARENT, DONOR)),
  limits = c(PARENT, DONOR),
  labels = star_labels,
  name   = "Genotype"
)

geno_chrom_plot <- function(ch, last = FALSE) {
  ggplot(df_geno %>% filter(chrom == ch)) +
    geom_rect(aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
    geno_fill +
    scale_x_continuous(labels = \(x) x / 1e6, expand = c(0, 0)) +
    scale_y_continuous(limits = c(1, length(line_order) + 1), expand = c(0, 0)) +
    labs(x = paste0("Chr ", ch, " position (Mb)"), y = NULL) +
    my_plot_theme +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin  = margin(t = 5, r = if (last) 5 else 1.5, b = 5, l = 1.5))
}

geno_plots <- map(seq_along(GENO_CHROM),
                  ~ geno_chrom_plot(GENO_CHROM[.x], last = .x == length(GENO_CHROM)))

# Width of each chromosome plot proportional to its length; together they
# take the same share of the row as the single genotype panel did (2 units)
chrom_len  <- map_dbl(GENO_CHROM, ~ max(df_geno$end[df_geno$chrom == .x]))
geno_width <- 2 * chrom_len / sum(chrom_len)

## Assemble panel C
# l = 14 on the embryo panel leaves room for the "C" tag on its left
embryo_p <- embryo_p + theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 14))
gravid_p <- gravid_p + theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 2))

panel_c_raw <- wrap_plots(c(list(embryo_p, gravid_p), geno_plots), nrow = 1) +
  plot_layout(widths = c(1, 1, geno_width), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

# wrap as one element so it gets a single tag ("C"), placed beside the
# top-left of the panels
panel_c <- wrap_elements(full = panel_c_raw) +
  theme(
    plot.margin = margin(t = 14, r = 5, b = 5, l = 5),
    plot.tag    = element_text(size = 10, face = "plain")
  )

# ==============================================================================
# 3. Patchwork Assembly
# ==============================================================================
layout_design <- "
  AB
  CC
"

supp_fig27 <- s27a + s27b + panel_c +
  plot_layout(design = layout_design, heights = c(1, 2)) +
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

supp_fig27   # view in RStudio

ggsave("supp_fig27.pdf",
       plot = supp_fig27,
       width = 180,
       height = 220*2/3,
       units = "mm",
       device = cairo_pdf)

# ==============================================================================
# 4. QTL Table Export
# ==============================================================================
EG4946_qtl <- getQTLTable(
  SNPset   = EG4946_data,
  method   = "Gprime",
  alpha    = 0.01,
  export   = TRUE,
  fileName = "EG4946_QTL_regions.csv"
)

# Named copies of the panel C data (used by the tables below)
EG4946_fitness     <- combined_long
EG4946_breakpoints <- df_geno
EG4946_summary     <- summary_df

# ==============================================================================
# 5. Panel C stats: descriptives, Dunnett vs N2, normality, Levene
# ==============================================================================

# Shapiro-Wilk that returns NA instead of erroring (n < 3, or all values identical)
safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

panel_c_stats <- function(fitness, donor, donor_label, line_order,
                          measures = c("Embryo", "Gravid")) {
  
  relabel <- function(x) ifelse(x == donor, donor_label, x)
  
  dunnett_full <- function(m, control, suffix) {
    d <- fitness %>% filter(measure == m) %>% droplevels()
    test <- DunnettTest(value ~ name, data = d, control = control)[[1]]
    tibble(measure = m,
           name    = sub(paste0("-", control, "$"), "", rownames(test)),
           diff    = test[, "diff"],
           lwr_ci  = test[, "lwr.ci"],
           upr_ci  = test[, "upr.ci"],
           p_adj   = test[, "pval"]) %>%
      mutate(signif = stars_from_p(p_adj),
             signif = ifelse(signif == "", "ns", signif)) %>%   # "ns" in the table (blank on the figure)
      rename_with(~ paste0(.x, "_vs_", suffix), c(diff, lwr_ci, upr_ci, p_adj, signif))
  }
  
  desc <- fitness %>%
    group_by(measure, name) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") %>%
    mutate(name = as.character(name))
  
  # Dunnett vs the N2 control only (same test as the stars in the figure)
  stats <- desc %>%
    left_join(map_dfr(measures, ~ dunnett_full(.x, donor, donor_label)),
              by = c("measure", "name")) %>%
    mutate(name = factor(relabel(name), levels = relabel(rev(line_order)))) %>%
    arrange(measure, name)
  
  normality <- map_dfr(measures, function(m) {
    d <- fitness %>% filter(measure == m) %>% droplevels()
    per_line <- d %>%
      group_by(name) %>%
      summarise(n         = n(),
                shapiro_W = unname(safe_shapiro(value)["W"]),
                shapiro_p = unname(safe_shapiro(value)["p"]),
                .groups = "drop") %>%
      mutate(name = relabel(as.character(name)))
    res <- safe_shapiro(residuals(aov(value ~ name, data = d)))
    bind_rows(per_line,
              tibble(name = "All lines (ANOVA residuals)", n = nrow(d),
                     shapiro_W = unname(res["W"]), shapiro_p = unname(res["p"]))) %>%
      mutate(measure = m, .before = 1)
  })
  
  levene <- map_dfr(measures, function(m) {
    d <- fitness %>% filter(measure == m) %>% droplevels()
    lev <- DescTools::LeveneTest(value ~ name, data = d)
    tibble(measure = m, F_value = lev[1, "F value"],
           df1 = lev[1, "Df"], df2 = lev[2, "Df"], levene_p = lev[1, "Pr(>F)"])
  })
  
  list(stats = stats, normality = normality, levene = levene)
}

# Relabel "Embryo"/"Gravid" -> "C: Average normalised embryo"/"C: % gravid"
# and call the column "panel", matching data_s27c_raw below.
# (panel_c_stats() rebuilds the tables from scratch, so rerunning is safe.)
relabel_measure <- function(df) {
  df %>%
    mutate(measure = unname(panel_labels_s27c[measure])) %>%
    rename(panel = measure)
}

# ---- Supplementary Figure 27C ----
supp27c <- panel_c_stats(EG4946_fitness, donor = DONOR,
                         donor_label = DONOR_LABEL, line_order = line_order)
supp27c <- lapply(supp27c, relabel_measure)

supp27c$stats
supp27c$normality
supp27c$levene

# ==============================================================================
# 6. Supplementary data
# ==============================================================================
# Panels A-B: QTL-seq SNP table and significant QTL regions
data_s27ab_snps <- EG4946_data
data_s27ab_qtl  <- EG4946_qtl

# Panel C: every replicate value (same layout as fig3e_raw / fig4e_raw / fig5de$raw)
data_s27c_raw <- EG4946_fitness %>%
  arrange(factor(measure, levels = c("Embryo", "Gravid")),
          factor(name, levels = rev(line_order))) %>%   # same order as the figure, top to bottom
  transmute(panel  = unname(panel_labels_s27c[measure]),
            strain = display_label(as.character(name)),  # JYIS8 shown as "N2", as on the figure
            value)

# Panel C: RIL genotype blocks (donor shown as "N2"; plotting column dropped)
data_s27c_breakpoints <- EG4946_breakpoints %>%
  mutate(sample = display_label(as.character(sample)),
         gt     = display_label(gt)) %>%
  dplyr::select(-index) %>%
  rename(strain = sample, genotype = gt)

data_s27c_raw
data_s27c_breakpoints

# Excel's row limit is 1,048,576; check the SNP table
nrow(data_s27ab_snps)

# ==============================================================================
# Tables for the supplementary workbooks
# ==============================================================================
#   Data:  data_s27ab_qtl, data_s27ab_snps (S27A-B),
#          data_s27c_raw, data_s27c_breakpoints (S27C)
#   Stats: supp27c$stats, supp27c$normality, supp27c$levene