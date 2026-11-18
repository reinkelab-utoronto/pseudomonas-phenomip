# ==============================================================================
# Bulk Segregant Analysis: all six populations
#
# Reads <POP>_allele_depths_columns.tsv for each population (made with
# bcftools query -H from <POP>_samples.vcf.gz) and, for each population, makes:
#
#   Part 1: per-replicate allele frequency, log10 odds ratio and delta AF plots
#           (LOESS + sliding window), infected replicate vs matching uninfected
#   Part 2: QTLseqr G' and QTL-seq (deltaSNP) analysis for each comparison:
#           <POP>_<comparison>.tsv, <POP>_<comparison>_results.csv, and the
#           G' and deltaSNP PDFs
#
# Run from the project folder (bac_repo_test). Inputs are read from
# BSA_results/<POP>/, and the QTLseqr files used by the figure scripts are
# written to the project folder itself:
#   <POP>_<comparison>.tsv and <POP>_<comparison>_results.csv
# (read by main_figure3-5.R, supp_fig20, supp_fig23, supp_fig27).
# The exploratory Part 1 plots and the QTLseqr G'/deltaSNP PDFs go to
# <PLOT_DIR>/<POP>/, with a summary of site counts in <PLOT_DIR>/site_counts.tsv.
#
# To change which populations run, edit RUN_POPS. To change comparisons,
# labels or colours, edit the design definitions in section 1.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(slider)
  library(QTLseqr)
})

# ── 0. Settings ───────────────────────────────────────────────────────────────
INPUT_DIR   <- "BSA_results"       # contains <POP>/<POP>_allele_depths_columns.tsv
RESULTS_DIR <- "."                 # .tsv / _results.csv read by the figure scripts
PLOT_DIR    <- "BSA_exploratory"   # Part 1 plots, G'/deltaSNP PDFs, site counts
RUN_POPS   <- c("JU2610", "ED3005", "JU2522", "WN2050", "JU2250", "EG4946")

# Optionally run only some populations from the command line, e.g.
#   Rscript BSA_all_populations.R JU2522 WN2050
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) RUN_POPS <- args

CHROMS        <- c("I", "II", "III", "IV", "V", "X")   # MtDNA excluded
MIN_SITE_DP   <- 10      # Part 1: minimum depth in every sample
SLIDE_BEFORE  <- 24      # Part 1: sliding window of 50 SNPs
SLIDE_AFTER   <- 25

# Part 2 (QTLseqr), same for every population
QTL_REF_AF        <- 0.20
QTL_MIN_TOTAL_DP  <- 0
QTL_MAX_TOTAL_DP  <- 400
QTL_MIN_SAMPLE_DP <- 40     # per bulk (replicates summed)
QTL_WINDOW        <- 1e6
QTL_POPSTRUC      <- "F9"
QTL_BULKSIZE      <- c(100, 100)
QTL_REPS          <- 10000
QTL_INTERVALS     <- c(95, 99)
QTL_Q             <- 0.01


# ── 1. Population designs ─────────────────────────────────────────────────────
# bulks:       bulk name -> sample columns (replicate 1, replicate 2)
# labels:      display name of each bulk in the odds ratio / delta AF legends
# colors:      one colour per replicate of each bulk
# comparisons: name -> c(Low bulk, High bulk) for QTLseqr
# part1:       comparisons plotted per replicate in Part 1 (vs uninfected)

design_myb114 <- list(
  bulks  = list(uni  = c("uni1",  "uni2"),
                low  = c("low1",  "low2"),
                full = c("full1", "full2")),
  labels = c(uni = "Uni", low = "Low", full = "Full"),
  colors = list(full = c("darkgreen", "green3"),
                low  = c("orange", "gold"),
                uni  = c("dodgerblue", "blue")),
  comparisons = list(uni_vs_full = c("uni", "full"),
                     uni_vs_low  = c("uni", "low"),
                     low_vs_full = c("low", "full")),
  part1 = c("uni_vs_full", "uni_vs_low")
)

design_pa14 <- list(
  bulks  = list(G6_uni  = c("G6_uni1", "G6_uni2"),
                G6_10pc = c("G6_10percentPA1", "G6_10percentPA2"),
                G8_uni  = c("G8_uni1", "G8_uni2"),
                G8_50pc = c("G8_50percentPA1", "G8_50percentPA2")),
  labels = c(G6_uni = "G6 Uni", G6_10pc = "G6 10%",
             G8_uni = "G8 Uni", G8_50pc = "G8 50%"),
  colors = list(G6_10pc = c("orange", "gold"),
                G8_50pc = c("darkgreen", "green3"),
                G6_uni  = c("dodgerblue", "blue"),
                G8_uni  = c("purple", "orchid")),
  comparisons = list(G6_uni_vs_low  = c("G6_uni", "G6_10pc"),
                     G8_uni_vs_half = c("G8_uni", "G8_50pc")),
  part1 = c("G6_uni_vs_low", "G8_uni_vs_half")
)

design_myb395 <- list(
  bulks  = list(uni = c("uni1", "uni2"),
                m   = c("m1", "m2")),
  labels = c(uni = "Uni", m = "M"),
  colors = list(m   = c("darkgreen", "green3"),
                uni = c("dodgerblue", "blue")),
  comparisons = list(uni_vs_m = c("uni", "m")),
  part1 = c("uni_vs_m")
)

DESIGNS <- list(JU2610 = design_myb114, ED3005 = design_myb114,
                JU2522 = design_pa14,   WN2050 = design_pa14,
                JU2250 = design_myb395, EG4946 = design_myb395)


# ── 2. Plot styling (same as the original script) ─────────────────────────────
part1_theme <- theme_bw() +
  theme(legend.position = "bottom",
        plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 16),
        strip.text   = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 16),
        legend.key.size = unit(0.8, "cm"))

my_plot_theme <- theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 20),
        axis.text    = element_text(size = 18),
        strip.text   = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18),
        legend.key.size = unit(1, "cm"))

fig_width  <- 6.3   # full printable A4 width
fig_height <- 4.5

save_part1 <- function(p, file) ggsave(file, p, width = 29.7, height = 21, units = "cm")


# ── 3. Read and clean one population ──────────────────────────────────────────
read_population <- function(pop) {
  file <- file.path(INPUT_DIR, pop, paste0(pop, "_allele_depths_columns.tsv"))
  if (!file.exists(file)) stop("Not found: ", file)
  d <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character")
  names(d) <- sub("^#", "", names(d))
  d$POS <- as.integer(d$POS)

  # Biallelic SNPs on the nuclear chromosomes only
  # (multi-allelic ALTs such as "A,G" have nchar > 1 and are removed here)
  d <- d[nchar(d$REF) == 1 & nchar(d$ALT) == 1 & d$CHROM %in% CHROMS, ]

  # Split every "ref,alt" sample column into <sample>_REF and <sample>_ALT
  samples <- setdiff(names(d), c("CHROM", "POS", "REF", "ALT"))
  for (s in samples) {
    d[[paste0(s, "_REF")]] <- as.integer(sub(",.*$", "", d[[s]]))
    d[[paste0(s, "_ALT")]] <- as.integer(sub("^.*,", "", d[[s]]))
    d[[s]] <- NULL
  }
  list(data = d, samples = samples, file = file)
}

check_design <- function(pop, samples, design) {
  expected <- unlist(design$bulks)
  missing  <- setdiff(expected, samples)
  if (length(missing) > 0) stop(pop, ": sample columns not found: ",
                                paste(missing, collapse = ", "))
}


# ── 4. Part 1: per-replicate AF, odds ratio and delta AF ──────────────────────
run_part1 <- function(pop, d, design, outdir) {
  samples <- unlist(design$bulks, use.names = FALSE)

  # Depth filter: every sample must have at least MIN_SITE_DP reads
  keep <- Reduce(`&`, lapply(samples, function(s)
    (d[[paste0(s, "_REF")]] + d[[paste0(s, "_ALT")]]) >= MIN_SITE_DP))
  d <- d[keep, ]

  # N2 (reference) allele frequency per sample
  for (s in samples) d[[paste0(s, "_AF")]] <- d[[paste0(s, "_REF")]] /
                                              (d[[paste0(s, "_REF")]] + d[[paste0(s, "_ALT")]])
  d$pos_Mb <- d$POS / 1e6
  d$chrom  <- factor(d$CHROM, levels = CHROMS)

  # Sample -> bulk / replicate lookup
  lookup <- bind_rows(lapply(names(design$bulks), function(b)
    data.frame(sample = design$bulks[[b]], bulk = b, rep = seq_along(design$bulks[[b]]))))
  lookup$label <- paste0(lookup$bulk, "_rep", lookup$rep)
  lookup$color <- mapply(function(b, r) design$colors[[b]][r], lookup$bulk, lookup$rep)
  af_colors <- setNames(lookup$color, lookup$label)

  # ── Allele frequency
  df_long <- d %>%
    select(chrom, pos_Mb, all_of(paste0(samples, "_AF"))) %>%
    pivot_longer(-c(chrom, pos_Mb), names_to = "sample", values_to = "allele_freq") %>%
    mutate(sample = sub("_AF$", "", sample)) %>%
    left_join(lookup[, c("sample", "label")], by = "sample")

  p <- ggplot(df_long, aes(pos_Mb, allele_freq, color = label)) +
    geom_smooth(se = FALSE) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(values = af_colors) +
    labs(x = "Genomic Position (Mb)", y = "N2 Allele Frequency", color = "Sample",
         title = "LOESS Smoothed Allele Frequency") +
    part1_theme
  save_part1(p, file.path(outdir, paste0(pop, "_allele_freq_loess.pdf")))

  df_long_smoothed <- df_long %>%
    arrange(chrom, pos_Mb) %>%
    group_by(chrom, label) %>%
    mutate(allele_freq_smooth = slide_dbl(allele_freq, mean, .before = SLIDE_BEFORE,
                                          .after = SLIDE_AFTER, .complete = TRUE)) %>%
    ungroup()

  p <- ggplot(df_long_smoothed, aes(pos_Mb, allele_freq_smooth, color = label)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_color_manual(values = af_colors) +
    labs(x = "Genomic Position (Mb)", y = "N2 Allele Frequency", color = "Sample",
         title = "Sliding Window Smoothed Allele Frequency") +
    part1_theme
  save_part1(p, file.path(outdir, paste0(pop, "_allele_freq_sliding_window.pdf")))
  rm(df_long, df_long_smoothed, p); invisible(gc())   # free memory before the next plots

  # ── Odds ratio and delta AF: infected replicate r vs uninfected replicate r
  clamp <- function(p) ifelse(p == 1, 1 - 1e-6, ifelse(p == 0, 1e-6, p))
  pairs <- list(); pair_colors <- c()
  for (comp in design$part1) {
    low  <- design$comparisons[[comp]][1]
    high <- design$comparisons[[comp]][2]
    for (r in seq_along(design$bulks[[high]])) {
      group <- paste(design$labels[[high]], paste0("Rep", r))
      pairs[[group]] <- c(design$bulks[[high]][r], design$bulks[[low]][r])
      pair_colors[group] <- design$colors[[high]][r]
    }
  }

  or_long <- bind_rows(lapply(names(pairs), function(g) {
    ph <- clamp(d[[paste0(pairs[[g]][1], "_AF")]])
    pl <- clamp(d[[paste0(pairs[[g]][2], "_AF")]])
    data.frame(chrom = d$chrom, pos_Mb = d$pos_Mb, group = g,
               log10_OR = log10((ph / (1 - ph)) / (pl / (1 - pl))),
               delta_AF = d[[paste0(pairs[[g]][1], "_AF")]] - d[[paste0(pairs[[g]][2], "_AF")]])
  }))

  or_smoothed <- or_long %>%
    arrange(chrom, pos_Mb) %>%
    group_by(chrom, group) %>%
    mutate(log10_OR_smooth = slide_dbl(log10_OR, mean, .before = SLIDE_BEFORE,
                                       .after = SLIDE_AFTER, .complete = TRUE),
           delta_AF_smooth = slide_dbl(delta_AF, mean, .before = SLIDE_BEFORE,
                                       .after = SLIDE_AFTER, .complete = TRUE)) %>%
    ungroup()

  pair_plot <- function(data, y, smooth, ylab, title) {
    p <- ggplot(data, aes(pos_Mb, .data[[y]], color = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
    p <- if (smooth) p + geom_smooth(se = FALSE) else p + geom_line(linewidth = 0.7)
    p + facet_wrap(~chrom, scales = "free_x") +
      scale_color_manual(values = pair_colors) +
      labs(x = "Genomic Position (Mb)", y = ylab, color = "Replicate", title = title) +
      part1_theme
  }

  save_part1(pair_plot(or_long, "log10_OR", TRUE, "log10 Odds Ratio (Infected / Control)",
                       "LOESS smoothed log10 Odds Ratio of N2 Allele Frequency"),
             file.path(outdir, paste0(pop, "_odds_ratio_loess.pdf")))
  save_part1(pair_plot(or_smoothed, "log10_OR_smooth", FALSE, "log10 Odds Ratio",
                       "Sliding Window Smoothed log10 Odds Ratio"),
             file.path(outdir, paste0(pop, "_odds_ratio_sliding_window.pdf")))
  invisible(gc())
  save_part1(pair_plot(or_long, "delta_AF", TRUE,
                       "Difference in Allele Frequency (Infected - Control)",
                       "LOESS Smoothed Difference in N2 Allele Frequency"),
             file.path(outdir, paste0(pop, "_delta_AF_loess.pdf")))
  save_part1(pair_plot(or_smoothed, "delta_AF_smooth", FALSE, "Difference in Allele Frequency",
                       "Sliding Window Smoothed Difference in Allele Frequency"),
             file.path(outdir, paste0(pop, "_delta_AF_sliding_window.pdf")))

  nrow(d)
}


# ── 5. Part 2: QTLseqr for each comparison ────────────────────────────────────
run_part2 <- function(pop, d, design, plotdir) {
  sum_bulk <- function(bulk, allele)
    Reduce(`+`, lapply(design$bulks[[bulk]], function(s) d[[paste0(s, "_", allele)]]))

  n_after <- c()
  for (comp in names(design$comparisons)) {
    message("  ", pop, " ", comp)
    low  <- design$comparisons[[comp]][1]
    high <- design$comparisons[[comp]][2]

    comp_df <- data.frame(CHROM = d$CHROM, POS = d$POS, REF = d$REF, ALT = d$ALT,
                          AD_REF.Low  = sum_bulk(low,  "REF"),
                          AD_ALT.Low  = sum_bulk(low,  "ALT"),
                          AD_REF.High = sum_bulk(high, "REF"),
                          AD_ALT.High = sum_bulk(high, "ALT"))
    tsv_file <- file.path(RESULTS_DIR, paste0(pop, "_", comp, ".tsv"))
    write.table(comp_df, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

    # One failing comparison should not stop the others
    n_after[comp] <- tryCatch({
      qtl <- importFromTable(file = tsv_file, highBulk = "High", lowBulk = "Low",
                             chromList = CHROMS, sep = "\t")
      qtl <- filterSNPs(SNPset = qtl, refAlleleFreq = QTL_REF_AF,
                        minTotalDepth = QTL_MIN_TOTAL_DP, maxTotalDepth = QTL_MAX_TOTAL_DP,
                        minSampleDepth = QTL_MIN_SAMPLE_DP)
      qtl <- runGprimeAnalysis(SNPset = qtl, windowSize = QTL_WINDOW,
                               outlierFilter = "deltaSNP")
      qtl <- runQTLseqAnalysis(SNPset = qtl, windowSize = QTL_WINDOW,
                               popStruc = QTL_POPSTRUC, bulkSize = QTL_BULKSIZE,
                               replications = QTL_REPS, intervals = QTL_INTERVALS)
      write.csv(qtl, file.path(RESULTS_DIR, paste0(pop, "_", comp, "_results.csv")),
                row.names = FALSE)

      pdf(file.path(plotdir, paste0(pop, "_", comp, "_Gprime.pdf")), fig_width, fig_height)
      print(plotQTLStats(SNPset = qtl, var = "Gprime", plotThreshold = TRUE, q = QTL_Q) +
              my_plot_theme)
      dev.off()

      pdf(file.path(plotdir, paste0(pop, "_", comp, "_deltaSNP.pdf")), fig_width, fig_height)
      print(plotQTLStats(SNPset = qtl, var = "deltaSNP", plotIntervals = TRUE) +
              my_plot_theme)
      dev.off()

      nrow(qtl)
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("  FAILED ", pop, " ", comp, ": ", conditionMessage(e))
      NA
    })
  }
  n_after
}


# ── 6. Run all populations ────────────────────────────────────────────────────
dir.create(PLOT_DIR, showWarnings = FALSE)
counts <- list()

for (pop in RUN_POPS) {
  message("\n=== ", pop, " ===")
  design <- DESIGNS[[pop]]
  outdir <- file.path(PLOT_DIR, pop)
  dir.create(outdir, showWarnings = FALSE)

  inp <- read_population(pop)
  check_design(pop, inp$samples, design)
  message("  ", nrow(inp$data), " biallelic SNPs on ", paste(CHROMS, collapse = ","),
          " (", basename(inp$file), ")")

  n_part1 <- run_part1(pop, inp$data, design, outdir); invisible(gc())
  n_part2 <- run_part2(pop, inp$data, design, outdir)

  counts[[pop]] <- data.frame(population = pop,
                              step = c("biallelic_SNPs", "part1_all_samples_DP10",
                                       paste0("QTLseqr_", names(n_part2))),
                              sites = c(nrow(inp$data), n_part1, unname(n_part2)))
}

site_counts <- bind_rows(counts)
# Append, so separate runs for different populations all end up in one table
counts_file <- file.path(PLOT_DIR, "site_counts.tsv")
if (file.exists(counts_file)) {
  old <- read.delim(counts_file, stringsAsFactors = FALSE)
  site_counts <- bind_rows(old[!old$population %in% RUN_POPS, ], site_counts)
}
write.table(site_counts, counts_file, sep = "\t", quote = FALSE, row.names = FALSE)
print(site_counts)
message("\nAll populations finished. Results in ", normalizePath(RESULTS_DIR),
        ", plots in ", PLOT_DIR, "/")
