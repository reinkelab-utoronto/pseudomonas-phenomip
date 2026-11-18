# ==============================================================================
# Write the four supplementary workbooks, each with a "Contents" sheet first
# ==============================================================================
# Run AFTER sourcing all main and supplementary figure scripts.
#
# Each workbook is described by a "manifest": one row per sheet, with
#   id          - data name shown in the Contents sheet (e.g. "Data S13")
#   figure      - figure/panel it belongs to (e.g. "Figure 3E")
#   description - short description for the Contents sheet
#   object      - the R object (or expression) written to that sheet
# The sheet name is built as "<id> (<figure>)", so the Contents sheet and the
# sheet tabs can never disagree. Edit a row here and both update on rerun.
# ==============================================================================
library(openxlsx)
library(dplyr)
library(tibble)

# Matrices (e.g. correlation matrices) keep their row names as a first column
as_sheet <- function(x) {
  if (is.matrix(x)) {
    x <- tibble::rownames_to_column(as.data.frame(x), "condition")
  }
  as.data.frame(x)
}

write_workbook_with_toc <- function(manifest, file, title, env = globalenv()) {
  
  # 1. Evaluate every object first, and stop with a clear list of any that fail
  values <- lapply(manifest$object, function(expr) {
    tryCatch(eval(parse(text = expr), envir = env), error = function(e) NULL)
  })
  failed <- manifest$object[vapply(values, is.null, logical(1))]
  if (length(failed) > 0) {
    stop("These objects don't exist or failed to evaluate (", file, "):\n  ",
         paste(failed, collapse = "\n  "),
         "\nSource the script that creates them, or fix the manifest.")
  }
  values <- lapply(values, as_sheet)
  
  # 2. Sheet names, checked against Excel's rules
  sheets <- paste0(manifest$id, " (", manifest$figure, ")")
  too_long <- sheets[nchar(sheets) > 31]
  if (length(too_long) > 0) stop("Sheet names over Excel's 31-character limit:\n  ",
                                 paste(too_long, collapse = "\n  "))
  if (anyDuplicated(sheets)) stop("Duplicate sheet names: ",
                                  paste(unique(sheets[duplicated(sheets)]), collapse = ", "))
  too_big <- sheets[vapply(values, nrow, integer(1)) >= 1048576]
  if (length(too_big) > 0) stop("Over Excel's row limit:\n  ", paste(too_big, collapse = "\n  "))
  
  # 3. Contents sheet
  wb <- createWorkbook()
  addWorksheet(wb, "Contents")
  
  toc <- tibble(
    `Data name`   = manifest$id,
    Figure        = manifest$figure,
    Description   = manifest$description,
    Rows          = vapply(values, nrow, integer(1)),
    Columns       = vapply(values, ncol, integer(1))
  )
  
  writeData(wb, "Contents", title, startRow = 1, startCol = 1)
  writeData(wb, "Contents",
            "Click a data name to go to its sheet. Each sheet is named '<Data name> (<Figure>)'.",
            startRow = 2, startCol = 1)
  addStyle(wb, "Contents", createStyle(textDecoration = "bold", fontSize = 13),
           rows = 1, cols = 1)
  
  header_row <- 4
  writeData(wb, "Contents", toc, startRow = header_row,
            headerStyle = createStyle(textDecoration = "bold",
                                      border = "bottom", fgFill = "#DDEBF7"))
  
  # Data names become links to their sheets
  for (i in seq_along(sheets)) {
    writeFormula(wb, "Contents", startCol = 1, startRow = header_row + i,
                 x = makeHyperlinkString(sheet = sheets[i], row = 1, col = 1,
                                         text = manifest$id[i]))
  }
  
  setColWidths(wb, "Contents", cols = 1:5, widths = c(18, 18, 70, 10, 10))
  addStyle(wb, "Contents", createStyle(wrapText = TRUE, valign = "top"),
           rows = header_row + seq_along(sheets), cols = 3, gridExpand = TRUE)
  freezePane(wb, "Contents", firstActiveRow = header_row + 1)
  
  # 4. Data sheets
  for (i in seq_along(sheets)) {
    addWorksheet(wb, sheets[i])
    writeData(wb, sheets[i], values[[i]],
              headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheets[i], firstRow = TRUE)
  }
  
  saveWorkbook(wb, file, overwrite = TRUE)
  message("Wrote ", file, " (", length(sheets), " data sheets + Contents)")
  invisible(toc)
}

# ==============================================================================
# 1. main_data.xlsx - data for the main figures
# ==============================================================================
main_data <- tribble(
  ~id,        ~figure,         ~description,                                                                  ~object,
  "Data S1",  "Table S1",      "All strains used in experiment, their corresponding pools, subpools and CaeNR sets",'read.csv("Table S1.csv")',
  "Data S2",  "Figure 1A",     "Sampling coordinates and pool of each wild isolate (world map)",             "coords",
  "Data S3",  "Figure 1B",     "Bacterial species tree (Newick format)",                                     "tree_df",
  "Data S4",  "Figure 1D",     "Replicate 1 and 2 z-scores per pool and condition",                          "total_correlation_df",
  "Data S5",  "Figure 1E",     "Mean z-scores from phenoMIP and mini-phenoMIP per strain and condition",     "combined_df",       
  "Data S6",  "Figure 2A",     "Hierarchical clustering dendrogram of strains (Newick format)",              "tree_2a_df",         
  "Data S7",  "Figure 2A",     "Average z-scores per strain and condition shown in the heatmap",             "heat_2a_df",
  "Data S8",  "Figure 2B",     "Average strain abundance per generation and condition of 5 candidate wild isolates", "four_strain_averaged",
  "Data S9",  "Figure 3A-B",   "QTL-seq significant regions (G', q = 0.01), JU2610",                         "ju2610_qtl",
  "Data S10", "Figure 3A-B",   "QTL-seq SNP-level statistics (G', delta SNP-index), JU2610",                 "ju2610_data",
  "Data S11", "Figure 3C-D",   "QTL-seq significant regions (G', q = 0.01), ED3005",                         "ed3005_qtl",
  "Data S12", "Figure 3C-D",   "QTL-seq SNP-level statistics (G', delta SNP-index), ED3005",                 "ed3005_data",
  "Data S13", "Figure 3E",     "Embryo counts and % gravid per replicate, JU2610 RILs, JU2610 and N2",       "fig3e_raw",
  "Data S14", "Figure 3E",     "Genotype blocks (breakpoints) of JU2610 RILs and parents",                   "JU2610_breakpoints",
  "Data S15", "Figure 4A-B",   "QTL-seq significant regions (G', q = 0.01), JU2522",                         "ju2522_qtl",
  "Data S16", "Figure 4A-B",   "QTL-seq SNP-level statistics (G', delta SNP-index), JU2522",                 "ju2522_data",
  "Data S17", "Figure 4C-D",   "QTL-seq SNP-level statistics (G', delta SNP-index), WN2050",                 "wn2050_data",
  "Data S18", "Figure 4E",     "Embryo counts and % gravid per replicate, JU2522 RILs, JU2522 and N2",       "fig4e_raw",
  "Data S19", "Figure 4E",     "Genotype blocks (breakpoints) of JU2522 RILs and parents",                   "df_geno_JU2522",     # CHECK: was labelled Figure 3E
  "Data S20", "Figure 5A-B",   "QTL-seq significant regions (G', q = 0.01), JU2250",                         "ju2250_qtl",
  "Data S21", "Figure 5A-B",   "QTL-seq SNP-level statistics (G', delta SNP-index), JU2250",                 "ju2250_data",
  "Data S22", "Figure 5C",     "Embryo counts and % gravid per replicate, JU2250 RILs, JU2250 and jyIs8",    "JU2250_fitness",
  "Data S23", "Figure 5C",     "Genotype blocks (breakpoints) of JU2250 RILs and parents",                   "JU2250_breakpoints",
  "Data S24", "Figure 5D-E",   "Embryo counts and % gravid per replicate, candidate wild isolates on MYb395","fig5de$raw",
  "Data S25", "Figure 6",      "% body area with FISH signal per replicate (mean of 20 worms)",              "fig6_raw"
)

# ==============================================================================
# 2. main_stats_data.xlsx - statistics for the main figures
# ==============================================================================
main_stats <- tribble(
  ~id,              ~figure,        ~description,                                                                       ~object,
  "Stats Data S1",  "Figure 1D",    "Correlation between replicate z-scores (R², Spearman rho) per method",            "r2_data_overall",
  "Stats Data S2",  "Figure 1D",    "Shapiro-Wilk normality tests (replicate z-scores and regression residuals)",      "normality_df_1d",
  "Stats Data S3",  "Figure 1E",    "Correlation between phenoMIP and mini-phenoMIP (R, R², rho, p-values)",           "r2_data_1e",
  "Stats Data S4",  "Figure 1E",    "Shapiro-Wilk normality tests (mean z-scores and regression residuals)",           "normality_df_1e",
  "Stats Data S5",  "Figure 2A",    "Shapiro-Wilk normality tests of average z-scores",                                "normality_2a",
  "Stats Data S6",  "Figure 3E",    "Descriptive statistics and Dunnett's test vs N2 (embryo, % gravid)",              "stats_3e",
  "Stats Data S7",  "Figure 3E",    "Shapiro-Wilk normality tests per line and on ANOVA residuals",                    "normality_3e",
  "Stats Data S8",  "Figure 4E",    "Descriptive statistics and Dunnett's test vs N2 (embryo, % gravid)",              "fig4e$stats",
  "Stats Data S9",  "Figure 4E",    "Shapiro-Wilk normality tests per line and on ANOVA residuals",                    "fig4e$normality",
  "Stats Data S10", "Figure 5C",    "Descriptive statistics and Dunnett's test vs jyIs8 (embryo, % gravid)",              "fig5c$stats",
  "Stats Data S11", "Figure 5C",    "Shapiro-Wilk normality tests per line and on ANOVA residuals",                    "fig5c$normality",
  "Stats Data S12", "Figure 5D-E",  "Descriptive statistics and Dunnett's test vs N2 (embryo, % gravid)",              "fig5de$stats",
  "Stats Data S13", "Figure 5D-E",  "Shapiro-Wilk normality tests per strain and on ANOVA residuals",                  "fig5de$normality",
  "Stats Data S14", "Figure 6",     "Descriptive statistics and Dunnett's test vs N2 (% FISH-positive area)",          "stats_fig6",
  "Stats Data S15", "Figure 6",     "Shapiro-Wilk normality tests per strain and on ANOVA residuals",                  "normality_fig6"
)

# ==============================================================================
# 3. supp_data.xlsx - data for the supplementary figures
# ==============================================================================
supp_data <- tribble(
  ~id,        ~figure,          ~description,                                                                      ~object,
  "Data S26", "Figure S1",      "Percentage of males per subpool during optimisation",                             "percent_males",              
  "Data S27", "Figure S2",      "MIP sequencing results per probe and sample",                                     "all_samples",
  "Data S28", "Figure S2",      "Pairwise Spearman correlations between MIP probes within each strain",            "all_cor_long",
  "Data S29", "Figure S3A",     "Replicate 1 and 2 z-scores per pathogenic condition",                             "total_correlation_df_forplotting",
  "Data S30", "Figure S3B",     "Control strain metrics per pool and pathogenic condition",                        "filtered_df_noMO",
  "Data S31", "Figure S4",      "Mini-phenoMIP replicate 1 and 2 z-scores per condition",                          "supp4_data",
  "Data S32", "Figure S5",      "PhenoMIP vs mini-phenoMIP mean z-scores per condition",                           "supp5_data",
  "Data S33", "Figure S6",      "Phylogenetic tree of strains (Newick format)",                                    "worm_species_newick",
  "Data S34", "Figure S6",      "Average z-scores of all strains shown in the heatmap",                            "combined_heat_df_allstrains",
  "Data S35", "Figure S7A",     "PCA loadings with strain metadata",                                               "pca_loadings_metadata",
  "Data S36", "Figure S7B",     "Broad-sense heritability (H²) per condition",                                     "h2_df_supp",
  "Data S37", "Figure S8A",     "Strain clustering dendrogram, AUC (Newick; Euclidean, complete linkage)",        "tree_s8a_strain_df",
  "Data S38", "Figure S8A",     "Condition clustering dendrogram, AUC (Newick; Euclidean, complete linkage)",     "tree_s8a_condition_df",
  "Data S39", "Figure S8B",     "Strain clustering dendrogram, regression (Newick; Euclidean, complete linkage)", "tree_s8b_strain_df",
  "Data S40", "Figure S8B",     "Condition clustering dendrogram, regression (Newick; Euclidean, complete linkage)", "tree_s8b_condition_df",
  "Data S41", "Figure S9",      "Average z-scores with sampling location of each strain",                          "total_location_df_clean",
  "Data S42", "Figure S10",     "Z-scores and genome divergence metrics per strain, PA14",                         "data_supp10",
  "Data S43", "Figure S11",     "Z-scores and genome divergence metrics per strain, PV",                           "data_supp11",
  "Data S44", "Figure S12",     "Z-scores and genome divergence metrics per strain, MYb114",                       "data_supp12",
  "Data S45", "Figure S13",     "Z-scores and genome divergence metrics per strain, MYb117",                       "data_supp13",
  "Data S46", "Figure S14",     "Z-scores and genome divergence metrics per strain, MYb395",                       "data_supp14",
  "Data S47", "Figure S15",     "Isotype phylogeny (CaeNDR 20231213), trimmed to strains shown (Newick)",          "tree_s15_df",
  "Data S48", "Figure S15",     "Average modified z-scores for all bacteria and microsporidia conditions",         "bac_mic_heatmap",
  "Data S49", "Figure S16A",    "Strain clustering dendrogram, AUC (Newick; Euclidean, complete linkage)",        "tree_s16a_strain_df",
  "Data S50", "Figure S16A",    "Condition clustering dendrogram, AUC (Newick; Euclidean, complete linkage)",     "tree_s16a_condition_df",
  "Data S51", "Figure S16B",    "Strain clustering dendrogram, regression (Newick; Euclidean, complete linkage)", "tree_s16b_strain_df",
  "Data S52", "Figure S16B",    "Condition clustering dendrogram, regression (Newick; Euclidean, complete linkage)", "tree_s16b_condition_df",
  "Data S53", "Figure S17A",    "Spearman correlation matrix between conditions, regression",                      "cor_matrix_slope",
  "Data S54", "Figure S17A",    "BH-adjusted p-values for the regression correlation matrix",                      "p_matrix_slope",
  "Data S55", "Figure S17B",    "Spearman correlation matrix between conditions, AUC",                             "cor_matrix_AUC",
  "Data S56", "Figure S17B",    "BH-adjusted p-values for the AUC correlation matrix",                             "p_matrix_AUC",
  "Data S57", "Figure S18",     "QTL intervals, peak markers and heritability for each trait",                     "qtl_s18",
  "Data S58", "Figure S18",     "Strains carrying the REF and ALT allele at each peak marker",                     "pxg_s18",
  "Data S59", "Figure S20A-B",  "QTL-seq significant regions, low dose, JU2610",                                   "ju2610_low_qtl",
  "Data S60", "Figure S20A-B",  "QTL-seq SNP-level statistics, low dose, JU2610",                                  "ju2610_low_data",
  "Data S61", "Figure S20C-D",  "QTL-seq significant regions, low dose, ED3005",                                   "ed3005_low_qtl",
  "Data S62", "Figure S20C-D",  "QTL-seq SNP-level statistics, low dose, ED3005",                                  "ed3005_low_data",
  "Data S63", "Figure S21A-B",  "BSA allele depths (REF,ALT reads) per sample, JU2610",                            "raw_s21_ju2610",
  "Data S64", "Figure S21C-D",  "BSA allele depths (REF,ALT reads) per sample, ED3005",                            "raw_s21_ed3005",
  "Data S65", "Figure S22",     "Embryo counts and % gravid per replicate of candidate wild isolates with variants in F11D5.5 gene", "data_s22",
  "Data S66", "Figure S23A-B",  "QTL-seq significant regions, low dose, JU2522",                                   "ju2522_low_qtl",
  "Data S67", "Figure S23A-B",  "QTL-seq SNP-level statistics, low dose, JU2522",                                  "ju2522_low_data",
  "Data S68", "Figure S23C-D",  "QTL-seq SNP-level statistics, low dose, WN2050",                                  "wn2050_low_data",
  "Data S69", "Figure S24A-B",  "BSA allele depths (REF,ALT reads) per sample, JU2522",                            "raw_s24_ju2522",
  "Data S70", "Figure S24C-D",  "BSA allele depths (REF,ALT reads) per sample, WN2050",                            "raw_s24_wn2050",
  "Data S71", "Figure S26",     "Embryo counts and % gravid per replicate, RILs on other bacteria (specificity test)", "raw_s26",
  "Data S72", "Figure S27A-B",  "QTL-seq SNP-level statistics, EG4946",                                            "EG4946_data",
  "Data S73", "Figure S27A-B",  "QTL-seq significant regions, EG4946",                                             "EG4946_qtl",
  "Data S74", "Figure S27C",    "Embryo counts and % gravid per replicate, EG4946 RILs, EG4946 and jyIs8",           "data_s27c_raw",
  "Data S75", "Figure S27C",    "Genotype blocks (breakpoints) of EG4946 RILs and parents",                        "data_s27c_breakpoints",
  "Data S76", "Figure S28A-B",    "BSA allele depths (REF,ALT reads) per sample, JU2250",                            "raw_s28_ju2250",
  "Data S77", "Figure S28C-D",  "BSA allele depths (REF,ALT reads) per sample, EG4946",                            "raw_s28_eg4946",
  "Data S78", "Figure S29",     "JU2250 and EG4946 RIL genotype blocks across all chromosomes",                    "s29_breakpoints",
  "Data S79", "Figure S31",     "Worm body area per replicate (from FISH images)",                                 "data_s31"
)

# ==============================================================================
# 4. supp_stats_data.xlsx - statistics for the supplementary figures
# ==============================================================================
supp_stats <- tribble(
  ~id,              ~figure,          ~description,                                                                     ~object,
  "Stats Data S16", "Figure S2",      "Shapiro-Wilk normality of SNV/WT read ratio per strain",                         "normality_ratio_strain", 
  "Stats Data S17", "Figure S3A",     "Replicate correlation per condition (R², Pearson r, Spearman rho, p-values)",   "stats_supp3a",
  "Stats Data S18", "Figure S3A",     "Shapiro-Wilk normality tests (replicate z-scores and regression residuals)",    "normality_supp3a",
  "Stats Data S19", "Figure S4",      "Mini-phenoMIP replicate correlation per condition (R², r, rho, p-values)",      "stats_supp4",
  "Stats Data S20", "Figure S4",      "Shapiro-Wilk normality tests (replicate z-scores and regression residuals)",    "normality_supp4",
  "Stats Data S21", "Figure S5",      "PhenoMIP vs mini-phenoMIP correlation per condition (R, R², rho, p-values)",    "stats_supp5",
  "Stats Data S22", "Figure S5",      "Shapiro-Wilk normality tests (mean z-scores and regression residuals)",         "normality_supp5",
  "Stats Data S23", "Figure S6", "Shapiro-Wilk normality of average z-scores (all strains and conditions pooled), per metric", "normality_S6",
  "Stats Data S24", "Figure S9",      "Dunn's post-hoc test between locations (two-sided, BH-adjusted)",               "dunn_s9",
  "Stats Data S25", "Figure S9",      "Kruskal-Wallis test across locations per method and condition",                 "kw_s9",
  "Stats Data S26", "Figure S9",      "Shapiro-Wilk normality per location (context only; tests are rank-based)",      "normality_s9",
  "Stats Data S27", "Figure S10",     "Correlation of z-scores with divergence metrics (R², r, rho, p-values), PA14",  "stats_supp10",
  "Stats Data S28", "Figure S10",     "Shapiro-Wilk normality of regression residuals, PA14",                          "normality_supp10",
  "Stats Data S29", "Figure S11",     "Correlation of z-scores with divergence metrics (R², r, rho, p-values), PV",    "stats_supp11",
  "Stats Data S30", "Figure S11",     "Shapiro-Wilk normality of regression residuals, PV",                            "normality_supp11",
  "Stats Data S31", "Figure S12",     "Correlation of z-scores with divergence metrics (R², r, rho, p-values), MYb114", "stats_supp12",
  "Stats Data S32", "Figure S12",     "Shapiro-Wilk normality of regression residuals, MYb114",                        "normality_supp12",
  "Stats Data S33", "Figure S13",     "Correlation of z-scores with divergence metrics (R², r, rho, p-values), MYb117", "stats_supp13",
  "Stats Data S34", "Figure S13",     "Shapiro-Wilk normality of regression residuals, MYb117",                        "normality_supp13",
  "Stats Data S35", "Figure S14",     "Correlation of z-scores with divergence metrics (R², r, rho, p-values), MYb395", "stats_supp14",
  "Stats Data S36", "Figure S14",     "Shapiro-Wilk normality of regression residuals, MYb395",                        "normality_supp14",
  "Stats Data S37", "Figure S17",     "Pairwise Spearman correlations between conditions (unadjusted and BH p)",       "stats_s17",
  "Stats Data S38", "Figure S18",     "GWAS summary: markers, significance threshold, QTL count, heritability",        "gwas_summary_s18",
  "Stats Data S39", "Figure S22",     "Descriptive statistics and Dunnett's test vs N2 (embryo, % gravid)",            "stats_s22",
  "Stats Data S40", "Figure S22",     "Shapiro-Wilk normality tests per strain and on ANOVA residuals",                "normality_s22",
  "Stats Data S41", "Figure S26A-D",  "Descriptive statistics and Dunnett's test vs N2 (embryo, % gravid)",            "stats_s26",                
  "Stats Data S42", "Figure S26A-D",  "Shapiro-Wilk normality tests per strain and on ANOVA residuals",                "normality_s26",            
  "Stats Data S43", "Figure S27C",    "Descriptive statistics and Dunnett's test vs jyIs8 (embryo, % gravid)",            "supp27c$stats",
  "Stats Data S44", "Figure S27C",    "Shapiro-Wilk normality tests per line and on ANOVA residuals",                  "supp27c$normality",
  "Stats Data S45", "Figure S31",     "Descriptive statistics and Dunnett's test vs N2 (body area)",                   "stats_s31",
  "Stats Data S46", "Figure S31",     "Shapiro-Wilk normality tests per strain and on ANOVA residuals",                "normality_s31"
)

# ==============================================================================
# 5. Write all four workbooks
# ==============================================================================
write_workbook_with_toc(main_data,  "main_data.xlsx",       "Supplementary data for the main figures")
write_workbook_with_toc(main_stats, "main_stats_data.xlsx", "Statistical tests for the main figures")
write_workbook_with_toc(supp_data,  "supp_data.xlsx",       "Supplementary data for the supplementary figures")
write_workbook_with_toc(supp_stats, "supp_stats_data.xlsx", "Statistical tests for the supplementary figures")