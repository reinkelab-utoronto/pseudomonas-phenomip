# Repository layout

This file lists every script and input file in the repository and where it lives. Everything else in
the folder after a run is output, written by the scripts. `.gitignore` keeps outputs out of the
repository: it ignores everything except the files listed here.

## Root (working directory for `pheno6_bac_v2.R`)

### Scripts (51)

| Group | Scripts |
|---|---|
| Setup | `00_install_packages.R` (run once) |
| Master | `pheno6_bac_v2.R` |
| BSA analysis | `BSA_all_populations.R` (run before the master; see below) |
| Abundance tables | `MIP.gc.abundance.R` (used by the setup scripts in `upstream/` and `minipheno/`) |
| Z-scores, PCA, clustering, heritability | `making_zscores.R`, `bac_PCA_all.R`, `bac_PCA_bycondition.R`, `bac_PCA_bycluster.R`, `bac_PCA_bycluster_mean.R`, `pheno6_controls_correlation.R`, `pheno6_bac_auc_heatmap_dendro.R`, `pheno6_bac_regression_heatmap_dendro.R`, `pheno6_bac_PC1_heatmap_dendro.R`, `auc_hclust.R`, `regression_hclust.R`, `auc_heritability_input.R`, `slope_heritability_input.R` |
| Maps and genetic distance | `worldmapv2.R`, `genetic_distance_zscorediff_auc.R`, `genetic_distance_zscorediff_regression.R`, `genetic_distance_PC1diff.R` |
| Main figures | `main_figure1.R` … `main_figure6.R` |
| Supplementary figures | `supp_fig1.R`, `supp_fig2.R`, `supp_fig3.R`, `supp_fig4_reordered.R`, `supp_fig5_reordered.R`, `supp_fig6_reodered.R`, `supp_fig7_reordered.R`, `supp_fig8_reordered.R`, `supp_fig9_reorderd.R`, `supp_fig15_reordered.R`, `supp_fig16_reordered.R`, `supp_fig17_reordered.R`, `supp_fig20_reordered.R`, `supp_fig21_reordered.R`, `supp_fig22.R`, `supp_fig23.R`, `supp_fig24.R`, `supp_fig25.R`, `supp_fig26_reordered.R`, `supp_fig27_reordered.R`, `supp_fig28_reordered.R`, `supp_fig29_reordered.R`, `supp_fig31.R` |
| Supplementary workbooks | `supp_data.R` |

Script names must stay exactly as they are (including `supp_fig6_reodered.R` and `supp_fig9_reorderd.R`),
because the master script calls them by name.

### Input files (36)

| Data | Files |
|---|---|
| Pooling plan and strain lists | `final_freezing_plan.csv`, `good_strains.csv`, `strains_analysed.csv`, `Table S1.csv` |
| Strain metadata and trees | `20231213_c_elegans_strain_data.csv`, `wild_isolate_variants.csv`, `WI.20231213.hard-filter.isotype.min4.tree`, `SpeciesTree_rooted copy.tree` |
| Microsporidia z-scores | `mic_all3data.csv` |
| % males | `pheno8_percent_males.csv` |
| RIL phenotypes and genotypes | `JU2610_RILs_M114_avg_embryos.csv`, `JU2610_RILs_M114_gravid.csv`, `JU2610_RILs_specificity_embryo.csv`, `JU2610_RILs_specificity_gravid.csv`, `JU2610_RILs_joint.filtered.tsv`, `JU2522_RILs_PA14_avg_embryos.csv`, `JU2522_RILs_PA14_gravid.csv`, `JU2522_RILs_specificity_embryos.csv`, `JU2522_RILs_specificity_gravid.csv`, `JU2522_RILs_joint.filtered.tsv`, `JU2250_RILs_avg_embryos.csv`, `JU2250_RILs_gravid.csv`, `JU2250_RILs_joint.filtered.tsv`, `EG4946_RILs_avg_embryos.csv`, `EG4946_RILs_gravid.csv`, `EG4946_RILs_joint.filtered.tsv` |
| Wild-isolate phenotypes | `JU2610_wild_isolates_canddiategenetest_avg_embryo.csv`, `JU2610_wild_isolates_canddiategenetest_gravid.csv`, `m395_wild_isolates_embryos.csv`, `m395_wild_isolates_gravid.csv` |
| FISH measurements | `M114_FISH_area.csv`, `M114_FISH_body.csv`, `M395_FISH_area.csv`, `M395_FISH_body.csv`, `PA14_FISH_area.csv`, `PA14_FISH_body.csv` |

## `BSA_results/`

| File | Role |
|---|---|
| `BSA_all.slurm` | cluster script (SLURM array, one task per population): reads → fastp → BWA-MEM → Picard MarkDuplicates → bcftools mpileup/call/filter → allele depths |
| `<POP>/<POP>_allele_depths.tsv` | CHROM, POS and per-sample allele depths (no header); read by `supp_fig21_reordered.R`, `supp_fig24.R` and `supp_fig28_reordered.R` |
| `<POP>/<POP>_allele_depths_columns.tsv` | the same with REF, ALT and a header of sample names; read by `BSA_all_populations.R` |
| `<POP>/<POP>_samples.txt`, `<POP>/<POP>_versions.txt` | sample order of the allele-depth columns, and the tool versions used |

`<POP>` is JU2610, ED3005, JU2522, WN2050, JU2250 or EG4946. The `*_BSA_results.tar.gz` archives from the cluster are not included.

`BSA_all_populations.R` (run from the project folder) reads these files and writes, to the project folder, the
QTLseqr tables used by the figure scripts (`<POP>_<comparison>.tsv` and `<POP>_<comparison>_results.csv`),
and to `BSA_exploratory/` the exploratory allele-frequency, odds-ratio, G′ and ΔSNP plots and `site_counts.tsv`.

## RIL folders (`JU2610_RILs/`, `JU2522_RILs/`, `JU2250_RILs/`, `EG4946_RILs/`)

| File | Role |
|---|---|
| `*.slurm` | cluster scripts: alignment and joint variant calling (bcftools mpileup/call, LOWQUAL flag) |
| `*.out`, `*.err` | the cluster logs of those runs |
| `*joint.filtered.vcf.gz` (+ `.csi`) | the filtered joint VCF of all RILs and both parents |
| `<POP>_vk_hmm_commands.sh` | the VCF-kit hmm command used to infer parental haplotypes |
| `*joint.filtered.tsv` | VCF-kit hmm output; a copy in the project folder, `<POP>_RILs_joint.filtered.tsv`, is read by the figure scripts |

Batch 1 (JU2250, EG4946) files use "NIL" in their names (e.g. `JU_NIL_A1`); these are the RILs described in the paper.

## `mip_results/`

The 129 phenoMIP 6.0 bacteria `*.results.txt` files: the G4 start populations and generations 5, 7 and
9 for every condition. Read by `upstream/remake_strain.average_removebadmips.R` and by `supp_fig2.R`.

## `upstream/` (working directory for stage 1)

| File | Role |
|---|---|
| `remake_strain.average_removebadmips.R` | `../mip_results/*.results.txt` → `pool1/2/3_strain.avg.badmipsremoved.txt` |
| `MIPs_graph_P1bac_261020.R`, `MIPs_graph_P2bac_261020.R`, `MIPs_graph_P3bac_261020.R` | strain averages → `PhenoMIP6.0_bac_P1/2/3/by.rep/growth.curve.data/GC.abundance.all.data.txt` (read by the master) |

## `minipheno/` (working directory for stage 2)

| File | Role |
|---|---|
| `mip_results/` | the 38 mini-phenoMIP `*.results.txt` files |
| `mini_pheno_remake_strain.average_removebadmips.R` | `mip_results/*.results.txt` → `strain.avg.renamed.txt` |
| `MIPs_graph_miniPheno_260613.R` | `strain.avg.renamed.txt` → `minipheno_graphing/by.rep/growth.curve.data/GC.abundance.all.data.txt` |
| `mini_pheno_analysis.R` | mini-phenoMIP z-scores; writes `minipheno_auc.csv` and `minipheno_regression.csv` to the root. Run by the master (with the root as working directory) before Figure 1, and again by `supp_fig4_reordered.R` |
| `correlation_allstatsAUC.R`, `correlation_allstatsregression.R` | phenoMIP vs mini-phenoMIP correlations; run by `main_figure1.R` and `supp_fig5_reordered.R` |

## `gwas/`

NemaScan outputs for the two traits with genome-wide significant QTL, MYb114 AUC z-score (LOCO model)
and MYb395 slope z-score (inbred model), downloaded from the CaeNDR portal and renamed with a trait
prefix where the portal gives a generic name:

| File (per trait) | Contents |
|---|---|
| `M114_avg_auc_zscore_loco.gwa`, `M395_avg_slope_zscore_inbred.gwa` | association results for every marker |
| `M114_avg_auc_zscore_loco_qtl.tsv`, `M395_avg_slope_zscore_inbred_qtl.tsv` | QTL interval and peak marker |
| `M114_genotype_matrix.tsv`, `M395_genotype_matrix.tsv` | strain genotypes (portal name: `genotype_matrix.tsv`) |
| `M114_independent_tests.txt`, `M395_independent_tests.txt` | number of independent tests (portal name: `independent_tests.txt`) |
| `M114_trait_heritability.tsv`, `M395_trait_heritability.tsv` | narrow-sense heritability (portal name: `trait_heritability.tsv`) |
| `NemaScan_Report_M114_avg_auc_zscore.html`, `NemaScan_Report_M395_avg_slope_zscore.html` | NemaScan reports; `supp_fig18.R` reads the strain phenotypes from these |

plus `supp_fig18.R`, which the master runs with `gwas/` as the working directory and which writes
`../supp_fig18.pdf`.

## `figures/`

Final PDFs of the figures as published: `Figure_1_Final.pdf` … `Figure_6_Final.pdf` and the code-generated
supplementary figures (`supp_fig*.pdf`). These are copies; the scripts write new versions to the root, which
are not tracked.

## Folders created by the scripts

`BSA_exploratory/`, `phenomip6.0_nemascan_input/`, `phenomip6.0_ranked_plots/`, `non_divergent_nemascan_input/` (the
NemaScan input files), `minipheno_growth_curves/`, `mini_pheno_regression_growth_curves/`, and the
`PhenoMIP6.0_bac_P*` and `minipheno_graphing` output folders. All are created automatically.

## How to run

| Stage | Working directory | Run |
|---|---|---|
| 1 | `upstream/` | `remake_strain.average_removebadmips.R`, then the three `MIPs_graph_P*bac_261020.R` |
| 2 | `minipheno/` | `mini_pheno_remake_strain.average_removebadmips.R`, then `MIPs_graph_miniPheno_260613.R` |
| 3 | root | `BSA_all_populations.R`, then `pheno6_bac_v2.R` |

Start each stage in a fresh R session.
