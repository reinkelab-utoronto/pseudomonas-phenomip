library(tidyverse)
library(dplyr)
library(reshape2)
library(tidyr)
library(e1071)
library(fitdistrplus)
library(data.table)

setwd("mip_results")

##Read files named xyz1111.csv, xyz2222.csv, etc.
filenames <- list.files(path=".",
                        pattern="*results.txt")

filenames_filtered <- filenames[!grepl("G1_|G2_|G3_|SM", filenames)]

# remove the suffix for list names
list_names <- sub("\\.results\\.txt$", "", filenames_filtered)

# read all files into a list
df_list <- setNames(
  lapply(filenames_filtered, function(f) {
    read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  }),
  list_names
)
setwd("..")

strains_analysed <- read.csv("strains_analysed.csv")



# Create the strain_to_pool lookup vector
strain_to_pool <- setNames(strains_analysed$pool, strains_analysed$strain)

# Apply to all dataframes in df_list
df_list <- lapply(df_list, function(df) {
  df$pool <- strain_to_pool[df$snv.strain]
  return(df)
})

df_list <- lapply(df_list, function(df) {
  df[!is.na(df$pool), ]
})

df_list$G4_P1BAC_START

# Split into P1_df, P2_df, P3_df based on names
P1_df <- df_list[grepl("P1", names(df_list))]
P2_df <- df_list[grepl("P2", names(df_list))]
P3_df <- df_list[grepl("P3", names(df_list))]


# Optional: Verify the split
cat("Summary of split:\n")
cat("Total dataframes:", length(df_list), "\n")
cat("P1 dataframes:", length(P1_df), "\n")
cat("P2 dataframes:", length(P2_df), "\n")
cat("P3 dataframes:", length(P3_df), "\n")

P1_df <- lapply(P1_df, function(df) {
  df$pool[df$pool == "control"] <- "pool1"
  df
})
P2_df <- lapply(P2_df, function(df) {
  df$pool[df$pool == "control"] <- "pool2"
  df
})
P3_df <- lapply(P3_df, function(df) {
  df$pool[df$pool == "control"] <- "pool3"
  df
})
# Loop through all dataframes in P1_df
for (i in seq_along(P1_df)) {
  P1_df[[i]] <- P1_df[[i]][P1_df[[i]]$pool == "pool1", ]
}

# Loop through all dataframes in P2_df
for (i in seq_along(P2_df)) {
  P2_df[[i]] <- P2_df[[i]][P2_df[[i]]$pool == "pool2", ]
}

# Loop through all dataframes in P3_df
for (i in seq_along(P3_df)) {
  P3_df[[i]] <- P3_df[[i]][P3_df[[i]]$pool == "pool3", ]
}

# --- 0. Combine all pools into one dataframe
all_samples <- imap_dfr(P1_df, ~ mutate(.x, sample = .y)) %>% mutate(pool = "P1") %>%
  bind_rows(
    imap_dfr(P2_df, ~ mutate(.x, sample = .y)) %>% mutate(pool = "P2")
  ) %>%
  bind_rows(
    imap_dfr(P3_df, ~ mutate(.x, sample = .y)) %>% mutate(pool = "P3")
  )

# --- 1. Create unique probe ID + parse sample metadata
all_samples <- all_samples %>%
  mutate(
    probe_id   = paste(snv.chr, snv.loc, sep = "_"),
    generation = str_extract(sample, "(?<=G)\\d+"),
    pool       = str_extract(sample, "(?<=P)\\d+"),        # extract from sample name
    replicate  = str_extract(sample, "(?<=P\\d)\\w+(?=_)"),
    condition  = str_extract(sample, "[^_]+$")
  )



# --- 2. Compute correlations per strain (now also per pool)
compute_strain_correlations <- function(df, method = "spearman") {
  
  df %>%
    split(list(.$snv.strain, .$pool), drop = TRUE) %>%   # split by strain AND pool
    lapply(function(strain_df) {
      
      strain_name <- unique(strain_df$snv.strain)
      pool_name   <- unique(strain_df$pool)
      
      # Pivot wide: samples as rows, probe_id as columns
      wide <- strain_df %>%
        dplyr::select(sample, probe_id, snv.wt.unique.ratio) %>%
        pivot_wider(names_from  = probe_id,
                    values_from = snv.wt.unique.ratio)
      
      mat <- wide %>% dplyr::select(-sample) %>% as.matrix()
      rownames(mat) <- wide$sample
      
      # --- Only remove probes that are ALL NA or have zero variance ---
      keep <- apply(mat, 2, function(x) {
        valid <- x[!is.na(x)]
        length(valid) >= 2 && var(valid) > 0
      })
      
      removed <- names(keep[!keep])
      if (length(removed) > 0) {
        message(sprintf("Pool %s Strain %s: removing uninformative probe(s): %s",
                        pool_name, strain_name, paste(removed, collapse = ", ")))
      }
      
      mat <- mat[, keep, drop = FALSE]
      
      # Skip only if fewer than 2 probes remain
      if (ncol(mat) < 2) {
        message(sprintf("Skipping Pool %s Strain %s: fewer than 2 informative probes remaining",
                        pool_name, strain_name))
        return(NULL)
      }
      
      # --- Correlation matrix ---
      cor_matrix <- cor(mat, method = method, use = "pairwise.complete.obs")
      
      # --- Long-format pairwise table ---
      cor_long <- cor_matrix %>%
        as.data.frame() %>%
        rownames_to_column("probe_1") %>%
        pivot_longer(-probe_1, names_to = "probe_2", values_to = "correlation") %>%
        filter(probe_1 < probe_2) %>%
        mutate(
          strain   = strain_name,
          pool     = pool_name,
          chr_1    = str_extract(probe_1, "^[^_]+"),
          chr_2    = str_extract(probe_2, "^[^_]+"),
          same_chr = chr_1 == chr_2
        ) %>%
        dplyr::select(pool, strain, probe_1, probe_2, chr_1, chr_2, same_chr, correlation) %>%
        arrange(desc(abs(correlation)))
      
      list(matrix = cor_matrix, long = cor_long)
    }) %>%
    Filter(Negate(is.null), .)
}

# --- 3. Run it
results <- compute_strain_correlations(all_samples, method = "spearman")

JU2001 <- all_samples %>%
  filter(snv.strain == "JU2001") %>%
  dplyr::select(MIP.name, snv.strain, snv.loc, snv.wt.unique.ratio, sample)

# --- 4. Combine all long-format tables into one master table
all_cor_long <- results %>%
  lapply(`[[`, "long") %>%
  bind_rows()

all_cor_long <- all_cor_long %>%
  mutate(r_squared = correlation^2)

# --- 5. Boxplot faceted by pool
# Formatting matched to main Figure 1 (theme_bw, base_size = 8,
# 8pt axis titles, 6pt axis text, 8pt bold strip text, 8pt bold plot title)
correlation_p <- ggplot(all_cor_long, aes(x = strain, y = correlation, fill = strain)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, colour = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "red") +
  facet_wrap(~ pool, ncol = 1, scales = "free") +
  labs(x = "Strain", y = "Spearman correlation") +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y     = element_text(size = 6),
    axis.title      = element_text(size = 8),
    strip.text      = element_text(size = 8, face = "bold"),
    plot.title      = element_text(size = 8, face = "bold")
  )

# Full-page figure, so dimensions kept at A4
ggsave("supp_fig2.pdf", plot = correlation_p, width = 180, height = 220, units = "mm")


# --- 6. Density plot coloured by pool
#ggplot(all_cor_long, aes(x = correlation, fill = pool, colour = pool)) +
#  geom_density(alpha = 0.3) +
#  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
#  labs(title = "Distribution of pairwise Spearman correlations by pool",
#       x = "Spearman correlation", y = "Density") +
#  theme_bw()
#
#all_cor_long <- all_cor_long %>%
#  mutate(pool = factor(pool))
#
#density_p <- ggplot(all_cor_long, aes(x = correlation, 
#                         colour = pool, 
#                         group = interaction(strain, pool))) +
#  geom_density(alpha = 0.3) +
#  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
#  labs(title = "Distribution of pairwise Spearman correlations by strain (coloured by pool)",
#       x = "Spearman correlation", y = "Density",
#       fill = "Pool", colour = "Pool") +
#  theme_bw() +
#  theme(legend.position = "bottom")
#
#ggsave("MIPs_density_removedbadmips.pdf", plot = density_p, height = 210, width = 297, units = "mm")


low_cor_strains <- all_cor_long %>%
  group_by(pool, strain) %>%
  summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_cor < 0.7)

# Count
nrow(low_cor_strains)

# View which ones
print(low_cor_strains)

write.csv(all_cor_long, "Pheno6_MIPs_correlation_removedbadmips.csv", row.names = F)


probe_table <- all_samples %>%
  dplyr::select(snv.strain, sample, probe_id, snv.wt.unique.ratio) %>%
  group_by(snv.strain) %>%
  group_modify(~ pivot_wider(.x, names_from = probe_id, values_from = snv.wt.unique.ratio)) %>%
  ungroup() %>%
  arrange(snv.strain, sample)

write.csv(probe_table, "Pheno6_probe_table_removedbadmips.csv", row.names = FALSE)

# Subset probe_table
low_cor_probe_table <- probe_table %>%
  filter(snv.strain %in% low_cor_strains$strain)
write.csv(low_cor_probe_table, "Pheno6_low_cor_probe_table_removedbadmips.csv", row.names = FALSE)


# Subset all_samples
low_cor_samples <- all_samples %>%
  filter(snv.strain %in% low_cor_strains$strain)
write.csv(low_cor_samples, "Pheno6_low_cor_results_removedbadmips.csv", row.names = FALSE)

#subset cor_long
low_cor_long <- all_cor_long %>%
  filter(strain %in% low_cor_strains$strain)

write.csv(low_cor_long, "Pheno6_low_cor_comparison_removedbadmips.csv", row.names = FALSE)



# Overall (all rows): Anderson-Darling, since Shapiro-Wilk is capped at 5000
library(nortest)
ratio <- na.omit(all_samples$snv.wt.unique.ratio)
ad_all <- ad.test(ratio)

normality_ratio_overall <- tibble(
  n     = length(ratio),
  test  = "Anderson-Darling",
  A     = unname(ad_all$statistic),
  p     = ad_all$p.value
)

normality_ratio_strain <- all_samples %>%
  group_by(pool, snv.strain) %>%
  summarise(
    n         = sum(!is.na(snv.wt.unique.ratio)),
    shapiro_W = unname(safe_shapiro(snv.wt.unique.ratio)["W"]),
    shapiro_p = unname(safe_shapiro(snv.wt.unique.ratio)["p"]),
    .groups = "drop"
  )
