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


# remove the suffix for list names
list_names <- sub("\\.results\\.txt$", "", filenames)

# read all files into a list
df_list <- setNames(
  lapply(filenames, function(f) {
    read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  }),
  list_names)

combined_df <- imap_dfr(df_list, ~ mutate(.x, sample = .y))

all_samples <- combined_df %>%
  mutate(
    probe_id   = paste(snv.chr, snv.loc, sep = "_"),
    generation = str_extract(sample, "(?<=G)\\d+"),
    replicate  = str_extract(sample, "(?<=R)\\d+"),
    condition  = str_extract(sample, "(?<=R\\d_).*(?=_S\\d+$)"),
    condition  = replace_na(condition, "START")
  )

#probes_to_remove <- c("III_10989744", "X_11299050", "III_4337091", "II_7590562",
#                      "III_7632363", "V_14343972", "V_9929788", "V_3365505", "X_5271327")
#
#all_samples <- all_samples %>%
#  filter(!probe_id %in% probes_to_remove)
  
avg_ratio_table <- all_samples %>%
  group_by(snv.strain, sample) %>%
  summarise(mean_snv_wt_ratio = mean(snv.wt.unique.ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = sample,
              values_from = mean_snv_wt_ratio) %>%
  arrange(snv.strain)

# Replace underscores with dots in column names
colnames(avg_ratio_table) <- gsub("_", ".", colnames(avg_ratio_table))
head(avg_ratio_table)


cols <- colnames(avg_ratio_table)

# G0 samples: replace the final .SXX with .START
cols <- gsub("^G0\\.(R\\d)\\.S\\d+$", "G0.\\1.START", cols)

# All other samples: remove the trailing .SXX
cols <- gsub("\\.S\\d+$", "", cols)

colnames(avg_ratio_table) <- cols

# Function to reformat a single column name
reformat_colname <- function(name) {
  parts <- strsplit(name, "\\.")[[1]]
  if (length(parts) != 3) return(name)  # leave snv.strain or unexpected names unchanged
  
  gen <- parts[1]   # e.g. "G5"
  rep <- parts[2]   # e.g. "P1A"
  cond <- parts[3]  # e.g. "LUA"
  
  # Extract generation number and zero-pad it, then prefix with E0
  gen_num <- sub("G", "", gen)
  new_gen <- paste0("E0", gen_num)  # G5 -> E05
  
  paste(cond, rep, new_gen, sep = ".")  # LUA.P1A.E05
}

# Apply to all column names, skipping snv.strain
colnames(avg_ratio_table) <- sapply(colnames(avg_ratio_table), function(name) {
  if (name == "snv.strain") name else reformat_colname(name)
})

avg_ratio_table <- avg_ratio_table %>%
  dplyr::select(-3)

colnames(avg_ratio_table)

colnames(avg_ratio_table) <- gsub("^START\\.R1", "MO.R00", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^UNI1\\.R1", "MO.R01", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^UNI1\\.R2", "MO.R02", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^UNI2\\.R1", "MO.R03", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^UNI2\\.R2", "MO.R04", colnames(avg_ratio_table))

colnames(avg_ratio_table) <- gsub("^PA\\.R1", "PA.R05", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^PA\\.R2", "PA.R06", colnames(avg_ratio_table))

colnames(avg_ratio_table) <- gsub("^PV\\.R1", "PV.R07", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^PV\\.R2", "PV.R08", colnames(avg_ratio_table))

colnames(avg_ratio_table) <- gsub("^M114\\.R1", "M114.R09", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^M114\\.R2", "M114.R10", colnames(avg_ratio_table))

colnames(avg_ratio_table) <- gsub("^M395\\.R1", "M395.R11", colnames(avg_ratio_table))
colnames(avg_ratio_table) <- gsub("^M395\\.R2", "M395.R12", colnames(avg_ratio_table))

colnames(avg_ratio_table)

sort_cols <- function(df) {
  # Separate snv.strain from the rest
  data_cols <- colnames(df)[colnames(df) != "snv.strain"]
  
  # Extract E and R numbers for sorting
  e_nums <- as.integer(sub(".*\\.E(\\d+)$", "\\1", data_cols))
  r_nums <- as.integer(sub(".*\\.R(\\d+)\\..*", "\\1", data_cols))
  
  # Sort by E number first, then R number
  sorted_cols <- data_cols[order(e_nums, r_nums)]
  
  df[, c("snv.strain", sorted_cols)]
}

avg_ratio_table <- sort_cols(avg_ratio_table)

colnames(avg_ratio_table)

setwd("..")

write.table(avg_ratio_table, "strain.avg.renamed.txt", sep = "\t", row.names = FALSE, quote = FALSE)
