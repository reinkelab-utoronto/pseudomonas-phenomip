library(tidyverse)
library(dplyr)
library(reshape2)
library(tidyr)
library(e1071)
library(fitdistrplus)
library(data.table)

setwd("../mip_results")

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
  list_names)

combined_df <- imap_dfr(df_list, ~ mutate(.x, sample = .y))

all_samples <- combined_df %>%
  mutate(
    probe_id   = paste(snv.chr, snv.loc, sep = "_"),
    generation = str_extract(sample, "(?<=G)\\d+"),
    pool       = str_extract(sample, "(?<=P)\\d+"),        # extract from sample name
    replicate  = str_extract(sample, "(?<=P\\d)\\w+(?=_)"),
    condition  = str_extract(sample, "[^_]+$")
  )

probes_to_remove <- c("III_10989744", "X_11299050", "III_4337091", "II_7590562",
                      "III_7632363", "V_14343972", "V_9929788", "V_3365505", "X_5271327")

all_samples <- all_samples %>%
  filter(!probe_id %in% probes_to_remove)
  
avg_ratio_table <- all_samples %>%
  group_by(snv.strain, sample) %>%
  summarise(mean_snv_wt_ratio = mean(snv.wt.unique.ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = sample,
              values_from = mean_snv_wt_ratio) %>%
  arrange(snv.strain)

# Replace underscores with dots in column names
colnames(avg_ratio_table) <- gsub("_", ".", colnames(avg_ratio_table))

# Replace NaN and NA with 0 (tibble-safe)
avg_ratio_table <- avg_ratio_table %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.na(.), 0, .)))

# Function to reformat a single column name
reformat_colname <- function(name) {
  parts <- strsplit(name, "\\.")[[1]]
  if (length(parts) != 3) return(name)  # leave snv.strain or unexpected names unchanged
  
  gen <- parts[1]   # e.g. "G5"
  rep <- parts[2]   # e.g. "P1A"
  cond <- parts[3]  # e.g. "M117"
  
  # Extract generation number and zero-pad it, then prefix with E0
  gen_num <- sub("G", "", gen)
  new_gen <- paste0("E0", gen_num)  # G5 -> E05
  
  paste(cond, rep, new_gen, sep = ".")  # M117.P1A.E05
}

# Apply to all column names, skipping snv.strain
colnames(avg_ratio_table) <- sapply(colnames(avg_ratio_table), function(name) {
  if (name == "snv.strain") name else reformat_colname(name)
})


# Helper to extract pool number from a column name
get_pool <- function(name) {
  if (name == "snv.strain") return(NA)
  rep <- strsplit(name, "\\.")[[1]][2]  # e.g. "P1A"
  sub("P(\\d).*", "\\1", rep)           # extract the digit after P
}

# Get pool assignment for each column
pools <- sapply(colnames(avg_ratio_table), get_pool)

# Split into 3 dataframes, always keeping snv.strain
pool1 <- avg_ratio_table[, is.na(pools) | pools == "1"]
pool2 <- avg_ratio_table[, is.na(pools) | pools == "2"]
pool3 <- avg_ratio_table[, is.na(pools) | pools == "3"]

# Pool 1
colnames(pool1) <- gsub("^START\\.P1BAC", "MO.R00", colnames(pool1))
# Start population: label it generation E00 (this was previously renamed by hand)
colnames(pool1) <- sub("^MO\\.R00\\.E04$", "MO.R00.E00", colnames(pool1))
colnames(pool1) <- gsub("^UN1BAC\\.P1A", "MO.R01", colnames(pool1))
colnames(pool1) <- gsub("^UN1BAC\\.P1B", "MO.R02", colnames(pool1))
colnames(pool1) <- gsub("^UN2BAC\\.P1A", "MO.R03", colnames(pool1))
colnames(pool1) <- gsub("^UN2BAC\\.P1B", "MO.R04", colnames(pool1))

# Pool 2
colnames(pool2) <- gsub("^START\\.P2BAC", "MO.R00", colnames(pool2))
# Start population: label it generation E00 (this was previously renamed by hand)
colnames(pool2) <- sub("^MO\\.R00\\.E04$", "MO.R00.E00", colnames(pool2))
colnames(pool2) <- gsub("^UN1BAC\\.P2A", "MO.R01", colnames(pool2))
colnames(pool2) <- gsub("^UN1BAC\\.P2B", "MO.R02", colnames(pool2))
colnames(pool2) <- gsub("^UN2BAC\\.P2A", "MO.R03", colnames(pool2))
colnames(pool2) <- gsub("^UN2BAC\\.P2B", "MO.R04", colnames(pool2))

# Pool 3
colnames(pool3) <- gsub("^START\\.P3BAC", "MO.R00", colnames(pool3))
# Start population: label it generation E00 (this was previously renamed by hand)
colnames(pool3) <- sub("^MO\\.R00\\.E04$", "MO.R00.E00", colnames(pool3))
colnames(pool3) <- gsub("^UN1BAC\\.P3A", "MO.R01", colnames(pool3))
colnames(pool3) <- gsub("^UN1BAC\\.P3B", "MO.R02", colnames(pool3))
colnames(pool3) <- gsub("^UN2BAC\\.P3A", "MO.R03", colnames(pool3))
colnames(pool3) <- gsub("^UN2BAC\\.P3B", "MO.R04", colnames(pool3))

# Pool 1
colnames(pool1) <- gsub("^PA14\\.P1A", "PA.R05", colnames(pool1))
colnames(pool1) <- gsub("^PA14\\.P1B", "PA.R06", colnames(pool1))
colnames(pool1) <- gsub("^PV\\.P1A",   "PV.R07", colnames(pool1))
colnames(pool1) <- gsub("^PV\\.P1B",   "PV.R08", colnames(pool1))
colnames(pool1) <- gsub("^M114\\.P1A", "M114.R09", colnames(pool1))
colnames(pool1) <- gsub("^M114\\.P1B", "M114.R10", colnames(pool1))
colnames(pool1) <- gsub("^M117\\.P1A", "M117.R11", colnames(pool1))
colnames(pool1) <- gsub("^M117\\.P1B", "M117.R12", colnames(pool1))
colnames(pool1) <- gsub("^M395\\.P1A", "M395.R13", colnames(pool1))
colnames(pool1) <- gsub("^M395\\.P1B", "M395.R14", colnames(pool1))

# Pool 2
colnames(pool2) <- gsub("^PA14\\.P2A", "PA.R05", colnames(pool2))
colnames(pool2) <- gsub("^PA14\\.P2B", "PA.R06", colnames(pool2))
colnames(pool2) <- gsub("^PV\\.P2A",   "PV.R07", colnames(pool2))
colnames(pool2) <- gsub("^PV\\.P2B",   "PV.R08", colnames(pool2))
colnames(pool2) <- gsub("^M114\\.P2A", "M114.R09", colnames(pool2))
colnames(pool2) <- gsub("^M114\\.P2B", "M114.R10", colnames(pool2))
colnames(pool2) <- gsub("^M117\\.P2A", "M117.R11", colnames(pool2))
colnames(pool2) <- gsub("^M117\\.P2B", "M117.R12", colnames(pool2))
colnames(pool2) <- gsub("^M395\\.P2A", "M395.R13", colnames(pool2))
colnames(pool2) <- gsub("^M395\\.P2B", "M395.R14", colnames(pool2))

# Pool 3
colnames(pool3) <- gsub("^PA14\\.P3A", "PA.R05", colnames(pool3))
colnames(pool3) <- gsub("^PA14\\.P3B", "PA.R06", colnames(pool3))
colnames(pool3) <- gsub("^PV\\.P3A",   "PV.R07", colnames(pool3))
colnames(pool3) <- gsub("^PV\\.P3B",   "PV.R08", colnames(pool3))
colnames(pool3) <- gsub("^M114\\.P3A", "M114.R09", colnames(pool3))
colnames(pool3) <- gsub("^M114\\.P3B", "M114.R10", colnames(pool3))
colnames(pool3) <- gsub("^M117\\.P3A", "M117.R11", colnames(pool3))
colnames(pool3) <- gsub("^M117\\.P3B", "M117.R12", colnames(pool3))
colnames(pool3) <- gsub("^M395\\.P3A", "M395.R13", colnames(pool3))
colnames(pool3) <- gsub("^M395\\.P3B", "M395.R14", colnames(pool3))


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

pool1 <- sort_cols(pool1)
pool2 <- sort_cols(pool2)
pool3 <- sort_cols(pool3)

setwd("../upstream")

write.table(pool1, "pool1_strain.avg.badmipsremoved.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pool2, "pool2_strain.avg.badmipsremoved.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pool3, "pool3_strain.avg.badmipsremoved.txt", sep = "\t", row.names = FALSE, quote = FALSE)
