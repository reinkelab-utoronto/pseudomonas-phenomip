# phenoMIP 6.0 pool 3: make GC.abundance.all.data.txt from the strain-average file.
# Run with the working directory set to the folder this script is in.
source("../MIP.gc.abundance.R")

data.file <- "pool3_strain.avg.badmipsremoved.txt"
file.exists(data.file)
output.directory <- "PhenoMIP6.0_bac_P3"
num.rep <- 14
num.exp <- 3

# Detection thresholds for this dataset (see README)
min.abundance.value <- 0.000449983
min.start.abundance.value <- 0.002249915

generate.gc.abundance(data.file,
                      output.directory,
                      num.rep,
                      num.exp,
                      min.abundance.value,
                      min.start.abundance.value)
