# phenoMIP 6.0 pool 2: make GC.abundance.all.data.txt from the strain-average file.
# Run with the working directory set to the folder this script is in.
source("../MIP.gc.abundance.R")

data.file <- "pool2_strain.avg.badmipsremoved.txt"
file.exists(data.file)
output.directory <- "PhenoMIP6.0_bac_P2"
num.rep <- 14
num.exp <- 3

# Detection thresholds for this dataset (see README)
min.abundance.value <- 6.98e-06
min.start.abundance.value <- 3.49e-05

generate.gc.abundance(data.file,
                      output.directory,
                      num.rep,
                      num.exp,
                      min.abundance.value,
                      min.start.abundance.value)
