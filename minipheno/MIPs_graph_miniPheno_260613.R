# mini-phenoMIP: make GC.abundance.all.data.txt from the strain-average file.
# Run with the working directory set to the folder this script is in.
source("../MIP.gc.abundance.R")

data.file <- "strain.avg.renamed.txt"
file.exists(data.file)
output.directory <- "minipheno_graphing"
num.rep <- 12
num.exp <- 3

# Detection thresholds for this dataset (see README)
min.abundance.value <- 0.0001
min.start.abundance.value <- 0.05

generate.gc.abundance(data.file,
                      output.directory,
                      num.rep,
                      num.exp,
                      min.abundance.value,
                      min.start.abundance.value)
