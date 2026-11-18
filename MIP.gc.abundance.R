# MIP.gc.abundance.R
#
# Produces only GC.abundance.all.data.txt (the growth-curve abundance table) from a
# strain-average file. This is a reduced version of the lab's MIP.graph.analysis.v14_6.R:
# the code that calculates this table is copied unchanged from that script, and every
# other output (heatmaps, growth-curve plots, normalized data, Excel files) is removed.
# The output is identical to the GC.abundance.all.data.txt that the full script writes.
#
# The two detection thresholds were global variables edited inside each copy of the
# original script. Here they are arguments, so one file serves every pool.
#
# Output: <output.directory>/by.rep/growth.curve.data/GC.abundance.all.data.txt
#         <output.directory>/run.info.txt (settings used and strains removed)

generate.gc.abundance = function(data.file, output.directory, num.rep, num.exp,
                                 min.abundance.value, min.start.abundance.value) {

  # From the original script's global variables
  gen.offset = -1

  output.directory = paste0(c("./", output.directory, "/"), collapse="")
  dir.create(output.directory, recursive=TRUE, showWarnings=FALSE)
  output.info = paste0(c(output.directory, "run.info.txt"), collapse = "")
  write("MIP.gc.abundance.R (GC table only, from MIP.graph.analysis.v14_6.R)", output.info, append=FALSE)
  write(paste0("Data file: ", data.file), output.info, append=TRUE)
  write(paste0("Minimum starting abundance: ", min.start.abundance.value), output.info, append=TRUE)
  write(paste0("Minimum abundance value: ", min.abundance.value), output.info, append=TRUE)
  write(paste0("Program start at ", Sys.time()), output.info, append=TRUE)

  # ---------- From generate.plots() ----------
  data.df.exp = read.table(data.file, header=TRUE, sep="\t", row.names=1)

  # Remove strains whose start abundance is below the start threshold
  trimmed.strains = which(data.df.exp[,1] < min.start.abundance.value)
  write(paste0("Total strains below start threshold: ", length(trimmed.strains)), output.info, append=TRUE)
  write(rownames(data.df.exp)[trimmed.strains], output.info, append=TRUE)
  if (length(trimmed.strains) > 0) {
    data.df.exp <- data.df.exp[-trimmed.strains,]
    write(paste0("Total strains remaining for analysis: ", nrow(data.df.exp)), output.info, append=TRUE)
  }

  # Raise values below the detection limit to the detection limit
  data.df.exp[data.df.exp < min.abundance.value] <- min.abundance.value

  col.info = unlist(strsplit(colnames(data.df.exp), split="[.]"))
  exp.names = unique(as.numeric(gsub('[A-Z]', "",unique(col.info[grep("E", col.info)])))) + gen.offset

  # Reorder columns from generation-major to replicate-major
  data.df.rep = data.df.exp[,0]
  data.df.rep = cbind(data.df.rep, data.df.exp[,1])
  data.df.rep.names = colnames(data.df.exp)[1]
  for (i in 2:(num.rep+1)) {
    curr.col = seq(i, ncol(data.df.exp), by = num.rep)
    data.df.rep = cbind(data.df.rep, data.df.exp[,curr.col])
    data.df.rep.names = c(data.df.rep.names, colnames(data.df.exp)[curr.col])
  }
  colnames(data.df.rep) = data.df.rep.names

  # ---------- From growth.curve(), as called by make.rep.plot.data() ----------
  series.data = data.df.rep
  start.col = 1

  # Only strains with a non-zero start abundance can have growth curves
  series.data = series.data[series.data[,start.col] > 0,]

  start.pop = series.data[,start.col]
  output.df = cbind(start.pop, row.names(series.data))
  colnames(output.df) = c("start.percent", "strain.name")
  rownames(output.df) = row.names(series.data)

  avg.df = 0
  for (i in 1:num.rep) {
    curr.start = 2+(i-1)*num.exp
    curr.end = curr.start + num.exp - 1
    curr.df = series.data[,curr.start:curr.end]
    curr.df = cbind(start.pop, curr.df)
    # log2 fold change relative to the start population
    curr.df = log2(curr.df / start.pop)
    output.df = cbind(output.df, curr.df)
    avg.df = avg.df + curr.df   # use.ctrl.avg is FALSE (spore.analysis = F)
  }
  avg.df = avg.df/num.rep
  avg.df[avg.df == 0] = 0.00001
  colnames(avg.df) = unlist(lapply(exp.names, function(x) paste(c(x, ".avg"), collapse="")))
  output.df = cbind(output.df, avg.df)

  growth.curve.dir = paste0(c(output.directory, "by.rep/growth.curve.data/"), collapse="")
  dir.create(growth.curve.dir, recursive=TRUE, showWarnings=FALSE)
  output.file = paste0(c(growth.curve.dir, "GC.abundance"), collapse="")
  write.table(output.df, file=paste0(c(output.file, "all.data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

  write(paste0("Program end at ", Sys.time()), output.info, append=TRUE)
  invisible(output.df)
}
