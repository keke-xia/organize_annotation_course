#!/usr/bin/env Rscript
# Initialize directories, run GENESPACE, and save the pangenome matrix as RDS.
suppressPackageStartupMessages({ library(GENESPACE) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: genespace.R <workingDirectory>")
wd <- args[1]


gpar <- init_genespace(
  wd = wd,
  path2mcscanx     = "/usr/local/bin/MCScanX",
  path2orthofinder = "/usr/local/bin/orthofinder"
)

out <- run_genespace(gpar, overwrite = TRUE)

pangenome <- query_pangenes(
  out, bed = NULL, refGenome = "TAIR10",
  transform = TRUE, showArrayMem = TRUE,
  showNSOrtho = TRUE, maxMem2Show = Inf
)

saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))