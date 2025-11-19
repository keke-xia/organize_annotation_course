#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GENESPACE))

# Get working directory from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: genespace.R <workingDirectory>")
wd <- args[1]

# Initialize GENESPACE parameters
gpar <- init_genespace(
  wd          = wd,
  path2mcscanx = "/usr/local/bin",  # <-- fixed: use directory, not /MCScanX
  verbose     = TRUE,
  nCores      = 20
)

# Run GENESPACE
out <- run_genespace(gpar, overwrite = TRUE)

# Build pangenome object
pangenome <- query_pangenes(
  out,
  bed          = NULL,
  refGenome    = "TAIR10",
  transform    = TRUE,
  showArrayMem = TRUE,
  showNSOrtho  = TRUE,
  maxMem2Show  = Inf
)

# Save pangenome matrix as RDS
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))