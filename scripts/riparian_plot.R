#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
})

# ----- get working directory from command line -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: riparian_plot.R <genespace_work_dir>")
}
wd <- args[1]

# ----- load or initialize GENESPACE parameters -----
# Try to reuse saved gsParam; if not found, re-init from wd
param_rds <- file.path(wd, "gsParam.rds")
param_rda <- file.path(wd, "gsParam.rda")

if (file.exists(param_rds)) {
  # gsParam stored as RDS
  gsParam <- readRDS(param_rds)
} else if (file.exists(param_rda)) {
  # gsParam stored as RDA (object name usually gsParam)
  load(param_rda)
  if (!exists("gsParam")) {
    stop("gsParam.rda loaded but object 'gsParam' not found.")
  }
} else {
  # Fall back to init_genespace: reuse existing Orthofinder/synteny results
  gsParam <- init_genespace(
    wd = wd,
    path2mcscanx = "/usr/local/bin",  # fixed path from course notes
    verbose = TRUE
  )
}

# ----- make riparian plot -----
# Use TAIR10 as reference; order genomes similar to course setting
out_pdf <- file.path(wd, "riparian_TAIR10_ref.pdf")

ripDat <- plot_riparian(
  gsParam   = gsParam,
  refGenome = "TAIR10",                     # reference genome
  genomeIDs = c("TAIR10", "Altai_5", "Ice_1", "Pa_1"),  # order in plot
  useOrder  = TRUE,                         # use gene rank order
  useRegions = TRUE,                        # aggregated syntenic regions
  minChrLen2plot = 100,                     # skip tiny scaffolds
  forceRecalcBlocks = TRUE,                 # recompute blocks using current params
  pdfFile = out_pdf                         # write to PDF
)

message("Riparian plot saved to: ", out_pdf)