#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(circlize)
})

# -------- options --------
option_list <- list(
  make_option("--gene_gff", type="character", help="Gene annotation GFF3 (e.g. filtered.genes.renamed.gff3)"),
  make_option("--te_gff",   type="character", help="TE annotation GFF3 from EDTA (e.g. assembly.mod.EDTA.TEanno.gff3)"),
  make_option("--fai",      type="character", help="FASTA index (.fai) of the assembly"),
  make_option("--outdir",   type="character", default=NULL, help="Output dir [default: <gene_gff_dir>/circos]"),
  make_option("--prefix",   type="character", default="circos_gene_TE", help="Output file prefix"),
  make_option("--scaf_top", type="integer", default=10, help="Top-K longest scaffolds to plot"),
  make_option("--win",      type="integer", default=100000, help="Window size for density")
)
opt <- parse_args(OptionParser(option_list=option_list))

# -------- validate inputs --------
if (is.null(opt$gene_gff) || is.null(opt$te_gff) || is.null(opt$fai)) {
  stop("Missing --gene_gff, --te_gff or --fai. Usage: --gene_gff <file.gff3> --te_gff <file.gff3> --fai <assembly.fai>")
}
if (!file.exists(opt$gene_gff)) stop(paste("Gene GFF not found:", opt$gene_gff))
if (!file.exists(opt$te_gff))   stop(paste("TE GFF not found:", opt$te_gff))
if (!file.exists(opt$fai))      stop(paste("FAI not found:", opt$fai))

if (is.null(opt$outdir)) opt$outdir <- file.path(dirname(opt$gene_gff), "circos")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -------- read FAI (scaffold lengths) --------
fai_dt <- fread(opt$fai, header = FALSE, sep = "\t",
                col.names = c("seqnames","length","off","line_blen","qual"))
fai_dt <- fai_dt %>% arrange(desc(length))
if (nrow(fai_dt) == 0) stop("Empty .fai file: ", opt$fai)

if (opt$scaf_top > 0 && opt$scaf_top < nrow(fai_dt)) {
  fai_dt <- fai_dt %>% slice_head(n = opt$scaf_top)
}
seq_order <- fai_dt$seqnames

# Ideogram data frame for circos
ideogram_df <- fai_dt %>%
  transmute(chr = seqnames, start = 0L, end = as.integer(length))

# -------- helper to read GFF (non-comment lines) --------
read_gff_no_comments <- function(path) {
  # Read GFF3 file, skip empty and comment lines, keep first 9 columns
  gff_lines <- readLines(path, warn = FALSE)
  gff_lines <- gff_lines[!grepl("^\\s*$", gff_lines)]
  gff_lines <- gff_lines[!grepl("^#", gff_lines)]
  if (length(gff_lines) == 0) stop("No non-comment lines in GFF: ", path)
  dt <- fread(text = gff_lines, sep = "\t", header = FALSE, fill = TRUE, quote = "")
  if (ncol(dt) < 5) stop("GFF has too few columns (<5): ", path)
  if (ncol(dt) > 9) dt <- dt[, 1:9]
  cn <- c("seqid","source","type","start","end","score","strand","phase","attr")
  setnames(dt, cn[1:ncol(dt)])
  if (!("attr" %in% names(dt))) dt$attr <- ""
  dt
}

# -------- read gene GFF --------
gene_dt <- read_gff_no_comments(opt$gene_gff)
gene_dt <- gene_dt[seqid %in% seq_order]

# Use 'gene' features as intervals for gene density; if missing, fallback to mRNA
if (!any(gene_dt$type == "gene")) {
  message("No 'gene' features found; using 'mRNA' as gene-level features.")
  gene_dt <- gene_dt[type == "mRNA"]
} else {
  gene_dt <- gene_dt[type == "gene"]
}

if (nrow(gene_dt) == 0) stop("No gene-level features on selected scaffolds.")

gene_df <- gene_dt %>%
  transmute(chr   = seqid,
            start = pmin(as.integer(start), as.integer(end)),
            end   = pmax(as.integer(start), as.integer(end)))

# Clip gene coordinates to scaffold boundaries
len_map <- setNames(ideogram_df$end, ideogram_df$chr)
gene_df <- gene_df %>%
  filter(chr %in% names(len_map)) %>%
  mutate(start = ifelse(start < 0, 0L, start),
         end   = pmin(end, len_map[chr]))

# -------- read TE GFF --------
te_dt <- read_gff_no_comments(opt$te_gff)
te_dt <- te_dt[seqid %in% seq_order]

if (nrow(te_dt) == 0) stop("No TE features on selected scaffolds.")

te_df <- te_dt %>%
  transmute(chr   = seqid,
            start = pmin(as.integer(start), as.integer(end)),
            end   = pmax(as.integer(start), as.integer(end)))

# Clip TE coordinates to scaffold boundaries
te_df <- te_df %>%
  filter(chr %in% names(len_map)) %>%
  mutate(start = ifelse(start < 0, 0L, start),
         end   = pmin(end, len_map[chr]))

# -------- plotting (PDF only) --------
pdf_file <- file.path(
  opt$outdir,
  sprintf("%s_gene_TE_density_top%d_win%d.pdf", opt$prefix, nrow(ideogram_df), opt$win)
)

# Single plotting function
plot_gene_te_density <- function() {
  pdf(pdf_file, width = 8, height = 8)
  on.exit(dev.off(), add = TRUE)

  # Reset circlize state
  circos.clear()

  # Circos parameters
  circos.par(
    gap.after = c(rep(2, nrow(ideogram_df) - 1), 8),
    start.degree = 90,
    track.margin = c(0.01, 0.01)
  )

  # Initialize genome ideogram
  circos.genomicInitialize(
    ideogram_df,
    sector.names = ideogram_df$chr,
    tickLabelsStartFromZero = TRUE
  )

  # Colors: one for genes, one for TEs
  col_gene <- "#0072B2"  # blue-ish
  col_te   <- "#D55E00"  # orange-ish

  # TE density track (outer)
  circos.genomicDensity(
    te_df,
    col = col_te,
    track.height = 0.15,
    window.size = opt$win,
    overlap = TRUE
  )

  # Gene density track (inner)
  circos.genomicDensity(
    gene_df,
    col = col_gene,
    track.height = 0.15,
    window.size = opt$win,
    overlap = TRUE
  )
}

plot_gene_te_density()

# -------- export summary (optional small table) --------
sum_tab <- data.frame(
  track = c("gene", "TE"),
  n_intervals = c(nrow(gene_df), nrow(te_df))
)
fwrite(
  sum_tab,
  file.path(opt$outdir, sprintf("%s_gene_TE_density_counts.tsv", opt$prefix)),
  sep = "\t"
)

message("Saved: ", pdf_file)
message("Saved: ", file.path(opt$outdir, sprintf("%s_gene_TE_density_counts.tsv", opt$prefix)))