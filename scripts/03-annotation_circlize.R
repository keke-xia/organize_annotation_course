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
  make_option("--gff", type="character", help="EDTA TE annotation GFF3 (e.g. assembly.mod.EDTA.TEanno.gff3)"),
  make_option("--fai", type="character", help="FASTA index (.fai) of the assembly"),
  make_option("--outdir", type="character", default=NULL, help="Output dir [default: <gff_dir>/circos]"),
  make_option("--prefix", type="character", default="circos_TE", help="Output file prefix"),
  make_option("--scaf_top", type="integer", default=10, help="Top-K longest scaffolds to plot"),
  make_option("--win", type="integer", default=100000, help="Window size for density"),
  make_option("--super", type="character", default="LTR/Gypsy,LTR/Copia,DNA/TIR,LINE,SINE",
              help="Comma separated superfamilies (prefix match), e.g. 'LTR/Gypsy,LTR/Copia'")
)
opt <- parse_args(OptionParser(option_list=option_list))

# -------- validate inputs --------
if (is.null(opt$gff) || is.null(opt$fai)) {
  stop("Missing --gff or --fai. Usage: --gff <file.gff3> --fai <assembly.fai>")
}
if (!file.exists(opt$gff)) stop(paste("GFF not found:", opt$gff))
if (!file.exists(opt$fai)) stop(paste("FAI not found:", opt$fai))

if (is.null(opt$outdir)) opt$outdir <- file.path(dirname(opt$gff), "circos")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

message("Selected tracks: ", opt$super)

# -------- read FAI (scaffold lengths) --------
fai_dt <- fread(opt$fai, header = FALSE, sep = "\t",
                col.names = c("seqnames","length","off","line_blen","qual"))
fai_dt <- fai_dt %>% arrange(desc(length))
if (nrow(fai_dt) == 0) stop("Empty .fai file: ", opt$fai)

if (opt$scaf_top > 0 && opt$scaf_top < nrow(fai_dt)) {
  fai_dt <- fai_dt %>% slice_head(n = opt$scaf_top)
}
seq_order <- fai_dt$seqnames

# ideogram for circos
ideogram_df <- fai_dt %>% transmute(chr = seqnames, start = 0L, end = as.integer(length))

# -------- read EDTA GFF3 (skip comments) --------
gff_lines <- readLines(opt$gff, warn = FALSE)
gff_lines <- gff_lines[!grepl("^\\s*$", gff_lines)]         # drop empty lines
gff_lines <- gff_lines[!grepl("^#", gff_lines)]             # drop comment lines

if (length(gff_lines) == 0) stop("No non-comment lines in GFF: ", opt$gff)

gff_dt <- fread(text = gff_lines, sep = "\t", header = FALSE, fill = TRUE, quote = "")
# standard GFF3 has 9 columns; keep at least 5, at most 9
if (ncol(gff_dt) < 5) stop("GFF has too few columns (<5) after skipping comments.")
if (ncol(gff_dt) > 9) gff_dt <- gff_dt[, 1:9]
cn <- c("seqid","source","type","start","end","score","strand","phase","attr")
setnames(gff_dt, cn[1:ncol(gff_dt)])
if (!("attr" %in% names(gff_dt))) gff_dt$attr <- ""

# 只保留要画的 scaffold
gff_dt <- gff_dt[seqid %in% seq_order]
if (nrow(gff_dt) == 0) stop("No GFF entries on the selected scaffolds (top ", length(seq_order), ").")

# -------- parse superfamily from attributes --------
parse_super <- function(a) {
  a <- as.character(a)
  # try "Classification=" like "Classification=LTR/Gypsy"
  m3 <- str_match(a, "(?i)Classification[=:]([^;\\s]+)")[,2]
  if (!is.na(m3)) return(m3)
  # try "Class=" and "Superfamily="
  m1 <- str_match(a, "(?i)Class[=:]([^;\\s]+)")[,2]
  m2 <- str_match(a, "(?i)Superfamily[=:]([^;\\s]+)")[,2]
  if (!is.na(m1) && !is.na(m2)) return(paste0(m1,"/",m2))
  if (!is.na(m1)) return(m1)
  if (!is.na(m2)) return(m2)
  # last resort: direct pattern like LTR/Gypsy
  m4 <- str_match(a, "([A-Za-z]+/[A-Za-z0-9_]+)")[,2]
  if (!is.na(m4)) return(m4)
  return("UNKNOWN")
}
sup <- vapply(gff_dt$attr, parse_super, FUN.VALUE = character(1))
sup <- gsub("%2F","/", sup, fixed = TRUE)
sup <- gsub("_INT$", "", sup)
sup <- toupper(sup)
gff_dt$superfamily <- sup

# filter requested superfamilies (prefix match)
sel_supers <- str_split(opt$super, ",")[[1]] %>% str_trim() %>% toupper()
keep <- Reduce(`|`, lapply(sel_supers, function(p) startsWith(gff_dt$superfamily, p)))
gff_sel <- gff_dt[keep, .(seqid, start = pmin(as.integer(start), as.integer(end)),
                          end = pmax(as.integer(start), as.integer(end)),
                          superfamily)]
if (nrow(gff_sel) == 0) stop("No entries match requested superfamilies: ", paste(sel_supers, collapse=", "))

# 确保位点在 scaffold 长度范围内
len_map <- setNames(ideogram_df$end, ideogram_df$chr)
gff_sel <- gff_sel[seqid %in% names(len_map)]
gff_sel$start[gff_sel$start < 0] <- 0L
gff_sel$end <- pmin(gff_sel$end, len_map[gff_sel$seqid])

# -------- plotting (PDF & PNG) --------
pdf_file <- file.path(opt$outdir, sprintf("%s_density_top%d_win%d.pdf", opt$prefix, nrow(ideogram_df), opt$win))
png_file <- sub("\\.pdf$", ".png", pdf_file)

plot_once <- function(out_dev = c("pdf","png")) {
  out_dev <- match.arg(out_dev)
  if (out_dev == "pdf") {
    pdf(pdf_file, width = 8, height = 8)
  } else {
    png(png_file, width = 1200, height = 1200, res = 150)
  }
  on.exit(dev.off(), add = TRUE)

  circos.clear()
  circos.par(gap.after = c(rep(2, nrow(ideogram_df)-1), 8),
             start.degree = 90, track.margin = c(0.01, 0.01))
  circos.genomicInitialize(ideogram_df, sector.names = ideogram_df$chr,
                           tickLabelsStartFromZero = TRUE)

  # color map for requested tracks (stable palette)
  track_names <- unique(sel_supers)
  palette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999")
  col_map <- setNames(palette[seq_along(track_names)], track_names)

  for (nm in track_names) {
    df_nm <- gff_sel %>%
      filter(startsWith(superfamily, nm)) %>%
      transmute(chr = seqid, start = as.integer(start), end = as.integer(end))
    if (nrow(df_nm) == 0) next
    circos.genomicDensity(df_nm, col = col_map[[nm]], track.height = 0.12,
                          window.size = opt$win, overlap = TRUE)
  }
}

# -------- plotting (PDF only) --------
pdf_file <- file.path(opt$outdir, sprintf("%s_density_top%d_win%d.pdf", opt$prefix, nrow(ideogram_df), opt$win))

plot_once <- function(out_dev = c("pdf")) {
  pdf(pdf_file, width = 8, height = 8)
  on.exit(dev.off(), add = TRUE)

  circos.clear()
  circos.par(gap.after = c(rep(2, nrow(ideogram_df)-1), 8),
             start.degree = 90, track.margin = c(0.01, 0.01))
  circos.genomicInitialize(ideogram_df, sector.names = ideogram_df$chr,
                           tickLabelsStartFromZero = TRUE)

  track_names <- unique(sel_supers)
  palette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999")
  col_map <- setNames(palette[seq_along(track_names)], track_names)

  for (nm in track_names) {
    df_nm <- gff_sel %>%
      filter(startsWith(superfamily, nm)) %>%
      transmute(chr = seqid, start = as.integer(start), end = as.integer(end))
    if (nrow(df_nm) == 0) next
    circos.genomicDensity(df_nm, col = col_map[[nm]], track.height = 0.12,
                          window.size = opt$win, overlap = TRUE)
  }
}

# only make pdf
plot_once("pdf")

# export a summary table
sum_tab <- gff_sel %>%
  mutate(fam = gsub("^([A-Z]+)/.*$", "\\1", superfamily)) %>%
  count(superfamily, name = "n_copies") %>%
  arrange(desc(n_copies))
fwrite(sum_tab, file.path(opt$outdir, sprintf("%s_density_counts.tsv", opt$prefix)), sep = "\t")

message("Saved: ", pdf_file)
message("Saved: ", file.path(opt$outdir, sprintf("%s_density_counts.tsv", opt$prefix)))