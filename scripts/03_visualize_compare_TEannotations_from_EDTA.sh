#!/bin/bash
#SBATCH --job-name=TE_sum_circos
#SBATCH --partition=pcmpg_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err

# ================= USER SETTINGS =================
# Project working directory (contains assembly fasta and EDTA results)
WORKDIR=/data/users/kxia/organize_annotation_course

# Input assembly (FASTA) and EDTA output directory
FASTA=$WORKDIR/assembly/assembly.fasta
EDTA_DIR=$WORKDIR/results/EDTA_annotation

# EDTA annotation products
SUM_FILE=$EDTA_DIR/assembly.fasta.mod.EDTA.TEanno.sum
GFF3=$EDTA_DIR/assembly.fasta.mod.EDTA.TEanno.gff3

# Output directory
OUTDIR=$WORKDIR/results/TE_summary_circos

# Circos (circlize) plotting parameters
TOP_SCAFFOLDS=6            # number of longest scaffolds as pseudo-chromosomes
WINDOW_SIZE=100000         # window size for density (bp)
R_PLOT_SCRIPT=/data/users/kxia/organize_annotation_course/scripts/03-annotation_circlize.R

# R runtime settings
USE_MODULE_R=true
R_MODULE=R/4.3.2-foss-2021a
R_USER_LIB=$WORKDIR/.Rlibs_R432
# =================================================

set -euo pipefail

mkdir -p "$OUTDIR" "$WORKDIR/logs"
cd "$OUTDIR"

echo "[I] Checking inputs..."
[[ -s "$GFF3" ]] || { echo "[E] Missing GFF3: $GFF3" >&2; exit 2; }
[[ -s "$FASTA" ]] || { echo "[E] Missing FASTA: $FASTA" >&2; exit 2; }

echo "[I] Ensuring FASTA index (.fai) exists..."
if [[ ! -s "${FASTA}.fai" ]]; then
  module load SAMtools || true
  if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$FASTA"
  else
    echo "[E] samtools not available. Please load it (module load SAMtools) or precompute $FASTA.fai" >&2
    exit 3
  fi
fi

# -------- Summarize TE content per superfamily --------
echo "[I] Summarizing TE content per superfamily..."

SUMMARY_TSV=TE_superfamily_summary.tsv
GENOME_LEN=$(awk '{s+=$2} END{print s+0}' "${FASTA}.fai")

if [[ -s "$SUM_FILE" ]]; then
  # Parse typical EDTA .sum lines; be case-insensitive.
  awk -v OFS='\t' '
    BEGIN{
      IGNORECASE=1
      pat="^(LTR/|TIR/|LINE|SINE|Helitron|PLE|Unknown)"
    }
    $0 ~ pat {
      name=$1; bp=""; pct=""
      for(i=1;i<=NF;i++){
        if($i ~ /^[0-9]+$/ && bp==""){ bp=$i }
        if($i ~ /\([0-9.]+%\)/){
          gsub(/[()]/,"",$i); split($i,a,"%"); pct=a[1]
        }
      }
      if(name!="" && bp!=""){
        if(pct==""){pct="NA"}
        print name, bp, pct
      }
    }
  ' "$SUM_FILE" > "$SUMMARY_TSV" || true
fi

# Fallback to GFF3-based calculation if no rows detected
if [[ ! -s "$SUMMARY_TSV" || ! -s "$SUM_FILE" || $(wc -l < "$SUMMARY_TSV") -lt 1 ]]; then
  echo "[W] Could not parse per-superfamily rows from .sum; computing from GFF3 instead."
  # Case-insensitive extraction of superfamily from attributes.
  awk -F'\t' -v OFS='\t' '
    BEGIN{ IGNORECASE=1 }
    function get_sf(attr,   sf){
      sf=""
      if (match(attr,/(^|;)classification=([^;]+)/,m)) sf=m[2]
      else if (match(attr,/(^|;)class=([^;]+)/,m)) sf=m[2]
      if (sf==""){
        if (match(attr,/(^|;)name=([^;]+)/,m)) sf=m[2]
      }
      gsub(/[, ].*$/,"",sf)
      return sf
    }
    $0 !~ /^#/ && $3!="long_terminal_repeat" && $3!="repeat_region" && $3!="target_site_duplication" {
      sf = get_sf($9)
      if (sf!=""){
        len = $5 - $4 + 1
        if (len>0) cov[sf]+=len
      }
    }
    END{
      for (k in cov){
        pct = (cov[k]>0 && '"$GENOME_LEN"'>0) ? (100.0*cov[k]/'"$GENOME_LEN"') : 0
        printf("%s\t%d\t%.4f\n", k, cov[k], pct)
      }
    }
  ' "$GFF3" | sort -k2,2nr > "$SUMMARY_TSV"
fi

echo "[I] Wrote per-superfamily summary: $SUMMARY_TSV"
echo "[I] Genome size (bp): $GENOME_LEN" > README.summary.txt
echo "[I] Columns: Superfamily\tCovered_bp\tPercent_of_genome" >> README.summary.txt

# -------- Prepare and run the circlize plotting --------
echo "[I] Preparing circlize plot..."
# Always (re)write the R script to ensure we use the robust version
mkdir -p "$(dirname "$R_PLOT_SCRIPT")"
cat > "$R_PLOT_SCRIPT" <<'RSCRIPT'
# circlize TE density plot (robust GFF3 parsing; no fread for GFF3)
# Inputs:
#   1) GFF3 path
#   2) FAI path (assembly.fasta.fai)
#   3) Output prefix
#   4) Top N scaffolds (integer)
#   5) Window size (bp)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(circlize)
  library(data.table)  # used only for FAI reading
})

# ---- Headless graphics setup (no X11 needed) ----
use_ragg <- requireNamespace("ragg", quietly = TRUE)
if (!use_ragg) {
  options(bitmapType = "cairo")  # fallback to cairo-based PNG
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript circos_TE_density.R <TEanno.gff3> <assembly.fai> <out_prefix> <top_scaffolds> <window_size>")
}
gff3        <- args[[1]]
fai         <- args[[2]]
out_prefix  <- args[[3]]
topN        <- as.integer(args[[4]])
win_size    <- as.integer(args[[5]])

# ---------- Helpers ----------
extract_superfamily <- function(attr) {
  # Case-insensitive extraction of classification/class/name, then trim tail after space/comma
  sf <- NA_character_
  m1 <- str_match(attr, "(?i)(^|;)classification=([^;]+)")
  if (!is.na(m1[1,3])) sf <- m1[1,3]
  if (is.na(sf)) {
    m2 <- str_match(attr, "(?i)(^|;)class=([^;]+)")
    if (!is.na(m2[1,3])) sf <- m2[1,3]
  }
  if (is.na(sf)) {
    m3 <- str_match(attr, "(?i)(^|;)name=([^;]+)")
    if (!is.na(m3[1,3])) sf <- m3[1,3]
  }
  if (!is.na(sf)) sf <- str_replace(sf, "[, ].*$", "")
  sf
}

# Strict 9-column GFF3 reader (no data.table::fread; ignore comments; force 9 columns)
read_gff3_strict9 <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(ln)]
  ln <- ln[!grepl("^#", ln)]  # drop "##..." headers and "###" separators
  if (length(ln) == 0) stop("No feature lines after removing comments.")

  sp <- strsplit(ln, "\t", fixed = TRUE)
  maxc <- max(vapply(sp, length, 1L))
  sp_pad <- lapply(sp, function(x) { length(x) <- maxc; x })
  g <- as.data.frame(do.call(rbind, sp_pad), stringsAsFactors = FALSE, optional = TRUE)

  if (ncol(g) < 9) g[, (ncol(g)+1):9] <- NA_character_
  if (ncol(g) > 9) g <- g[, 1:9, drop = FALSE]

  colnames(g) <- c("chr","src","type","start","end","score","strand","phase","attr")

  suppressWarnings({
    g$start <- as.integer(g$start)
    g$end   <- as.integer(g$end)
  })
  g <- g[!is.na(g$start) & !is.na(g$end) & g$end >= g$start, , drop = FALSE]

  # Drop sub-features we don’t want to count separately
  exclude_feats <- c("long_terminal_repeat","repeat_region","target_site_duplication")
  g <- g[!(g$type %in% exclude_feats), , drop = FALSE]

  # Compute length and Superfamily
  g$len <- as.integer(g$end - g$start + 1L)
  g$Superfamily <- vapply(g$attr, extract_superfamily, FUN.VALUE = character(1))
  g <- g[!is.na(g$Superfamily) & g$len > 0, , drop = FALSE]
  g
}

# ---------- Read FAI (scaffold lengths) ----------
ideogram <- data.table::fread(fai, header = FALSE, sep = "\t",
                              col.names = c("chr","len","off","line_bases","line_width"))
ideogram <- ideogram %>%
  arrange(desc(len)) %>%
  dplyr::slice(seq_len(min(n(), topN))) %>%
  transmute(chr, start = 0L, end = as.integer(len))
if (nrow(ideogram) == 0) stop("No scaffolds found in FAI.")

# ---------- Read GFF3 ----------
gff <- read_gff3_strict9(gff3)
gff <- dplyr::filter(gff, chr %in% ideogram$chr)

# ---- Clip coordinates to scaffold bounds and sanitize ----
len_map <- setNames(as.integer(ideogram$end), ideogram$chr)

gff$start[is.na(gff$start)] <- 0L
gff$end[is.na(gff$end)]     <- 0L
gff$start <- pmax(0L, gff$start)
gff$end <- mapply(function(ch, e) {
  L <- len_map[[ch]]
  if (!is.na(L)) min(as.integer(e), L) else as.integer(e)
}, gff$chr, gff$end)
gff <- dplyr::filter(gff, end >= start)
gff$len <- as.integer(gff$end - gff$start + 1L)
gff <- dplyr::filter(gff, len > 0)

message("Parsed GFF3 (strict9): ", nrow(gff), " features on ", length(unique(gff$chr)), " scaffolds.")

# ---------- Select top superfamilies by covered bp ----------
top_sf <- gff %>%
  group_by(Superfamily) %>%
  summarise(bp = sum(len), .groups = "drop") %>%
  arrange(desc(bp)) %>%
  dplyr::slice_head(n = 5) %>%
  pull(Superfamily)
if (length(top_sf) == 0) stop("No superfamilies detected for plotting.")
message("Selected superfamilies: ", paste(top_sf, collapse = ", "))

# Build list for circos.genomicDensity
sf_list <- lapply(top_sf, function(sf) {
  gff %>% filter(Superfamily == sf) %>% select(chr, start, end)
})
names(sf_list) <- top_sf

# ---------- Plot (PDF) ----------
pdf(paste0(out_prefix, "_TEdensity.pdf"), width = 9, height = 9)
circos.clear()
circos.par("track.height" = 0.12, start.degree = 90, gap.degree = 4)
circos.genomicInitialize(ideogram, sector.names = ideogram$chr, tickLabelsStartFromZero = TRUE)
for (sf in names(sf_list)) {
  circos.genomicDensity(sf_list[[sf]], col = NA, border = NA,
                        window.size = win_size, track.height = 0.12)
  circos.text(CELL_META$xcenter, CELL_META$ycenter + 0.25, sf, cex = 0.5,
              facing = "bending.inside", niceFacing = TRUE,
              sector.index = CELL_META$sector.index,
              track.index = get.current.track.index())
}
dev.off()

# ---------- Plot (PNG; headless) ----------
png_file <- paste0(out_prefix, "_TEdensity.png")
if (use_ragg) {
  ragg::agg_png(png_file, width = 2400, height = 2400, units = "px", res = 300)
} else {
  grDevices::png(png_file, width = 2400, height = 2400, res = 300, type = "cairo")
}
circos.clear()
circos.par("track.height" = 0.12, start.degree = 90, gap.degree = 4)
circos.genomicInitialize(ideogram, sector.names = ideogram$chr, tickLabelsStartFromZero = TRUE)
for (sf in names(sf_list)) {
  circos.genomicDensity(sf_list[[sf]], col = NA, border = NA,
                        window.size = win_size, track.height = 0.12)
  circos.text(CELL_META$xcenter, CELL_META$ycenter + 0.25, sf, cex = 0.5,
              facing = "bending.inside", niceFacing = TRUE,
              sector.index = CELL_META$sector.index,
              track.index = get.current.track.index())
}
dev.off()

cat("Done. Figures written to: ", paste0(out_prefix, "_TEdensity.[pdf|png]"), "\n")
RSCRIPT

# -------- R runtime bootstrap --------
if [[ "$USE_MODULE_R" == "true" ]]; then
  module purge
  module load "$R_MODULE"
  export R_LIBS_USER="$R_USER_LIB"
  mkdir -p "$R_LIBS_USER"
  # Ensure required packages (including ragg for headless PNG if available)
  Rscript -e 'pkgs <- c("data.table","dplyr","stringr","circlize","ragg");
              ip <- rownames(installed.packages());
              need <- setdiff(pkgs, ip);
              if (length(need)) install.packages(need, repos="https://cloud.r-project.org");
              for (p in pkgs) cat(p,":", if (p %in% ip) as.character(utils::packageVersion(p)) else "installed now", "\n")'
else
  module load Anaconda3
  eval "$(conda shell.bash hook)"
  conda activate "${CONDA_ENV_NAME:?Set CONDA_ENV_NAME or switch to module R}"
fi

# -------- Run the circlize plotting --------
OUT_PREFIX=$OUTDIR/circos_${TOP_SCAFFOLDS}scaff_win${WINDOW_SIZE}
echo "[I] Running circlize plot..."
Rscript "$R_PLOT_SCRIPT" "$GFF3" "${FASTA}.fai" "$OUT_PREFIX" "$TOP_SCAFFOLDS" "$WINDOW_SIZE"

echo "[I] All done."
echo "[I] Summary table: $SUMMARY_TSV"
echo "[I] Figures: ${OUT_PREFIX}_TEdensity.pdf / .png"