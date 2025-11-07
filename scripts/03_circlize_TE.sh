#!/bin/bash
#SBATCH --job-name=circlize_TE
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err
set -euo pipefail

# ===== Load R environment =====
module load R/4.3.2-foss-2021a
export R_LIBS_USER="/data/users/kxia/Rlibs/4.3"
mkdir -p "$R_LIBS_USER"

# Install required R packages into user lib (idempotent)
Rscript -e 'repos <- c(CRAN="https://cloud.r-project.org");
            pkgs <- c("optparse","data.table","dplyr","stringr","circlize");
            inst <- rownames(installed.packages(lib.loc=Sys.getenv("R_LIBS_USER")));
            need <- setdiff(pkgs, inst);
            if(length(need)) install.packages(need, repos=repos, lib=Sys.getenv("R_LIBS_USER"), Ncpus=2)'

# ===== Project paths (adjusted to your layout) =====
WORKDIR="/data/users/kxia/organize_annotation_course"
EDTA_DIR="$WORKDIR/results/EDTA_annotation"

# GFF3 from EDTA all-TE annotation
GFF="$EDTA_DIR/assembly.fasta.mod.EDTA.TEanno.gff3"

# Assembly FASTA and index (.fai) under your assembly folder
FASTA="$WORKDIR/assembly/assembly.fasta"
FAI="${FASTA}.fai"

# Output directory for circos/circlize figures
CIRCDIR="$WORKDIR/results/TE_summary_circos"
mkdir -p "$CIRCDIR"

# ===== Ensure FASTA index exists (samtools faidx) =====
if [[ ! -s "$FAI" ]]; then
  module spider SAMtools >/dev/null 2>&1 || true
  module load SAMtools || true
  if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$FASTA"
  else
    echo "WARNING: samtools not available; please create $FAI manually: samtools faidx $FASTA"
    exit 1
  fi
fi

echo ">>> GFF : $GFF"
echo ">>> FASTA: $FASTA"
echo ">>> FAI : $FAI"
echo ">>> OUT : $CIRCDIR"

# ===== Run the circlize plotting script =====
# NOTE:
#   --gff expects the EDTA TEanno.gff3
#   --fai expects the FASTA index produced by `samtools faidx`
#   --scaf_top is how many longest scaffolds (pseudo-chromosomes) to show
#   --win is the density window size in bp
#   --super is a comma-separated list of superfamilies to plot (optional, depends on your R script)
Rscript "$WORKDIR/scripts/03-annotation_circlize.R" \
  --gff "$GFF" \
  --fai "$FAI" \
  --outdir "$CIRCDIR" \
  --prefix "assembly" \
  --scaf_top 10 \
  --win 100000 \
  --super "LTR/Gypsy,LTR/Copia,DNA/TIR,LINE,SINE"

echo ">>> Done."