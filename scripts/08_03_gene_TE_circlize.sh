#!/bin/bash
#SBATCH --job-name=circlize_gene_TE
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

# ===== Project paths =====
WORKDIR="/data/users/kxia/organize_annotation_course"

# Gene annotation (MAKER final filtered genes)
GENE_GFF="$WORKDIR/gene_annotation/final/filtered.genes.renamed.gff3"

# TE annotation (EDTA)
TE_GFF="$WORKDIR/results/EDTA_annotation/assembly.fasta.mod.EDTA.TEanno.gff3"

# Assembly FASTA and index (.fai)
FASTA="$WORKDIR/assembly/assembly.fasta"
FAI="${FASTA}.fai"

# Output directory for circos/circlize figures
CIRCDIR="$WORKDIR/results/Gene_TE_circos"
mkdir -p "$CIRCDIR"

# ===== Ensure FASTA index exists (samtools faidx) =====
if [[ ! -s "$FAI" ]]; then
  module spider SAMtools >/dev/null 2>&1 || true
  module load SAMtools || true
  if command -v samtools >/dev-null 2>&1; then
    samtools faidx "$FASTA"
  else
    echo "WARNING: samtools not available; please create $FAI manually: samtools faidx $FASTA"
    exit 1
  fi
fi

echo ">>> Gene GFF : $GENE_GFF"
echo ">>> TE GFF   : $TE_GFF"
echo ">>> FASTA    : $FASTA"
echo ">>> FAI      : $FAI"
echo ">>> OUTDIR   : $CIRCDIR"

# ===== Run the circlize plotting script =====
Rscript "$WORKDIR/scripts/08-gene_TE_circlize.R" \
  --gene_gff "$GENE_GFF" \
  --te_gff "$TE_GFF" \
  --fai "$FAI" \
  --outdir "$CIRCDIR" \
  --prefix "assembly_gene_TE" \
  --scaf_top 10 \
  --win 100000

echo ">>> Done."