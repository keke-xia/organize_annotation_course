#!/bin/bash
#SBATCH --job-name=TE_annotation
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=2-00:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err

WORKDIR=/data/users/kxia/organize_annotation_course
INPUT_FASTA=$WORKDIR/assembly/assembly.fasta
OUTDIR=$WORKDIR/results/EDTA_annotation

CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif

CDS=/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated

# ================== start==================

mkdir -p "$OUTDIR" "$WORKDIR/logs"
cd "$OUTDIR" || exit 1

#clean environment variables for Apptainer
export BASH_ENV=/dev/null
export ENV=/dev/null

#temp directory for EDTA
export TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR"

# Run EDTA
apptainer exec \
  --cleanenv \
  --bind "$WORKDIR:$WORKDIR","$OUTDIR:$OUTDIR","/data:/data" \
  "$CONTAINER" \
  bash -lc "
    set -euo pipefail
    cd '$OUTDIR'
    EDTA.pl \
      --genome '$INPUT_FASTA' \
      --species others \
      --step all \
      --sensitive 1 \
      --cds '$CDS' \
      --anno 1 \
      --force 1 \
      --threads $SLURM_CPUS_PER_TASK
  "

echo '✅ EDTA done, saved at:'
echo "   $OUTDIR"