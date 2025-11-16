#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=InterProScan
#SBATCH --output=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.err

set -euo pipefail

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/kxia/organize_annotation_course/gene_annotation/final"
INP="${WORKDIR}/assembly.all.maker.proteins.fasta.renamed.fasta"

mkdir -p "$(dirname "$WORKDIR")/logs"

# Run InterProScan inside container (Pfam only)
apptainer exec \
  --bind $COURSEDIR/data/interproscan-5.70-102.0/data:/opt/interproscan/data \
  --bind $WORKDIR \
  --bind $COURSEDIR \
  --bind $SCRATCH:/temp \
  $COURSEDIR/containers/interproscan_latest.sif \
  /opt/interproscan/interproscan.sh \
  -appl pfam --disable-precalc -f TSV \
  --goterms --iprlookup --seqtype p \
  -i "$INP" -o "${WORKDIR}/output.iprscan"

echo "[OK] InterProScan -> ${WORKDIR}/output.iprscan"