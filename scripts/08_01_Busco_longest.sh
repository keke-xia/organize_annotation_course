#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=busco_longest
#SBATCH --output=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.err
set -euo pipefail

module load BUSCO/5.4.2-foss-2021a

WORKDIR=/data/users/kxia/organize_annotation_course/gene_annotation/final
INP=$WORKDIR/proteins.renamed.longest.fasta
LINEAGE=brassicales_odb10          # or embryophyta_odb10
OUTNAME=busco_longest_output

busco -i "$INP" -l "$LINEAGE" -o "$OUTNAME" -m proteins
echo "[✓] BUSCO finished -> $WORKDIR/$OUTNAME"

