#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=GFF_update_iprscan
#SBATCH --output=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.err
set -euo pipefail
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"
WORKDIR=/data/users/kxia/organize_annotation_course/gene_annotation/final

cd "$WORKDIR"
gff="/data/users/kxia/organize_annotation_course/gene_annotation/final/assembly.all.maker.noseq.gff.renamed.gff"

# Merge iprscan annotations into GFF
"$MAKERBIN/ipr_update_gff" "$gff" output.iprscan > "${gff}.iprscan.gff"
echo "[OK] Updated GFF -> ${gff}.iprscan.gff"

#calculate AED
perl "$MAKERBIN/AED_cdf_generator.pl" -b 0.025 "$gff" > assembly.all.maker.renamed.gff.AED.txt
echo "[OK] AED -> assembly.all.maker.renamed.gff.AED.txt"


QP="${gff}.iprscan.gff"

perl "$MAKERBIN/quality_filter.pl" -s "$QP" > "${gff}_iprscan_quality_filtered.gff"
echo "[OK] Filtered -> ${gff}_iprscan_quality_filtered.gff"

IN="${gff}_iprscan_quality_filtered.gff"

# Keep key gene features only
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" \
  "$IN" > filtered.genes.renamed.gff3

# Quick check of feature types
cut -f3 filtered.genes.renamed.gff3 | sort | uniq

echo "[OK] Feature-filtered -> filtered.genes.renamed.gff3"