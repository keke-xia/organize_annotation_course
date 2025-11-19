#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --job-name=UniProt_blastp
#SBATCH --output=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/gene_annotation/logs/%x_%j.err
set -euo pipefail

# Load BLAST+
module load BLAST+/2.15.0-gompi-2021a

# Paths
WORK=/data/users/kxia/organize_annotation_course
FINAL=$WORK/gene_annotation/final
QUERY=$FINAL/proteins.renamed.filtered.fasta

# ---------------- UniProt BLAST ----------------
DB=/data/courses/assembly-annotation-course/CDS_annotation/data/uniprot/uniprot_viridiplantae_reviewed.fa
OUT1=$FINAL/uniprot_vs_proteins.tsv

# Run blastp against UniProt (makeblastdb is already done)
blastp -query "$QUERY" -db "$DB" \
  -num_threads 10 -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "$OUT1"

# Keep only best hit per query (current criterion: smallest bitscore column #12)
# NOTE: if you prefer lowest evalue, change to: -k11,11g
sort -k1,1 -k12,12g "$OUT1" | sort -u -k1,1 --merge > "${OUT1}.besthits"

# ---------------- Map UniProt functions back to GFF / FASTA ----------------
MAKERBIN=/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin
UNIPROT=/data/courses/assembly-annotation-course/CDS_annotation/data/uniprot/uniprot_viridiplantae_reviewed.fa
BEST=$FINAL/uniprot_vs_proteins.tsv.besthits

PROT=$FINAL/proteins.renamed.filtered.fasta
GFF=$FINAL/filtered.genes.renamed.gff3

# Make copies that will be overwritten with functional annotations
cp -f "$PROT" "$FINAL/proteins.renamed.filtered.fasta.Uniprot"
cp -f "$GFF"  "$FINAL/filtered.genes.renamed.gff3.Uniprot.gff3"

# Write functional annotations into fasta/gff
"$MAKERBIN/maker_functional_fasta" "$UNIPROT" "$BEST" "$PROT" > "$FINAL/proteins.renamed.filtered.fasta.Uniprot"
"$MAKERBIN/maker_functional_gff"   "$UNIPROT" "$BEST" "$GFF"  > "$FINAL/filtered.genes.renamed.gff3.Uniprot.gff3"

echo "[OK] UniProt annotations written to:"
echo " - $FINAL/proteins.renamed.filtered.fasta.Uniprot"
echo " - $FINAL/filtered.genes.renamed.gff3.Uniprot.gff3"

# ---------------- TAIR10 BLAST ----------------
# Now, get the best blast hit with Arabidopsis thaliana TAIR10 representative gene models:

# Use TAIR10 peptide DB here (this was commented before, causing the bug)
DB=/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_pep_20110103_representative_gene_model
OUT2=$FINAL/tair10_vs_proteins.tsv

# blastp to TAIR10
blastp -query "$QUERY" -db "$DB" \
  -num_threads 10 -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "$OUT2"

# Keep best hit per query (same criterion as above to stay consistent)
sort -k1,1 -k12,12g "$OUT2" | sort -u -k1,1 --merge > "${OUT2}.besthits"

echo "[OK] TAIR besthits -> ${OUT2}.besthits"