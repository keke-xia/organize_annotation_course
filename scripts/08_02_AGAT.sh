conda activate agat
# --- Paths ---
WORKDIR=/data/users/kxia/organize_annotation_course/gene_annotation/final
cd "$WORKDIR"

GFF=filtered.genes.renamed.gff3
OUT=annotation.stat

# --- Run AGAT statistics ---
agat_sp_statistics.pl -i "$GFF" -o "$OUT"

echo "[✓] AGAT statistics -> $WORKDIR/$OUT"