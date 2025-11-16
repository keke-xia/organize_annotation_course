#run in terminal
WORK=/data/users/kxia/organize_annotation_course
FINAL=$WORK/gene_annotation/final
ACC=Pa-1                             # <<< set your accession code
GS=$WORK/genespace_work

mkdir -p "$GS/bed" "$GS/peptide"

# 1) Extract gene features to BED (0-based start)
grep -P "\tgene\t" "$FINAL/filtered.genes.renamed.gff3" > "$GS/temp_genes.gff3"
awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
  "$GS/temp_genes.gff3" > "$GS/bed/${ACC}.bed"
rm -f "$GS/temp_genes.gff3"

# 2) Prepare peptide fasta (headers must match gene names: strip "-R*")
awk '/^>/{sub(/-R.*/,"",$0)}1' "$FINAL/proteins.renamed.longest.fasta" > "$GS/peptide/${ACC}.fa"

# 3) Copy TAIR10 reference files
cp /data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10.bed "$GS/bed/TAIR10.bed"
cp /data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10.fa  "$GS/peptide/TAIR10.fa"


# Prepare GENESPACE inputs (BED + peptide FASTA) from Lian_et_al GFFs
# - BED is built from GFF 'gene' features (0-based start).
# - Peptide FASTA is copied from the course 'protein' folder if present.
#   (Headers in the provided *.fa already match gene names.)


# --- Paths (edit if needed) ---
COURSE_BASE=/data/courses/assembly-annotation-course/CDS_annotation/data/Lian_et_al
GFF_DIR=/data/courses/assembly-annotation-course/CDS_annotation/data/Lian_et_al/gene_gff/selected/
PROT_DIR=$COURSE_BASE/protein/selected

# Your GENESPACE working directory
GS=/data/users/kxia/organize_annotation_course/genespace_work
mkdir -p "$GS/bed" "$GS/peptide"

# --- Process each accession ---
#acc="Altai-5"   # e.g., Altai-5 from Altai-5.EVM.v3.5.ann.protein_coding_genes.gff
acc="Ice-1"
gff="$GFF_DIR/${acc}.EVM.v3.5.ann.protein_coding_genes.gff"
echo "[I] Processing accession: $acc"


  # 1) Extract gene features to BED (0-based start; 4th column = gene ID from attributes)
  #    Keep only lines where the third column is 'gene'.
  grep -P "\tgene\t" "$gff" > "$GS/bed/.tmp_${acc}.gff3"
  awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
      "$GS/bed/.tmp_${acc}.gff3" > "$GS/bed/${acc}.bed"
  rm -f "$GS/bed/.tmp_${acc}.gff3"

  # 2) Peptide FASTA: prefer course-provided file if it exists
  if [[ -s "$PROT_DIR/${acc}.protein.faa" ]]; then
    cp -f "$PROT_DIR/${acc}.protein.faa" "$GS/peptide/${acc}.fa"
  else
    echo "[W] Missing peptide FASTA for $acc at: $PROT_DIR/${acc}.protein.faa"
    echo "    Please check the course 'protein' folder or provide a corresponding ${acc}.protein.faa"
  fi

