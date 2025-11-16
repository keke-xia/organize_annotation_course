#run in terminal
module load UCSC-Utils/448-foss-2021a

WORKDIR=/data/users/kxia/organize_annotation_course/gene_annotation/final
cd "$WORKDIR"

# Make ID list from filtered GFF
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' > list.txt

# Filter FASTA by ID list
faSomeRecords assembly.all.maker.transcripts.fasta.renamed.fasta list.txt \
  transcripts.renamed.filtered.fasta

faSomeRecords assembly.all.maker.proteins.fasta.renamed.fasta list.txt \
  proteins.renamed.filtered.fasta

echo "[OK] Wrote transcripts.renamed.filtered.fasta and proteins.renamed.filtered.fasta"