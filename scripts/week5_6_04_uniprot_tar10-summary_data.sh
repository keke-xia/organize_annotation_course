# ---- paths (adapt if needed) ----
FINAL=/data/users/kxia/organize_annotation_course/gene_annotation/final
GFF=$FINAL/filtered.genes.renamed.gff3
UNI_BEST=$FINAL/uniprot_vs_proteins.tsv.besthits
TAIR_BEST=$FINAL/tair10_vs_proteins.tsv.besthits
PROT_FASTA=$FINAL/proteins.renamed.filtered.fasta   # BLAST 使用的 query

# 1) total number of proteins (number of query sequences)
total_prot=$(grep -c "^>" "$PROT_FASTA")

# 2) number of proteins with UniProt best hit (unique query IDs in UniProt besthits)
uniprot_with_hit=$(cut -f1 "$UNI_BEST" | sort -u | wc -l)

# 3) number of proteins with TAIR10 best hit (unique query IDs in TAIR10 besthits)
tair_with_hit=$(cut -f1 "$TAIR_BEST" | sort -u | wc -l)

# 4) proteins without hits
uniprot_without_hit=$(( total_prot - uniprot_with_hit ))
tair_without_hit=$(( total_prot - tair_with_hit ))

# ---- print summary ----
echo "Total proteins (queries in BLAST): $total_prot"
echo "Proteins with UniProt hits:        $uniprot_with_hit"
echo "Proteins without UniProt hits:     $uniprot_without_hit"
echo "Proteins with TAIR10 hits:         $tair_with_hit"
echo "Proteins without TAIR10 hits:      $tair_without_hit"