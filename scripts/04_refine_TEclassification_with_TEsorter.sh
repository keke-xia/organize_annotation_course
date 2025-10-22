#!/bin/bash
#SBATCH --job-name=TEsorter_refine
#SBATCH --partition=pcmpg_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err

# ================= USER SETTINGS =================
# Project workdir (contains results/EDTA_annotation)
WORKDIR=/data/users/kxia/organize_annotation_course

# EDTA outputs
EDTA_DIR=$WORKDIR/results/EDTA_annotation
TELIB=$EDTA_DIR/assembly.fasta.mod.EDTA.TElib.fa

# Apptainer image for TEsorter
TESORTER_SIF=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
TESORTER_DB=rexdb-plant

# Output directory for this step
OUTDIR=$WORKDIR/results/TEsorter_refine
# =================================================

set -euo pipefail

mkdir -p "$OUTDIR"/{logs,tmp}
cd "$OUTDIR"

echo "[I] Checking inputs..."
[[ -s "$TELIB" ]] || { echo "[E] Missing EDTA TE library: $TELIB" >&2; exit 2; }
[[ -s "$TESORTER_SIF" ]] || { echo "[E] Missing TEsorter SIF: $TESORTER_SIF" >&2; exit 2; }

# ---------------------------
# Step 1: Extract Copia/Gypsy sequences with seqkit
# ---------------------------
echo "[I] Extracting Copia and Gypsy sequences from EDTA TElib..."

# Check seqkit availability; try to load a module name if missing
if ! command -v seqkit >/dev/null 2>&1; then
  module load SeqKit 2>/dev/null || true
fi
command -v seqkit >/dev/null 2>&1 || { echo "[E] seqkit not available. Please load it on this cluster." >&2; exit 3; }

# NOTE:
#  - EDTA library headers generally contain the classification token.
#  - We use regex (-r) to match case-insensitively: Copia / Gypsy
seqkit grep -r -i -p "Copia" "$TELIB" > Copia_sequences.fa
seqkit grep -r -i -p "Gypsy" "$TELIB" > Gypsy_sequences.fa

# Basic stats for sanity check
echo "[I] Copia count:  $(grep -c '^>' Copia_sequences.fa || true)"
echo "[I] Gypsy count:  $(grep -c '^>' Gypsy_sequences.fa || true)"

# ---------------------------
# Step 2: Run TEsorter on Copia and Gypsy (in Apptainer)
# ---------------------------
echo "[I] Running TEsorter (database: $TESORTER_DB) ..."

# Use CPUs from SLURM if present
P=${SLURM_CPUS_PER_TASK:-4}

# Copia
if [[ -s Copia_sequences.fa ]]; then
  echo "[I] TEsorter on Copia_sequences.fa ..."
  apptainer exec --bind "$WORKDIR" "$TESORTER_SIF" \
    TEsorter Copia_sequences.fa -db "$TESORTER_DB" -p "$P"
else
  echo "[W] Copia_sequences.fa is empty; skipping."
fi

# Gypsy
if [[ -s Gypsy_sequences.fa ]]; then
  echo "[I] TEsorter on Gypsy_sequences.fa ..."
  apptainer exec --bind "$WORKDIR" "$TESORTER_SIF" \
    TEsorter Gypsy_sequences.fa -db "$TESORTER_DB" -p "$P"
else
  echo "[W] Gypsy_sequences.fa is empty; skipping."
fi

# (Optional) You can also classify the entire TElib in one go for completeness:
# apptainer exec --bind "$WORKDIR" "$TESORTER_SIF" \
#   TEsorter "$TELIB" -db "$TESORTER_DB" -p "$P"

# ---------------------------
# Step 3: Collect key outputs and build merged tables
# ---------------------------
echo "[I] Collecting outputs and summarizing..."

# Expected key files (created in current OUTDIR)
COPIA_CLS=Copia_sequences.fa.${TESORTER_DB}.cls.tsv
GYPSY_CLS=Gypsy_sequences.fa.${TESORTER_DB}.cls.tsv

# Copy/rename to a standard location (keep as-is since we already run in OUTDIR)
[[ -s "$COPIA_CLS" ]] && echo "[I] Found $COPIA_CLS"
[[ -s "$GYPSY_CLS" ]] && echo "[I] Found $GYPSY_CLS"

# Merge cls tables (keep header once; handle missing files gracefully)
MERGED_CLS=TEsorter_merged_cls.tsv
{
  # print header: if exists in Copia; else from Gypsy; else a generic header
  if [[ -s "$COPIA_CLS" ]]; then
    head -n1 "$COPIA_CLS"
  elif [[ -s "$GYPSY_CLS" ]]; then
    head -n1 "$GYPSY_CLS"
  else
    echo -e "TE\tOrder\tSuperfamily\tClade\tEvalue\tCoverage\tIdentity\tHMM_Length\tQuery_Length"
  fi

  # append data rows (skip header)
  [[ -s "$COPIA_CLS" ]] && tail -n +2 "$COPIA_CLS"
  [[ -s "$GYPSY_CLS" ]] && tail -n +2 "$GYPSY_CLS"
} > "$MERGED_CLS"

# Quick summaries (per Superfamily and per Clade)
awk -F'\t' 'NR>1 && $3!="" {cnt[$3]++} END{print "Superfamily\tCount"; for(k in cnt) printf("%s\t%d\n",k,cnt[k]) | "sort -k2,2nr"}' "$MERGED_CLS" > summary_superfamily.tsv || true
awk -F'\t' 'NR>1 && $4!="" {cnt[$4]++} END{print "Clade\tCount"; for(k in cnt) printf("%s\t%d\n",k,cnt[k]) | "sort -k2,2nr"}' "$MERGED_CLS" > summary_clade.tsv || true

echo "[I] Wrote:"
echo "    - $MERGED_CLS"
echo "    - summary_superfamily.tsv"
echo "    - summary_clade.tsv"

# ---------------------------
# Step 4: Tips for integrating back to annotations
# ---------------------------
cat > HOWTO_integrate_back.txt <<'EOF'
How to integrate TEsorter clade calls back into your TE annotations (GFF3):

1) Map keys:
   - TEsorter 'TE' column corresponds to fasta record names (from Copia/Gypsy_sequences.fa).
   - EDTA GFF3 attributes typically include Name=... (family) and classification=... (e.g., LTR/Gypsy).
   - If your GFF3 'Name' matches the fasta headers (or derivable by stripping suffixes), you can join on that.

2) Suggested workflow in R:
   - Read the GFF3 (ignore comments), extract 'Name' and 'classification' from the attributes.
   - Read TEsorter_merged_cls.tsv, keep columns: TE, Superfamily, Clade.
   - Normalize keys:
       * remove trailing tokens like "_INT" in TEsorter TE if needed
       * or map by your library header convention
   - Left-join GFF3 <- TEsorter classes by the normalized key.
   - Write a new GFF3 or a table with an extra attribute, e.g., "Clade=Tat" appended to column 9.

3) Caveats:
   - Not all records will have a clade (NA). Keep them with Superfamily only.
   - Keep homology vs structural provenance if needed (method= in GFF3).
   - When in doubt, keep original GFF3 and write a separate table (Name,Superfamily,Clade).

EOF

echo "[I] Done. Outputs are in: $OUTDIR"