#!/bin/bash
#SBATCH --job-name=LTR_identity
#SBATCH --partition=pcmpg_el8 
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err


# ================= USER-DEFINED VARIABLES =================
# Working directory (should contain results/EDTA_annotation)
WORKDIR=/data/users/kxia/organize_annotation_course

# TEsorter container path
TESORTER_SIF=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif

# Provided R plotting script (do not rename)
PLOT_R=/data/users/kxia/organize_annotation_course/scripts/02-full_length_LTRs_identity.R

# EDTA result directory (based on course default output structure)
EDTA_DIR=$WORKDIR/results/EDTA_annotation
RAW_DIR=$EDTA_DIR/assembly.fasta.mod.EDTA.raw

# Key input files (according to course instructions)
GFF_INTact=$RAW_DIR/assembly.fasta.mod.LTR.intact.raw.gff3
LTR_RAW_FA=$RAW_DIR/assembly.fasta.mod.LTR.raw.fa

# Output directory for this step
OUTDIR=$WORKDIR/results/full_length_LTR_identity
# ==========================================================

set -euo pipefail

mkdir -p "$OUTDIR"/{plots,logs,tmp}
cd "$OUTDIR"

echo "==> Checking input files"
if [[ ! -s "$GFF_INTact" ]]; then
  echo "GFF file not found: $GFF_INTact" >&2; exit 2
fi
if [[ ! -s "$LTR_RAW_FA" ]]; then
  echo "Raw LTR fasta file not found: $LTR_RAW_FA" >&2; exit 2
fi

echo "==> Step 1: Extract LTR percent identity from intact LTR-RT annotations (optional preview table)"
# Generate a simple verification table (Name\tIdentity).
# The R script will also parse these values internally; this table is just for quick inspection.
awk -F'\t' '
  $3!="long_terminal_repeat" && $3!="repeat_region" && $3!="target_site_duplication" {
    split($9,a,";");
    name=""; id="";
    for(i in a){
      if(a[i] ~ /^Name=/){ sub(/^Name=/,"",a[i]); name=a[i] }
      if(a[i] ~ /^ltr_identity=/){ sub(/^ltr_identity=/,"",a[i]); id=a[i] }
    }
    if(name!="" && id!=""){ print name"\t"id }
  }' "$GFF_INTact" | sort -u > tmp/Name_Identity.tsv

echo "==> Step 2: Run TEsorter for clade classification (rexdb-plant)"
# Run within OUTDIR so output *.cls.tsv is created here
apptainer exec --bind "$WORKDIR" "$TESORTER_SIF" \
  TEsorter "$LTR_RAW_FA" -db rexdb-plant

# The output prefix matches the input; the result is: assembly.fasta.mod.LTR.raw.fa.rexdb-plant.cls.tsv
CLS_TSV=assembly.fasta.mod.LTR.raw.fa.rexdb-plant.cls.tsv
if [[ ! -s "$CLS_TSV" ]]; then
  echo "TEsorter classification file not found: $CLS_TSV" >&2; exit 3
fi

echo "==> Step 3: Prepare expected filenames for the course R script (create symbolic links)"
# The provided R script expects exactly these two filenames:
# genomic.fna.mod.LTR.intact.raw.gff3
# genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv
ln -sf "$GFF_INTact" genomic.fna.mod.LTR.intact.raw.gff3
ln -sf "$CLS_TSV"   genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv

echo "==> Step 4: Run the R plotting script (outputs will be stored under plots/)"
module load R
Rscript "$PLOT_R"

# Move generated figures to plots/
shopt -s nullglob
for f in 01_LTR_Copia_Gypsy_cladelevel.png 01_LTR_Copia_Gypsy_cladelevel.pdf; do
  [[ -s "$f" ]] && mv -f "$f" plots/
done

echo "==> All steps completed successfully. Final figures are located in: $OUTDIR/plots/"