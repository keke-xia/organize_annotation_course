#!/bin/bash
#SBATCH --job-name=TE_dynamics
#SBATCH --partition=pcmpg_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err

# ================= USER SETTINGS =================
WORKDIR=/data/users/kxia/organize_annotation_course
EDTA_DIR=$WORKDIR/results/EDTA_annotation
RM_OUT=/data/users/kxia/organize_annotation_course/results/EDTA_annotation/assembly.fasta.mod.EDTA.anno/assembly.fasta.mod.out
OUTDIR=$WORKDIR/results/TE_dynamics

# parseRM.pl script location (downloaded or provided)
PARSER=/data/users/kxia/organize_annotation_course/scripts/05-parseRM.pl

# Provided course R script for TE landscape
PLOT_R=/data/users/kxia/organize_annotation_course/scripts/06-plot_div.R

# Substitution rate for Brassicaceae (Kagale et al. 2014)
RATE=8.22e-9
# =================================================

set -euo pipefail
mkdir -p "$OUTDIR" "$WORKDIR/logs"
cd "$OUTDIR"

echo "[I] Checking inputs..."
[[ -s "$RM_OUT" ]] || { echo "[E] Missing RepeatMasker output: $RM_OUT" >&2; exit 2; }
[[ -s "$PARSER" ]] || { echo "[E] Missing parseRM.pl: $PARSER" >&2; exit 3; }
[[ -s "$PLOT_R" ]] || { echo "[E] Missing plot_div.R: $PLOT_R" >&2; exit 4; }

# ---------------------------
# Step 1: Parse RepeatMasker output
# ---------------------------
echo "[I] Running parseRM.pl to compute corrected divergence..."

module purge
module load BioPerl/1.7.8-GCCcore-10.3.0

perl "$PARSER" -i "$RM_OUT" -l 50,1 -v > parseRM.log 2>&1

# Expected output: a file named "$RM_OUT.divsum"
# parseRM.pl v5.x may output Rclass.tab instead of .divsum
DIVSUM="${RM_OUT}.divsum"
ALT_RCLASS="${RM_OUT}.landscape.Div.Rclass.tab"

if [[ ! -s "$DIVSUM" && -s "$ALT_RCLASS" ]]; then
  DIVSUM="$ALT_RCLASS"
fi

if [[ ! -s "$DIVSUM" ]]; then
  echo "[E] Expected output not found: $DIVSUM or $ALT_RCLASS" >&2
  echo "    Check parseRM.log for details (parseRM.pl may have changed output format)." >&2
  exit 5
fi

cp "$DIVSUM" "$OUTDIR/"
echo "[I] parseRM.pl finished. File copied to $OUTDIR."

# ---------------------------
# Step 2: Generate TE landscape plot
# ---------------------------
echo "[I] Generating TE landscape plot with R..."

module load R/4.3.2-foss-2021a
export R_LIBS_USER="$WORKDIR/.Rlibs_R432"
mkdir -p "$R_LIBS_USER"

Rscript "$PLOT_R" "$DIVSUM"

# Move figures to the output folder
mv -f *.pdf *.png "$OUTDIR" 2>/dev/null || true

# ---------------------------
# Step 3: Estimate TE insertion ages
# ---------------------------
echo "[I] Estimating insertion ages based on substitution rate r=$RATE"

OUT_TSV="$OUTDIR/TE_insertion_age.tsv"

awk -v rate="$RATE" '
BEGIN{
  OFS="\t"
  print "Superfamily","Corrected_divergence_K","Estimated_age_years"
}
NR>1 && $1!~"^#" {
  K=$2/100  # convert percentage divergence to proportion
  T=K/(2*rate)
  print $1, K, T
}' "$DIVSUM" > "$OUT_TSV"

echo "[I] Wrote estimated insertion ages to $OUT_TSV"

echo "[✓] TE dynamics workflow completed."
echo "    Key outputs:"
echo "      - ${DIVSUM}"
echo "      - TE landscape plots (PDF/PNG)"
echo "      - $OUT_TSV"