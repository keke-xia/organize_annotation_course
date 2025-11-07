#!/bin/bash
#SBATCH --job-name=TE_parseRM
#SBATCH --partition=pcmpg_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err
set -euo pipefail

# ===== User settings =====
WORKDIR=/data/users/kxia/organize_annotation_course
EDTA_DIR=$WORKDIR/results/EDTA_annotation
RM_OUT=$EDTA_DIR/assembly.fasta.mod.EDTA.anno/assembly.fasta.mod.out
OUTDIR=$WORKDIR/results/TE_dynamics
PARSER=$WORKDIR/scripts/05-parseRM.pl

mkdir -p "$OUTDIR" "$WORKDIR/logs"
cd "$OUTDIR"

echo "[I] Checking inputs..."
[[ -s "$RM_OUT" ]]  || { echo "[E] Missing RepeatMasker output: $RM_OUT" >&2; exit 2; }
[[ -s "$PARSER" ]]  || { echo "[E] Missing parseRM.pl: $PARSER" >&2; exit 3; }

# BioPerl is required for parseRM.pl
module purge
module load BioPerl/1.7.8-GCCcore-10.3.0

# Run parseRM.pl to create landscape tables
# -l 50,1: minimum length and bin size; adjust if needed
perl "$PARSER" -i "$RM_OUT" -l 50,1 -v > "$OUTDIR/parseRM.log" 2>&1 || {
  echo "[E] parseRM.pl failed. See $OUTDIR/parseRM.log"; exit 4; }

# Decide which landscape table to use (Rname preferred)
RNAME_TAB="${RM_OUT}.landscape.Div.Rname.tab"
RCLASS_TAB="${RM_OUT}.landscape.Div.Rclass.tab"
CHOSEN_TAB=""
if [[ -s "$RNAME_TAB" ]]; then
  CHOSEN_TAB="$RNAME_TAB"
elif [[ -s "$RCLASS_TAB" ]]; then
  CHOSEN_TAB="$RCLASS_TAB"
else
  echo "[E] No landscape tables found after parseRM."; exit 5
fi

cp -f "$CHOSEN_TAB" "$OUTDIR/"
echo "[✓] parseRM done. Landscape table copied to $OUTDIR/"
echo "    Using: $CHOSEN_TAB"

# === Extra: make .divsum and .summary for plot_div.R ===
REPEATMASKER_DIR="/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"
EDTA_ANN="$EDTA_DIR/assembly.fasta.mod.EDTA.anno"

[[ -s "$EDTA_ANN/assembly.fasta.mod.align" ]] || { echo "[E] Missing .align: $EDTA_ANN/assembly.fasta.mod.align"; exit 6; }

"$REPEATMASKER_DIR"/util/calcDivergenceFromAlign.pl -s \
  "$EDTA_ANN/assembly.fasta.mod.align" \
  > "$EDTA_ANN/assembly.fasta.mod.out.divsum"

"$REPEATMASKER_DIR"/util/buildSummary.pl \
  "$EDTA_ANN/assembly.fasta.mod.out" \
  > "$EDTA_ANN/assembly.fasta.mod.out.summary"

echo "[✓] Produced:"
ls -l "$EDTA_ANN"/assembly.fasta.mod.out.{divsum,summary}