#!/bin/bash
#SBATCH --job-name=TE_landscape
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err
set -euo pipefail

# --- Load R ---
module load R/4.3.2-foss-2021a

# === Paths (edit only if your layout differs) ===
cd /data/users/kxia/organize_annotation_course/results/EDTA_annotation
WORKDIR="/data/users/kxia/organize_annotation_course"
# EDTA annotation dir containing assembly.fasta.mod.out (+ .divsum/.summary)
EDTA_ANN="$WORKDIR/results/EDTA_annotation/assembly.fasta.mod.EDTA.anno"
# 'BASE' is the path prefix without extension; plot_div.R will read <BASE>.divsum and <BASE>.summary
BASE="$EDTA_ANN/assembly.fasta.mod.out"
# Output dir for plots
PLOTDIR="$WORKDIR/results/TE_dynamics/Plots"
# Course plotting script
PLOT_R="/data/users/kxia/organize_annotation_course/scripts/06-plot_div.R"
# ================================================

mkdir -p "$PLOTDIR"

# Sanity check: required inputs
[[ -s "${BASE}.divsum" && -s "${BASE}.summary" ]] || {
  echo "[E] Missing ${BASE}.divsum or ${BASE}.summary"; exit 2; }

# Run plotting script provided by the course
# -i takes base path (no extension); -o is output directory
Rscript "$PLOT_R" -i "$BASE" -o "$PLOTDIR"

echo "[✓] TE landscape plotting finished -> $PLOTDIR"