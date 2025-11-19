#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=Riparian_plot
#SBATCH --output=/data/users/kxia/organize_annotation_course/genespace_work/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/genespace_work/logs/%x_%j.err

set -euo pipefail

COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
WD=/data/users/kxia/organize_annotation_course/genespace_work

apptainer exec \
  --env LC_ALL=C.UTF-8,LANG=C.UTF-8 \
  --bind "$COURSEDIR" \
  --bind /data/users/kxia/organize_annotation_course \
  "$COURSEDIR/containers/genespace_latest.sif" \
  Rscript /data/users/kxia/organize_annotation_course/scripts/riparian_plot.R "$WD"