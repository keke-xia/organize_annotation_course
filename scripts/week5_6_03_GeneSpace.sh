#!/bin/bash
#SBATCH -p pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --job-name=Run_GeneSpace
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
  Rscript /data/users/kxia/organize_annotation_course/scripts/genespace.R "$WD"