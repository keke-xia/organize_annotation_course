#!/bin/bash
#SBATCH --job-name=MakerControl
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/%x_%j.err

set -euo pipefail

WORKDIR=/data/users/kxia/organize_annotation_course/gene_annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Generate control files with the container
# maker -CTL creates maker_opts.ctl, maker_bopts.ctl, maker_exe.ctl in current dir
apptainer exec --bind "$WORKDIR" \
  $COURSEDIR/containers/MAKER_3.01.03.sif \
  maker -CTL