#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --mem=64G
#SBATCH -p pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=Maker
#SBATCH --output=/data/users/kxia/organize_annotation_course/logs/Maker_%j.out
#SBATCH --error=/data/users/kxia/organize_annotation_course/logs/Maker_%j.err

set -euo pipefail

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/kxia/organize_annotation_course/gene_annotation"
REPEATMASKER_DIR="$COURSEDIR/softwares/RepeatMasker"

module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a
export PATH="$PATH:$REPEATMASKER_DIR"

cd "$WORKDIR"

mpiexec --oversubscribe -n 50 apptainer exec \
  --bind "$SCRATCH":/TMP \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$COURSEDIR/containers/MAKER_3.01.03.sif" \
  maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl