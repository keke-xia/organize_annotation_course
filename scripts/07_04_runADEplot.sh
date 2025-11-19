#!/bin/bash
#SBATCH --job-name=AED_plot
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/AED_plot_%j.out
#SBATCH --error=logs/AED_plot_%j.err

# Load R module
module load R/4.3.2

# Go to the folder containing the AED input & where you want the PDF saved
cd /data/users/kxia/organize_annotation_course/gene_annotation/final

# Run the R script
Rscript /data/users/kxia/organize_annotation_course/scripts/07-ADE_plot.R