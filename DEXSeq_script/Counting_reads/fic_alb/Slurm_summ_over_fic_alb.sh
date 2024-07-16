#!/bin/bash
#SBATCH --job-name=summarize_overlaps_fic_alb
#SBATCH --output=summarize_overlaps_fic_alb.out
#SBATCH --error=summarize_overlaps_fic_alb.err
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=05:00:00


/beegfs/data/soft/R-4.3.1/bin/Rscript  summarize_overlaps_fic_alb.R
