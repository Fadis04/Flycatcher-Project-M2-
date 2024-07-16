#!/bin/bash
#SBATCH --job-name=summarize_overlaps_fic_hyp
#SBATCH --output=summarize_overlaps_fic_hyp.out
#SBATCH --error=summarize_overlaps_fic_hyp.err
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=05:00:00


/beegfs/data/soft/R-4.3.1/bin/Rscript  summarize_overlaps_fic_hyp.R
