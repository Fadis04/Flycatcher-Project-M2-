#!/bin/bash
#SBATCH --job-name=summarize_overlaps_f1_hyb
#SBATCH --output=summarize_overlaps_f1_hyb.out
#SBATCH --error=summarize_overlaps_f1_hyb.err
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=110G
#SBATCH --time=01-00:00:00


/beegfs/data/soft/R-4.3.1/bin/Rscript  summarize_overlaps_f1_hyb.R
