#!/bin/bash
#SBATCH --job-name=SRR9617545_Alignment
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --output=/beegfs/home/fslimi/scripts/Sc_mapping_readsSRR9617545_%j_Alig.out
#SBATCH --error=/beegfs/home/fslimi/scripts/Sc_mapping_readsSRR9617545_%j_Alig.err

# Run STAR alignment with gunzip command for decompression
STAR --outFileNamePrefix "/beegfs/data/fslimi/fic_alb_bam/SRR9617545_"      --quantMode GeneCounts      --readFilesCommand "gunzip -c"      --outSAMtype BAM SortedByCoordinate      --runThreadN 1      --genomeDir "/beegfs/data/fslimi/genome_index/"      --readFilesIn "/beegfs/data/fslimi/RNAseq/fic_alb/SRR9617545/SRR9617545_1.fastq.gz" "/beegfs/data/fslimi/RNAseq/fic_alb/SRR9617545/SRR9617545_2.fastq.gz"      --outFilterMultimapNmax 1      --outFilterMismatchNmax 3      --twopassMode Basic      --outSAMstrandField intronMotif
