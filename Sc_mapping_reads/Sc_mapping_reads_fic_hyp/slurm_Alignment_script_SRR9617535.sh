#!/bin/bash
#SBATCH --job-name=SRR9617535_Alignment
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --output=/beegfs/home/fslimi/scripts/Sc_mapping_reads_fic_hyp/SRR9617535_%j_Alig.out
#SBATCH --error=/beegfs/home/fslimi/scripts/Sc_mapping_reads_fic_hyp/SRR9617535_%j_Alig.err

# Run STAR alignment with gunzip command for decompression
STAR --outFileNamePrefix "/beegfs/data/fslimi/fic_hyp_bam/SRR9617535_"      --quantMode GeneCounts      --readFilesCommand "gunzip -c"      --outSAMtype BAM SortedByCoordinate      --runThreadN 1      --genomeDir "/beegfs/data/fslimi/genome_index/"      --readFilesIn "/beegfs/data/fslimi/RNAseq/fic_hyp/SRR9617535/SRR9617535_1.fastq.gz" "/beegfs/data/fslimi/RNAseq/fic_hyp/SRR9617535/SRR9617535_2.fastq.gz"      --outFilterMultimapNmax 1      --outFilterMismatchNmax 3      --twopassMode Basic      --outSAMstrandField intronMotif
