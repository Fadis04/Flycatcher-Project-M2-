#!/bin/bash
#SBATCH --job-name=SRR9617580_Alignment
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --output=/beegfs/home/fslimi/scripts/mapped_reads_f1_hyb_output-err/SRR9617580_%j_Alig.out
#SBATCH --error=/beegfs/home/fslimi/scripts/mapped_reads_f1_hyb_output-err/SRR9617580_%j_Alig.err

# Run STAR alignment with gunzip command for decompression
STAR --outFileNamePrefix "/beegfs/data/fslimi/f1_hyb_bam/SRR9617580_"      --quantMode GeneCounts      --readFilesCommand "gunzip -c"      --outSAMtype BAM SortedByCoordinate      --runThreadN 1      --genomeDir "/beegfs/data/fslimi/genome_index/"      --readFilesIn "/beegfs/data/fslimi/RNAseq/hybrids/SRR9617580/SRR9617580_1.fastq.gz" "/beegfs/data/fslimi/RNAseq/hybrids/SRR9617580/SRR9617580_2.fastq.gz"      --outFilterMultimapNmax 1      --outFilterMismatchNmax 3      --twopassMode Basic      --outSAMstrandField intronMotif
