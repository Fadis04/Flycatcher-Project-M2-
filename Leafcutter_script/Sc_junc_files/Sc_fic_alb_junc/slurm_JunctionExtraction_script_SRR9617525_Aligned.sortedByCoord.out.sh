#!/bin/bash
#SBATCH --job-name=SRR9617525_Aligned.sortedByCoord.out.bam_JunctionExtraction
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --output=/beegfs/home/fslimi/scripts/Sc_fic_alb_junc/SRR9617525_Aligned.sortedByCoord.out.bam_%j_Junc.out
#SBATCH --error=/beegfs/home/fslimi/scripts/Sc_fic_alb_junc/SRR9617525_Aligned.sortedByCoord.out.bam_%j_Junc.err

# Define the directory containing the BAM file
bam_dir="/beegfs/data/fslimi/fic_bam/fic_alb_bam/"

# Define the output directory for junction files
output_dir="/beegfs/data/fslimi/fic_junc/fic_alb_junc/"

# Define the name of the BAM file to test
bamfile="SRR9617525_Aligned.sortedByCoord.out.bam"

# Move to the directory containing the BAM files
cd "$bam_dir" || exit

# Print a message indicating which file is being converted
echo "Converting ${bamfile} to ${bamfile}.junc"

# Index the BAM file using samtools
samtools index "${bamfile}"

# Extract junctions using regtools
regtools junctions extract -s XS -a 8 -m 50 -M 500000 "${bamfile}" -o "${output_dir}${bamfile}.junc"

# Add the path to the junction file to the text file
echo "${output_dir}${bamfile}.junc" >> "${bam_dir}juncfiles.txt"

echo "Junction file for ${bamfile} created successfully."
