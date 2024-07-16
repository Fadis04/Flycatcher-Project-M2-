import pandas as pd
import numpy as np

# Define paths to input and output files
input_file = '/home/sfadi/Bureau/final_read_counts_all_samples_fixed.txt'
output_file = '/home/sfadi/Bureau/normalized_read_counts.txt'

# Read the fixed file into a pandas DataFrame
df = pd.read_csv(input_file, sep='\t', index_col=0)

# Step 1: Normalize by sequencing depth (to Reads Per Million, RPM)
# Calculate total reads per sample
total_reads = df.sum()

# Calculate RPM for each exon
rpm_df = df.div(total_reads) * 1e6  # multiply by 10^6 for RPM

# Step 2: Normalize by expression levels
# Calculate total expression (sum of all exons) for each gene in each sample
total_gene_expression = df.groupby(df.index.str.split(':').str[0]).sum()

# Normalize each exon count by the total expression of the gene it belongs to
normalized_df = df.div(total_gene_expression.loc[df.index.str.split(':').str[0]].values)

# Write the normalized data to a new file
normalized_df.to_csv(output_file, sep='\t')

print(f"Normalization complete. Normalized data saved to {output_file}")
