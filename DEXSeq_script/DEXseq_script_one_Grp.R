# Load necessary libraries
library(DEXSeq)
library(SummarizedExperiment)

# Read the summarized experiment object
se <- readRDS("/home/sfadi/Bureau/summarized_experiment.rds")

# Read the group file
group_file <- read.table("/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_f1_Testis/alb_x_f1_Testis.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(group_file) <- c("Sample", "Species")

# Match the samples in the group file with the colnames of the SummarizedExperiment object
matching_indices <- match(group_file$Sample, colnames(se))
se <- se[, matching_indices]

# Assign the condition and libType to the colData of the SummarizedExperiment object
colData(se)$condition <- factor(group_file$Species)
colData(se)$libType <- factor(rep("single-end", ncol(se)))  # Adjust this if you have different library types

# Create the DEXSeqDataSet object
dxd <- DEXSeqDataSetFromSE(se, design = ~ sample + exon + condition:exon)

# Print the DEXSeqDataSet object to confirm
print(dxd)

# Save the DEXSeqDataSet object
saveRDS(dxd, "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_f1_Testis/dxd.rds")

# Confirm saving
print("DEXSeqDataSet object saved successfully.")

# Normalisation
dxd = estimateSizeFactors( dxd )

# Dispersion estimation
dxd = estimateDispersions( dxd )

# Plotting the per-exon dispersion estimates versus the mean normalised count
plotDispEsts( dxd )

# Testing for differential exon usage
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
 dxr1 = DEXSeqResults( dxd )

