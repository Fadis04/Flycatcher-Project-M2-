# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define file paths
dexseq_genes_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_hyp_Testis/unique_genes.txt"
non_significant_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_hyp_Testis/non_significant.txt"
padj_na_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_hyp_Testis/padj_na.txt"
significant_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_hyp_Testis/significant.txt"

# Read the files
dexseq_genes <- readLines(dexseq_genes_file)
non_significant <- read.table(non_significant_file, header = TRUE, stringsAsFactors = FALSE)
padj_na <- read.table(padj_na_file, header = TRUE, stringsAsFactors = FALSE)
significant <- read.table(significant_file, header = TRUE, stringsAsFactors = FALSE)

# Function to count occurrences of unique genes in a given data frame
count_genes <- function(df, gene_list) {
  df_filtered <- df %>% filter(genes %in% gene_list)
  return(length(unique(df_filtered$genes)))
}

# Count the occurrences in each file
count_non_significant <- count_genes(non_significant, dexseq_genes)
count_padj_na <- count_genes(padj_na, dexseq_genes)
count_significant <- count_genes(significant, dexseq_genes)

# Create a data frame for plotting
plot_data <- data.frame(
  category = c("Non-significant (Leafcutter)", "Padj NA (Leafcutter)", "Significant (Leafcutter)"),
  count = c(count_non_significant, count_padj_na, count_significant)
)

# Print the plot data to debug
print("Plot data:")
print(plot_data)

# Plot the data
ggplot(plot_data, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Abundance of DEXSeq-Identified Unique Genes in Leafcutter Results",
       x = "Category",
       y = "Number of Unique Genes Detected") +
  theme(legend.position = "none")
  
  
  
  # Load necessary libraries
library(ggplot2)
library(dplyr)

# Define file paths
dexseq_significant_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_hyp_Testis/significant_exons.txt"
dexseq_na_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_hyp_Testis/na_exons.txt"
dexseq_non_significant_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_hyp_Testis/non_significant_exons.txt"
leafcutter_unique_genes_file <- "/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_hyp_Testis/unique_genes.txt"

# Read DEXSeq files
dexseq_significant <- read.table(dexseq_significant_file, header = TRUE, stringsAsFactors = FALSE)
dexseq_na <- read.table(dexseq_na_file, header = TRUE, stringsAsFactors = FALSE)
dexseq_non_significant <- read.table(dexseq_non_significant_file, header = TRUE, stringsAsFactors = FALSE)

# Read unique genes from Leafcutter
leafcutter_unique_genes <- readLines(leafcutter_unique_genes_file)

# Function to filter and count genes
count_genes_in_files <- function(file_data, unique_genes) {
  file_genes <- unique(file_data$groupID)
  overlap_genes <- intersect(file_genes, unique_genes)
  return(length(overlap_genes))
}

# Count genes in each DEXSeq file
count_significant <- count_genes_in_files(dexseq_significant, leafcutter_unique_genes)
count_na <- count_genes_in_files(dexseq_na, leafcutter_unique_genes)
count_non_significant <- count_genes_in_files(dexseq_non_significant, leafcutter_unique_genes)

# Create plot data
plot_data <- data.frame(
  category = c("Significant (DEXSeq)", "Padj NA (DEXSeq)", "Non-significant (DEXSeq)"),
  count = c(count_significant, count_na, count_non_significant)
)

# Plotting
ggplot(plot_data, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Abundance of Unique Genes Detected by Leafcutter in DEXSeq Results",
       x = "Category",
       y = "Number of Unique Genes Detected") +
  theme(legend.position = "none")


