import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define the file paths
file_paths = {
    'Brain': {
        'PiedxCollared': '/home/sfadi/Bureau/DSA_all_species_Brain/Leafcutter_output/alb_x_hyp_Brain/unique_significant_genes_clusters.txt',
        'PiedxF1':  '/home/sfadi/Bureau/DSA_all_species_Brain/Leafcutter_output/hyp_x_f1_Brain/unique_significant_genes_clusters.txt',
        'CollaredxF1': '/home/sfadi/Bureau/DSA_all_species_Brain/Leafcutter_output/alb_x_f1_Brain/unique_significant_genes_clusters.txt'
    },
    'Heart': {
        'PiedxCollared': '/home/sfadi/Bureau/DSA_all_species_Heart/Leafcutter_output/alb_x_hyp_Heart/unique_significant_genes_clusters.txt',
        'PiedxF1': '/home/sfadi/Bureau/DSA_all_species_Heart/Leafcutter_output/hyp_x_f1_Heart/unique_significant_genes_clusters.txt',
        'CollaredxF1':  '/home/sfadi/Bureau/DSA_all_species_Heart/Leafcutter_output/alb_x_f1_Heart/unique_significant_genes_clusters.txt',
    },
    'Kidney': {
        'PiedxCollared': '/home/sfadi/Bureau/DSA_all_species_Kidney/Leafcutter_output/alb_x_hyp_Kidney/unique_significant_genes_clusters.txt',
        'PiedxF1':'/home/sfadi/Bureau/DSA_all_species_Kidney/Leafcutter_output/hyp_x_f1_Kidney/unique_significant_genes_clusters.txt',
        'CollaredxF1':  '/home/sfadi/Bureau/DSA_all_species_Kidney/Leafcutter_output/alb_x_f1_Kidney/unique_significant_genes_clusters.txt',
    },
    'Liver': {
       'PiedxCollared': '/home/sfadi/Bureau/DSA_all_species_Liver/Leafcutter_output/alb_x_hyp_Liver/unique_significant_genes_clusters.txt',
        'PiedxF1': '/home/sfadi/Bureau/DSA_all_species_Liver/Leafcutter_output/hyp_x_f1_Liver/unique_significant_genes_clusters.txt',
        'CollaredxF1':'/home/sfadi/Bureau/DSA_all_species_Liver/Leafcutter_output/alb_x_f1_Liver/unique_significant_genes_clusters.txt',
    },
    'Testis': {
        'PiedxCollared': '/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_hyp_Testis/unique_significant_genes_clusters.txt',
        'PiedxF1': '/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/hyp_x_f1_Testis/unique_significant_genes_clusters.txt',
        'CollaredxF1': '/home/sfadi/Bureau/DSA_all_species_Testis/Leafcutter_output/alb_x_f1_Testis/unique_significant_genes_clusters.txt',
    }
}

# Function to read the genes from a file
def read_genes(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df['genes'].unique()

# Read genes data
genes_data = {tissue: {group: read_genes(file_paths[tissue][group]) for group in file_paths[tissue]} for tissue in file_paths}

# Count the number of affected genes and find common genes
counts = []
for tissue, groups in genes_data.items():
    tissue_counts = {'Tissue': tissue}
    for group, genes in groups.items():
        tissue_counts[group] = len(genes)
    counts.append(tissue_counts)

# Create a DataFrame from the counts
df_counts = pd.DataFrame(counts)

# Plot the data
fig, ax = plt.subplots(figsize=(14, 8))
bar_width = 0.2
index = range(len(df_counts))

bar1 = plt.bar(index, df_counts['PiedxCollared'], bar_width, label='PiedxCollared', color='blue')
bar2 = plt.bar([i + bar_width for i in index], df_counts['PiedxF1'], bar_width, label='PiedxF1', color='orange')
bar3 = plt.bar([i + 2 * bar_width for i in index], df_counts['CollaredxF1'], bar_width, label='CollaredxF1', color='green')

plt.xlabel('Tissue', fontsize=20)  # Increase font size for x-axis label
plt.ylabel('Number of Affected Genes', fontsize=16)  # Adjust font size for y-axis label if needed
plt.title('Number of Affected Genes Across Tissues for Each Group', fontsize=24)  # Increase font size for title
plt.xticks([i + 1.5 * bar_width for i in index], df_counts['Tissue'], fontsize=16)  # Increase font size for x-axis ticks

# Add counts on top of the bars
def add_labels(bars, fontsize=12):
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{int(height)}', 
                 ha='center', va='bottom', fontsize=fontsize, weight='bold')

add_labels(bar1, fontsize=18)
add_labels(bar2, fontsize=18)
add_labels(bar3, fontsize=18)

# Move legend to top left, make it larger
plt.legend(loc='upper left', fontsize=20)

plt.tight_layout()
plt.show()




