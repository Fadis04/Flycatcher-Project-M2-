import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt

# Use a non-interactive backend
import matplotlib
matplotlib.use('Agg')

# Load metadata (groups_file_PiedxCollared.txt)
metadata_file = "/home/sfadi/Bureau/PCA_Plot_PC/groups_file_PiedxCollared.txt"
metadata = pd.read_csv(metadata_file, sep='\t')

# Load proportions data
proportions_file = "/home/sfadi/Bureau/normalized_read_counts.txt"
proportions_data = pd.read_csv(proportions_file, sep='\t', index_col=0)

# Transpose proportions data
proportions_data = proportions_data.T
proportions_data.columns = proportions_data.columns.str.strip()

# Split the concatenated column into separate columns
metadata[['Sample', 'Species', 'Tissue']] = metadata['Sample Species Tissue'].str.split(expand=True)
metadata.drop(columns=['Sample Species Tissue'], inplace=True)

# Now merge metadata with proportions data
merged_data = pd.merge(metadata, proportions_data, left_on='Sample', right_index=True, how='inner')

# Extract numeric columns for PCA
pca_data = merged_data.drop(['Sample', 'Species', 'Tissue'], axis=1)

# Check for NaN values before imputation
print("Number of NaN values before imputation:", pca_data.isna().sum().sum())

# Impute missing values with the mean of each column
imputer = SimpleImputer(strategy='mean')
pca_data_imputed = imputer.fit_transform(pca_data)

# Check for NaN values after imputation
print("Number of NaN values after imputation:", pd.DataFrame(pca_data_imputed).isna().sum().sum())

# Standardize the data
scaler = StandardScaler()
pca_data_scaled = scaler.fit_transform(pca_data_imputed)

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pca_data_scaled)
explained_variance = pca.explained_variance_ratio_ * 100  # Percentage of variance explained by each component

# Create DataFrame from PCA results
df_pca_result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
df_pca_result['Species'] = merged_data['Species'].replace({
    'Ficedula_hypoleuca': 'Pied',  # Replace Ficedula_hypoleuca with Pied
    'Ficedula_albicollis': 'Collared'  # Replace Ficedula_albicollis with Collared
})
df_pca_result['Tissue'] = merged_data['Tissue']

# Plot PCA results
plt.figure(figsize=(14, 10))  # Increased figure size for better presentation

# Define colors for Pied and Collared
colors = {'Pied': '#0072BD', 'Collared': '#D95319'}

# Define marker shapes for 5 tissues
markers = {'Brain': 'o', 'Heart': 's', 'Kidney': 'D', 'Liver': '^', 'Testis': 'v'}

for tissue, marker in markers.items():
    for specie in df_pca_result['Species'].unique():
        subset = df_pca_result[(df_pca_result['Tissue'] == tissue) & (df_pca_result['Species'] == specie)]
        plt.scatter(subset['PC1'],
                    subset['PC2'],
                    color=colors[specie],
                    marker=marker,
                    label=f'{tissue} - {specie}',
                    alpha=0.7,
                    s=300)  # Further increased marker size

# Add legend for tissue shapes
handles = []
for tissue, marker in markers.items():
    handles.append(plt.scatter([], [], marker=marker, color='gray', label=tissue, alpha=0.7, s=300))

# Add legend for species colors
legend_colors = {'Pied': '#0072BD', 'Collared': '#D95319'}
for specie, color in legend_colors.items():
    handles.append(plt.scatter([], [], marker='o', color=color, label=specie, alpha=0.7, s=300))

plt.title('Principal Component Analysis of Exon Usage in Flycatcher Species Across Different Tissues', fontsize=20)
plt.xlabel(f'PC1 ({explained_variance[0]:.2f}% variance)', fontsize=18)
plt.ylabel(f'PC2 ({explained_variance[1]:.2f}% variance)', fontsize=18)

# Adjust legend position and labels
plt.legend(handles=handles, loc='upper right', title='Species and Tissue', fontsize=14, title_fontsize=16)
plt.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Set x and y axis ticks to units of 100
plt.xticks(range(int(min(df_pca_result['PC1'])//100)*100, int(max(df_pca_result['PC1'])//100)*100 + 100, 100))
plt.yticks(range(int(min(df_pca_result['PC2'])//100)*100, int(max(df_pca_result['PC2'])//100)*100 + 100, 100))

# Save the plot to a file
plt.savefig('/home/sfadi/Bureau/pca_plotDEXSEq_Normalized.png')

