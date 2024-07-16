import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.transforms as transforms

# File paths
proportions_file = "/home/sfadi/Bureau/tmpmerge_proportions.txt"
metadata_file = "/home/sfadi/Bureau/PCA_plot_Tissues/liver_grp_file"

# Load metadata and proportions data
metadata = pd.read_csv(metadata_file, delim_whitespace=True)
proportions_data = pd.read_csv(proportions_file, sep='\t')

# Clean 'Sample' column in metadata
metadata['Sample'] = metadata['Sample'].str.strip()

# Filter metadata for liver tissue
metadata_liver = metadata[metadata['Tissue'] == 'Liver']

# Filter metadata for specific species and count
species_count = {
    'Ficedula_hypoleuca': 5,
    'Ficedula_albicollis': 5,
    'F1_hybrid': 3
}

filtered_metadata = pd.DataFrame(columns=metadata_liver.columns)

for species, count in species_count.items():
    species_subset = metadata_liver[metadata_liver['Species'] == species].head(count)
    filtered_metadata = pd.concat([filtered_metadata, species_subset])

# Merge metadata with proportions data based on 'Sample'
proportions_data.columns = proportions_data.columns.str.strip()
merged_data = proportions_data.loc[:, proportions_data.columns.isin(filtered_metadata['Sample'])]

# Reindex the metadata to match the merged proportions data columns
filtered_metadata = filtered_metadata.set_index('Sample').reindex(merged_data.columns).reset_index()

# Check if merged data is empty
if merged_data.empty:
    print("No matching data found after merging metadata and proportions data. Please check the input files.")
else:
    # Transpose merged_data to have samples as rows and features as columns
    merged_data = merged_data.T

    # Standardize the data
    scaler = StandardScaler()
    pca_data_scaled = scaler.fit_transform(merged_data)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(pca_data_scaled)
    explained_variance = pca.explained_variance_ratio_ * 100  # Percentage of variance explained by each component

    # Create DataFrame from PCA results
    df_pca_result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    df_pca_result['Species'] = filtered_metadata['Species'].replace({
        'Ficedula_hypoleuca': 'Pied',
        'Ficedula_albicollis': 'Collared',
        'F1_hybrid': 'Hybrid'
    }).values

    # Define colors for species
    colors = {
        'Pied': '#0072BD',
        'Collared': '#D95319',
        'Hybrid': 'goldenrod'
    }

    # Plot PCA results
plt.figure(figsize=(14, 10))

# Plotting each species with their respective colors
for species, color in colors.items():
    subset = df_pca_result[df_pca_result['Species'] == species]
    if subset.empty:
        print(f"No data for species: {species}")
        continue
    plt.scatter(subset['PC1'], subset['PC2'], color=color, label=species, alpha=0.7, s=1000)  # Increased marker size

# Function to plot ellipses with thicker lines
def plot_ellipse(mean, cov, ax, n_std=1.0, facecolor='none', edgecolor='black', alpha=0.3, linewidth=2.0, **kwargs):
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor,
                      edgecolor=edgecolor, alpha=alpha, linewidth=linewidth, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mean[0]
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mean[1]
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Plot ellipses for all species
for species, color in colors.items():
    subset = df_pca_result[df_pca_result['Species'] == species]
    if not subset.empty:
        mean = subset[['PC1', 'PC2']].mean().values
        cov = np.cov(subset[['PC1', 'PC2']], rowvar=False)
        plot_ellipse(mean, cov, plt.gca(), edgecolor=color, alpha=0.3)

plt.title('PCA Plot: Intron Usage Variation in Liver Tissue', fontsize=20)
plt.xlabel(f'PC1 ({explained_variance[0]:.2f}% variance)', fontsize=18)
plt.ylabel(f'PC2 ({explained_variance[1]:.2f}% variance)', fontsize=18)

plt.legend()
plt.grid(False)  # Turn off grid to ensure no lines
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.show()

