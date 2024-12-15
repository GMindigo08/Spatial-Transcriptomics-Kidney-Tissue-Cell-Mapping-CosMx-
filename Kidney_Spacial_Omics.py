# Spatial Omics - Kidney CosMx
# Michael Gage



#----------------------------------------------- Load in Libraries/Data ------------------------------------------------

# Loading libraries

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import title


# Load AnnData object (assumes multiple samples)
adata = sc.read("C:/Users/micha/OneDrive/Documents/Bioinformatics Masters/INDEPENDENT PROJECTS/Spacial Transcriptomics -"
                " Kidney/kidney_data.h5ad")



#---------------------------------------------Subset Sample of Interest-------------------------------------------------

# Subset the data for sample of interest (Healthy - HK3039, Diseased - HK2844)
sample_id = "HK3039"
if "sample" in adata.obs.columns and sample_id in adata.obs["sample"].unique():
    adata_subset = adata[adata.obs["sample"] == sample_id].copy()
    print(f"adata_subset created with shape: {adata_subset.shape}")
else:
    raise ValueError(f"Sample ID {sample_id} not found in adata.obs['sample'] or 'sample' column missing!")



# --------------------------------- Plot Cell Type Clusters and Associated Marker Genes ---------------------------------

# Plot a basic UMAP based on cell types
if "X_umap" in adata_subset.obsm:
    sc.pl.umap(adata_subset, color="cellType_CosMx_1", title= f"Sample {sample_id}", legend_loc="on data")


# Validate cell-type clustering with marker genes (to ensure clusters are biologically meaningful)
sc.tl.rank_genes_groups(adata_subset, groupby="cellType_CosMx_1", method="wilcoxon")

# Plot top 10 marker genes for each cell type
sc.pl.rank_genes_groups(adata_subset, n_genes=10, sharey=False)



# ------------------------------------- Plot Spatial Expression for Cell Type ---------------------------------------

# Specify the desired cell types
#Cell Type List: "TAL", "PT", "Immune", "Endothelium", "PC", "Injured TAL", "Fibro", "VSMC/Mes", "IC/DCT", "iPT", "PEC","DCT", "Podo"

cell_types = ["TAL", "PT", "Immune", "Endothelium", "PC", "Injured TAL", "Fibro", "VSMC/Mes", "IC/DCT", "iPT", "PEC","DCT", "Podo"]

# Verify spatial coordinates
if "spatial" in adata_subset.obsm:
    print("Spatial keys in adata_subset.obsm:", adata_subset.obsm.keys())
    print("Shape of spatial coordinates:", adata_subset.obsm["spatial"].shape)
else:
    raise KeyError("'spatial' key not found in adata_subset.obsm!")

# Ensure spatial data is in the correct format
if not isinstance(adata_subset.obsm["spatial"], np.ndarray):
    adata_subset.obsm["spatial"] = np.array(adata_subset.obsm["spatial"])
    print("Reformatted spatial data shape:", adata_subset.obsm["spatial"].shape)

# Validate that the cell type column exists
if "cellType_CosMx_1" not in adata_subset.obs.columns:
    raise KeyError("'cellType_CosMx_1' column not found in adata_subset.obs!")

# Extract cell type labels
cell_type_labels = adata_subset.obs["cellType_CosMx_1"].astype(str)
print(f"Unique cell types: {cell_type_labels.unique()}")

# Create a color mapping for specified cell types
unique_cell_types = list(np.unique(cell_type_labels))
specified_colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))  # Discrete colormap
cell_type_to_color = {cell_type: specified_colors[i] for i, cell_type in enumerate(cell_types)}
cell_type_to_color["Other"] = "gray"  # Assign gray for points not in specified cell types

# Assign each point a color based on its cell type
point_colors = [
    cell_type_to_color[cell_type] if cell_type in cell_types else "gray"
    for cell_type in cell_type_labels
]

# Plot the tissue section
plt.figure(figsize=(10, 8))
plt.scatter(
    adata_subset.obsm["spatial"][:, 0],  # X coordinates
    adata_subset.obsm["spatial"][:, 1],  # Y coordinates
    c=point_colors,  # Colors based on cell type
    s=0.1,  # Small point size
    alpha=0.8
)

# Add legend for specified cell types
for cell_type, color in cell_type_to_color.items():
    if cell_type != "Other":  # Skip "Other" from the legend
        plt.scatter([], [], c=color, label=cell_type, s=5, alpha=0.8)

plt.legend(loc="upper right", fontsize="small", title="Cell Type", markerscale=5)
plt.title(f"Spatial Expression by Specified Cell Types, Sample: {sample_id}")
plt.xlabel("Spatial1")
plt.ylabel("Spatial2")
plt.tight_layout()

# Show plot
plt.show()




# ------------------------------------- Plot Spatial Expression for Marker Genes ---------------------------------------


# Define marker genes
marker_genes = ["HAVCR1"]

if "sample" in adata.obs.columns and sample_id in adata.obs["sample"].unique():
    adata_subset = adata[adata.obs["sample"] == sample_id].copy()
    print(f"adata_subset created with shape: {adata_subset.shape}")
else:
    raise ValueError(f"Sample ID {sample_id} not found in adata.obs['sample'] or 'sample' column missing!")

# Verify spatial coordinates
if "spatial" in adata_subset.obsm:
    print("Spatial keys in adata_subset.obsm:", adata_subset.obsm.keys())
    print("Shape of spatial coordinates:", adata_subset.obsm["spatial"].shape)
else:
    raise KeyError("'spatial' key not found in adata_subset.obsm!")

# Ensure spatial data is in the correct format
if not isinstance(adata_subset.obsm["spatial"], np.ndarray):
    adata_subset.obsm["spatial"] = np.array(adata_subset.obsm["spatial"])
    print("Reformatted spatial data shape:", adata_subset.obsm["spatial"].shape)

# Initialize a categorical array to store the label of the "most expressed gene" for each point
expression_labels = np.full(adata_subset.n_obs, "None")  # Default "None" for no gene expression

# Assign a gene label to each cell based on maximum expression
for gene in marker_genes:
    if gene in adata_subset.var_names:
        expression_values = adata_subset[:, gene].X.toarray().flatten()  # Convert sparse matrix to dense array
        mask = expression_values > 0  # Find cells where the gene is expressed
        expression_labels[mask] = gene  # Assign the gene's name to these cells
    else:
        print(f"Marker gene {gene} not found in adata_subset.var_names!")

# Create a unique color for each gene and set "None" to gray
unique_labels = list(np.unique(expression_labels))  # All unique labels, including "None"
colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))  # Discrete colormap
label_to_color = {label: colors[i] for i, label in enumerate(unique_labels)}  # Map each label to a color
label_to_color["None"] = "gray"  # Explicitly set "None" to gray

# Plot the tissue section with unique colors for each gene
plt.figure(figsize=(10, 8))

# Plot "None" points first (background layer)
none_mask = expression_labels == "None"
plt.scatter(
    adata_subset.obsm["spatial"][none_mask, 0],  # X coordinates
    adata_subset.obsm["spatial"][none_mask, 1],  # Y coordinates
    c="gray",  # Gray color for "None"
    s=0.1,  # Small point size
    label="None",
    alpha=0.5  # Slight transparency for "None"
)

# Overlay marker gene points (foreground layers)
for label, color in label_to_color.items():
    if label != "None":  # Skip "None" for this loop
        mask = expression_labels == label
        plt.scatter(
            adata_subset.obsm["spatial"][mask, 0],  # X coordinates
            adata_subset.obsm["spatial"][mask, 1],  # Y coordinates
            c=[color],  # Color for this gene
            s=0.1,  # Small point size
            label=label,
            alpha=0.8  # Full opacity for marker genes
        )

# Add legend and labels
plt.legend(loc="upper right", fontsize="small", title="Genes", markerscale=5)
plt.title(f"{marker_genes} Spatial Expression (Sample: {sample_id})")
plt.xlabel("Spatial1")
plt.ylabel("Spatial2")
plt.tight_layout()

# Show plot
plt.show()

