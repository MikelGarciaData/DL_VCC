import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from pathlib import Path

from eda_toolkit import create_subset_h5ad_by_size

current_path = Path(__file__).parent
data_folder = (current_path.parent / 'data' / 'subset_5000cells_2000genes.h5ad').resolve()


# subset = create_subset_h5ad_by_size(
#     input_path=data_folder,
#     output_path="subset_5000cells_2000genes.h5ad",
#     n_cells=5000,
#     n_genes=2000
# )

adata = ad.read_h5ad("subset_5000cells_2000genes.h5ad")

# Compute total counts per cell
adata.obs["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()
print(adata.obs["total_counts"])

plt.hist(adata.obs["total_counts"], bins=50)
plt.xlabel("Total counts per cell")
plt.ylabel("Number of cells")
plt.title("Distribution of total counts per cell")
plt.show()


# Compute total expression per gene
adata.var["total_counts"] = np.array(adata.X.sum(axis=0)).flatten()

plt.hist(adata.var["total_counts"], bins=50)
plt.xlabel("Total counts per gene")
plt.ylabel("Number of genes")
plt.title("Distribution of total counts per gene")
plt.show()


# Number of cells per perturbation
print(adata.obs["target_gene"].value_counts().head(10))

# Number of batches
print(adata.obs["batch"].value_counts())