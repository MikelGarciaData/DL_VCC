import scanpy as sc
import os
import matplotlib.pyplot as plt

# 1. SETUP
# Define the seed for reproducibility
SEED = 42 

selected_genes = [
    "DNMT1", "INSIG1", "PHF10", "BRD9", "KLF10", "STAT6", "XRCC4", "MED13L", 
    "DAXX", "FDPS", "MAP3K7", "CASP2", "NDUFB4", "CREG1", "KAT2A", "USP22", 
    "IKBKG", "MED24", "OXA1L", "SALL4", "SIN3B", "STX4", "RRM1", "SV2A", 
    "SRC", "NCK2", "UBE3C", "CAMSAP2", "SMARCA5", "CASP3"
]

save_dir = "results_final"
os.makedirs(save_dir, exist_ok=True)

# 2. LOAD DATA
adata_path = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/adata_all_scvi_50.h5ad"
adata = sc.read_h5ad(adata_path)

# 3. SUBSET
subset_list = selected_genes + ['non-targeting']
mask = adata.obs['target_gene'].isin(subset_list)
adata_subset = adata[mask].copy()

# 4. RE-COMPUTE UMAP (With Seed)
print("Re-computing UMAP for the subset...")

# Neighbors (Deterministic usually, but good to be safe)
if 'X_scvi' in adata_subset.obsm:
    sc.pp.neighbors(adata_subset, use_rep='X_scvi', n_neighbors=15, random_state=SEED)
else:
    if 'X_pca' not in adata_subset.obsm:
        sc.pp.pca(adata_subset)
    sc.pp.neighbors(adata_subset, n_neighbors=15, random_state=SEED)

# UMAP (Stochastic -> Needs Seed)
sc.tl.umap(adata_subset, random_state=SEED)  # <--- SEED ADDED HERE

# 5. PLOT & SAVE
print("Plotting...")
plt.figure(figsize=(10, 8))

sc.pl.umap(
    adata_subset, 
    color='target_gene', 
    title='Selected Genes vs Non-Targeting',
    frameon=False,
    size=50,
    legend_loc='on data', 
    show=False
)

save_path = os.path.join(save_dir, "subset_genes_umap_fixed.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved to: {save_path}")