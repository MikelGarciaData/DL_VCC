import scanpy as sc
import os
import matplotlib.pyplot as plt
import anndata as ad

# ==============================
# 1. SETUP
# ==============================
SAVE_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8"
target_genes = ['SMARCA5', 'KAT2A', 'DNMT1', 'FDPS', 'SIN3B']

palette = {
    'Start (t=0)': '#d3d3d3',  # Light Grey
    'SMARCA5': 'tab:blue',
    'KAT2A':   'tab:orange',
    'DNMT1':   'tab:green',
    'FDPS':    'tab:red',
    'SIN3B':   'tab:purple'
}

# ==============================
# 2. LOAD & PREPARE DATA
# ==============================

# --- A. Process t=0 (Control) ---
path_t0 = os.path.join(SAVE_DIR, "predictions_t_0.00.h5ad")
print(f"Loading Start data: {path_t0}")
adata_t0 = sc.read_h5ad(path_t0)
adata_t0 = adata_t0[adata_t0.obs['condition'].isin(target_genes)].copy()
adata_t0.obs['plot_label'] = 'Start (t=0)'

# --- B. Process t=1 (Endpoint) ---
path_t1 = os.path.join(SAVE_DIR, "predictions_t_1.00.h5ad")
print(f"Loading End data: {path_t1}")
adata_t1 = sc.read_h5ad(path_t1)
adata_t1 = adata_t1[adata_t1.obs['condition'].isin(target_genes)].copy()
adata_t1.obs['plot_label'] = adata_t1.obs['condition']

# --- C. Combine ---
print("Concatenating t=0 and t=1...")
adata_combined = ad.concat([adata_t0, adata_t1], join='outer')

# *** FIX: Make observation names unique so sorting works ***
adata_combined.obs_names_make_unique()

# ==============================
# 3. COMPUTE UMAP
# ==============================
print("Computing Neighbors and UMAP...")
sc.pp.neighbors(adata_combined, use_rep='X', n_neighbors=30)
sc.tl.umap(adata_combined)

# ==============================
# 4. PLOT
# ==============================
# Sort so 'Start (t=0)' is plotted first (background) and colored genes on top
adata_combined.obs['z_order'] = adata_combined.obs['plot_label'].apply(lambda x: 0 if x == 'Start (t=0)' else 1)

# This line caused the error before; now it will work because names are unique
adata_combined = adata_combined[adata_combined.obs.sort_values('z_order').index]

# Create Plot
fig, ax = plt.subplots(figsize=(8, 8))

sc.pl.umap(
    adata_combined,
    color="plot_label",
    palette=palette,
    size=50,
    ax=ax,
    show=False,
    title="Gene Trajectory: Start (t=0) vs End (t=1)",
    frameon=False
)

ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")

# ==============================
# 5. SAVE
# ==============================
save_path = os.path.join(SAVE_DIR, "trajectory_t0_vs_t1_genes.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSaved comparison plot to: {save_path}")