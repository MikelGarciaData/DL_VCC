import scanpy as sc
import os
import matplotlib.pyplot as plt

# ==============================
# 1. SETUP
# ==============================
SAVE_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_128_4"
time_points = [0.00, 0.25, 0.50, 0.75, 1.00]

# Create figure with 5 subplots
fig, axes = plt.subplots(1, 5, figsize=(25, 5))

print("Generating 5 side-by-side plots using sc.pl.umap...")

# ==============================
# 2. LOOP & PLOT
# ==============================
for i, t in enumerate(time_points):
    fname = f"predictions_t_{t:.2f}.h5ad"
    path = os.path.join(SAVE_DIR, fname)
    
    # Get the current axis (subplot)
    ax = axes[i]
    
    if os.path.exists(path):
        print(f"Processing t={t}...")
        
        # 1. Load Data
        adata_t = sc.read_h5ad(path)
        
        # 2. Compute UMAP (Independent for each timepoint)
        sc.pp.neighbors(adata_t, use_rep='X', n_neighbors=30)
        sc.tl.umap(adata_t)
        
        # 3. Create 'sample' column if it doesn't exist
        # (Mapping 'condition' (Gene Name) to 'sample' so your snippet works)
        if 'sample' not in adata_t.obs:
            adata_t.obs['sample'] = adata_t.obs['condition']

        # 4. Plot using Scanpy
        # We pass 'ax=ax' to put it in the specific panel
        # We pass 'show=False' so it doesn't pop up immediately
        sc.pl.umap(
            adata_t,
            color="sample",
            # Setting a smaller point size to prevent overlap
            size=2,
            ax=ax,
            show=False, 
            title=f"Time: {t}",
            # Optional: Turn off legend for first 4 plots to save space
            legend_loc='right margin' if i == 4 else None 
        )

    else:
        print(f"Warning: File {fname} not found.")
        ax.text(0.5, 0.5, "Missing File", ha='center')
        ax.axis('off')

# ==============================
# 3. SAVE
# ==============================
plt.tight_layout()
save_path = os.path.join(SAVE_DIR, "trajectory_scanpy_umap.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSaved plot to: {save_path}")