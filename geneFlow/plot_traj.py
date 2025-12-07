import scanpy as sc
import os
import matplotlib.pyplot as plt

# ==============================
# 1. SETUP
# ==============================
SAVE_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8"
time_points = [0.00, 0.25, 0.50, 0.75, 1.00]

# Create figure with 5 subplots
fig, axes = plt.subplots(1, 5, figsize=(25, 5))

print("Generating UMAP plots...")

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
        
        # Load Data
        adata_t = sc.read_h5ad(path)
        
        # Compute UMAP (Independent for each timepoint)
        sc.pp.neighbors(adata_t, use_rep='X', n_neighbors=30)
        sc.tl.umap(adata_t, random_state=42)
    

        # Plot
        sc.pl.umap(
            adata_t,
            color="condition",
            size=2,
            ax=ax,
            show=False, 
            title=f"Time: {t}",
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