
# ============================
# 1) LOAD THE DATA
# ============================
adata = sc.read_h5ad(adata_path)
print("Loaded:", adata)


# ============================
# 3) SUBSAMPLING
# ============================

df = adata.obs.copy()

# PARAMETERS YOU CONTROL
n_perts = 10          # number of perturbations to sample
cells_per_pert = 40   # number of cells per perturbation
n_controls = 10       # number of control cells

# list all perturbations except controls
all_perts = df["target_gene"].unique().tolist()
all_perts = [p for p in all_perts if p != "non-targeting"]

# randomly select N perturbations
chosen_perts = np.random.choice(all_perts, size=n_perts, replace=False)
print("Chosen perturbations:", chosen_perts)

cells_keep = []

# ---- sample controls ----
control_cells = df.index[df["target_gene"] == "non-targeting"]
control_sub = np.random.choice(control_cells, size=n_controls, replace=False)
cells_keep.extend(control_sub)

# ---- sample K cells from each perturbation ----
for p in chosen_perts:
    ids = df.index[df["target_gene"] == p]
    if len(ids) < cells_per_pert:
        print(f"Skipping {p}: only {len(ids)} cells available.")
        continue
    sampled = np.random.choice(ids, size=cells_per_pert, replace=False)
    cells_keep.extend(sampled)

# Final dataset
adata_sub = adata[cells_keep].copy()
print("Final subset:", adata_sub)
print("Total cells:", adata_sub.n_obs)
