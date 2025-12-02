import scanpy as sc
import scvi
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import os

# =========================================================
# CONFIGURATION
# =========================================================
INPUT_PATH = "data/adata_Training.h5ad"
TRAIN_OUT  = "data/mini_train_random.h5ad"
TEST_OUT   = "data/mini_test_random.h5ad"

# --- SAMPLING SETTINGS ---
N_RANDOM_PERTS = 20     # How many different perturbation genes to pick randomly
CELLS_PER_PERT = 100    # How many cells to keep per perturbation
N_CONTROLS     = 1000   # How many control cells to keep (Need enough for x0)
# -------------------------

# --- MODEL SETTINGS ---
N_LATENT = 30
BATCH_KEY = "batch"
GENE_COL = "target_gene"
CONTROL_LABEL = "non-targeting"

def main():
    print(f"1. Loading {INPUT_PATH}...")
    # Using backed='r' is safer for memory, but for subsetting we load fully
    adata = sc.read_h5ad(INPUT_PATH)
    
    os.makedirs(os.path.dirname(TRAIN_OUT), exist_ok=True)

    # =====================================================
    # STEP 1: RANDOM SUBSAMPLING
    # =====================================================
    print("2. Selecting Random Subset...")
    df = adata.obs
    cells_to_keep = []

    # --- A. Sample Controls ---
    # We need a good base of control cells (x0)
    control_indices = df.index[df[GENE_COL] == CONTROL_LABEL].tolist()
    if len(control_indices) > N_CONTROLS:
        sampled_controls = np.random.choice(control_indices, size=N_CONTROLS, replace=False)
    else:
        sampled_controls = control_indices
    cells_to_keep.extend(sampled_controls)
    print(f"   Selected {len(sampled_controls)} Control cells.")

    # --- B. Identify Valid Perturbations ---
    # 1. Count cells per gene
    gene_counts = df[GENE_COL].value_counts()
    
    # 2. Filter: Must be a perturbation AND have enough cells
    valid_genes = []
    for gene, count in gene_counts.items():
        if gene != CONTROL_LABEL and count >= CELLS_PER_PERT:
            valid_genes.append(gene)
            
    print(f"   Found {len(valid_genes)} genes with at least {CELLS_PER_PERT} cells.")

    # --- C. Randomly Select Genes ---
    # Randomly pick N genes from the valid list
    if len(valid_genes) < N_RANDOM_PERTS:
        print(f"   Warning: Only {len(valid_genes)} valid genes found. Using all of them.")
        chosen_genes = valid_genes
    else:
        chosen_genes = np.random.choice(valid_genes, size=N_RANDOM_PERTS, replace=False)
    
    print(f"   Randomly Selected Genes: {chosen_genes}")

    # --- D. Sample Cells from Chosen Genes ---
    for gene in chosen_genes:
        gene_indices = df.index[df[GENE_COL] == gene].tolist()
        # Randomly sample specific cells for this gene
        sampled_cells = np.random.choice(gene_indices, size=CELLS_PER_PERT, replace=False)
        cells_to_keep.extend(sampled_cells)

    # --- E. Create Subset AnnData ---
    adata_sub = adata[cells_to_keep].copy()
    print(f"   Final Mini-Dataset Shape: {adata_sub.shape}")

    # =====================================================
    # STEP 2: SCVI EMBEDDING
    # =====================================================
    
    print("3. Running scVI Embedding...")
    
    # 1. Setup Layers
    adata_sub.layers["counts"] = adata_sub.X.copy()
    
    # 2. Preprocessing (Normalize & Log1p & HVG)
    # We re-calculate HVG specifically for this subset to ensure good embedding
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    sc.pp.highly_variable_genes(
        adata_sub, 
        n_top_genes=1000, 
        batch_key=BATCH_KEY, 
        subset=True
    )
    
    # 3. Train scVI
    scvi.model.SCVI.setup_anndata(adata_sub, layer="counts", batch_key=BATCH_KEY)
    model = scvi.model.SCVI(adata_sub, n_latent=N_LATENT)
    
    print("   Training scVI (Fast Mode)...")
    model.train(max_epochs=40) 
    
    adata_sub.obsm["X_scvi"] = model.get_latent_representation()
    print(f"   Embedding 'X_scvi' created with shape {adata_sub.obsm['X_scvi'].shape}")

    # =====================================================
    # STEP 3: TRAIN / TEST SPLIT
    # =====================================================
    print("4. Splitting Data (Gene-based)...")
    
    # Split the CHOSEN RANDOM GENES into train/test sets
    # We do NOT split cells randomly; we hold out specific genes.
    if len(chosen_genes) > 1:
        train_genes, test_genes = train_test_split(chosen_genes, test_size=0.2, random_state=42)
    else:
        train_genes = chosen_genes
        test_genes = []

    # Create Masks
    # Train set = All Train Genes + All Controls
    train_mask = adata_sub.obs[GENE_COL].isin(train_genes) | (adata_sub.obs[GENE_COL] == CONTROL_LABEL)
    
    # Test set = All Test Genes + All Controls
    test_mask = adata_sub.obs[GENE_COL].isin(test_genes) | (adata_sub.obs[GENE_COL] == CONTROL_LABEL)
    
    adata_train = adata_sub[train_mask].copy()
    adata_test = adata_sub[test_mask].copy()

    # =====================================================
    # STEP 4: SAVE
    # =====================================================
    print(f"5. Saving files...")
    adata_train.write(TRAIN_OUT)
    adata_test.write(TEST_OUT)
    
    # Save the gene list so you know what is in the test set
    if len(test_genes) > 0:
        np.savetxt("mini_test_genes.txt", test_genes, fmt='%s')

    print("-" * 40)
    print("MINI DATASET CREATION COMPLETE")
    print(f"Train File: {TRAIN_OUT} ({adata_train.shape})")
    print(f"Test File:  {TEST_OUT} ({adata_test.shape})")
    print("-" * 40)

if __name__ == "__main__":
    main()