import scanpy as sc
import numpy as np
import pandas as pd
import os
import torch
import scvi
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt

# ==============================
# 0. FIX FOR PYTORCH 2.6+ ERROR
# ==============================
# This block forces torch.load to accept older scVI files
# regardless of the new security defaults.
original_torch_load = torch.load

def permissive_torch_load(*args, **kwargs):
    if 'weights_only' not in kwargs:
        kwargs['weights_only'] = False
    return original_torch_load(*args, **kwargs)

torch.load = permissive_torch_load
print("Applied fix for PyTorch 2.6+ weights_only error.")

# ==============================
# 1. SETUP
# ==============================
GENE_LIST = [
    "DNMT1", "INSIG1", "PHF10", "BRD9", "KLF10", "STAT6", "XRCC4", "MED13L", 
    "DAXX", "FDPS", "MAP3K7", "CASP2", "NDUFB4", "CREG1", "KAT2A", "USP22", 
    "IKBKG", "MED24", "OXA1L", "SALL4", "SIN3B", "STX4", "RRM1", "SV2A", 
    "SRC", "NCK2", "UBE3C", "CAMSAP2", "SMARCA5", "CASP3"
]

# Paths
REAL_DATA_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/adata_all_scvi_50.h5ad"
PRED_DATA_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8/predictions_t_1.00.h5ad"
MODEL_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/scvi_model_dir_50"
SAVE_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8"
os.makedirs(SAVE_DIR, exist_ok=True)
os.makedirs(SAVE_DIR, exist_ok=True)

# ==============================
# 2. LOAD DATA & MODEL
# ==============================
print("Loading datasets...")
adata_real = sc.read_h5ad(REAL_DATA_PATH)
adata_pred = sc.read_h5ad(PRED_DATA_PATH)

# Filter predictions for the final timepoint (t=1.0)
if 'time' in adata_pred.obs:
    print("Filtering predictions for t=1.0...")
    adata_pred = adata_pred[adata_pred.obs['time'] == 1.0].copy()

print(f"Loading scVI Model from {MODEL_DIR}...")
# This should now work thanks to the fix above
model = scvi.model.SCVI.load(MODEL_DIR, adata=adata_real)
print("SUCCESS: Model loaded!")

# ==============================
# 3. CORE DECODER FUNCTION
# ==============================
def decode_and_calculate_lfc(model, z_perturbed, z_control):
    """
    Decodes latent vectors to gene expression and calculates LFC.
    FIX: Now handles scVI models trained with Batch Correction.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.module.to(device)
    model.module.eval()
    
    # Standardize library size (10k counts)
    lib_size = 10e4 
    
    with torch.no_grad():
        # --- 1. Decode PERTURBED Vectors ---
        z_p = torch.tensor(z_perturbed, dtype=torch.float32).to(device)
        lib_p = torch.full((z_p.shape[0], 1), lib_size).to(device)
        
        # FIX: Create dummy batch indices (Batch 0) for the decoder
        # This prevents the "cat not provided" error
        batch_p = torch.zeros((z_p.shape[0], 1), device=device)
        
        # Pass batch_index=batch_p
        out_p = model.module.generative(z=z_p, library=lib_p, batch_index=batch_p)
        expr_p = out_p['px_rate'].cpu().numpy()

        # --- 2. Decode CONTROL Vectors ---
        z_c = torch.tensor(z_control, dtype=torch.float32).to(device)
        lib_c = torch.full((z_c.shape[0], 1), lib_size).to(device)
        
        # Fix for control too
        batch_c = torch.zeros((z_c.shape[0], 1), device=device)
        
        out_c = model.module.generative(z=z_c, library=lib_c, batch_index=batch_c)
        expr_c = out_c['px_rate'].cpu().numpy()

    # 3. Calculate Log Fold Change (LFC)
    mean_log_p = np.log2(expr_p + 1).mean(axis=0)
    mean_log_c = np.log2(expr_c + 1).mean(axis=0)
    
    return mean_log_p - mean_log_c

# ==============================
# 4. PREPARE REAL CONTROL BASELINE
# ==============================
print("Extracting Real Control Latents...")
ctrl_mask = adata_real.obs['target_gene'] == 'non-targeting'
z_control_real = adata_real.obsm['X_scvi'][ctrl_mask]

# ==============================
# 5. MAIN LOOP
# ==============================
metrics = []
print(f"Comparing DGE for {len(GENE_LIST)} genes...")

for gene in tqdm(GENE_LIST):
    if gene not in adata_pred.obs['condition'].values:
        print(f"Skipping {gene} (not in predictions)")
        continue

    # --- PART A: REAL DATA DGE (Using Standard API) ---
    try:
        de_res = model.differential_expression(
            groupby="target_gene", 
            group1=gene, 
            group2="non-targeting", 
            mode="change" 
        )
        real_lfc = de_res['lfc_mean'].sort_index()
    except Exception as e:
        print(f"Could not calculate Real DGE for {gene}: {e}")
        continue

    # --- PART B: PREDICTED DATA DGE (Using Manual Decoding) ---
    pred_mask = adata_pred.obs['condition'] == gene
    z_pred = adata_pred.X[pred_mask]
    
    # Decode vs Real Control Latents
    pred_lfc_values = decode_and_calculate_lfc(model, z_pred, z_control_real)
    pred_lfc = pd.Series(pred_lfc_values, index=adata_real.var_names).sort_index()

    # --- PART C: COMPARE ---
    # Correlation across ALL genes
    slope, intercept, r_val, p_val, std_err = stats.linregress(real_lfc, pred_lfc)
    r2_all = r_val ** 2
    
    # Correlation across TOP 50 Significant Real Genes
    top_50_genes = real_lfc.abs().nlargest(50).index
    slope_top, intercept_top, r_val_top, p_val_top, std_err_top = stats.linregress(
        real_lfc[top_50_genes], 
        pred_lfc[top_50_genes]
    )
    r2_top = r_val_top ** 2

    metrics.append({
        'gene': gene,
        'R2_All_Genes': r2_all,
        'R2_Top50_Genes': r2_top
    })
    
    # Save a scatter plot for the first few genes
    if len(metrics) <= 3:
        plt.figure(figsize=(6, 6))
        # Find genes in common to plot
        common = real_lfc.index.intersection(pred_lfc.index)
        
        plt.scatter(real_lfc, pred_lfc, color='lightgrey', s=5, alpha=0.5, label='All Genes')
        plt.scatter(real_lfc[top_50_genes], pred_lfc[top_50_genes], color='red', s=15, label='Top 50 Real')
        
        plt.plot([-3, 3], [-3, 3], 'k--', alpha=0.5) # Identity line
        plt.xlabel("Real LFC (scVI Bayesian)")
        plt.ylabel("Predicted LFC (Decoded)")
        plt.title(f"{gene}: DGE Correlation (R2={r2_all:.2f})")
        plt.legend()
        plt.savefig(os.path.join(SAVE_DIR, f"scatter_dge_{gene}.png"))
        plt.close()

# ==============================
# 6. SAVE SUMMARY
# ==============================
df_results = pd.DataFrame(metrics)
csv_path = os.path.join(SAVE_DIR, "dge_validation_metrics.csv")
df_results.to_csv(csv_path, index=False)

print("\n--- Summary ---")
print(df_results[['gene', 'R2_All_Genes', 'R2_Top50_Genes']].head())
print(f"Full results saved to: {csv_path}")
print(f"Average R2 (Top 50 Genes): {df_results['R2_Top50_Genes'].mean():.3f}")