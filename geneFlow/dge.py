import scanpy as sc
import numpy as np
import pandas as pd
import os
import torch
import scvi
from tqdm import tqdm
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================
# 1. SETUP
# ==============================
GENE_LIST = [
    "DNMT1", "INSIG1", "PHF10", "BRD9", "KLF10", "STAT6", "XRCC4", "MED13L", 
    "DAXX", "FDPS", "MAP3K7", "CASP2", "NDUFB4", "CREG1", "KAT2A", "USP22", 
    "IKBKG", "MED24", "OXA1L", "SALL4", "SIN3B", "STX4", "RRM1", "SV2A", 
    "SRC", "NCK2", "UBE3C", "CAMSAP2", "SMARCA5", "CASP3"
]

REAL_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/adata_all_scvi_50.h5ad"
PRED_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8/predictions_t_1.00.h5ad"
MODEL_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/scvi_model_dir_50"
SAVE_DIR = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/results_experiment_200_epoch_256_8"

# Threshold: We consider the Top 50 real genes as "Truly Differentially Expressed"
TOP_N = 50 

os.makedirs(SAVE_DIR, exist_ok=True)

# ==============================
# 2. LOAD & PREPARE
# ==============================
adata_real = sc.read_h5ad(REAL_PATH)
adata_pred = sc.read_h5ad(PRED_PATH)
if 'time' in adata_pred.obs:
    adata_pred = adata_pred[adata_pred.obs['time'] == 1.0].copy()

model = scvi.model.SCVI.load(MODEL_DIR, adata=adata_real)

# Extract Real Control Latents for baseline comparison
ctrl_mask = adata_real.obs['target_gene'] == 'non-targeting'
z_control_real = adata_real.obsm['X_scvi'][ctrl_mask]

# ==============================
# 3. HELPER: LFC CALCULATION
# ==============================
def get_model_lfc(model, z_perturbed, z_control):
    """Decodes latents and calculates LFC using model decoder."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.module.to(device)
    model.module.eval()
    lib_size = 10e4
    
    with torch.no_grad():
        z_p = torch.tensor(z_perturbed, dtype=torch.float32).to(device)
        out_p = model.module.generative(z=z_p, library=torch.full((z_p.shape[0],1), lib_size).to(device), batch_index=None)
        
        z_c = torch.tensor(z_control, dtype=torch.float32).to(device)
        out_c = model.module.generative(z=z_c, library=torch.full((z_c.shape[0],1), lib_size).to(device), batch_index=None)
        
    mean_log_p = np.log2(out_p['px_rate'].cpu().numpy() + 1).mean(axis=0)
    mean_log_c = np.log2(out_c['px_rate'].cpu().numpy() + 1).mean(axis=0)
    return mean_log_p - mean_log_c

# ==============================
# 4. LOOP: COLLECT SCORES
# ==============================
all_true_labels = []   # 1 if gene is in Top 50 Real, 0 otherwise
all_pred_scores = []   # The absolute LFC predicted by model
predicted_top_genes = [] # For frequency analysis

print(f"Evaluating Classification Metrics (AUROC) on {len(GENE_LIST)} perturbations...")

for gene in tqdm(GENE_LIST):
    if gene not in adata_pred.obs['condition'].values: continue

    # --- A. GROUND TRUTH (Real Data) ---
    try:
        # Get Real LFC from scVI
        de_res = model.differential_expression(groupby="target_gene", group1=gene, group2="non-targeting", mode="change")
        real_lfc = de_res['lfc_mean'].sort_index()
    except:
        continue

    # Define "True Positives": The Top 50 genes with highest ABSOLUTE change
    # Create a binary mask: 1 = Changed, 0 = Didn't Change
    real_abs = real_lfc.abs()
    top_genes_mask = real_abs >= real_abs.nlargest(TOP_N).min()
    true_labels = top_genes_mask.astype(int).values
    
    # --- B. PREDICTION (Model Data) ---
    z_pred = adata_pred.X[adata_pred.obs['condition'] == gene]
    pred_lfc_raw = get_model_lfc(model, z_pred, z_control_real)
    pred_lfc_series = pd.Series(pred_lfc_raw, index=adata_real.var_names).sort_index()
    
    # We use Absolute Predicted LFC as the "Score" (Confidence that it changed)
    pred_scores = pred_lfc_series.abs().values
    
    # --- C. STORE DATA ---
    all_true_labels.extend(true_labels)
    all_pred_scores.extend(pred_scores)
    
    # Store the names of the Top 50 Predicted genes for artifact detection
    top_pred_names = pred_lfc_series.abs().nlargest(TOP_N).index.tolist()
    predicted_top_genes.extend(top_pred_names)

# ==============================
# 5. METRICS: AUROC & CONFUSION
# ==============================
print("\nCalculating Metrics...")

# 1. AUROC (Global)
auroc = roc_auc_score(all_true_labels, all_pred_scores)
fpr, tpr, thresholds = roc_curve(all_true_labels, all_pred_scores)

print(f"Global AUROC: {auroc:.3f}")
# (0.5 = Random Guessing, 1.0 = Perfect, >0.7 = Good)

# 2. Confusion Matrix (at Top 50 Cutoff)
# To make a confusion matrix, we must binarize the predictions.
# We will take the top N predicted scores as "1", rest "0"
# Since we did this per gene in the loop, let's approximate by thresholding the scores
# Note: For strict confusion matrix, we usually do it per perturbation, but global is fine for summary.
binary_preds = np.zeros_like(all_pred_scores)
# This is a simplification; ideally we threshold per perturbation.
# But looking at AUROC is usually more robust.

# ==============================
# 6. ARTIFACT ANALYSIS (Common DGEs)
# ==============================
# Count how often each gene appears in the "Predicted Top 50" across all perturbations
from collections import Counter
counts = Counter(predicted_top_genes)
common_df = pd.DataFrame(counts.most_common(20), columns=['Gene', 'Frequency'])
common_df['Percentage'] = (common_df['Frequency'] / len(GENE_LIST)) * 100

print("\n--- Most Frequent Predicted DE Genes (Possible Artifacts) ---")
print(common_df.head(10))
# If a gene has 100%, it means the model ALWAYS predicts it changes, regardless of perturbation.

# ==============================
# 7. VISUALIZATION
# ==============================
plt.figure(figsize=(12, 5))

# Plot 1: ROC Curve
plt.subplot(1, 2, 1)
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {auroc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Global AUROC (Ability to detect DE genes)')
plt.legend(loc="lower right")

# Plot 2: Artifacts (Bar Chart)
plt.subplot(1, 2, 2)
sns.barplot(x='Percentage', y='Gene', data=common_df, palette='viridis')
plt.title(f"Most Common Predicted DE Genes\n(Frequency across {len(GENE_LIST)} perturbations)")
plt.xlabel("Frequency (%)")
plt.axvline(50, color='r', linestyle='--', label='50% Threshold') # Warning line

plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "auroc_and_artifacts.png"), dpi=300)
plt.close()

# Save Artifact CSV
common_df.to_csv(os.path.join(SAVE_DIR, "artifact_genes.csv"), index=False)