import torch
import torch.nn as nn
import numpy as np
import scanpy as sc
import pandas as pd
import math
import argparse
import json
import os

# ==========================================
# 1. Architecture (Must match your training)
# ==========================================

class SinusoidalPosEmb(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.dim = dim
    def forward(self, x):
        device = x.device
        half_dim = self.dim // 2
        emb = math.log(10000) / (half_dim - 1)
        emb = torch.exp(torch.arange(half_dim, device=device) * -emb)
        emb = x[:, None] * emb[None, :]
        emb = torch.cat((emb.sin(), emb.cos()), dim=-1)
        return emb

class AdaLN(nn.Module):
    def __init__(self, dim, cond_dim):
        super().__init__()
        self.norm = nn.LayerNorm(dim, elementwise_affine=False)
        self.proj = nn.Sequential(nn.SiLU(), nn.Linear(cond_dim, dim * 2))
        nn.init.zeros_(self.proj[1].weight)
        nn.init.zeros_(self.proj[1].bias)
    def forward(self, x, cond):
        scale_shift = self.proj(cond).unsqueeze(1)
        scale, shift = scale_shift.chunk(2, dim=-1)
        return self.norm(x) * (1 + scale) + shift

class TransformerBlock(nn.Module):
    def __init__(self, dim, num_heads, cond_dim, mlp_ratio=4):
        super().__init__()
        self.attn_norm = AdaLN(dim, cond_dim)
        self.attn = nn.MultiheadAttention(dim, num_heads, batch_first=True)
        self.mlp_norm = AdaLN(dim, cond_dim)
        self.mlp = nn.Sequential(nn.Linear(dim, dim * mlp_ratio), nn.GELU(), nn.Linear(dim * mlp_ratio, dim))
    def forward(self, x, cond):
        x_norm = self.attn_norm(x, cond)
        attn_out, _ = self.attn(x_norm, x_norm, x_norm)
        x = x + attn_out
        x_norm = self.mlp_norm(x, cond)
        x = x + self.mlp(x_norm)
        return x

class VectorTransformer(nn.Module):
    def __init__(self, input_dim, cond_vec_dim, hidden_dim=256, num_layers=6, num_heads=4):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        self.time_mlp = nn.Sequential(
            SinusoidalPosEmb(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.gene_vec_proj = nn.Linear(cond_vec_dim, hidden_dim)
        self.blocks = nn.ModuleList([
            TransformerBlock(hidden_dim, num_heads, cond_dim=hidden_dim)
            for _ in range(num_layers)
        ])
        self.final_norm = AdaLN(hidden_dim, hidden_dim)
        self.head = nn.Linear(hidden_dim, input_dim)

    def forward(self, x, t, y_vec):
        t_emb = self.time_mlp(t.squeeze(-1))
        y_emb = self.gene_vec_proj(y_vec)
        cond = t_emb + y_emb
        h = self.input_proj(x).unsqueeze(1)
        for block in self.blocks:
            h = block(h, cond)
        h = self.final_norm(h, cond)
        return self.head(h).squeeze(1)

# ==========================================
# 2. ODE Solver (Prediction Logic)
# ==========================================

@torch.no_grad()
def ode_solve(model, x0, gene_vec, steps=20, device='cuda'):
    """
    Moves cells from Control (x0) to Perturbed state over 'steps'.
    """
    model.eval()
    x = x0.clone().to(device)
    gene_vec = gene_vec.to(device)
    
    # Broadcast gene vector to match batch size
    if gene_vec.dim() == 1:
        gene_vec = gene_vec.unsqueeze(0).expand(x.shape[0], -1)
    
    dt = 1.0 / steps
    for i in range(steps):
        t_val = i / steps
        t = torch.full((x.shape[0], 1), t_val, device=device)
        v = model(x, t, gene_vec)
        x = x + v * dt
    return x.cpu().numpy()

# ==========================================
# 3. Helpers
# ==========================================

def load_gene_vectors(file_path):
    print("Loading gene vectors...")
    gene_map = {}
    with open(file_path, 'r') as f:
        # Check header
        pos = f.tell()
        line = f.readline().split()
        if not (len(line) == 2 and line[0].isdigit()):
            f.seek(pos)
        
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue
            try:
                gene_map[parts[0].upper()] = torch.tensor(np.array(parts[1:], dtype=np.float32))
            except: continue
    return gene_map

# ==========================================
# 4. Main Script
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Multi-Gene Perturbation Predictor")
    
    # Required Paths
    parser.add_argument("--checkpoint_dir", type=str, required=True, help="Folder containing best_model.pt and config.json")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to .h5ad containing Control cells")
    parser.add_argument("--gene_vec_path", type=str, required=True, help="Path to gene_vectors.txt")
    parser.add_argument("--output_path", type=str, default="predictions.h5ad", help="Where to save the result")
    
    # Inputs
    parser.add_argument("--genes", type=str, required=True, help="Comma-separated list of genes to perturb (e.g. TP53,MYC)")
    parser.add_argument("--n_cells", type=int, default=500, help="How many cells to predict per gene")
    
    # Data Config (UPDATED DEFAULTS)
    parser.add_argument("--emb_key", type=str, default="X_scvi", help="Key in .obsm where embeddings are stored")
    parser.add_argument("--control_label", type=str, default="non-targeting", help="Label in gene_col that identifies controls")
    parser.add_argument("--gene_col", type=str, default="target_gene", help="Column name for gene labels")
    
    # Model Architecture (MUST MATCH TRAINING)
    parser.add_argument("--hidden_dim", type=int, default=128)
    parser.add_argument("--num_layers", type=int, default=4)
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--steps", type=int, default=20, help="Integration steps")

    args = parser.parse_args()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Running on {device}")

    # 1. Load Model
    config_path = f"{args.checkpoint_dir}/config.json"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config not found at {config_path}")
        
    with open(config_path, 'r') as f:
        config = json.load(f)

    model = VectorTransformer(
        input_dim=config['input_dim'],
        cond_vec_dim=config['cond_vec_dim'],
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        num_heads=args.num_heads
    ).to(device)
    
    model_path = f"{args.checkpoint_dir}/best_model.pt"
    print(f"Loading weights from {model_path}...")
    model.load_state_dict(torch.load(model_path, map_location=device))
    
    # 2. Load Data & Vectors
    print("Loading data...")
    adata = sc.read_h5ad(args.adata_path)
    gene_map = load_gene_vectors(args.gene_vec_path)
    
    # 3. Extract Control Pool
    if args.gene_col in adata.obs:
        ctrl_mask = (adata.obs[args.gene_col] == args.control_label)
        X_control_pool = adata.obsm[args.emb_key][ctrl_mask]
    else:
        print(f"Warning: Column '{args.gene_col}' not found. Using ALL cells as control pool.")
        X_control_pool = adata.obsm[args.emb_key]

    n_controls = len(X_control_pool)
    if n_controls == 0:
        raise ValueError(f"No control cells found with label '{args.control_label}' in column '{args.gene_col}'")
    
    print(f"Control pool size: {n_controls} cells")

    # 4. Prediction Loop
    gene_list = [g.strip() for g in args.genes.split(',')]
    print(f"Generating {args.n_cells} cells for each of: {gene_list}")
    
    output_adatas = []
    
    for gene in gene_list:
        gene_u = gene.upper()
        
        # A. Check Vector
        if gene_u not in gene_map:
            print(f"Skipping {gene}: No vector found in {args.gene_vec_path}")
            continue
            
        # B. Sample Control Cells (Randomly)
        # If we want more cells than exist in control, we must replace=True
        replace = args.n_cells > n_controls
        indices = np.random.choice(n_controls, args.n_cells, replace=replace)
        
        x0_batch = torch.tensor(X_control_pool[indices], dtype=torch.float32)
        y_vec = gene_map[gene_u]
        
        # C. Predict
        x_pred = ode_solve(model, x0_batch, y_vec, steps=args.steps, device=device)
        
        # D. Store in Temporary AnnData
        # We create an AnnData for this specific gene's predictions
        tmp_adata = sc.AnnData(shape=(args.n_cells, 1)) # Dummy shape
        tmp_adata.obsm[args.emb_key] = x_pred
        
        # Use the requested column name (target_gene)
        tmp_adata.obs[args.gene_col] = gene         
        tmp_adata.obs['origin_type'] = 'predicted'
        
        output_adatas.append(tmp_adata)
        print(f"  -> Generated {args.n_cells} cells for {gene}")

    # 5. Merge and Save
    if len(output_adatas) > 0:
        final_adata = sc.concat(output_adatas)
        final_adata.write(args.output_path)
        print(f"\nSaved combined predictions to: {args.output_path}")
        print(f"Total cells predicted: {final_adata.n_obs}")
        print(f"Data stored in .obsm['{args.emb_key}']")
    else:
        print("No predictions were generated (check gene names).")

if __name__ == "__main__":
    main()