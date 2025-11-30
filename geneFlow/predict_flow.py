import torch
import torch.nn as nn
import numpy as np
import scanpy as sc
import math
import argparse
import json
import os

# ==========================================
# 1. Model Architecture (Must match training)
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
        self.proj = nn.Sequential(
            nn.SiLU(),
            nn.Linear(cond_dim, dim * 2)
        )
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
        self.mlp = nn.Sequential(
            nn.Linear(dim, dim * mlp_ratio),
            nn.GELU(),
            nn.Linear(dim * mlp_ratio, dim)
        )

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
# 2. ODE Solver (Euler Method)
# ==========================================

def simulate_flow(model, x_start, gene_vec, steps=100, device='cpu'):
    """
    Simulates the flow from t=0 to t=1 using Euler integration.
    x_start: (Batch, Dim) - Control cells
    gene_vec: (Dim,) - The vector for the perturbation to simulate
    """
    model.eval()
    
    # Batch the gene vector to match x_start
    batch_size = x_start.shape[0]
    y_vec_batch = gene_vec.unsqueeze(0).repeat(batch_size, 1).to(device)
    
    x_t = x_start.clone().to(device)
    dt = 1.0 / steps
    
    with torch.no_grad():
        for i in range(steps):
            t_val = i * dt
            t_batch = torch.full((batch_size, 1), t_val, device=device)
            
            # Predict velocity (dx/dt)
            velocity = model(x_t, t_batch, y_vec_batch)
            
            # Update state: x_{t+1} = x_t + v * dt
            x_t = x_t + velocity * dt
            
    return x_t

# ==========================================
# 3. Main Logic
# ==========================================

def load_target_gene_vector(file_path, target_gene):
    """Parses text file to find specific gene vector"""
    target_gene = target_gene.upper()
    with open(file_path, 'r') as f:
        # Skip header logic similar to training
        first_line = f.readline().strip().split()
        if not (len(first_line) == 2 and first_line[0].isdigit()):
            f.seek(0)
            
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue
            if parts[0].upper() == target_gene:
                return torch.tensor(np.array(parts[1:], dtype=np.float32))
    
    raise ValueError(f"Gene {target_gene} not found in {file_path}")

def main():
    parser = argparse.ArgumentParser(description="Predict with OT-CFM")
    parser.add_argument("--checkpoint", type=str, required=True, help="Path to .pt model file")
    parser.add_argument("--config", type=str, required=True, help="Path to .json config file")
    parser.add_argument("--gene_vec_path", type=str, required=True, help="Path to gene2vec txt file")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to .h5ad with control cells")
    parser.add_argument("--control_label", type=str, default="non-targeting", help="Label for control cells")
    parser.add_argument("--target_gene", type=str, required=True, help="Gene symbol to simulate (e.g., TP53)")
    parser.add_argument("--num_cells", type=int, default=100, help="Number of control cells to simulate")
    parser.add_argument("--steps", type=int, default=100, help="Number of ODE steps (higher = more accurate)")
    parser.add_argument("--output_path", type=str, default="prediction.npy", help="Where to save result")
    
    # Model Architecture args (must match training!)
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--num_layers", type=int, default=6)
    parser.add_argument("--num_heads", type=int, default=4)

    args = parser.parse_args()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # 1. Load Config & Model
    print(f"Loading config from {args.config}...")
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    print(f"Loading weights from {args.checkpoint}...")
    model = VectorTransformer(
        input_dim=config['input_dim'],
        cond_vec_dim=config['cond_vec_dim'],
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        num_heads=args.num_heads
    ).to(device)
    
    model.load_state_dict(torch.load(args.checkpoint, map_location=device))
    
    # 2. Get Gene Vector
    print(f"Looking up vector for {args.target_gene}...")
    target_vec = load_target_gene_vector(args.gene_vec_path, args.target_gene)
    print(f"Found vector of dim {len(target_vec)}")

    # 3. Get Control Cells
    print(f"Loading control cells from {args.adata_path}...")
    adata = sc.read_h5ad(args.adata_path)
    
    # Assuming embeddings are in X_pca, change if you saved config differently
    # Ideally, save the emb_key in config too, but for now we assume X_pca or check adata
    emb_key = "X_pca" 
    if emb_key not in adata.obsm:
         # Fallback try
         emb_key = list(adata.obsm.keys())[0]
    
    print(f"Using embedding: {emb_key}")
    X = adata.obsm[emb_key]
    genes = adata.obs["target_gene"].values # Adjust column name if needed
    
    ctrl_mask = (genes == args.control_label)
    all_ctrls = X[ctrl_mask]
    
    if len(all_ctrls) == 0:
        raise ValueError(f"No control cells found with label {args.control_label}")
        
    # Sample random controls
    indices = np.random.choice(len(all_ctrls), size=min(args.num_cells, len(all_ctrls)), replace=False)
    x_control = torch.tensor(all_ctrls[indices], dtype=torch.float32)
    print(f"Simulating {len(x_control)} cells...")

    # 4. Run Simulation (ODE Solver)
    print("Starting ODE simulation...")
    predicted_cells = simulate_flow(
        model, 
        x_control, 
        target_vec, 
        steps=args.steps, 
        device=device
    )
    
    # 5. Save Results
    result_np = predicted_cells.cpu().numpy()
    np.save(args.output_path, result_np)
    print(f"Saved predictions to {args.output_path}")
    print(f"Shape: {result_np.shape}")

if __name__ == "__main__":
    main()