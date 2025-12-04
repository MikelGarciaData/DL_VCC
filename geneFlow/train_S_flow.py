import torch
import torch.nn as nn
import numpy as np
import scanpy as sc
import math
import argparse
import os
import json
from tqdm import tqdm
from sklearn.model_selection import train_test_split

# ==========================================
# 1. Neural Network Architecture (Same as before)
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
# 2. Data Loading Helpers
# ==========================================

def load_gene_vectors(file_path):
    print(f"Loading gene vectors from {file_path}...")
    gene_map = {}
    detected_dim = None

    with open(file_path, 'r') as f:
        # Check header
        first_line = f.readline().strip().split()
        if len(first_line) == 2 and first_line[0].isdigit():
            pass # It was a header
        else:
            f.seek(0)

        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue
            
            gene = parts[0].upper()
            try:
                vec = np.array(parts[1:], dtype=np.float32)
                if detected_dim is None:
                    detected_dim = len(vec)
                    print(f"Detected gene vector dimension: {detected_dim}")
                
                if len(vec) == detected_dim:
                    gene_map[gene] = torch.tensor(vec)
            except ValueError:
                continue
                
    print(f"Loaded {len(gene_map)} gene vectors.")
    return gene_map, detected_dim

# ==========================================
# 3. Training Logic
# ==========================================

class FlowMatchingEngine:
    def __init__(self, model, device):
        self.model = model
        self.device = device
        self.optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
        self.criterion = nn.MSELoss()

    def compute_loss(self, x_control, x_pert, gene_vecs):
        batch_size = x_control.shape[0]
        t = torch.rand(batch_size, 1, device=self.device)
        
        # OT Path
        x_t = (1 - t) * x_control + t * x_pert
        u_target = x_pert - x_control
        u_pred = self.model(x_t, t, gene_vecs)
        
        return self.criterion(u_pred, u_target)

def get_batch(control_tensor, pert_dict, gene_map, active_genes, batch_size, device):
    batch_genes = np.random.choice(active_genes, batch_size)
    x_pert_list = []
    gene_vec_list = []
    
    for gene in batch_genes:
        pert_cells = pert_dict[gene]
        idx = torch.randint(0, len(pert_cells), (1,))
        x_pert_list.append(pert_cells[idx])
        gene_vec_list.append(gene_map[gene])
        
    x_pert = torch.cat(x_pert_list, dim=0).to(device)
    gene_vecs = torch.stack(gene_vec_list, dim=0).to(device)
    
    idx_ctrl = torch.randint(0, len(control_tensor), (batch_size,))
    x_control = control_tensor[idx_ctrl].to(device)
    
    return x_control, x_pert, gene_vecs

# ==========================================
# 4. Main Execution (Simplified)
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Train OT-CFM (Simple Split)")
    parser.add_argument("--adata_path", type=str, required=True)
    parser.add_argument("--gene_vec_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, default="./checkpoints")
    parser.add_argument("--emb_key", type=str, default="X_pca")
    parser.add_argument("--gene_col", type=str, default="target_gene")
    parser.add_argument("--control_label", type=str, default="non-targeting")
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--val_split", type=float, default=0.1, help="Fraction of genes to hold out for validation")
    
    # Architecture
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--num_layers", type=int, default=6)
    parser.add_argument("--num_heads", type=int, default=4)

    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Running on device: {device}")
    
    os.makedirs(args.save_dir, exist_ok=True)

    # 1. Load Data
    print(f"Reading AnnData from {args.adata_path}...")
    adata = sc.read_h5ad(args.adata_path)
    
    X_emb = adata.obsm[args.emb_key]
    target_genes = adata.obs[args.gene_col].values
    input_dim = X_emb.shape[1]
    print(f"Detected Cell Embedding Dimension: {input_dim}")

    # 2. Load Gene Vectors
    gene_map, cond_vec_dim = load_gene_vectors(args.gene_vec_path)
    
    # 3. Filter & Prepare
    available_genes = set(gene_map.keys())
    data_genes = set(np.unique(target_genes))
    valid_pert_genes = list((available_genes & data_genes) - {args.control_label})
    print(f"Found {len(valid_pert_genes)} valid perturbation genes.")

    ctrl_mask = (target_genes == args.control_label)
    x_control_all = torch.tensor(X_emb[ctrl_mask], dtype=torch.float32)
    print(f"Number of control cells: {len(x_control_all)}")
    
    pert_dict_all = {}
    for g in valid_pert_genes:
        mask = (target_genes == g)
        pert_dict_all[g] = torch.tensor(X_emb[mask], dtype=torch.float32)

    # 4. Simple Train/Val Split (No Cross Validation loop)
    # We hold out random genes to check generalization, but only train once.
    train_genes, val_genes = train_test_split(valid_pert_genes, test_size=args.val_split, random_state=42)
    print(f"Training on {len(train_genes)} genes. Validating on {len(val_genes)} genes.")

    # 5. Initialize Model
    model = VectorTransformer(
        input_dim=input_dim,
        cond_vec_dim=cond_vec_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        num_heads=args.num_heads
    ).to(device)
    
    engine = FlowMatchingEngine(model, device)
    
    best_val_loss = float('inf')
    val_loss = 0.0

    # 6. Training Loop
    pbar = tqdm(range(args.epochs), desc="Training")
    for epoch in pbar:
        
        # Train Step
        x0, x1, y_vec = get_batch(
            x_control_all, pert_dict_all, gene_map, train_genes, args.batch_size, device
        )
        
        engine.optimizer.zero_grad()
        loss = engine.compute_loss(x0, x1, y_vec)
        loss.backward()
        engine.optimizer.step()
        train_loss = loss.item()
        
        # Validation Step (Every 50 epochs)
        if epoch % 50 == 0:
            with torch.no_grad():
                x0_v, x1_v, y_vec_v = get_batch(
                    x_control_all, pert_dict_all, gene_map, val_genes, args.batch_size, device
                )
                v_loss = engine.compute_loss(x0_v, x1_v, y_vec_v)
                val_loss = v_loss.item()
            
            # Save Best Model
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(model.state_dict(), f"{args.save_dir}/best_model.pt")
        
        pbar.set_postfix({'T_Loss': f"{train_loss:.4f}", 'V_Loss': f"{val_loss:.4f}"})

    # Save Final Config
    config = {
        'input_dim': input_dim,
        'cond_vec_dim': cond_vec_dim,
        'train_genes': train_genes,
        'val_genes': val_genes
    }
    with open(f"{args.save_dir}/config.json", 'w') as f:
        json.dump(config, f)
        
    print("\nTraining Complete. Best model saved as 'best_model.pt'")

if __name__ == "__main__":
    main()