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
from torch.utils.data import Dataset, DataLoader
from collections import Counter

# ==========================================
# 1. Neural Network Architecture
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
        scale_shift = self.proj(cond)
        scale, shift = scale_shift.chunk(2, dim=-1)
        return self.norm(x) * (1 + scale) + shift

class ResBlock(nn.Module):
    def __init__(self, dim, cond_dim, dropout=0.0):
        super().__init__()
        self.cond_norm = AdaLN(dim, cond_dim)
        self.mlp = nn.Sequential(
            nn.SiLU(),
            nn.Linear(dim, dim),
            nn.Dropout(dropout),
            nn.SiLU(),
            nn.Linear(dim, dim),
            nn.Dropout(dropout)
        )

    def forward(self, x, cond):
        x_norm = self.cond_norm(x, cond)
        return x + self.mlp(x_norm)

class VectorFlowNet(nn.Module):
    def __init__(self, input_dim, cond_vec_dim, hidden_dim=256, num_layers=6):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        self.time_mlp = nn.Sequential(
            SinusoidalPosEmb(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.gene_vec_proj = nn.Linear(cond_vec_dim, hidden_dim)
        
        combined_cond_dim = hidden_dim * 2

        self.blocks = nn.ModuleList([
            ResBlock(hidden_dim, cond_dim=combined_cond_dim)
            for _ in range(num_layers)
        ])
        
        self.final_norm = AdaLN(hidden_dim, combined_cond_dim)
        self.head = nn.Linear(hidden_dim, input_dim)

    def forward(self, x, t, y_vec):
        t_emb = self.time_mlp(t.squeeze(-1))
        y_emb = self.gene_vec_proj(y_vec)
        cond = torch.cat([t_emb, y_emb], dim=-1)
        
        h = self.input_proj(x)
        for block in self.blocks:
            h = block(h, cond)
        h = self.final_norm(h, cond)
        return self.head(h)

# ==========================================
# 2. Helpers: Loading, Solving, Dataset
# ==========================================

def load_gene_vectors(file_path):
    print(f"Loading gene vectors from {file_path}...")
    gene_map = {}
    detected_dim = None

    with open(file_path, 'r') as f:
        first_line = f.readline().strip().split()
        if len(first_line) == 2 and first_line[0].isdigit():
            pass 
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
                if len(vec) == detected_dim:
                    gene_map[gene] = torch.tensor(vec)
            except ValueError:
                continue
    print(f"Loaded {len(gene_map)} gene vectors.")
    return gene_map, detected_dim

@torch.no_grad()
def ode_solve_trajectory(model, x_control, gene_vec, steps=20):
    """
    Integrates from t=0 to t=1, saving snapshots at t=0, 0.25, 0.5, 0.75, 1.0.
    """
    # Ensure steps is a multiple of 4 so we hit exact quarter points
    if steps % 4 != 0:
        steps = ((steps // 4) + 1) * 4
    
    dt = 1.0 / steps
    x_t = x_control.clone()
    batch_size = x_control.shape[0]
    device = x_control.device

    # We want snapshots at these exact time fractions
    save_points = {0.0, 0.25, 0.5, 0.75, 1.0}
    snapshots = {0.0: x_t.clone()} # Save t=0

    for i in range(steps):
        t_val = i / steps
        t = torch.full((batch_size, 1), t_val, device=device)
        
        velocity = model(x_t, t, gene_vec)
        x_t = x_t + velocity * dt
        
        # Calculate the time at the NEXT step
        next_t = (i + 1) / steps
        
        # If the next step is one of our target points, save it
        # (Using a small epsilon for float comparison safety)
        if any(abs(next_t - target) < 1e-6 for target in save_points):
            # Round the key to nice float (e.g. 0.25)
            clean_time = round(next_t, 2)
            snapshots[clean_time] = x_t.clone()
            
    return snapshots

class SingleCellPairDataset(Dataset):
    def __init__(self, x_control, pert_dict, gene_map, active_genes):
        self.x_control = x_control
        self.pert_cells = []
        self.pert_conds = []
        
        print(f"Preparing dataset for {len(active_genes)} genes...")
        for gene in active_genes:
            if gene not in pert_dict: continue
            cells = pert_dict[gene]
            vec = gene_map[gene]
            
            self.pert_cells.append(cells)
            n_cells = cells.shape[0]
            conds = vec.unsqueeze(0).repeat(n_cells, 1) 
            self.pert_conds.append(conds)
            
        if len(self.pert_cells) > 0:
            self.pert_cells = torch.cat(self.pert_cells, dim=0)
            self.pert_conds = torch.cat(self.pert_conds, dim=0)
        else:
            self.pert_cells = torch.empty(0)
            self.pert_conds = torch.empty(0)

    def __len__(self):
        return self.pert_cells.shape[0]

    def __getitem__(self, idx):
        x1 = self.pert_cells[idx]
        y_vec = self.pert_conds[idx]
        rand_idx = torch.randint(0, len(self.x_control), (1,)).item()
        x0 = self.x_control[rand_idx]
        return x0, x1, y_vec

class FlowMatchingEngine:
    def __init__(self, model, device, lr=1e-4):
        self.model = model
        self.device = device
        self.optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        self.criterion = nn.MSELoss()

    def compute_loss(self, x_control, x_pert, gene_vecs):
        batch_size = x_control.shape[0]
        t = torch.rand(batch_size, 1, device=self.device)
        x_t = (1 - t) * x_control + t * x_pert
        u_target = x_pert - x_control
        u_pred = self.model(x_t, t, gene_vecs)
        return self.criterion(u_pred, u_target)

# ==========================================
# 3. Main Execution
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Train OT-CFM (ResNet-MLP)")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to .h5ad file")
    parser.add_argument("--gene_vec_path", type=str, required=True, help="Path to gene embeddings")
    parser.add_argument("--save_dir", type=str, default="./checkpoints")
    parser.add_argument("--emb_key", type=str, default="X_scvi", help="Key in adata.obsm")
    parser.add_argument("--gene_col", type=str, default="target_gene", help="Key in adata.obs")
    parser.add_argument("--control_label", type=str, default="non-targeting", help="Label for control cells")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--val_split", type=float, default=0.2)
    parser.add_argument("--hidden_dim", type=int, default=128)
    parser.add_argument("--num_layers", type=int, default=6)
    parser.add_argument("--pred_steps", type=int, default=20, help="Number of ODE solver steps for prediction should be divisible by 4")

    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Running on device: {device}")
    os.makedirs(args.save_dir, exist_ok=True)

    # 1. Load Data
    print(f"Reading AnnData from {args.adata_path}...")
    try:
        adata = sc.read_h5ad(args.adata_path)
    except Exception as e:
        print(f"Error loading AnnData: {e}")
        return

    X_emb = adata.obsm[args.emb_key]
    target_genes = adata.obs[args.gene_col].values.astype(str)
    gene_counts = Counter(target_genes)
    input_dim = X_emb.shape[1]
    
    # 2. Load Gene Vectors
    gene_map, cond_vec_dim = load_gene_vectors(args.gene_vec_path)
    if cond_vec_dim is None: return

    # 3. Filter Valid Genes
    available_genes = set(gene_map.keys())
    data_genes = set(np.unique(target_genes))
    valid_pert_genes = list((available_genes & data_genes) - {args.control_label})
    
    # 4. Prepare Tensors
    print("Preparing tensors...")
    ctrl_mask = (target_genes == args.control_label)
    x_control_all = torch.tensor(X_emb[ctrl_mask], dtype=torch.float32)
    print(f"Total control cells available: {len(x_control_all)}")
    
    pert_dict_all = {}
    for g in valid_pert_genes:
        mask = (target_genes == g)
        if np.sum(mask) > 0:
            pert_dict_all[g] = torch.tensor(X_emb[mask], dtype=torch.float32)

    # 5. Split and Create Datasets
    train_genes, val_genes = train_test_split(valid_pert_genes, test_size=args.val_split, random_state=42)
    
    train_ds = SingleCellPairDataset(x_control_all, pert_dict_all, gene_map, train_genes)
    val_ds = SingleCellPairDataset(x_control_all, pert_dict_all, gene_map, val_genes)
    
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, num_workers=0)

    # 6. Initialize Model
    model = VectorFlowNet(
        input_dim=input_dim,
        cond_vec_dim=cond_vec_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers
    ).to(device)
    
    engine = FlowMatchingEngine(model, device)
    
    patience = 10
    counter = 0
    best_val_loss = float('inf')
    
    # --- NEW: History Dictionary ---
    loss_history = {
        'train_loss': [],
        'val_loss': []
    }

    # 7. Training Loop
    print("\nStarting Training...")
    for epoch in range(args.epochs):
        model.train()
        train_loss_sum = 0
        batches = 0
        
        pbar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{args.epochs}")
        
        for x0, x1, y_vec in pbar:
            x0, x1, y_vec = x0.to(device), x1.to(device), y_vec.to(device)
            engine.optimizer.zero_grad()
            loss = engine.compute_loss(x0, x1, y_vec)
            loss.backward()
            engine.optimizer.step()
            train_loss_sum += loss.item()
            batches += 1
            pbar.set_postfix({'Loss': f"{loss.item():.4f}"})
            
        avg_train_loss = train_loss_sum / batches if batches > 0 else 0
        
        # Validation
        model.eval()
        val_loss_sum = 0
        val_batches = 0
        with torch.no_grad():
            for x0_v, x1_v, y_vec_v in val_loader:
                x0_v, x1_v, y_vec_v = x0_v.to(device), x1_v.to(device), y_vec_v.to(device)
                v_loss = engine.compute_loss(x0_v, x1_v, y_vec_v)
                val_loss_sum += v_loss.item()
                val_batches += 1
        
        avg_val_loss = val_loss_sum / val_batches if val_batches > 0 else 0
        
        # --- NEW: Record History ---
        loss_history['train_loss'].append(avg_train_loss)
        loss_history['val_loss'].append(avg_val_loss)
        
        print(f"Epoch {epoch+1}: Train Loss: {avg_train_loss:.5f} | Val Loss: {avg_val_loss:.5f}")

        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            torch.save(model.state_dict(), f"{args.save_dir}/best_model.pt")
            counter = 0
            print("  -> New best model saved.")
        else:
            counter += 1
            print(f"  -> No improvement. Patience: {counter}/{patience}")
            if counter >= patience:
                print("\nEarly stopping triggered!")
                break
    
    # 8. Save Configuration & History
    config = {
        'input_dim': int(input_dim),
        'cond_vec_dim': int(cond_vec_dim),
        'train_genes': train_genes,
        'val_genes': val_genes
    }
    with open(f"{args.save_dir}/config.json", 'w') as f:
        json.dump(config, f, indent=4)
        
    # --- NEW: Save History ---
    with open(f"{args.save_dir}/loss_history.json", 'w') as f:
        json.dump(loss_history, f, indent=4)
        
    print("Configuration and loss history saved.")

    # ==========================================
    # 9. PREDICTION PHASE (5 Separate Files)
    # ==========================================
    print("\n==========================================")
    print("Starting Trajectory Prediction (Saving 5 separate files)...")
    
    model.load_state_dict(torch.load(f"{args.save_dir}/best_model.pt"))
    model.eval()
    
    # Define the time points we want to save
    time_points = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    # Create a dictionary to hold lists for EACH time point separately
    # Structure: { 0.0: {'embs': [], 'genes': []}, 0.25: {...}, ... }
    data_by_time = {t: {'embs': [], 'genes': []} for t in time_points}
    
    print(f"Simulating {len(val_genes)} genes...")

    for gene in tqdm(val_genes, desc="Predicting Trajectories"):
        if gene not in gene_map: continue
        
        n_required = gene_counts[gene]
        if n_required == 0: continue
        
        # 1. Prepare Inputs
        indices = torch.randint(0, len(x_control_all), (n_required,))
        x_curr = x_control_all[indices].to(device)
        
        g_vec = gene_map[gene].to(device)
        g_vec_batch = g_vec.unsqueeze(0).repeat(n_required, 1)
        
        # 2. Run Solver (Returns dict: {0.0: tensor, 0.25: tensor...})
        snapshots = ode_solve_trajectory(model, x_curr, g_vec_batch, steps=args.steps)
        
        # 3. Sort data into the correct buckets
        for t in time_points:
            if t in snapshots:
                # Store the embedding
                data_by_time[t]['embs'].append(snapshots[t].cpu().numpy())
                # Store the gene label
                data_by_time[t]['genes'].extend([gene] * n_required)

    # 4. Save 5 separate files
    print("\nSaving files...")
    for t in time_points:
        container = data_by_time[t]
        
        if len(container['embs']) > 0:
            # Concatenate all genes for this specific time point
            X_full = np.concatenate(container['embs'], axis=0)
            
            # Create AnnData
            adata_t = sc.AnnData(X=X_full)
            adata_t.obs['condition'] = container['genes']
            adata_t.obs['type'] = 'predicted'
            adata_t.obs['time_point'] = t  # Metadata tag
            
            # Format filename (e.g., predictions_t_0.25.h5ad)
            fname = f"predictions_t_{t:.2f}.h5ad"
            save_path = os.path.join(args.save_dir, fname)
            
            adata_t.write(save_path)
            print(f"Saved: {save_path} | Shape: {adata_t.shape}")
        else:
            print(f"Warning: No data found for t={t}")

if __name__ == "__main__":
    main()