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

# ==========================================
# 1. Neural Network Architecture (ResNet-MLP)
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
    """
    Adaptive Layer Normalization.
    Modulates the layer norm statistics based on the condition (time + gene).
    """
    def __init__(self, dim, cond_dim):
        super().__init__()
        self.norm = nn.LayerNorm(dim, elementwise_affine=False)
        self.proj = nn.Sequential(
            nn.SiLU(),
            nn.Linear(cond_dim, dim * 2) 
        )
        # Initialize scaling parameters to zero for stability
        nn.init.zeros_(self.proj[1].weight)
        nn.init.zeros_(self.proj[1].bias)

    def forward(self, x, cond):
        # cond -> [scale, shift]
        scale_shift = self.proj(cond)
        scale, shift = scale_shift.chunk(2, dim=-1)
        return self.norm(x) * (1 + scale) + shift

class ResBlock(nn.Module):
    """
    A simple Residual Block with Adaptive Normalization.
    Structure: x -> AdaLN -> SiLU -> Linear -> SiLU -> Linear -> + x
    """
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
        # Pre-Norm architecture
        x_norm = self.cond_norm(x, cond)
        return x + self.mlp(x_norm)

class VectorFlowNet(nn.Module):
    def __init__(self, input_dim, cond_vec_dim, hidden_dim=256, num_layers=6):
        super().__init__()
        
        # 1. Input Projection
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        
        # 2. Time Embedding (Sinusoidal)
        self.time_mlp = nn.Sequential(
            SinusoidalPosEmb(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        
        # 3. Gene Vector Projection
        self.gene_vec_proj = nn.Linear(cond_vec_dim, hidden_dim)
        
        # --- CALCULATION FOR NEW COND DIMENSION ---
        # We will concatenate Time (hidden_dim) and Gene (hidden_dim)
        combined_cond_dim = hidden_dim * 2

        # 4. The Backbone
        self.blocks = nn.ModuleList([
            # Pass the combined dimension to the blocks
            ResBlock(hidden_dim, cond_dim=combined_cond_dim)
            for _ in range(num_layers)
        ])
        
        # 5. Output Head
        self.final_norm = AdaLN(hidden_dim, hidden_dim)
        self.head = nn.Linear(hidden_dim, input_dim)

    def forward(self, x, t, y_vec):
        """
        x: (batch, input_dim) -> Cell Embedding
        t: (batch, 1) -> Time [0,1]
        y_vec: (batch, cond_vec_dim) -> Gene Embedding
        """
        t_emb = self.time_mlp(t.squeeze(-1))   # (batch, hidden_dim)
        y_emb = self.gene_vec_proj(y_vec)      # (batch, hidden_dim)
        
        # --- MODIFIED: CONCATENATION ---
        # Concatenate along the feature dimension (dim=-1)
        cond = torch.cat([t_emb, y_emb], dim=-1) # (batch, 2 * hidden_dim)
        
        h = self.input_proj(x)
        for block in self.blocks:
            h = block(h, cond)
        h = self.final_norm(h, cond)
        
        return self.head(h)

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


class SingleCellPairDataset(Dataset):
    def __init__(self, x_control, pert_dict, gene_map, active_genes):
        """
        x_control: Tensor of all control cells (N_ctrl, dim)
        pert_dict: Dict {gene_name: Tensor(N_cells, dim)}
        gene_map: Dict {gene_name: Tensor(vec_dim)}
        active_genes: List of genes to include in this dataset
        """
        self.x_control = x_control
        
        # Flatten the perturbation data for easy indexing
        self.pert_cells = []
        self.pert_conds = []
        
        print(f"Preparing dataset for {len(active_genes)} genes...")
        for gene in active_genes:
            if gene not in pert_dict: continue
            
            cells = pert_dict[gene]
            vec = gene_map[gene]
            
            # Add cells
            self.pert_cells.append(cells)
            
            # Repeat gene vector for N cells
            n_cells = cells.shape[0]
            conds = vec.unsqueeze(0).repeat(n_cells, 1) 
            self.pert_conds.append(conds)
            
        if len(self.pert_cells) > 0:
            self.pert_cells = torch.cat(self.pert_cells, dim=0)
            self.pert_conds = torch.cat(self.pert_conds, dim=0)
        else:
            print("Warning: No cells found for the provided genes.")
            self.pert_cells = torch.empty(0)
            self.pert_conds = torch.empty(0)

    def __len__(self):
        return self.pert_cells.shape[0]

    def __getitem__(self, idx):
        # 1. Target (Perturbed Cell)
        x1 = self.pert_cells[idx]
        
        # 2. Condition (Gene Vector)
        y_vec = self.pert_conds[idx]
        
        # 3. Source (Random Control Cell) - Independent Coupling
        rand_idx = torch.randint(0, len(self.x_control), (1,)).item()
        x0 = self.x_control[rand_idx]
        
        return x0, x1, y_vec

# ==========================================
# 3. Training Engine
# ==========================================

class FlowMatchingEngine:
    def __init__(self, model, device, lr=1e-4):
        self.model = model
        self.device = device
        self.optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        self.criterion = nn.MSELoss()

    def compute_loss(self, x_control, x_pert, gene_vecs):
        batch_size = x_control.shape[0]
        t = torch.rand(batch_size, 1, device=self.device)
        
        # OT Path: Linear Interpolation
        x_t = (1 - t) * x_control + t * x_pert
        
        # Target Velocity: Difference vector
        u_target = x_pert - x_control
        
        # Predict Velocity
        u_pred = self.model(x_t, t, gene_vecs)
        
        return self.criterion(u_pred, u_target)

# ==========================================
# 4. Main Execution
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Train OT-CFM (ResNet-MLP)")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to .h5ad file")
    parser.add_argument("--gene_vec_path", type=str, required=True, help="Path to gene embeddings (txt/vec)")
    parser.add_argument("--save_dir", type=str, default="./checkpoints")
    parser.add_argument("--emb_key", type=str, default="X_scvi", help="Key in adata.obsm")
    parser.add_argument("--gene_col", type=str, default="target_gene", help="Key in adata.obs for gene names")
    parser.add_argument("--control_label", type=str, default="non-targeting", help="Label for control cells")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--val_split", type=float, default=0.2)
    parser.add_argument("--hidden_dim", type=int, default=128)
    parser.add_argument("--num_layers", type=int, default=6)

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
    input_dim = X_emb.shape[1]
    print(f"Detected Cell Embedding Dimension: {input_dim}")

    # 2. Load Gene Vectors
    gene_map, cond_vec_dim = load_gene_vectors(args.gene_vec_path)
    if cond_vec_dim is None:
        print("Failed to load gene vectors.")
        return

    # 3. Filter Valid Genes
    available_genes = set(gene_map.keys())
    data_genes = set(np.unique(target_genes))
    # Intersection of genes in data and genes in vector file, excluding control
    valid_pert_genes = list((available_genes & data_genes) - {args.control_label})
    print(f"Found {len(valid_pert_genes)} valid perturbation genes.")

    if len(valid_pert_genes) == 0:
        print("No valid genes found. Check column names and gene vector file.")
        return

    # 4. Prepare Tensors
    print("Preparing tensors...")
    ctrl_mask = (target_genes == args.control_label)
    x_control_all = torch.tensor(X_emb[ctrl_mask], dtype=torch.float32)
    print(f"Number of control cells: {len(x_control_all)}")
    
    pert_dict_all = {}
    for g in valid_pert_genes:
        mask = (target_genes == g)
        if np.sum(mask) > 0:
            pert_dict_all[g] = torch.tensor(X_emb[mask], dtype=torch.float32)

    # 5. Split and Create Datasets
    train_genes, val_genes = train_test_split(valid_pert_genes, test_size=args.val_split, random_state=42)
    print(f"Train genes: {len(train_genes)} | Val genes: {len(val_genes)}")
    
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
    
    best_val_loss = float('inf')
    # 1. Define Early Stopping parameters
    patience = 10         # Stop if no improvement for 10 epochs
    counter = 0           # Tracks how many bad epochs in a row
    best_val_loss = float('inf')

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
        print(f"Epoch {epoch+1}: Train Loss: {avg_train_loss:.5f} | Val Loss: {avg_val_loss:.5f}")

        # Check for improvement
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            torch.save(model.state_dict(), f"{args.save_dir}/best_model.pt")
            counter = 0  # Reset counter because we improved!
            print("  -> New best model saved.")
        else:
            counter += 1 # No improvement, increment counter
            print(f"  -> No improvement. Patience: {counter}/{patience}")
            
            # 3. Stop if patience is exceeded
            if counter >= patience:
                print("\nEarly stopping triggered!")
                break
    
    # Save Config
    config = {
        'input_dim': input_dim,
        'cond_vec_dim': cond_vec_dim,
        'train_genes': train_genes,
        'val_genes': val_genes
    }
    with open(f"{args.save_dir}/config.json", 'w') as f:
        json.dump(config, f)
        
    print("\nTraining Complete. Best model saved.")

if __name__ == "__main__":
    main()