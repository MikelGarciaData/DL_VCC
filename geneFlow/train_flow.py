import torch
import torch.nn as nn
import numpy as np
import scanpy as sc
import math
import argparse
import os
import json
from tqdm import tqdm
from sklearn.model_selection import KFold

# ==========================================
# 1. Neural Network Architecture (Dynamic)
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
    """
    Loads gene vectors. Auto-detects dimension based on the first valid line.
    Assumes format: GENE_SYMBOL value1 value2 ...
    """
    print(f"Loading gene vectors from {file_path}...")
    gene_map = {}
    detected_dim = None

    with open(file_path, 'r') as f:
        # Try to skip header if it exists (heuristic: 2 numbers)
        first_line = f.readline().strip().split()
        if len(first_line) == 2 and first_line[0].isdigit():
            pass # It was a header
        else:
            f.seek(0) # Reset if no header

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
# 3. Training & Validation Logic
# ==========================================

class FlowMatchingEngine:
    def __init__(self, model, device):
        self.model = model
        self.device = device
        self.optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
        self.criterion = nn.MSELoss()

    def compute_loss(self, x_control, x_pert, gene_vecs):
        """
        Computes conditional flow matching loss.
        x_control: (B, Dim) - Source
        x_pert: (B, Dim) - Target
        gene_vecs: (B, Cond_Dim) - Conditioning
        """
        batch_size = x_control.shape[0]
        t = torch.rand(batch_size, 1, device=self.device)
        
        # OT Path Interpolation
        # x_t = (1-t)x0 + t*x1
        x_t = (1 - t) * x_control + t * x_pert
        
        # Target Velocity
        # u = x1 - x0
        u_target = x_pert - x_control
        
        # Prediction
        u_pred = self.model(x_t, t, gene_vecs)
        
        return self.criterion(u_pred, u_target)

    def run_epoch(self, loader, is_train=True):
        if is_train:
            self.model.train()
        else:
            self.model.eval()
            
        total_loss = 0
        
        # We assume loader yields (x_control_batch, x_pert_batch, gene_vec_batch)
        # But since our data structure is custom, we iterate manually in the main loop
        # This function is a placeholder for standard loader logic if you used DataLoader
        return 0

# ==========================================
# 4. Main HPC Execution
# ==========================================

def get_batch(control_tensor, pert_dict, gene_map, active_genes, batch_size, device):
    """
    Custom sampler that balances batches across active perturbation genes.
    """
    # 1. Sample Gene Labels
    batch_genes = np.random.choice(active_genes, batch_size)
    
    x_pert_list = []
    gene_vec_list = []
    
    for gene in batch_genes:
        # Get a random cell from this gene's perturbation group
        pert_cells = pert_dict[gene]
        idx = torch.randint(0, len(pert_cells), (1,))
        x_pert_list.append(pert_cells[idx])
        
        # Get the gene vector
        gene_vec_list.append(gene_map[gene])
        
    x_pert = torch.cat(x_pert_list, dim=0).to(device)
    gene_vecs = torch.stack(gene_vec_list, dim=0).to(device)
    
    # 2. Sample Control Cells (Source)
    # We sample randomly from controls. 
    # OT-CFM learns the coupling, but sampling random controls is standard for training.
    idx_ctrl = torch.randint(0, len(control_tensor), (batch_size,))
    x_control = control_tensor[idx_ctrl].to(device)
    
    return x_control, x_pert, gene_vecs

def main():
    parser = argparse.ArgumentParser(description="Train OT-CFM on HPC")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to .h5ad file")
    parser.add_argument("--gene_vec_path", type=str, required=True, help="Path to gene2vec txt file")
    parser.add_argument("--save_dir", type=str, default="./checkpoints", help="Directory to save weights")
    parser.add_argument("--emb_key", type=str, default="X_pca", help="Key in obsm to use (e.g., X_pca, X_scvi)")
    parser.add_argument("--gene_col", type=str, default="target_gene", help="Column name for perturbation labels")
    parser.add_argument("--control_label", type=str, default="non-targeting", help="Label for control cells")
    parser.add_argument("--k_folds", type=int, default=5, help="Number of cross-validation folds")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of training steps/epochs")
    parser.add_argument("--batch_size", type=int, default=128, help="Batch size")
    
    # Model Architecture Arguments
    parser.add_argument("--hidden_dim", type=int, default=256, help="Hidden dimension size of the Transformer")
    parser.add_argument("--num_layers", type=int, default=6, help="Number of Transformer blocks")
    parser.add_argument("--num_heads", type=int, default=4, help="Number of attention heads")

    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Running on device: {device}")
    
    os.makedirs(args.save_dir, exist_ok=True)

    # 1. Load Data
    print(f"Reading AnnData from {args.adata_path}...")
    adata = sc.read_h5ad(args.adata_path)
    
    # 2. Extract Embeddings
    if args.emb_key not in adata.obsm:
        raise ValueError(f"Embedding key {args.emb_key} not found in adata.obsm")
    
    X_emb = adata.obsm[args.emb_key]
    target_genes = adata.obs[args.gene_col].values
    
    input_dim = X_emb.shape[1]
    print(f"Detected Cell Embedding Dimension: {input_dim}")

    # 3. Load Gene Vectors
    gene_map, cond_vec_dim = load_gene_vectors(args.gene_vec_path)
    
    # 4. Filter Data
    # Only keep perturbations that we have a vector for
    available_genes = set(gene_map.keys())
    data_genes = set(np.unique(target_genes))
    
    # Genes that are in data AND have a vector (excluding control)
    valid_pert_genes = list((available_genes & data_genes) - {args.control_label})
    valid_pert_genes.sort() # Ensure deterministic order for KFold
    print(f"Found {len(valid_pert_genes)} valid perturbation genes with vectors.")

    # 5. Prepare Tensors
    # Controls
    ctrl_mask = (target_genes == args.control_label)
    x_control_all = torch.tensor(X_emb[ctrl_mask], dtype=torch.float32)
    print(f"Number of control cells: {len(x_control_all)}")
    
    # Perturbations Dict (Load to CPU memory first, move to GPU in batch)
    pert_dict_all = {}
    for g in valid_pert_genes:
        mask = (target_genes == g)
        pert_dict_all[g] = torch.tensor(X_emb[mask], dtype=torch.float32)

    # 6. Cross Validation Loop
    kf = KFold(n_splits=args.k_folds, shuffle=True, random_state=42)
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(valid_pert_genes)):
        print(f"\n=== Starting Fold {fold+1}/{args.k_folds} ===")
        
        # Split genes
        train_genes = [valid_pert_genes[i] for i in train_idx]
        val_genes = [valid_pert_genes[i] for i in val_idx]
        
        print(f"Train Genes: {len(train_genes)} | Val Genes: {len(val_genes)}")
        
        # Initialize Model
        model = VectorTransformer(
            input_dim=input_dim,
            cond_vec_dim=cond_vec_dim, # Dynamic from file
            hidden_dim=args.hidden_dim,
            num_layers=args.num_layers,
            num_heads=args.num_heads
        ).to(device)
        
        engine = FlowMatchingEngine(model, device)
        
        # Save validation loss history
        history = {'train_loss': [], 'val_loss': []}
        best_val_loss = float('inf')
        
        # Training Loop
        pbar = tqdm(range(args.epochs), desc=f"Fold {fold+1}")
        for epoch in pbar:
            
            # --- TRAIN STEP ---
            # Sample batch from TRAIN genes
            x0, x1, y_vec = get_batch(
                x_control_all, pert_dict_all, gene_map, train_genes, args.batch_size, device
            )
            
            engine.optimizer.zero_grad()
            loss = engine.compute_loss(x0, x1, y_vec)
            loss.backward()
            engine.optimizer.step()
            train_loss = loss.item()
            
            # --- VALIDATION STEP (Every 50 epochs) ---
            val_loss = 0.0
            if epoch % 50 == 0:
                with torch.no_grad():
                    # Sample batch from VALIDATION genes (Unseen during training)
                    x0_val, x1_val, y_vec_val = get_batch(
                        x_control_all, pert_dict_all, gene_map, val_genes, args.batch_size, device
                    )
                    v_loss = engine.compute_loss(x0_val, x1_val, y_vec_val)
                    val_loss = v_loss.item()
                
                # Checkpointing based on validation
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    torch.save(model.state_dict(), f"{args.save_dir}/best_model_fold_{fold}.pt")
                    
                history['val_loss'].append(val_loss)
            
            history['train_loss'].append(train_loss)
            pbar.set_postfix({'T_Loss': f"{train_loss:.4f}", 'V_Loss': f"{val_loss:.4f}"})

        # Save Final Model config for this fold
        config = {
            'input_dim': input_dim,
            'cond_vec_dim': cond_vec_dim,
            'train_genes': train_genes,
            'val_genes': val_genes
        }
        with open(f"{args.save_dir}/config_fold_{fold}.json", 'w') as f:
            json.dump(config, f)
            
    print("\nTraining Complete. Models saved.")

if __name__ == "__main__":
    main()