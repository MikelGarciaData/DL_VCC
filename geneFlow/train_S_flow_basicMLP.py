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
# 1. Neural Network Architecture (MLP Version)
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


class MLPVectorField(nn.Module):
    """
    Minimal MLP alternative to the previous VectorTransformer.

    Learns u(x_t, t, gene_vec) = MLP([x_t, cond]),
    where cond = time_mlp(t) + gene_vec_proj(gene_vec).
    """
    def __init__(self, input_dim, cond_vec_dim, hidden_dim=256, num_layers=4):
        super().__init__()

        # Time and gene conditioning -> hidden_dim
        self.time_mlp = nn.Sequential(
            SinusoidalPosEmb(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.gene_vec_proj = nn.Linear(cond_vec_dim, hidden_dim)

        # MLP over [x_t, cond]
        mlp_layers = []
        in_dim = input_dim + hidden_dim  # concat x_t and cond

        # num_layers here controls depth of the MLP
        # we'll use (num_layers - 1) hidden layers + final output layer
        for _ in range(max(num_layers - 1, 1)):
            mlp_layers.append(nn.Linear(in_dim, hidden_dim))
            mlp_layers.append(nn.GELU())
            in_dim = hidden_dim

        # Final layer back to input_dim (velocity in embedding space)
        mlp_layers.append(nn.Linear(in_dim, input_dim))

        self.mlp = nn.Sequential(*mlp_layers)

    def forward(self, x, t, y_vec):
        """
        x:      (batch, input_dim)        - current state x_t
        t:      (batch, 1)                - scalar time in [0,1]
        y_vec:  (batch, cond_vec_dim)     - gene vector

        returns:
            u_pred: (batch, input_dim)    - predicted velocity
        """
        t_emb = self.time_mlp(t.squeeze(-1))        # (batch, hidden_dim)
        y_emb = self.gene_vec_proj(y_vec)           # (batch, hidden_dim)
        cond = t_emb + y_emb                        # (batch, hidden_dim)

        z = torch.cat([x, cond], dim=-1)            # (batch, input_dim + hidden_dim)
        u_pred = self.mlp(z)                        # (batch, input_dim)
        return u_pred

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
            # It was a header like "N D"
            pass
        else:
            # No header: rewind
            f.seek(0)

        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue

            gene = parts[0].upper()
            try:
                vec = np.array(parts[1:], dtype=np.float32)
                if detected_dim is None:
                    detected_dim = len(vec)
                    print(f"Detected gene vector dimension: {detected_dim}")

                if len(vec) == detected_dim:
                    gene_map[gene] = torch.tensor(vec)
            except ValueError:
                # skip malformed lines
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
# 4. Main Execution
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Train OT-CFM (Simple Split, MLP version)")
    parser.add_argument("--adata_path", type=str, required=True)
    parser.add_argument("--gene_vec_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, default="./checkpoints")
    parser.add_argument("--emb_key", type=str, default="X_pca")
    parser.add_argument("--gene_col", type=str, default="target_gene")
    parser.add_argument("--control_label", type=str, default="non-targeting")
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--val_split", type=float, default=0.1,
                        help="Fraction of genes to hold out for validation")

    # Architecture
    parser.add_argument("--hidden_dim", type=int, default=256,
                        help="Hidden dimension for conditioning and MLP")
    parser.add_argument("--num_layers", type=int, default=6,
                        help="Number of layers in the MLP (depth)")

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
    train_genes, val_genes = train_test_split(
        valid_pert_genes,
        test_size=args.val_split,
        random_state=42
    )
    print(f"Training on {len(train_genes)} genes. Validating on {len(val_genes)} genes.")

    # 5. Initialize Model (MLP-based vector field)
    model = MLPVectorField(
        input_dim=input_dim,
        cond_vec_dim=cond_vec_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers
    ).to(device)

    engine = FlowMatchingEngine(model, device)

    best_val_loss = float('inf')
    val_loss = 0.0

    # 6. Training Loop
    pbar = tqdm(range(args.epochs), desc="Training")
    for epoch in pbar:

        # Train step
        x0, x1, y_vec = get_batch(
            x_control_all, pert_dict_all, gene_map, train_genes, args.batch_size, device
        )

        engine.optimizer.zero_grad()
        loss = engine.compute_loss(x0, x1, y_vec)
        loss.backward()
        engine.optimizer.step()
        train_loss = loss.item()

        # Validation Step (Every 50 epochs)
        if epoch % 50 == 0 and len(val_genes) > 0:
            with torch.no_grad():
                x0_v, x1_v, y_vec_v = get_batch(
                    x_control_all, pert_dict_all, gene_map, val_genes, args.batch_size, device
                )
                v_loss = engine.compute_loss(x0_v, x1_v, y_vec_v)
                val_loss = v_loss.item()

            # Save Best Model
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(model.state_dict(), os.path.join(args.save_dir, "best_model.pt"))

        pbar.set_postfix({'T_Loss': f"{train_loss:.4f}", 'V_Loss': f"{val_loss:.4f}"})

    # Save Final Config
    config = {
        'input_dim': input_dim,
        'cond_vec_dim': cond_vec_dim,
        'train_genes': train_genes,
        'val_genes': val_genes
    }
    with open(os.path.join(args.save_dir, "config.json"), 'w') as f:
        json.dump(config, f)

    print("\nTraining Complete. Best model saved as 'best_model.pt'")

if __name__ == "__main__":
    main()
