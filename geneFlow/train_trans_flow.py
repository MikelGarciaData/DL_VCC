# save as: train_transformer.py
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
# 1. MODEL ARCHITECTURE (DiT)
# ==========================================

class TimestepEmbedder(nn.Module):
    def __init__(self, hidden_dim, frequency_embedding_size=256):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(frequency_embedding_size, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.frequency_embedding_size = frequency_embedding_size

    @staticmethod
    def timestep_embedding(t, dim, max_period=10000):
        half = dim // 2
        freqs = torch.exp(
            -math.log(max_period) * torch.arange(start=0, end=half, dtype=torch.float32) / half
        ).to(device=t.device)
        args = t[:, None].float() * freqs[None]
        embedding = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
        if dim % 2:
            embedding = torch.cat([embedding, torch.zeros_like(embedding[:, :1])], dim=-1)
        return embedding

    def forward(self, t):
        t_freq = self.timestep_embedding(t, self.frequency_embedding_size)
        return self.mlp(t_freq)

class DiTBlock(nn.Module):
    def __init__(self, hidden_dim, num_heads, mlp_ratio=4.0):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_dim, elementwise_affine=False, eps=1e-6)
        self.attn = nn.MultiheadAttention(hidden_dim, num_heads=num_heads, batch_first=True)
        self.norm2 = nn.LayerNorm(hidden_dim, elementwise_affine=False, eps=1e-6)
        
        mlp_hidden_dim = int(hidden_dim * mlp_ratio)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_dim, mlp_hidden_dim),
            nn.GELU(),
            nn.Linear(mlp_hidden_dim, hidden_dim)
        )
        
        # adaLN modulation: regresses scale/shift parameters from the condition
        self.adaLN_modulation = nn.Sequential(
            nn.SiLU(),
            nn.Linear(hidden_dim, 6 * hidden_dim, bias=True)
        )

    def forward(self, x, c):
        shift_msa, scale_msa, gate_msa, shift_mlp, scale_mlp, gate_mlp = self.adaLN_modulation(c).chunk(6, dim=1)
        
        x_norm1 = self.norm1(x) * (1 + scale_msa.unsqueeze(1)) + shift_msa.unsqueeze(1)
        attn_out, _ = self.attn(x_norm1, x_norm1, x_norm1)
        x = x + gate_msa.unsqueeze(1) * attn_out
        
        x_norm2 = self.norm2(x) * (1 + scale_mlp.unsqueeze(1)) + shift_mlp.unsqueeze(1)
        mlp_out = self.mlp(x_norm2)
        x = x + gate_mlp.unsqueeze(1) * mlp_out
        return x

class TransformerFlowNet(nn.Module):
    def __init__(self, input_dim, cond_vec_dim, hidden_dim=256, num_layers=4, num_heads=8, num_tokens=8):
        super().__init__()
        self.num_tokens = num_tokens
        self.hidden_dim = hidden_dim
        
        # Project Cell Vector -> Sequence of Tokens
        self.x_embedder = nn.Linear(input_dim, hidden_dim * num_tokens)
        
        # Embed Condition (Time + Gene)
        self.t_embedder = TimestepEmbedder(hidden_dim)
        self.y_embedder = nn.Linear(cond_vec_dim, hidden_dim)
        
        # Transformer Backbone
        self.blocks = nn.ModuleList([
            DiTBlock(hidden_dim, num_heads) for _ in range(num_layers)
        ])
        
        # Final Decoder
        self.final_norm = nn.LayerNorm(hidden_dim, elementwise_affine=False, eps=1e-6)
        self.final_adaLN = nn.Sequential(
            nn.SiLU(),
            nn.Linear(hidden_dim, 2 * hidden_dim, bias=True)
        )
        self.output_proj = nn.Linear(hidden_dim * num_tokens, input_dim)
        self.initialize_weights()

    def initialize_weights(self):
        # Zero-init the last AdaLN layers so the model starts as Identity function
        for block in self.blocks:
            nn.init.constant_(block.adaLN_modulation[-1].weight, 0)
            nn.init.constant_(block.adaLN_modulation[-1].bias, 0)
        nn.init.constant_(self.final_adaLN[-1].weight, 0)
        nn.init.constant_(self.final_adaLN[-1].bias, 0)
        nn.init.constant_(self.output_proj.weight, 0)
        nn.init.constant_(self.output_proj.bias, 0)

    def forward(self, x, t, y_vec):
        # x: (Batch, Input_Dim) -> (Batch, Num_Tokens, Hidden_Dim)
        x = self.x_embedder(x).view(-1, self.num_tokens, self.hidden_dim)
        
        t_emb = self.t_embedder(t.squeeze(-1))
        y_emb = self.y_embedder(y_vec)
        c = t_emb + y_emb # Condition
        
        for block in self.blocks:
            x = block(x, c)
            
        shift, scale = self.final_adaLN(c).chunk(2, dim=1)
        x = self.final_norm(x) * (1 + scale.unsqueeze(1)) + shift.unsqueeze(1)
        
        x = x.flatten(1)
        return self.output_proj(x)

# ==========================================
# 2. DATA UTILS
# ==========================================

def load_gene_vectors(file_path):
    print(f"Loading gene vectors from {file_path}...")
    gene_map = {}
    detected_dim = None

    with open(file_path, 'r') as f:
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
    return gene_map, detected_dim

class SingleCellPairDataset(Dataset):
    def __init__(self, x_control, pert_dict, gene_map, active_genes):
        self.x_control = x_control
        self.pert_cells = []
        self.pert_conds = []
        
        print(f"Building dataset for {len(active_genes)} genes...")
        for gene in active_genes:
            if gene not in pert_dict: continue
            cells = pert_dict[gene]
            vec = gene_map[gene]
            
            # Add cells
            self.pert_cells.append(cells)
            # Add gene condition (repeated for each cell)
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
        x1 = self.pert_cells[idx]         # Target (Perturbed)
        y_vec = self.pert_conds[idx]      # Condition (Gene)
        
        # Random Control Cell coupling
        rand_idx = torch.randint(0, len(self.x_control), (1,)).item()
        x0 = self.x_control[rand_idx]
        return x0, x1, y_vec

# ==========================================
# 3. MAIN TRAINING LOOP
# ==========================================

def main():
    # HARDCODED PATHS (You can change these or use argparse)
    ADATA_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/helpers/data_sub/adata_all_scvi.h5ad"
    GENE_VEC_PATH = "/dtu/blackhole/1a/187738/deep/DL_VCC/geneFlow/149_genePT.txt" # Created by prepare_data.py
    SAVE_DIR = "./checkpoints_transformer"
    
    # PARAMETERS
    BATCH_SIZE = 64  # Lower if out of memory
    EPOCHS = 50
    LR = 1e-4
    EMB_KEY = "X_scvi"
    GENE_COL = "target_gene"
    CONTROL_LABEL = "non-targeting"
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')
    print(f"Using device: {device}")
    os.makedirs(SAVE_DIR, exist_ok=True)

    # 1. Load Data
    adata = sc.read_h5ad(ADATA_PATH)
    X_emb = adata.obsm[EMB_KEY]
    target_genes = adata.obs[GENE_COL].values.astype(str)
    input_dim = X_emb.shape[1]
    
    # 2. Load Gene Vectors
    gene_map, cond_vec_dim = load_gene_vectors(GENE_VEC_PATH)
    print(f"Cell Dim: {input_dim}, Gene Vector Dim: {cond_vec_dim}")

    # 3. Prepare Tensors
    ctrl_mask = (target_genes == CONTROL_LABEL)
    x_control_all = torch.tensor(X_emb[ctrl_mask], dtype=torch.float32)
    
    valid_pert_genes = list(gene_map.keys())
    pert_dict_all = {}
    for g in valid_pert_genes:
        # Case insensitive match for safety
        mask = np.char.upper(target_genes) == g
        if np.sum(mask) > 0:
            pert_dict_all[g] = torch.tensor(X_emb[mask], dtype=torch.float32)

    # 4. Split Genes (Zero-Shot setting)
    train_genes, val_genes = train_test_split(valid_pert_genes, test_size=0.1, random_state=42)
    
    train_ds = SingleCellPairDataset(x_control_all, pert_dict_all, gene_map, train_genes)
    val_ds = SingleCellPairDataset(x_control_all, pert_dict_all, gene_map, val_genes)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)

    # 5. Initialize Model
    model = TransformerFlowNet(
        input_dim=input_dim,
        cond_vec_dim=cond_vec_dim,
        hidden_dim=256,
        num_layers=4,
        num_heads=8,
        num_tokens=8
    ).to(device)
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)
    criterion = nn.MSELoss()

    print("\nStarting Training...")
    best_loss = float('inf')
    
    for epoch in range(EPOCHS):
        model.train()
        train_loss = 0
        pbar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{EPOCHS}")
        
        for x0, x1, y_vec in pbar:
            x0, x1, y_vec = x0.to(device), x1.to(device), y_vec.to(device)
            
            # Flow Matching Loss Calculation
            t = torch.rand(x0.shape[0], 1, device=device)
            x_t = (1 - t) * x0 + t * x1    # Path
            u_target = x1 - x0             # Target Velocity
            
            u_pred = model(x_t, t, y_vec)  # Prediction
            
            loss = criterion(u_pred, u_target)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
            pbar.set_postfix({'Loss': f"{loss.item():.4f}"})
            
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for x0, x1, y_vec in val_loader:
                x0, x1, y_vec = x0.to(device), x1.to(device), y_vec.to(device)
                t = torch.rand(x0.shape[0], 1, device=device)
                x_t = (1 - t) * x0 + t * x1
                u_target = x1 - x0
                u_pred = model(x_t, t, y_vec)
                val_loss += criterion(u_pred, u_target).item()
        
        avg_val_loss = val_loss / len(val_loader)
        print(f"Val Loss: {avg_val_loss:.5f}")
        
        if avg_val_loss < best_loss:
            best_loss = avg_val_loss
            torch.save(model.state_dict(), f"{SAVE_DIR}/best_model.pt")

    print("Training Complete.")

if __name__ == "__main__":
    main()