import os
from pathlib import Path
import h5py
import scanpy as sc
import polars as pl
import numpy as np
from scipy.sparse import issparse

import torch
from torch.utils.data import DataLoader, Dataset

from PairDataset import AnnDataPairDataset
from FlowMatchingTransformer import FlowMatchingTransformer  
from FlowMatchingWrapper import FlowMatchingWrapper

from anndata import AnnData
import scipy.sparse as sp

current_path = Path(__file__).parent
data_folder = (current_path.parent / 'data' / 'adata_hvg2000.h5ad').resolve()

print(f"loading data from {data_folder}")
adata_hvg = sc.read_h5ad(data_folder)

X = adata_hvg.X.toarray() if issparse(adata_hvg.X) else adata_hvg.X
X = np.log1p(X)
X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-6)
adata_hvg.X = X

base_dataset = AnnDataPairDataset(
    adata_hvg,
    key="X",
    p_std=1.0,
)

dataset = FlowMatchingWrapper(base_dataset)
loader = torch.utils.data.DataLoader(dataset, batch_size=16, shuffle=True, pin_memory=True) # pin_memory only if accelerator

model = FlowMatchingTransformer(
    dim=2000,
    hidden_dim=256,
    n_heads=4,
    n_layers=2
)


optim = torch.optim.Adam(model.parameters(), lr=1e-4)

loss_history = []

for epoch in range(100):
    for xt, t, vt in loader:

        v_pred = model(xt, t)
        loss = torch.nn.functional.mse_loss(v_pred, vt)

        optim.zero_grad()
        loss.backward()
        optim.step()

    # save loss
    loss_history.append(loss.item())

    print(f"epoch {epoch}: loss = {loss.item():.4f}")

    # === checkpoint every 25 epochs ===
    if (epoch + 1) % 25 == 0:

        ckpt_path = f"transformerFM_epoch_{epoch+1}.pt"
        torch.save({
            "epoch": epoch + 1,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optim.state_dict(),
            "loss": loss.item(),
        }, ckpt_path)

        # === pickle with loss history ===
        pkl_path = "loss_history.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(
                {
                    "epochs": list(range(1, len(loss_history) + 1)),
                    "loss": loss_history
                },
                f
            )

        print(f"[checkpoint] model saved as {ckpt_path}")
        print(f"[checkpoint] loss saved as {pkl_path}")

# ==== plot loss ====
plt.figure(figsize=(8, 5))
plt.plot(range(1, len(loss_history) + 1), loss_history)
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Flow Matching Training Loss")
plt.grid(True)
plt.show()