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
from FlowMatchingNet import FlowMatchingNet
from FlowMatchingWrapper import FlowMatchingWrapper

# Pero si las definiste en el mismo script o arriba, necesitas:
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
    key="X",        # use hvg
    p_std=1.0,
)

dataset = FlowMatchingWrapper(base_dataset)
loader = torch.utils.data.DataLoader(dataset, batch_size=1024, shuffle=True)

model = FlowMatchingNet(dim=2000)
optim = torch.optim.Adam(model.parameters(), lr=1e-4)

for epoch in range(2):
    for xt, t, vt in loader:
        # xt, t, vt = xt.cuda(), t.cuda(), vt.cuda()

        v_pred = model(xt, t)
        loss = torch.nn.functional.mse_loss(v_pred, vt)

        optim.zero_grad()
        loss.backward()
        optim.step()

    print(f"epoch {epoch}: loss = {loss.item():.4f}")

