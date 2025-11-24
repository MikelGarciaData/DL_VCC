# flow_matching_pytorch.py
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
from typing import Optional

# -------------------------
# Model: MLP takes (x, t) -> speed
# -------------------------
class FlowMatchingNet(nn.Module):
    def __init__(self, dim: int, hidden_dim: int = 512, n_hidden: int = 2):
        """
        dim: dimension of embeddings (D)
        hidden_dim: hidden layer's size
        n_hidden: nr hidden layers (other than input/output layer)
        """
        super().__init__()
        # priject t to a dimension of y and concatenate x
        self.t_proj = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU()
        )
        layers = []
        in_dim = dim + hidden_dim
        for _ in range(n_hidden):
            layers.append(nn.Linear(in_dim, hidden_dim))
            layers.append(nn.SiLU())
            in_dim = hidden_dim
        layers.append(nn.Linear(hidden_dim, dim))  # output = speed (dim)
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor, t: torch.Tensor):
        """
        x: (B, D)
        t: (B, 1) values in [0,1]
        """
        t_emb = self.t_proj(t)          # (B, hidden_dim)
        inp = torch.cat([x, t_emb], dim=1)
        v = self.net(inp)               # (B, D)
        return v
