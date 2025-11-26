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
        super().__init__()

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

        layers.append(nn.Linear(hidden_dim, dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor, t: torch.Tensor):
        # Accept shapes (B,), (B,1), (B,1,1)
        t = t.reshape(-1, 1)

        t_emb = self.t_proj(t)     # (B, hidden_dim)

        inp = torch.cat([x, t_emb], dim=1)  # (B, D + hidden_dim)
        v = self.net(inp)
        return v
