import torch
import torch.nn as nn

class StateEmbeddingAE(nn.Module):
    def __init__(self, n_genes, embed_dim=64):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_genes, 512),
            nn.ReLU(),
            nn.Linear(512, embed_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(embed_dim, 512),
            nn.ReLU(),
            nn.Linear(512, n_genes),
            nn.Softplus()  # for only positive vals
        )

    def forward(self, x):
        z = self.encoder(x)
        recon = self.decoder(z)
        return z, recon
