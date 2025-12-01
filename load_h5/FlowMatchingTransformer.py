import torch
import torch.nn as nn

class FlowMatchingTransformer(nn.Module):
    """
    Simple Transformer encoder for flow matching.
    Treats gene expression vector (dim,) as a sequence of tokens.
    """
    def __init__(self, dim: int, hidden_dim: int = 256, n_heads: int = 4, n_layers: int = 2):
        super().__init__()

        self.dim = dim
        self.token_proj = nn.Linear(1, hidden_dim)   # map scalar gene to embedding

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_dim,
            nhead=n_heads,
            dim_feedforward=hidden_dim * 4,
            batch_first=True
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        # Time embedding
        self.t_proj = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )

        # Output head → predict velocity
        self.out = nn.Sequential(
            nn.Linear(hidden_dim + hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, dim)
        )

    def forward(self, x, t):
        """
        x: (B, dim)
        t: (B, 1)
        """
        B, D = x.shape
        
        # convert to sequence (B, D, 1)
        seq = x.unsqueeze(-1)

        # (B, D, hidden_dim)
        seq = self.token_proj(seq)

        # Transformer encoding
        h = self.encoder(seq)     # (B, D, hidden_dim)

        # Pool (mean)
        h = h.mean(dim=1)         # (B, hidden_dim)

        # Time embedding
        t_emb = self.t_proj(t.reshape(-1, 1))

        # concat features + t
        z = torch.cat([h, t_emb], dim=1)

        return self.out(z)
