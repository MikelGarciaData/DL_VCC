import numpy as np
import torch
from torch.utils.data import Dataset
from anndata import AnnData
import scipy.sparse as sp


class AnnDataPairDataset(Dataset):
    """
    Dataset for Flow Matching that uses AnnData as obejctive distribution q(x),
    and a isotropic normal distribution as prior distribution p(x).
    """

    def __init__(
        self,
        adata: AnnData,
        key: str = "X",            # Use next matrices: "X", "obsm:SE", etc.
        p_std: float = 1.0,
        rng: np.random.RandomState = None,
    ):
        """
        Parameters
        ----------
        adata : AnnData
            Data loaded with scanpy
        key : str
            "X" → adata.X
            "obsm:SE" → adata.obsm["SE"]
            "layers:raw" → adata.layers["raw"]
        p_std : float
            std of prior distribution p = N(0, p_std^2 I)
        """
        self.p_std = float(p_std)
        self.rng = rng or np.random.RandomState(0)

        # --- select matrix ---
        if key == "X":
            X = adata.X
        elif key.startswith("obsm:"):
            name = key.split(":")[1]
            X = adata.obsm[name]
        elif key.startswith("layers:"):
            name = key.split(":")[1]
            X = adata.layers[name]
        else:
            raise ValueError(f"key '{key}' no reconocido")

        # --- sparse -> dense ---
        if sp.issparse(X):
            X = X.toarray()

        # cast to float32 numpy
        self.X_q = np.asarray(X, dtype=np.float32)
        self.N, self.D = self.X_q.shape

    def __len__(self):
        return self.N

    def __getitem__(self, idx):
        # sample x1 ~ q (empirical data)
        i = self.rng.randint(0, self.N)
        x1 = self.X_q[i]

        # sample x0 ~ p (isotropic gaussian)
        x0 = self.rng.normal(loc=0.0, scale=self.p_std, size=self.D).astype(np.float32)

        return x0, x1

