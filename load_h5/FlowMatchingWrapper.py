import torch
from torch.utils.data import Dataset

class FlowMatchingWrapper(Dataset):
    def __init__(self, base_dataset):
        self.base = base_dataset

    def __len__(self):
        return len(self.base)

    def __getitem__(self, idx):
        x0, x1 = self.base[idx]

        # Convertir a tensores
        x0 = torch.from_numpy(x0)
        x1 = torch.from_numpy(x1)

        # samplear t uniformemente
        t = torch.rand(1)

        # generar punto interpolado
        xt = (1 - t) * x0 + t * x1

        # velocidad real
        vt = x1 - x0

        return xt, t.unsqueeze(0), vt
