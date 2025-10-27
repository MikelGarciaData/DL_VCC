import os
from pathlib import Path
import h5py
import scanpy as sc
import polars as pl
import numpy as np
from scipy.sparse import issparse

current_path = Path(__file__).parent
data_folder = (current_path.parent / 'data' / 'adata_Training.h5ad').resolve()
# data_folder = os.path.abspath(data_folder)

print(f"loading data from {data_folder}")
# Cargar el archivo .h5ad
adata = sc.read_h5ad(data_folder)

# Ejemplo: mostrar resumen
print(adata)

# print(adata.X)     # expression matrix (genes x cells)
# print(adata.obs)      # cell metadata (DataFrame)
# print(adata.var)       # gene metadata (DataFrame)
# print(adata.obsm)       # embeddings (e.g., PCA, UMAP, etc.)

for key in adata.__dict__.keys():
    print(f"- {key}")

# Sparse matrix in Scipy
print(type(adata._X))
print(adata._X.shape)

# Dataframe
print(type(adata.obs))
print(type(adata.var))

# Observation/variable-level matrices 
# Aligned mapping object
print(type(adata.obsm)) # UMAP embedding
# .obsm matrices must length equal to the number of observations as .n_obs


# Extraer nombres de genes y células
print("Extract gene and cell names...")

genes = adata.var_names.to_list()
cells = adata.obs_names.to_list()

# Obtener la matriz de expresión
X = adata.X

print("Creating dataframe of cell and genes...")

if issparse(X):
    X = X.tocoo()  # formato coordenado: más fácil para reconstruir
    df = pl.DataFrame({
        "CELL": [cells[i] for i in X.row],
        "GEN": [genes[j] for j in X.col],
        "count": X.data
    })
else:
    # Si no es sparse (denso)
    df = pl.DataFrame({
        "GEN": np.repeat(genes, len(cells)),
        "CELL": np.tile(cells, len(genes)),
        "count": X.T.flatten()  # transpuesta para coincidir con orden (gene, cell)
    })

print(df.head())
print(df.shape)
