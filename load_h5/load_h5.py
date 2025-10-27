import os
import h5py
import scanpy as sc

datapath = "data/adata_Training.h5ad"

print(f"loading data from {datapath}")
# Cargar el archivo .h5ad
adata = sc.read_h5ad(datapath)

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