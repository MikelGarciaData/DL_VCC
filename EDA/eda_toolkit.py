import anndata as ad
import pandas as pd
import numpy as np

def create_subset_h5ad(
    input_path: str,
    output_path: str,
    cell_filter: pd.Series = None,
    gene_filter: pd.Series = None
):
    """
    Create a subset of an AnnData (.h5ad) file and save it to disk.

    Parameters
    ----------
    input_path : str
        Path to the input .h5ad file containing the full dataset.
    output_path : str
        Path to save the subset .h5ad file.
    cell_filter : pandas.Series, optional
        Boolean mask (length = n_obs) to filter cells.
        Example: adata.obs['target_gene'] == 'TP53'
    gene_filter : pandas.Series, optional
        Boolean mask (length = n_vars) to filter genes.
        Example: adata.var_names.isin(['TP53', 'GAPDH'])
    
    Returns
    -------
    ad.AnnData
        The subset AnnData object (also written to `output_path`).
    
    Example
    -------
    >>> subset = create_subset_h5ad(
    ...     input_path="adata_Training.h5ad",
    ...     output_path="adata_TP53_subset.h5ad",
    ...     cell_filter=adata.obs['target_gene'] == 'TP53'
    ... )
    """

    # Load the full AnnData object
    print(f"Loading AnnData from {input_path}...")
    adata = ad.read_h5ad(input_path)

    # Apply filters if provided
    if cell_filter is not None:
        print(f"Filtering cells: {cell_filter.sum()} cells retained.")
        adata = adata[cell_filter, :].copy()

    if gene_filter is not None:
        print(f"Filtering genes: {gene_filter.sum()} genes retained.")
        adata = adata[:, gene_filter].copy()

    # Save the subset
    adata.write_h5ad(output_path)
    print(f"Subset saved to {output_path}.")

    return adata


def create_subset_h5ad_by_size(
    input_path: str,
    output_path: str,
    n_cells: int = None,
    n_genes: int = None,
    random_state: int = 42
):
    """
    Create a subset of an AnnData (.h5ad) file by specifying 
    the number of cells (rows) and/or genes (columns).

    Parameters
    ----------
    input_path : str
        Path to the input .h5ad file containing the full dataset.
    output_path : str
        Path to save the subset .h5ad file.
    n_cells : int, optional
        Number of cells (rows) to include in the subset. 
        If None, all cells are kept.
    n_genes : int, optional
        Number of genes (columns) to include in the subset. 
        If None, all genes are kept.
    random_state : int, default=42
        Random seed for reproducibility when sampling.

    Returns
    -------
    ad.AnnData
        The subset AnnData object (also written to `output_path`).

    Example
    -------
    >>> subset = create_subset_h5ad_by_size(
    ...     input_path="adata_Training.h5ad",
    ...     output_path="adata_subset_10k_5k.h5ad",
    ...     n_cells=10000,
    ...     n_genes=5000
    ... )
    """

    print(f"Loading AnnData from {input_path}...")
    adata = ad.read_h5ad(input_path)

    # Randomly select a subset of cells
    if n_cells is not None:
        n_cells = min(n_cells, adata.n_obs)
        np.random.seed(random_state)
        cell_indices = np.random.choice(adata.n_obs, n_cells, replace=False)
        adata = adata[cell_indices, :].copy()
        print(f"Selected {n_cells} cells.")

    # Randomly select a subset of genes
    if n_genes is not None:
        n_genes = min(n_genes, adata.n_vars)
        np.random.seed(random_state + 1)
        gene_indices = np.random.choice(adata.n_vars, n_genes, replace=False)
        adata = adata[:, gene_indices].copy()
        print(f"Selected {n_genes} genes.")

    # Save the subset to disk
    adata.write_h5ad(output_path)
    print(f"Subset saved to {output_path}.")

    return adata