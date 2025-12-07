import scvi
import scanpy as sc
import torch
import os

def run_scvi_embedding(file_path, output_path, n_latent=50):
    """
    Loads AnnData, trains scVI with 50 latent dimensions, and saves the result.
    """
    
    #Load your AnnData
    print(f"Loading data from {file_path}...")
    adata = sc.read_h5ad(file_path)

    # Preprocessing
    if 'counts' not in adata.layers:
        print("Warning: 'counts' layer not found. Copying .X to layers['counts'].")
        adata.layers['counts'] = adata.X.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata 

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=9500,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="batch",
    )

    # Setup AnnData for scVI
    print("Setting up AnnData for scVI...")
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", 
        batch_key="batch"
    )

    #Initialize the model
    print(f"Initializing scVI model with n_latent={n_latent}...")
    model = scvi.model.SCVI(adata, n_latent=n_latent)

    # Train the model
    print("Training model...")
    model.train(
        max_epochs=None,  
        accelerator="gpu", 
        devices=1           
    )

    print("Extracting latent representation...")
    latent_representation = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent_representation
    print(f"Saving results to {output_path}...")
    adata.write_h5ad(output_path)
    
    #SAVE THE MODEL WEIGHTS
    model_dir = os.path.join(os.path.dirname(output_path), "scvi_model_dir_50")
    
    print(f"Saving model weights to {model_dir}...")
    model.save(model_dir, overwrite=True)
    print("Done! Model and data saved.")
    print("Done! Embedding stored in adata.obsm['X_scvi']")

if __name__ == "__main__":
    if torch.cuda.is_available():
        print(f"GPU detected: {torch.cuda.get_device_name(0)}")
    else:
        print("WARNING: No GPU detected.")

    INPUT_FILE = "/dtu/blackhole/1a/187738/vcc_data/adata_Training.h5ad"
    folder = "data_sub"    
    os.makedirs(folder, exist_ok=True)
    OUTPUT_FILE = os.path.join(folder, "adata_all_scvi_50.h5ad")
    
    run_scvi_embedding(INPUT_FILE, OUTPUT_FILE)