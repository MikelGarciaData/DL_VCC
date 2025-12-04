import scvi
import scanpy as sc
import torch
import os

def run_scvi_embedding(file_path, output_path, n_latent=100):
    """
    Loads AnnData, trains scVI with 100 latent dimensions, and saves the result.
    """
    
    # 1. Load your AnnData
    print(f"Loading data from {file_path}...")
    adata = sc.read_h5ad(file_path)

    # Basic preprocessing if not already done (scVI expects raw counts usually)
    # Ideally, adata.layers['counts'] should contain raw integer counts
    if 'counts' not in adata.layers:
        print("Warning: 'counts' layer not found. Copying .X to layers['counts'].")
        # Ensure we preserve raw counts
        adata.layers['counts'] = adata.X.copy()

    # 2. Setup AnnData for scVI
    # batch_key: Change this to the column name in adata.obs that corresponds to your batch/sample
    # layer: scVI is best trained on raw counts
    print("Setting up AnnData for scVI...")
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", 
        batch_key=None  # Set this to "sample" or "batch" if you have batch effects!
    )

    # 3. Initialize the model
    # n_latent=100: This sets the embedding dimension to 100 as requested
    print(f"Initializing scVI model with n_latent={n_latent}...")
    model = scvi.model.SCVI(adata, n_latent=n_latent)

    # 4. Train the model
    # scVI automatically detects GPU. explicit strict check can be added via accelerator="gpu"
    print("Training model...")
    model.train(
        max_epochs=None,    # Default heuristic is usually good
        accelerator="gpu",  # Force GPU usage
        devices=1           # Use 1 GPU
    )

    # 5. Extract the 100-dimension embedding
    print("Extracting latent representation...")
    latent_representation = model.get_latent_representation()

    # Store it in the AnnData object
    adata.obsm["X_scvi"] = latent_representation

    # 6. Save the results
    print(f"Saving results to {output_path}...")
    adata.write_h5ad(output_path)
    
    # Optional: Save just the model for later use
    # model.save("my_scvi_model/")
    
    print("Done! Embedding stored in adata.obsm['X_scVI']")

if __name__ == "__main__":
    # Check for GPU availability
    if torch.cuda.is_available():
        print(f"GPU detected: {torch.cuda.get_device_name(0)}")
    else:
        print("WARNING: No GPU detected. Training will be slow.")

    # Replace these with your actual file paths
    INPUT_FILE = "../../adata_subb.h5ad"
    folder = "data_sub"    
    os.makedirs(folder, exist_ok=True)
    OUTPUT_FILE = os.path.join(folder, "adata_subb_scvi.h5ad")
    
    run_scvi_embedding(INPUT_FILE, OUTPUT_FILE)