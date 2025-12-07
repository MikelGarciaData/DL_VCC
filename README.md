# GeneFlow: Conditional Flow Matching for Single-Cell Perturbation

## Jupyter Notebook to run and understand the model: [geneFlow.ipynb](https://github.com/MikelGarciaData/DL_VCC/blob/main/geneFlow.ipynb)

This repository **GeneFlow**, a generative model designed to predict the effects of genetic perturbations on single cells.

The working and theory behind the model can be found here: [GeneFlow Paper](https://github.com/MikelGarciaData/DL_VCC/blob/main/Deep_Learning_project.pdf)

GeneFlow aims for **Zero-Shot Generalization**: predicting the transcriptomic phenotype of a cell after perturbing a gene the model has *never seen before*, provided it has access to a vector representation of that gene (e.g., from GenePT).

<img width="427" height="315" alt="Screenshot 2025-12-07 at 22 15 14" src="https://github.com/user-attachments/assets/f7ab6f60-ebb8-4644-be66-075011d8fe5c" />

The white are the control cells and C1 C2 are the perturbed states. We want to predict the perturbed states.
## 1. The Algorithm

GeneFlow uses the ideas of **Flow Matching**. It learns a continuous **Velocity Field** that pushes cells directly from a source distribution (Control) to a target distribution (Perturbed) along a straight-line path.

### Training Procedure

The model minimizes the difference between the predicted velocity and the true direction between control and perturbed states.

<img width="568" height="634" alt="Screenshot 2025-12-07 at 21 54 54" src="https://github.com/user-attachments/assets/76295b9b-6ca6-4d9e-870b-199da7f1418f" />

## 2. Model Architecture

The core model (`VectorFlowNet`) is a **Conditional ResNet** that conditions the entire network on the specific gene perturbation being applied.

### Key Components

1.  **`SinusoidalPosEmb`**: Encodes the continuous time scalar $t \in [0, 1]$ into a high-dimensional vector. This tells the network "where" in the interpolation process the cell currently is.
2.  **`AdaLN` (Adaptive Layer Normalization)**: This is the critical conditioning mechanism. Instead of simply concatenating the gene vector, `AdaLN` uses the gene embedding $y$ to predict the `scale` (gamma) and `shift` (beta) parameters of the normalization layers.
    * *Effect:* This effectively "reprograms" the network's behavior based on *which* gene is being perturbed, allowing generalization to unseen gene vectors.
3.  **`ResBlock`**: Standard residual blocks (MLP + Nonlinearity + Dropout) that process the cell state while being modulated by the `AdaLN` layers.

## 3. Training Strategy (Zero-Shot)

To ensure the model isn't just memorizing data, we use a **Leave-Gene-Out** strategy:

1.  **Gene Splitting:** The script identifies all unique perturbation genes. It randomly sets aside a percentage (e.g., 20%) as a **Validation Set**.
2.  **Separation:** During training, the model **never sees** cells perturbed by the validation genes.
3.  **Checkpointing:** The model weights (`best_model.pt`) are saved only when the loss on these *unseen* genes improves.

## 4. How to Run

### Step 1: Prepare Data
You need two files:
1.  **`data.h5ad`**: A Scanpy object containing your single-cell data.
    * `.obsm['X_scvi']`: Latent embeddings of the cells (or PCA).
    * `.obs['target_gene']`: Column indicating the perturbed gene (use 'non-targeting' for controls).
2.  **`gene_vectors.txt`**: A text file mapping gene names to vectors.
    * Format: `GENE_NAME 0.123 0.456 ...`

### Step 2: Execution (Training & Prediction)

`train_m_flow.py` **Arguments**



| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--adata_path` | str | Required | Path to the input `.h5ad` file containing the single-cell data. |
| `--gene_vec_path` | str | Required | Path to the file containing gene embeddings. |
| `--save_dir` | str | `./checkpoints` | Directory where model checkpoints, configuration, loss history, and predictions will be saved. |
| `--emb_key` | str | `X_scvi` | Key in `adata.obsm` that contains the embeddings to be used as input. |
| `--gene_col` | str | `target_gene` | Column in `adata.obs` specifying the gene associated with each cell. |
| `--control_label` | str | `non-targeting` | Label identifying control cells in `adata.obs[gene_col]`. |
| `--epochs` | int | `50` | Number of training epochs. |
| `--batch_size` | int | `128` | Number of samples per training batch. |
| `--val_split` | float | `0.2` | Fraction of perturbed genes to use for validation. |
| `--hidden_dim` | int | `128` | Hidden dimension size for the neural network. |
| `--num_layers` | int | `6` | Number of residual blocks in the VectorFlowNet model. |
| `--seed` | int | `42` | Random seed for reproducibility. |
| `--pred_steps` | int | `20` | Number of ODE solver steps during trajectory prediction (should be divisible by 4). |

---

```bash
python main.py \
  --adata_path "./data/perturbation_data.h5ad" \
  --gene_vec_path "./data/gene_vectors.txt" \
  --save_dir "./checkpoints_geneflow" \
  --emb_key "X_scvi" \
  --control_label "non-targeting" \
  --gene_col "target_gene" \
  --epochs 100 \
  --batch_size 128 \
  --hidden_dim 256 \
  --num_layers 6 \
  --pred_steps 20
```


