# Optimal Transport Conditional Flow Matching (OT-CFM) for Single-Cell Perturbation

This repository contains a PyTorch implementation of **Optimal Transport Conditional Flow Matching (OT-CFM)** designed to predict the effects of genetic perturbations on single cells.

This model is capable of **Zero-Shot Generalization**: it can predict the transcriptomic phenotype of a cell after perturbing a gene the model has *never seen before*, provided it has access to a vector representation of that gene (e.g., from Gene2Vec or an LLM).

## 1. The Algorithm: What is OT-CFM?
<img width="2178" height="828" alt="image" src="https://github.com/user-attachments/assets/387fe811-81cb-4a8f-9e62-059e15e822f1" />


### The Core Concept

Traditional generative models (like GANs or VAEs) try to map noise to data. **Flow Matching** is different: it learns a continuous **Velocity Field** (a vector field) that pushes cells from a source distribution (Control/Healthy) to a target distribution (Perturbed/Disease) over a virtual time interval $t \in [0, 1]$. 

In this code, we use **Optimal Transport (OT)** flow matching. This assumes the most efficient path between a control cell $x_0$ and a perturbed cell $x_1$ is a straight line.

### The "Flow" in the Code

Mathematically, the probability path is defined as a linear interpolation:

$$
x_t = (1 - t)x_0 + t x_1
$$

The **velocity** (how $x_t$ changes with time) is simply the derivative with respect to $t$:

$$
u_t(x_t) = \frac{d}{dt}x_t = x_1 - x_0
$$

### The Objective

The neural network $v_\theta(x_t, t, y)$ tries to predict this velocity. The loss function in `train_flow.py` minimizes the squared difference between the predicted velocity and the true straight-line direction:

$$
\mathcal{L} = || v_\theta(x_t, t, gene\_vector) - (x_1 - x_0) ||^2
$$

## 2. Model Architecture

The core model is a **Conditional Transformer**. It does not treat the cell as a sequence of tokens (like text), but rather uses the Transformer blocks for their powerful conditioning mechanisms.

### Key Components

#### 1. `SinusoidalPosEmb`

Encodes the continuous time scalar $t \in [0, 1]$ into a high-dimensional vector. This allows the network to understand "where" in the flow process the cell currently is (e.g., "we are 50% transformed").

#### 2. `AdaLN` (Adaptive Layer Normalization)

This is the **conditioning mechanism**. Instead of just concatenating the gene vector to the input, we use it to modulate the network layers.

* **Input:** The Gene Vector (200-dim) + Time Embedding.

* **Action:** It predicts `scale` (gamma) and `shift` (beta) parameters.

* **Result:** It effectively "reprograms" the normalization layers of the Transformer based on which gene is being perturbed. This is what allows the model to generalize to new genes. 

#### 3. `VectorTransformer` (The Main Class)

* **Input:** A cell embedding (e.g., PCA coordinates).

* **Process:** It projects the cell to a hidden dimension, passes it through Transformer blocks modulated by AdaLN, and outputs a vector of the same size as the input.

* **Output:** The **Velocity Vector**. If you add this vector to the current cell state, you move closer to the perturbed state.

## 3.0 Cross-Validation Strategy

Standard random splits are invalid for this task because they would leak gene information. We use **Leave-Gene-Group-Out Cross-Validation**:

1. **Group K-Fold:** The script identifies all unique perturbation genes (e.g., TP53, KRAS, MYC...).

2. **Splitting:** It splits these *genes* into 5 folds.

3. **Zero-Shot Test:** In Fold 1, the model might train on TP53 and KRAS but is validated on MYC.

4. **Goal:** This proves the model isn't memorizing "Cell A turns into Cell B", but is learning "Gene Vector V causes a transformation in Direction D".

## 3.1 Training Strategy (Simple Split)

We use a **Train/Validation Split** strategy to ensure the model generalizes to unseen genes without the computational cost of K-Fold Cross-Validation.

1. **Gene Splitting:** The script identifies all unique perturbation genes (e.g., TP53, KRAS, MYC...). It randomly sets aside 10% of these genes as a **Validation Set**.

2. **Zero-Shot Validation:** During training, the model **never sees** cells perturbed by the validation genes. Every 50 epochs, we pause and ask the model to predict the velocity for these unseen genes.

3. **Checkpointing:** We save the model weights (`best_model.pt`) only when the validation loss improves. This prevents overfitting and ensures the final model is the one that generalized best.


## 4. How to Run on HPC

### Step 1: Prepare Data

You need two files:

1. **`data.h5ad`**: Scanpy object with `.obsm['X_pca']` (or `X_scvi`) and a column `target_gene`. 

2. **`gene_vectors.txt`**: Text file with format `GENE 0.123 0.456 ...`.

### Step 2: Create Submission Script (`run.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=OT-CFM
#SBATCH --gres=gpu:1           # Request 1 GPU
#SBATCH --mem=32G              # Memory
#SBATCH --time=24:00:00        # Time limit

# Load environment
source activate my_env

# Run Training
python train_flow.py \
  --adata_path "/path/to/your_data.h5ad" \
  --gene_vec_path "/path/to/gene_vectors.txt" \
  --save_dir "./checkpoints_v1" \
  --emb_key "X_pca" \
  --control_label "non-targeting" \
  --epochs 2000 \
  --k_folds 5 \
  --hidden_dim 512 \
  --num_layers 8 \
  --num_heads 8
```

## 5. Prediction

To predict the effect of a gene (e.g., "EGFR") on control cells:

```bash
python predict_flow.py \
  --checkpoint "./checkpoints_v1/best_model_fold_0.pt" \
  --config "./checkpoints_v1/config_fold_0.json" \
  --gene_vec_path "/path/to/gene_vectors.txt" \
  --adata_path "/path/to/your_data.h5ad" \
  --target_gene "EGFR" \
  --num_cells 500 \
  --output_path "prediction_EGFR.npy" \
  --hidden_dim 512 \
  --num_layers 8 \
  --num_heads 8
```

```bash
python predict_ot_cfm.py \
  --checkpoint "./checkpoints_simple/best_model.pt" \
  --config "./checkpoints_simple/config.json" \
  --gene_vec_path "esm2_gene_vectors_small.txt" \
  --adata_path "/path/to/your_data.h5ad" \
  --target_gene "EGFR" \
  --num_cells 500 \
  --output_path "prediction_EGFR.npy" \
  --hidden_dim 512 \
  --num_layers 8 \
  --num_heads 8
```

**Note:** You must use the exact same architecture flags (`hidden_dim`, `num_layers`, `num_heads`) as you did in training.

### Predicted Data Format

The output is a **Numpy Array (`.npy`)** containing the predicted latent embeddings.

* **Shape:** `(500, 50)` (assuming 500 cells and 50 PCA components). 

* **Usage:** Load this into Python, creates a new AnnData object, and visualize with UMAP to see how the predicted cells shift away from the controls.
