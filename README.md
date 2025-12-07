# geneFlow

## Jupyter Notebook to run the model: [geneFlow.ipynb](https://github.com/MikelGarciaData/DL_VCC/blob/main/geneFlow.ipynb)


# GeneFlow: Conditional Flow Matching for Single-Cell Perturbation

This repository contains the PyTorch implementation of **GeneFlow**, a generative model designed to predict the effects of genetic perturbations on single cells.

GeneFlow aims for **Zero-Shot Generalization**: predicting the transcriptomic phenotype of a cell after perturbing a gene the model has *never seen before*, provided it has access to a vector representation of that gene (e.g., from GenePT).

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

The script handles both training and the subsequent trajectory generation in one go.

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

Deep Learning Virtual Cell Challenge

[Environment Setup](environment.yml)

raw data should be put here:
[raw data](data/raw_data/)

a viewer for the h5 file can be found here:
[viewer](https://myhdf5.hdfgroup.org/)

We are gonna use the cell embedding model from the STATE paper (named SE model) to embed the gene counts data into a latent space. The model can be found here:
https://huggingface.co/arcinstitute/SE-600M

Then we will use flow matching diffusion model for training. It will train on the embedded data and will train to predict how the cell state changes, from control cells (target_gene = "non-targeting") to pertubed state (the genes with a target gene)

<img alt="image" src="https://github.com/user-attachments/assets/363caaf8-58b2-4ffa-9a05-e2d0fbeaa545" />


<img alt="image" src="https://github.com/user-attachments/assets/b4f967fb-9cfa-41f8-8755-ba0b36cb508b" />

The white are the control cells and C1 C2 are the perturbed states.

Once we get the pertubed state in the latent space we will decode it back the gene counts. And then run benchmarks studies. The model will be trained on 150 pertubations roughly 145,000 cells have them. One cell has one genetic pertubation. And we 38,000 control cells with no pertubation.

<img alt="cee75766-089d-40b9-922b-ce853601472e" src="https://github.com/user-attachments/assets/4d9c2e08-977d-4244-b4d4-b38069980ce5" />

We can see that the control cells have a similar gene expression (the blob in the middle) but also most of the genetic pertubations have very low effect on the over all gene expression. There are some that have higher effect, they are in the outside clusters.

Once the flow model learns to change the state of the cell, we will validate it with new data that wasn't in the training set. New pertubations, same cell type.

Learning stuff: 

[what diffusion models are](https://youtu.be/iv-5mZ_9CPY?si=yb6IEbxzJG7K6-al)


[flow matching 1](https://youtu.be/7cMzfkWFWhI?si=trWp7UBoSivH-uwf)

[flow matching 2](https://www.youtube.com/watch?v=7NNxK3CqaDk)


[Flow matching for single cell data (Marginal flow matching - maybe use this model)](https://youtu.be/I6zCSrs60eA?si=8xoLRO9o8wIlpBDg)

You can read more about the SE model in this paper: 
https://github.com/MikelGarciaData/DL_VCC/blob/main/docs/state_paper.pdf


