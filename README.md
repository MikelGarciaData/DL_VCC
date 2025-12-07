# GeneFlow: Conditional Flow Matching for Single-Cell Perturbation

## Jupyter Notebook to run and understand the model: [geneFlow.ipynb](https://github.com/MikelGarciaData/DL_VCC/blob/main/geneFlow.ipynb)

## Main file for GeneFlow model: [train_m_flow.py](https://github.com/MikelGarciaData/DL_VCC/blob/main/geneFlow/train_m_flow.py)

This repository of **GeneFlow**, a generative model designed to predict the effects of genetic perturbations on single cells.

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

## Predictions:

For predictions we picked up random control cells, and predicted the effects of the genes from the validation set on those control cells.

We ploted the trajectory of the predictions:

<img width="7471" height="1466" alt="trajectory_scanpy_umap_256" src="https://github.com/user-attachments/assets/8c26284a-aecc-463e-9109-e20d80e9d524" />

We can see that **SMARCA5** forms a cluster seperate from the other perturbations. Most of the perturbations don't change the cell state by much hence they are clustered together. Also it must have been difficult for our model to learn the perturbations of new genes as there we only 119 genes in the training set and there might have not been similar genes. But it managed to push the cells perturbed with **SMARCA5** towards a similar state. There must have been a gene similar to SMARCA5 in the training set, and the model learnt its effects.

In the future we will train a decoder that can decode the latent space to get the real gene expressions of the cells. This will allow us to evaluate our model rigorously for true biological predictions. 


