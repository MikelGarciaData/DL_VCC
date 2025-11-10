# DL_VCC
Deep Learning Virtual Cell Challenge

[Environment Setup](environment.yml)

raw data should be put here:
[raw data](data/raw_data/)

a viewer for the h5 file can be found here:
[viewer](https://myhdf5.hdfgroup.org/)

We are gonna use the cell embedding model from the STATE paper (named SE model) to embed the gene counts data into a latent space. The model can be found here:
https://huggingface.co/arcinstitute/SE-600M

Then we will use flow matching diffusion model for training. It will train on the embedded data and will train to predict how the cell state changes, from control cells (target_gene = "non-targeting") to pertubed state (the genes with a target gene)

<img width="1280" height="306" alt="image" src="https://github.com/user-attachments/assets/b4f967fb-9cfa-41f8-8755-ba0b36cb508b" />

The white are the control cells and C1 C2 are the perturbed states.

Once we get the pertubed state in the latent space we will decode it back the gene counts. And then run benchmarks studies. The model will be trained on 150 pertubations roughly 145,000 cells have them. One cell has one genetic pertubation. And we 38,000 control cells with no pertubation.

<img width="515" height="409" alt="cee75766-089d-40b9-922b-ce853601472e" src="https://github.com/user-attachments/assets/4d9c2e08-977d-4244-b4d4-b38069980ce5" />

We can see that the control cells have a similar gene expression (the blob in the middle) but also most of the genetic pertubations have very low effect on the over all gene expression. There are some that have higher effect, they are in the outside clusters.

Once the flow model learns to change the state of the cell, we will validate it with new data that wasn't in the training set. New pertubations, same cell type.

Learning stuff: 

[what diffusion models are](https://youtu.be/iv-5mZ_9CPY?si=yb6IEbxzJG7K6-al)


[flow matching 1](https://youtu.be/7cMzfkWFWhI?si=trWp7UBoSivH-uwf)

[flow matching 2]([https://youtu.be/7cMzfkWFWhI?si=trWp7UBoSivH-uwf](https://www.youtube.com/watch?v=7NNxK3CqaDk))


[Flow matching for single cell data (Marginal flow matching - maybe use this model)](https://youtu.be/I6zCSrs60eA?si=8xoLRO9o8wIlpBDg)

You can read more about the SE model in this paper: 
https://github.com/MikelGarciaData/DL_VCC/blob/main/docs/2025.06.26.661135v2.full.pdf



