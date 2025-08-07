import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from scipy import io

# Diretório base onde estão os arquivos
base_dir = "Anndata"

matrix_path = os.path.join(base_dir, "matrix", "matrix.mtx")
genes_path = os.path.join(base_dir, "gene_names", "gene_names.csv")
meta_path = os.path.join(base_dir, "metadata", "metadata.csv")
pca_path = os.path.join(base_dir, "pca", "pca.csv")
harmony_path = os.path.join(base_dir, "harmony", "harmony.csv")


X = io.mmread(matrix_path).transpose().tocsr()


metadata = pd.read_csv(meta_path, index_col=0)


with open(genes_path, 'r') as f:
    gene_names = f.read().splitlines()


adata = anndata.AnnData(X=X, obs=metadata)
adata.var_names = gene_names
adata.obs_names = metadata.index


pca = pd.read_csv(pca_path, index_col=0)
pca = pca.loc[adata.obs_names]
adata.obsm["X_pca"] = pca.to_numpy()


harmony = pd.read_csv(harmony_path, index_col=0)
harmony = harmony.loc[adata.obs_names]
adata.obsm["X_harmony"] = harmony.to_numpy()


adata.obsm["X_umap"] = np.vstack([adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy()]).T


adata.layers["counts"] = adata.X.copy()

sc.pl.umap(adata, color=["seurat_clusters"])

adata.obsm["X_harmony"][0, :5]
adata.obsm["X_umap"][:5]
adata.obsm["X_pca"][:5]
adata.obsm["counts"][0, :5]