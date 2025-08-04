import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from scipy import io

# Diretório principal onde estão os dados exportados do Seurat
base_dir = "Anndata"

# Subpastas organizadas
matrix_dir = os.path.join(base_dir, "matrix")
pca_dir = os.path.join(base_dir, "pca")
gene_dir = os.path.join(base_dir, "gene_names")
meta_dir = os.path.join(base_dir, "metadata")


# Detectar automaticamente os nomes das amostras
samples = sorted([
    f.replace("_matrix.mtx", "")
    for f in os.listdir(matrix_dir)
    if f.endswith("_matrix.mtx")
])

print("Amostras encontradas:", samples)

# Lista para guardar os objetos AnnData
adatas = []

# Loop por cada amostra
for sample in samples:
    print(f"Processando {sample}...")

    # Caminhos dos arquivos
    matrix_path = os.path.join(matrix_dir, f"{sample}_matrix.mtx")
    genes_path = os.path.join(gene_dir, f"{sample}_gene_names.csv")
    meta_path = os.path.join(meta_dir, f"{sample}_metadata.csv")
    pca_path = os.path.join(pca_dir, f"{sample}_pca.csv")

    # Carregar matriz de contagem
    X = io.mmread(matrix_path).transpose().tocsr()

    # Carregar nomes dos genes
    gene_names = pd.read_csv(genes_path)["gene"].tolist()

    # Carregar metadados
    metadata = pd.read_csv(meta_path, index_col=0)

    # Criar objeto AnnData
    adata = anndata.AnnData(X=X, obs=metadata)
    adata.var_names = gene_names
    adata.obs_names = metadata.index

    # Carregar PCA
    pca = pd.read_csv(pca_path, index_col=0)
    adata.obsm["X_pca"] = pca.loc[adata.obs_names].to_numpy()

    # Carregar UMAP (já está nos metadados)
    adata.obsm["X_umap"] = np.vstack([
        adata.obs["UMAP_1"].to_numpy(),
        adata.obs["UMAP_2"].to_numpy()
    ]).T

    # Adicionar nome da amostra
    adata.obs["amostra"] = sample

    # Salvar individualmente (opcional)
    adata.write(os.path.join(base_dir, f"{sample}.h5ad"))

    # Adicionar à lista
    adatas.append(adata)

# 🔗 Juntar todas as amostras em um único objeto
adata_all = adatas[0].concatenate(*adatas[1:], batch_key="amostra")

# 🔹 Visualizar UMAP colorido por amostra
sc.pl.umap(adata_all, color=["amostra", "seurat_clusters"])

# 🔹 Salvar o conjunto integrado
adata_all.write(os.path.join(base_dir, "dados_integrados.h5ad"))



# Caminho para o arquivo
adata = sc.read("Anndata/dados_integrados.h5ad")

# Carregar metadados extras
extra_meta = pd.read_csv("metadados_extra.csv", index_col=0)
meta_amostra = pd.read_csv("metadados_por_amostra.csv")

# Verificar se os barcodes batem
print(set(extra_meta.index).issubset(set(adata.obs_names)))

# Adicionar ao AnnData
adata.obs = adata.obs.join(extra_meta)

# Mapear cada coluna
for col in ["paciente", "condição"]:
    mapping = dict(zip(meta_amostra["amostra"], meta_amostra[col]))
    adata.obs[col] = adata.obs["amostra"].map(mapping)
