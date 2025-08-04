import scanpy as sc
import harmonypy as hm

# Carregar dados
adata = sc.read("Anndata/dados_integrados.h5ad")

# Pré-processamento
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, subset=True)
sc.pp.scale(adata)
sc.tl.pca(adata)

# Rodar Harmony
ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, 'amostra')

# Substituir PCA por Harmony
adata.obsm["X_harmony"] = ho.Z_corr.T

# UMAP e clustering
sc.pp.neighbors(adata, use_rep="X_harmony")
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Visualizar
sc.pl.umap(adata, color=["amostra", "leiden"])


import scvi

# Preparar dados
scvi.model.SCVI.setup_anndata(adata, batch_key="amostra")

# Treinar modelo
model = scvi.model.SCVI(adata)
model.train()

# Obter embeddings
adata.obsm["X_scVI"] = model.get_latents()

# UMAP e clustering
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Visualizar
sc.pl.umap(adata, color=["amostra", "leiden"])


# Assumindo que adata.obs["labels"] tem algumas anotações
scvi.model.SCANVI.setup_anndata(adata, batch_key="amostra", labels_key="labels")

scanvi_model = scvi.model.SCANVI(adata, unlabeled_category="unknown")
scanvi_model.train()

# Embeddings
adata.obsm["X_scanVI"] = scanvi_model.get_latents()

# UMAP e clustering
sc.pp.neighbors(adata, use_rep="X_scanVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata, color=["amostra", "leiden", "labels"])


from scib.metrics import nmi, ari, silhouette, graph_connectivity

# Exemplo com Harmony
nmi_score = nmi(adata, group_key="cell_type", label_key="leiden")
ari_score = ari(adata, group_key="cell_type", label_key="leiden")
sil_score = silhouette(adata, embed="X_harmony", group_key="amostra")
gc_score = graph_connectivity(adata, label_key="cell_type")

print(f"NMI: {nmi_score:.3f}, ARI: {ari_score:.3f}, Silhouette: {sil_score:.3f}, GraphConn: {gc_score:.3f}")


import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

# Carregar seu objeto AnnData
import scanpy as sc
adata = sc.read("Anndata/dados_integrados.h5ad")

# Verifique os embeddings disponíveis
print("Embeddings:", adata.obsm.keys())


# Verifique os rótulos para avaliação
print("Labels disponíveis:", adata.obs.columns)

label_key = "leiden"  # ou "cell_type"

# Função para graph connectivity
def graph_connectivity(embedding, labels, k=10):
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(embedding)
    _, indices = nbrs.kneighbors(embedding)
    scores = []
    for i in range(len(embedding)):
        neighbors = indices[i][1:]  # exclui a própria célula
        same_label = sum(labels[j] == labels[i] for j in neighbors)
        scores.append(same_label / k)
    return np.mean(scores)

# Avaliar cada embedding
embeddings_to_test = ["X_pca", "X_harmony", "X_scVI", "X_scanVI"]
results = []

for emb in embeddings_to_test:
    print(f"Avaliando {emb}...")
    X = adata.obsm[emb]
    labels = adata.obs[label_key].values

    sil_score = silhouette_score(X, labels)
    gc_score = graph_connectivity(X, labels)

    results.append({
        "Embedding": emb,
        "Silhouette Score": round(sil_score, 3),
        "Graph Connectivity": round(gc_score, 3)
    })

# Mostrar resultados
df_results = pd.DataFrame(results)
print(df_results)
