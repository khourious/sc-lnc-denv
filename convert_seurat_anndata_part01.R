library(Seurat)
library(Matrix)
library(tidyverse)


setwd("~/denv/")
rds_files <- list.files("integration", pattern = "*.rds", full.names = TRUE)

dir.create("Anndata", showWarnings = FALSE)
dir.create("Anndata/matrix", showWarnings = FALSE)
dir.create("Anndata/pca", showWarnings = FALSE)
dir.create("Anndata/gene_names", showWarnings = FALSE)
dir.create("Anndata/metadata", showWarnings = FALSE)
dir.create("Anndata/harmony", showWarnings = FALSE)


# --- Matriz de contagem
counts_matrix <- GetAssayData(combined, assay = "RNA", slot = "counts")
writeMM(counts_matrix, file = paste0("Anndata/matrix/", "matrix.mtx"))

# --- Nomes dos genes
write.table(data.frame('gene'=rownames(counts_matrix)),
            file = paste0("Anndata/gene_names/", "gene_names.csv"), 
            quote=F, row.names=F,col.names=f)

# --- Embeddings do Harmony
harmony_embeddings <- combined@reductions$harmony@cell.embeddings
write.csv(harmony_embeddings,
            file = paste0("Anndata/harmony", "harmony.csv"),
            quote = FALSE, row.names = TRUE)


# --- Embeddings do UMAP (opcional)
combined$barcode <- colnames(combined)

combined$UMAP_1 <- combined@reductions$umap@cell.embeddings[,1]
combined$UMAP_2 <- combined@reductions$umap@cell.embeddings[,2]


# --- Metadados
write.csv(combined@meta.data,
          file = paste0("Anndata/metadata/", "metadata.csv"),
          quote = FALSE, row.names = TRUE)
