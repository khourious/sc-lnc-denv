library(Seurat)
library(Matrix)
library(tidyverse)

# Caminho da pasta com os arquivos .rds
rds_files <- list.files("qc_doublets", pattern = "*_singlets.rds", full.names = TRUE)

if (!dir.exists("Anndata/matrix")) {
  dir.create("Anndata/matrix")
}

if (!dir.exists("Anndata/pca")) {
  dir.create("Anndata/pca")
}

if (!dir.exists("Anndata/gene_names")) {
  dir.create("Anndata/gene_names")
}

if (!dir.exists("Anndata/metadata")) {
  dir.create("Anndata/metadata")
}


for (rds_path in rds_files) {
  sample_name <- tools::file_path_sans_ext(basename(rds_path))
  message("Processando: ", sample_name)
  
  cells <- readRDS(rds_path)
  
  # Exportar matriz de contagem
  counts_matrix <- GetAssayData(NML, assay = "RNA", slot = "counts")
  writeMM(counts_matrix, file = paste0("Anndata/matrix", sample_name, "_matrix.mtx"))
  
  # Exportar nomes dos genes
  write.table(data.frame(gene = rownames(counts_matrix)),
              file = paste0("Anndata/gene_names", sample_name, "_gene_names.csv"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Exportar PCA
  write.csv(NML@reductions$pca@cell.embeddings,
            file = paste0("Anndata/pca", sample_name, "_pca.csv"),
            quote = FALSE, row.names = TRUE)
  
  # Adicionar UMAP e barcode
  cells$barcode <- colnames(NML)
  cells$UMAP_1 <- Embeddings(NML, "umap")[, 1]
  cells$UMAP_2 <- Embeddings(NML, "umap")[, 2]
  
  # Exportar metadata
  write.csv(cells@meta.data,
            file = paste0("Anndata/metadata", sample_name, "_metadata.csv"),
            quote = FALSE, row.names = TRUE)
}
