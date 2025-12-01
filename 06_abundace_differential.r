library(dplyr)
library(Matrix)
#BiocManager::install("DESeq2")
#BiocManager::install("EnhancedVolcano")
library(ggplot2)
library(ggrepel)
library(plyr)
library(biomaRt)
library(Seurat)
library(NMF)
library(miloR)
library(Matrix)
library(SingleCellExperiment)
install.packages("dplyr")
BiocManager::install("Biobase")
BiocManager::install("biomaRt")
install.packages("NMF")
BiocManager::install("miloR")
BiocManager::install("Seurat")

# --- Pipeline ---
setwd("~/denv")
seurat_integrado <- readRDS("RDS/sct_harmony_merged.rds")
integrated_sct_harmony <- seurat_integrado

# Anotação de genes
genes <- rownames(integrated_sct_harmony)

# Usar biomaRt para obter anotações
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Obter anotações relevantes
annot <- getBM(attributes = c("external_gene_name","ensembl_gene_id","gene_biotype","chromosome_name","start_position"),
               filters = "external_gene_name", values = genes, mart = mart)

# Remover duplicatas mantendo a primeira ocorrência
annot_df <- annot[!duplicated(annot$external_gene_name), ]

# Definir nomes das linhas como nomes dos genes
rownames(annot_df) <- annot_df$external_gene_name

# Filtrar para manter apenas os genes presentes no objeto Seurat
#verificar o gene_biotype
nc_genes_all <- rownames(annot_df)[grepl("lncRNA|miRNA|misc_RNA|lnc|linc|antisense|processed_pseudogene|snoRNA|snRNA|transcribed_processed_pseudogene|transcribed_unitary_pseudogene|unprocessed_pseudogene|vault_RNA", 
                                          annot_df$gene_biotype, ignore.case = TRUE)]

lnc_genes_all <- rownames(annot_df)[grepl("lncRNA|miRNA", 
                                         annot_df$gene_biotype, ignore.case = TRUE)]

# Genes codificadores de proteínas
pc_genes_all  <- rownames(annot_df)[annot_df$gene_biotype == "protein_coding"]

# Genes imunoglobulinas
ig_genes_all  <- rownames(annot_df)[grepl("^IG", rownames(annot_df))]

# salvar no objeto para referência
integrated_sct_harmony@misc$annot_df <- annot_df
integrated_sct_harmony@misc$lncRNA_all <- lnc_genes_all
integrated_sct_harmony@misc$nc_genes_all <- nc_genes_all
integrated_sct_harmony@misc$pc_genes_all <- pc_genes_all
integrated_sct_harmony@misc$ig_genes_all <- ig_genes_all

# Detectar os HVGs incluindo lncRNAs

# Definir o ensaio padrão para SCT
DefaultAssay(integrated_sct_harmony) <- "SCT"

# HVGs padrão
hvg <- VariableFeatures(integrated_sct_harmony)

# lncRNAs detectáveis no objeto
lnc_detect <- intersect(lnc_genes_all, rownames(integrated_sct_harmony))

# opcional: filtrar lncRNAs por detecção mínima por célula (ex.: >= 0.5% das células)
pct_detect <- rowSums(GetAssayData(integrated_sct_harmony, layer = "data")[lnc_detect, ] > 0) / ncol(integrated_sct_harmony)
lnc_keep <- lnc_detect[pct_detect >= 0.005] # ajuste 0.005 = 0.5%
hvg2 <- unique(c(hvg, lnc_keep))

# PCA sobre HVGs expandidos
integrated_sct_harmony <- RunPCA(integrated_sct_harmony, features = hvg2, verbose = FALSE)
ElbowPlot(integrated_sct_harmony, ndims = 50)

#NMF para identificar tópicos
# matriz genes x células com dados log-normalizados (SCT data)
mat <- as.matrix(GetAssayData(integrated_sct_harmony, slot = "data")[hvg2, ])

# rodar NMF; testar ranks 6 a 10 e seeds diferentes
rank <- 8
nmf_res <- nmf(mat, rank = rank, nrun = 5, seed = 123)
W <- basis(nmf_res)  # genes x topics
H <- coef(nmf_res)   # topics x cells

# adicionar scores de topic ao meta
for(i in seq_len(nrow(H))) {
  integrated_sct_harmony[[paste0("topic", i)]] <- H[i, ]
}

# extrair top genes por topic
top_genes_by_topic <- apply(W, 2, function(x) names(sort(x, decreasing = TRUE))[1:30])

_________________________________________________________________________________________
# Subset por amostras pediátricas e condições DF/DHF
meta <- integrated_sct_harmony@meta.data

# identificar sample_ids pediátricos
child_samples <- unique(meta$orig.ident[meta$age == "child"])

# subset por sample_id e por condição DF ou DHF
cells_child <- colnames(integrated_sct_harmony)[
  integrated_sct_harmony$orig.ident %in% child_samples &
  integrated_sct_harmony$dengue_classification %in% c("DF", "DHF", "control") &
  integrated_sct_harmony$group %in% c("acute")
]

children_obj <- subset(integrated_sct_harmony, cells = cells_child)
_________________________________________________________________________________________
# Analise com miloR
# converter Seurat para SingleCellExperiment
sce_children <- as.SingleCellExperiment(children_obj, assay = "SCT", layer = "data")

# Construir Milo a partir do objeto subsetado
milo_children <- Milo(sce_children)

# ajustar k conforme número total de células pediátricas
#n_cells_total <- ncol(children_obj)
#k_use <- ifelse(n_cells_total < 5000, 20,
#                ifelse(n_cells_total < 20000, 15, 20))

k_use <- 20

# construir grafo e nhoods
milo_children <- buildGraph(milo_children, k = k_use, d = 30)
milo_children <- makeNhoods(milo_children, prop = 0.1, k = k_use, d = 30)

# Contar células por nhood usando o meta do subset
milo_children <- countCells(  milo_children,
                              meta.data = children_obj@meta.data,
                              samples   = "orig.ident" 
)

# Preparar design corretamente a partir do meta do subset
meta_sub <- children_obj@meta.data

# garantir que condition seja factor com níveis DF e DHF
meta_sub$condition <- factor(meta_sub$dengue_classification, levels = c("DF", "DHF"))
design.df <- meta_sub %>%
  distinct(orig.ident, dengue_classification) %>%
  mutate(condition = factor(dengue_classification, levels = c("DF","DHF")))

# converter para data.frame e definir rownames
design.df <- as.data.frame(design.df)
rownames(design.df) <- design.df$orig.ident
design.df$orig.ident <- NULL   # remover coluna redundante

# se preselected_flag não existir, crie com FALSE
#if(!"preselected_flag" %in% colnames(meta_sub)) meta_sub$preselected_flag <- FALSE

#design_df <- data.frame(condition = meta_sub$condition,
#                        preselected = meta_sub$preselected_flag)

# Testar nhoods controlando preselected_flag
#da_res_children <- testNhoods(milo_children, design = ~ condition + preselected, fdr.weighting = TRUE)
da_res <- testNhoods(
  milo_children,
  design.df = design.df,
  design    = ~ condition,
  fdr.weighting = "neighbour-distance"
)

#Construção de grafo de nhoods 
milo_children <- buildNhoodGraph(milo_children)

# Visualizar resultados
plotNhoodGraphDA(milo_children, da_res, layout = "UMAP", alpha = 0.3)
plotNhoodGraphDA(milo_children, da_res, layout = "UMAP", alpha = 0.3, significance = 0.05)

#associnar informações
annotate_nhoods_by_celltype <- function(milo_obj, meta_df, celltype_col = "cell_type",
                                        min_cells_in_nhood = 5, tie_method = c("first","random")) {
  tie_method <- match.arg(tie_method)
  
  # 1) Matriz célula × nhood (dgCMatrix) e alinhamento de células
  nhood_mat <- milo_obj@nhoods
  stopifnot(inherits(nhood_mat, "dgCMatrix"))
  cells <- rownames(nhood_mat)
  if (is.null(cells)) stop("nhoods matrix must have rownames = cell barcodes")
  if (!all(cells %in% rownames(meta_df))) {
    stop("Some cells in milo_obj@nhoods are missing from meta_df rownames")
  }
  meta_df <- meta_df[cells, , drop = FALSE]
  
  # 2) Vetor de tipos celulares
  cell_types <- meta_df[[celltype_col]]
  if (is.null(cell_types)) stop(sprintf("Column '%s' not found in meta_df", celltype_col))
  cell_types <- as.character(cell_types)
  
  # 3) Matriz indicador células × tipos (esparsa)
  #    Cria um fator ordenado para tipos e usa modelo matriz para indicador 0/1
  ct_factor <- factor(cell_types)
  ct_levels <- levels(ct_factor)
  # Modelo matriz sem intercepto: uma coluna por nível
  X_ct <- sparse.model.matrix(~ ct_factor - 1)  # n_cells × n_types
  colnames(X_ct) <- ct_levels
  rownames(X_ct) <- cells
  
  # 4) Contagens por tipo × nhood: (n_types × n_cells) %*% (n_cells × n_nhoods)
  #    Resultado: n_types × n_nhoods
  type_counts <- t(X_ct) %*% nhood_mat
  
  # 5) Tamanho do nhood (número de células) e proporções
  nhood_sizes <- colSums(nhood_mat)
  # Evitar divisão por zero
  safe_sizes <- ifelse(nhood_sizes == 0, 1, nhood_sizes)
  type_props <- sweep(type_counts, 2, safe_sizes, "/")
  
  # 6) Rótulo por maioria com resolução de empates
  #    Maior contagem; em empate: "first" usa a primeira ordem, "random" sorteia
  max_idx <- apply(type_counts, 2, function(col) {
    if (all(col == 0)) return(NA_integer_)
    m <- max(col)
    w <- which(col == m)
    if (length(w) == 1) w else if (tie_method == "first") w[1] else sample(w, 1)
  })
  majority_label <- ifelse(is.na(max_idx), NA_character_, ct_levels[max_idx])
  majority_prop <- mapply(function(j, i) {
    if (is.na(i)) return(NA_real_)
    as.numeric(type_props[i, j])
  }, seq_len(ncol(type_counts)), max_idx)
  
  # 7) Data.frame de anotação por nhood
  nhood_anno <- data.frame(
    Nhood = seq_len(ncol(nhood_mat)),
    nhood_size = as.integer(nhood_sizes),
    majority_cell_type = majority_label,
    majority_prop = majority_prop,
    small_nhood_flag = nhood_sizes < min_cells_in_nhood,
    stringsAsFactors = FALSE
  )
  
  nhood_anno
}


# Construir anotação por nhood
nhood_anno <- annotate_nhoods_by_celltype(
  milo_obj      = milo_children,
  meta_df       = children_obj@meta.data,
  celltype_col  = "cell_type",
  min_cells_in_nhood = 5,
  tie_method    = "first"   # ou "random"
)

# Juntar com resultados de DA
da_res_anno <- merge(da_res, nhood_anno, by = "Nhood", all.x = TRUE)

# Condição pela direção do logFC
da_res_anno$condition <- ifelse(da_res_anno$logFC > 0, "DHF",
                                ifelse(da_res_anno$logFC < 0, "DF", NA))

# Tabela final (opcional: use SpatialFDR se preferir correção espacial)
sig_table <- subset(
  da_res_anno,
  FDR < 0.05,
  select = c("Nhood", "logFC", "FDR", "majority_cell_type", "majority_prop",
             "nhood_size", "small_nhood_flag", "condition")
)

sig_spatial_table <- subset(
  da_res_anno,
  SpatialFDR < 0.05, 
  select = c("Nhood", "logFC", "FDR", "majority_cell_type", "majority_prop",
             "nhood_size", "small_nhood_flag", "condition")
)


# Salvar tudo
saveRDS(milo_children, file = "RDS/milo_children.rds")

write.csv(da_res, file =  "miloR/da_res.csv", row.names = FALSE)
write.csv(sig_spatial_table, file = "miloR/sig_spatial_table.csv", row.names = FALSE)
write.csv(sig_table, file = "miloR/sig_table.csv", row.names = FALSE)
