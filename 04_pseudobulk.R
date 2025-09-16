library(Matrix)
library(dplyr)
library(DESeq2)
BiocManager::install("DESeq2")

seurat_integrado <- readRDS("RDS/sct_harmony_merged.rds")

table(seurat_integrado$SCT_snn_res.0.2)

cluster_names <- c(
  "0" = "NK cell",
  "1" = "T Cell CD4+",
  "2" = "B cell",
  "3" = "Plasmablast",
  "4" = "Monocyte",
  "5" = "T Cell CD4+",
  "6" = "T Cell CD8+",
  "7" = "T Cell CD8+",
  "8" = "T Cell CD8+",
  "9" = "Plasmablast",
  "10" = "B cell",
  "11" = "Dendritic Cell",
  "12" = "Plasma Cell",
  "13" = "Platelets",
  "14" = "Red Cell"
)


seurat_integrado$cell_type <- mapvalues(
  seurat_integrado$SCT_snn_res.0.2,
  from = names(cluster_names),
  to = cluster_names
)

DimPlot(seurat_integrado, group.by = "cell_type", label = TRUE) + ggtitle("Anotação Manual dos Clusters")

write.csv(seurat_integrado@meta.data, "metadata_cluster_anotado.csv")

# --- 1. Extrair metadados

meta <- seurat_integrado@meta.data
meta$cell_id <- rownames(meta)

# --- 2. Escolher tipo celular
meta_mono <- meta[meta$cell_type == "Monocyte", ]


# --- 3. Agrupar por amostra
cells_by_sample <- split(meta_mono$cell_id, meta_mono$orig.ident)


# --- 4. Somar contagens por amostra

DefaultAssay(seurat_integrado) <- "SCT"
counts <- GetAssayData(seurat_integrado, layer = "counts")

# Soma por amostra
pseudobulk_counts <- sapply(cells_by_sample, function(cells) {
  rowSums(counts[, cells, drop = FALSE])
})
table(meta_mono$dengue_classification)
table(meta_mono$orig.ident, meta_mono$dengue_classification)


# --- 5. DeSeq2

# Coldata: Metadados por amostra
amostras <- colnames(pseudobulk_counts)

# Filtrar metadados para essas amostras
sample_meta <- meta_mono %>%
  dplyr::select(orig.ident, dengue_classification)

sample_meta <- unique(sample_meta)

# Objeto deseq2
dds <- DESeqDataSetFromMatrix(pseudobulk_counts,
                              colData = sample_meta,
                              design = ~ dengue_classification)


# Executar Deseq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("dengue_classification", "DF", "DHF"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

library(ggplot2)

res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: DF vs DHF (Monocytes)", x = "log2 Fold Change", y = "-log10 FDR")

library(pheatmap)

top_genes <- res_df %>% arrange(padj) %>% head(30) %>% pull(gene)
heatmap_data <- pseudobulk_counts[top_genes, ]

pheatmap(log2(heatmap_data + 1), cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = sample_meta["dengue_classification"],
         main = "Top 30 Genes - Monocytes")

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = res_df$gene,
  mart = mart
)

res_annot <- merge(res_df, annot, by.x = "gene", by.y = "external_gene_name")

coding_genes <- res_annot[res_annot$gene_biotype == "protein_coding", ]
lncRNAs <- res_annot[grepl("lncRNA", res_annot$gene_biotype), ]
