library(dplyr)
library(Matrix)
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(plyr)
library(biomaRt)
library(Seurat)

# --- Pipeline ---
setwd("~/denv")
seurat_integrado <- readRDS("RDS/sct_harmony_merged.rds")
table(seurat_integrado$SCT_snn_res.0.2)
seurat_integrado$SCT_snn_res.0.2

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

saveRDS(seurat_integrado, file = "RDS/sct_harmony_merged.rds")

# --- 1. Extrair metadados

meta <- seurat_integrado@meta.data
meta$cell_id <- rownames(meta)

# --- 2. Counts
DefaultAssay(seurat_integrado) <- "SCT"
counts <- GetAssayData(seurat_integrado, layer = "counts")

# --- 3. Escolher tipo celular: monócitos
meta_mono <- meta[meta$cell_type == "Monocyte", ]

# --- 4. Tabela para o DESeq2 ---

# Agrupar células por amostra
cells_by_sample <- split(meta_mono$cell_id, meta_mono$orig.ident)

# Somar contagens por amostra
pseudobulk_counts <- sapply(cells_by_sample, function(cells) {
  rowSums(counts[, cells, drop = FALSE])
})

# Construir metadados por amostra
sample_meta <- meta_mono %>%
  filter(orig.ident %in% colnames(pseudobulk_counts)) %>%
  distinct(orig.ident, dengue_classification) %>%
  as.data.frame()


# Garantir que as colunas de pseudobulk_counts e sample_meta estejam na mesma ordem
rownames(sample_meta) <- sample_meta$orig.ident #definir rownames
sample_meta <- sample_meta[colnames(pseudobulk_counts), , drop = FALSE] #reordenar sample_meta pelo colnames de pseudobulk_counts

# --- 5. DESeq2 ---

# Criar objeto
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = sample_meta,
  design = ~ dengue_classification
)

# executar Deseq2
dds <- DESeq(dds)

# Extrair comparações
# comparar DF vs DHF
res <- results(dds, contrast = c("dengue_classification", "DF", "DHF")) #comparar DF vs DHF
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# comparar DF vs control
resdf <- results(dds, contrast = c("dengue_classification", "DF", "control")) #comparar DF vs control
resdf_df <- as.data.frame(resdf)
resdf_df$gene <- rownames(resdf_df)

# comparar DF vs control
resdhf <- results(dds, contrast = c("dengue_classification", "DHF", "control")) #comparar DF vs control
resdhf_df <- as.data.frame(resdhf)
resdhf_df$gene <- rownames(resdhf_df)

# --- 6. Filtrar genes significativos ---
# --- dengue: DF vs DHF ---
res_df$significant <- "NS" 
res_df$significant[res_df$padj <= 0.05 & abs(res_df$log2FoldChange) >= 1]  <- "Up"
res_df$significant[res_df$padj <= 0.05 & res_df$log2FoldChange <= -1] <- "Down"

# --- DF vs control ---
resdf_df$significant <- "NS" 
resdf_df$significant[resdf_df$padj <= 0.05 & abs(resdf_df$log2FoldChange) >= 1]  <- "Up"
resdf_df$significant[resdf_df$padj <= 0.05 & resdf_df$log2FoldChange <= -1] <- "Down"

# --- DHF vs control ---
resdhf_df$significant <- "NS" 
resdhf_df$significant[resdhf_df$padj <= 0.05 & abs(resdhf_df$log2FoldChange) >= 1]  <- "Up"
resdhf_df$significant[resdhf_df$padj <= 0.05 & resdhf_df$log2FoldChange <= -1] <- "Down"

# --- Salvar resultados---
dir.create("DEG_results/monocytes/plots", recursive = TRUE)
write.csv(res_df, "DEG_results/monocytes/DF_DHF.csv", row.names = FALSE)
write.csv(resdf_df, "DEG_results/monocytes/DF_control.csv", row.names = FALSE)
write.csv(resdhf_df, "DEG_results/monocytes/DhF_control.csv", row.names = FALSE)

# --- Visualizações ---
# Volcano plot
keyvals <- ifelse(res_df$log2FoldChange < -1 & res_df$padj <= 0.05 , 'royalblue',
                  ifelse(res_df$log2FoldChange > 1 & res_df$padj <= 0.05, 'red', 'gray'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up Regulated'
names(keyvals)[keyvals == 'gray'] <- 'No Significance'
names(keyvals)[keyvals == 'royalblue'] <- 'Down Regulated'
png("DEG_results/monocytes/plots/volcano_DF_DHF.png", width=3000, height=2500, units="px", pointsize=20, bg="white",res=300)
EnhancedVolcano(res_df, 
                lab=rownames(res_df),title = 'DF vs DHF (Monocytes)', 
                subtitle = ".", x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=1, pointSize = 3.0, labSize = 4.0, 
                colCustom = keyvals, xlim = c(-3, 3), ylim = c(0,5))  + theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                )
dev.off()



keyvals <- ifelse(resdf_df$log2FoldChange < -2 & resdf_df$padj <= 0.05 , 'royalblue',
                  ifelse(resdf_df$log2FoldChange > 2 & resdf_df$padj <= 0.05, 'red', 'gray'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up Regulated'
names(keyvals)[keyvals == 'gray'] <- 'No Significance'
names(keyvals)[keyvals == 'royalblue'] <- 'Down Regulated'
png("DEG_results/monocytes/plots/volcano_DF_control.png", width=3000, height=2500, units="px", pointsize=20, bg="white",res=300)
EnhancedVolcano(resdf_df, 
                lab=rownames(resdf_df),title = 'DF vs Control (Monocytes)', 
                subtitle = ".", x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 4.0, 
                colCustom = keyvals, xlim = c(-8, 8), ylim = c(0,15))  + theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                )
dev.off()


keyvals <- ifelse(resdhf_df$log2FoldChange < -2 & resdhf_df$padj <= 0.05 , 'royalblue',
                  ifelse(resdhf_df$log2FoldChange > 2 & resdhf_df$padj <= 0.05, 'red', 'gray'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up Regulated'
names(keyvals)[keyvals == 'gray'] <- 'No Significance'
names(keyvals)[keyvals == 'royalblue'] <- 'Down Regulated'
png("DEG_results/monocytes/plots/volcano_DHF_control.png", width=3000, height=2500, units="px", pointsize=20, bg="white",res=300)
EnhancedVolcano(resdhf_df, 
                lab=rownames(resdhf_df),title = 'DHF vs Control (Monocytes)', 
                subtitle = ".", x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 4.0, 
                colCustom = keyvals , xlim = c(-8, 8), ylim = c(0,15)) + theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                )
dev.off()
