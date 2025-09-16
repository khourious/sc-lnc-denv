library(plyr)
library(ggplot2)
library(biomaRt)
library(Seurat)


seurat_integrado <- readRDS("caminho/para/seu_arquivo.rds")
table(seurat_integrado$seurat_clusters)

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
  seurat_integrado$seurat_clusters,
  from = names(cluster_names),
  to = cluster_names
)

DimPlot(seurat_integrado, group.by = "cell_type", label = TRUE) + ggtitle("Anotação Manual dos Clusters")

write.csv(seurat_integrado@meta.data, "metadata_cluster_anotado.csv")


mono <- subset(seurat_integrado, subset = cell_type == "Monocyte")

# DENV Fever vs DENV Hemorrhagic Fever
markers_mono <- FindMarkers(
  mono,
  ident.1 = "DF",
  ident.2 = "DHF",
  group.by = "dengue_classification",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# DENV Fever vs Controle
markers_mono_DF <- FindMarkers(
  mono,
  ident.1 = "DF",
  ident.2 = "control",
  group.by = "dengue_classification",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# DENV Hemorrhagic Fever vs Controle
markers_mono_DHF <- FindMarkers(
  mono,
  ident.1 = "DHF",
  ident.2 = "control",
  group.by = "dengue_classification",
  min.pct = 0.25,
  logfc.threshold = 0.25
)


# Top 10 genes
head(markers_mono_DF , 10)

# Plot de violino
VlnPlot(seurat_integrado, features = c("CD3D", "MS4A1"), group.by = "cell_type")

# Heatmap
top_genes <- rownames(markers_mono_DF)[1:20]
DoHeatmap(seurat_integrado, features = top_genes, group.by = "cell_type")


# ViolinPlot
markers <- markers_mono_DF
markers$gene <- rownames(markers)
markers$log10_pval <- -log10(markers$p_val_adj)

ggplot(markers, aes(x = avg_log2FC, y = log10_pval)) +
  geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(
    title = "MOnocytes: DENV Fever vs Controle",
    x = "Log2 Fold Change",
    y = "-log10(p-valor ajustado)"
  )


# loop para os tipos celulares

# Inicializa o mart do Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

cell_types <- unique(seurat_integrado$cell_type)
marker_list <- list()

for (ct in cell_types) {
  subset_ct <- subset(seurat_integrado, subset = cell_type == ct)
  
  # Comparação 1: DENV Fever vs Controle
  markers_fever <- FindMarkers(
    subset_ct,
    ident.1 = "DF",
    ident.2 = "control",
    group.by = "dengue_classification",
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  markers_fever$gene <- rownames(markers_fever)
  markers_fever$cell_type <- ct
  markers_fever$comparison <- "DF_vs_control"
  
  # Anotação com biomaRt
  gene_annot <- getBM(
    attributes = c("external_gene_name", "gene_biotype"),
    filters = "external_gene_name",
    values = markers_fever$gene,
    mart = ensembl
  )
  
  gene_annot$category <- ifelse(
    gene_annot$gene_biotype == "protein_coding", "coding",
    ifelse(grepl("lncRNA", gene_annot$gene_biotype), "lncRNA", "other")
  )
  
  markers_fever_annotated <- merge(
    markers_fever,
    gene_annot,
    by.x = "gene",
    by.y = "external_gene_name",
    all.x = TRUE
  )
  
  write.csv(
    markers_fever_annotated,
    file = paste0("markers_", ct, "_DF_vs_control_annotated.csv"),
    row.names = FALSE
  )
  
  # Comparação 2: DENV Hemorrhagic Fever vs Controle
  markers_hemo <- FindMarkers(
    subset_ct,
    ident.1 = "DHF",
    ident.2 = "control",
    group.by = "dengue_classification",
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  markers_hemo$gene <- rownames(markers_hemo)
  markers_hemo$cell_type <- ct
  markers_hemo$comparison <- "DHF_vs_control"
  
  gene_annot_hemo <- getBM(
    attributes = c("external_gene_name", "gene_biotype"),
    filters = "external_gene_name",
    values = markers_hemo$gene,
    mart = ensembl
  )
  
  gene_annot_hemo$category <- ifelse(
    gene_annot_hemo$gene_biotype == "protein_coding", "coding",
    ifelse(grepl("lncRNA", gene_annot_hemo$gene_biotype), "lncRNA", "other")
  )
  
  markers_hemo_annotated <- merge(
    markers_hemo,
    gene_annot_hemo,
    by.x = "gene",
    by.y = "external_gene_name",
    all.x = TRUE
  )
  
  write.csv(
    markers_hemo_annotated,
    file = paste0("markers_", ct, "_DHF_vs_control_annotated.csv"),
    row.names = FALSE
  )
  
  # Armazena na lista (opcional)
  marker_list[[paste0(ct, "_DF")]] <- markers_fever_annotated
  marker_list[[paste0(ct, "_DHF")]] <- markers_hemo_annotated
}

library(dplyr)
library(ggplot2)

markers_annotated %>%
  group_by(cluster, category) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = cluster, y = freq, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Proporção de genes coding vs lncRNA por cluster", y = "Proporção", x = "Cluster") +
  theme_minimal()




library(Seurat)
library(biomaRt)
library(ggplot2)
library(dplyr)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
cell_types <- unique(seurat_integrado$cell_type)

for (ct in cell_types) {
  subset_ct <- subset(seurat_integrado, subset = cell_type == ct)
  
  for (comp in c("DF", "DHF")) {
    comparison_label <- paste0(comp, "_vs_control")
    
    # Expressão diferencial
    markers <- FindMarkers(
      subset_ct,
      ident.1 = comp,
      ident.2 = "control",
      group.by = "dengue_classification",
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    markers$gene <- rownames(markers)
    markers$cell_type <- ct
    markers$comparison <- comparison_label
    
    # Anotação com biomaRt
    gene_annot <- getBM(
      attributes = c("external_gene_name", "gene_biotype"),
      filters = "external_gene_name",
      values = markers$gene,
      mart = ensembl
    )
    
    gene_annot$category <- ifelse(
      gene_annot$gene_biotype == "protein_coding", "coding",
      ifelse(grepl("lncRNA", gene_annot$gene_biotype), "lncRNA", "other")
    )
    
    markers_annotated <- merge(
      markers,
      gene_annot,
      by.x = "gene",
      by.y = "external_gene_name",
      all.x = TRUE
    )
    
    # Salva CSV
    write.csv(
      markers_annotated,
      file = paste0("markers_", ct, "_", comparison_label, "_annotated.csv"),
      row.names = FALSE
    )
    
    # Volcano Plot
    markers_annotated$log10_pval <- -log10(markers_annotated$p_val_adj)
    markers_annotated$significant <- with(markers_annotated, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
    
    volcano <- ggplot(markers_annotated, aes(x = avg_log2FC, y = log10_pval)) +
      geom_point(aes(color = significant)) +
      scale_color_manual(values = c("gray", "red")) +
      theme_minimal() +
      labs(
        title = paste("Volcano Plot:", ct, comparison_label),
        x = "Log2 Fold Change",
        y = "-log10(p-valor ajustado)"
      )
    
    ggsave(
      filename = paste0("volcano_", ct, "_", comparison_label, ".png"),
      plot = volcano,
      width = 7, height = 5
    )
    
    # Heatmap dos top 20 genes
    top_genes <- markers_annotated %>%
      arrange(p_val_adj) %>%
      slice_head(n = 20) %>%
      pull(gene)
    
    heatmap <- DoHeatmap(
      subset_ct,
      features = top_genes,
      group.by = "dengue_classification"
    ) + ggtitle(paste("Heatmap:", ct, comparison_label))
    
    ggsave(
      filename = paste0("heatmap_", ct, "_", comparison_label, ".png"),
      plot = heatmap,
      width = 8, height = 6
    )
  }
}

