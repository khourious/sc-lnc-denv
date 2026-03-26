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
seurat_integrado <- readRDS("sct_harmony_merged_novo.rds")
seurat_integrado <- integrated_sct_harmony
table(seurat_integrado$SCT_snn_res.0.3)
seurat_integrado$SCT_snn_res.0.3

sub <- subset(seurat_integrado, idents = 0)
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.06)
sub <- RunUMAP(sub, dims = 1:6)
DimPlot(sub, label = TRUE)

FeaturePlot(sub, features = c("STAT4","PRKCH","FYN","TOX","IFNG", "CD3D", "CD4")) # tcd efector
FeaturePlot(sub, features = c("HBB","HBA1","HBA2")) # Red Cells

Idents(sub) <- sub$SCT_snn_res.0.3
sub <- SCTransform(sub, assay = "RNA", verbose = TRUE)
sub <- PrepSCTFindMarkers(sub, assay = "SCT", verbose = TRUE)
markers_cluster0<- FindMarkers(sub, ident.1 = 0, min.pct = 0.25)
head(markers_cluster0, n = 20)


sub1 <- subset(seurat_integrado, idents = 4)
sub1 <- FindVariableFeatures(sub1)
sub1 <- ScaleData(sub1)
sub1 <- RunPCA(sub1)
sub1 <- FindNeighbors(sub1, dims = 1:5)
sub1 <- FindClusters(sub1, resolution = 0.07)
sub1 <- RunUMAP(sub1, dims = 1:6)
DimPlot(sub1, label = TRUE)

FeaturePlot(sub1, features = c("PPBP","GP9","ITGA2B","TUBB1")) # plaletes
FeaturePlot(sub1, features = c("STAT4","PRKCH","FYN","TOX","IFNG"))

Idents(sub) <- sub$SCT_snn_res.0.3
sub <- SCTransform(sub, assay = "RNA", verbose = TRUE)
sub <- PrepSCTFindMarkers(sub, assay = "SCT", verbose = TRUE)
markers_cluster0<- FindMarkers(sub, ident.1 = 0, min.pct = 0.25)
head(markers_cluster0, n = 20)



# Verifique os clusters atuaisIdents(seurat_integrado) <- seurat_integrado$SCT_snn_res.0.3

# Exclui o cluster 14
seurat_integrado <- subset(seurat_integrado, idents = setdiff(levels(Idents(seurat_integrado)), "14"))
seurat_integrado <- subset(seurat_integrado, idents = setdiff(levels(Idents(seurat_integrado)), "16"))

sub$refined_clusters <- Idents(sub)
seurat_integrado$refined_clusters <- NA
seurat_integrado$refined_clusters[Cells(sub)] <- sub$refined_clusters

# Confirme que foi removido
table(Idents(seurat_integrado))


DimPlot(
  seurat_integrado,
  reduction = "umap",
  label = TRUE,        # coloca o número dentro do cluster
  label.size = 6,      # tamanho da fonte
  repel = TRUE         # evita sobreposição
)

# --- identificando os clusters ----

FeaturePlot(seurat_integrado, features = c("GNLY","NKG7","KLRD1","PRF1")) # NK cluster 2
FeaturePlot(seurat_integrado, features = c("CD3D","CD3E","CD4","IL7R")) # T CD4 cluster 0 e 13
FeaturePlot(seurat_integrado, features = c("CD3D","CD3E","CD8A","CD8B","GZMB", "GZMK")) # T CD8 cluster 3
FeaturePlot(seurat_integrado, features = c("MS4A1","CD79A","CD79B")) # B cell cluster 1
FeaturePlot(seurat_integrado, features = c("MZB1","JCHAIN","IGKC","IGHG1")) # Plasmablast cluster 6 7 11 10 9
FeaturePlot(seurat_integrado, features = c("CD14","LYZ","S100A8","S100A9")) # Monocyte classical cluster 5
FeaturePlot(seurat_integrado, features = c("FCGR3A","MS4A7")) # Monocyte non classical cluster 8
FeaturePlot(seurat_integrado, features = c("FCER1A","CLEC10A","LILRA4","IL3RA","HLA-DRA")) # Dendritic cluster 15
FeaturePlot(seurat_integrado, features = c("CD27","BANK1","HLA-DRA")) # b cell memory
FeaturePlot(seurat_integrado, features = c("PPBP","GP9","ITGA2B","TUBB1")) # Platelets
FeaturePlot(seurat_integrado, features = c("HBB","HBA1","HBA2")) # Red Cells

genes_markers_effector <- list(
  TCell    = c("CD3E"), #cluster 0  2  3  4  13
  T_naive  = c("IL7R","CCR7","SELL", "CD62L"), # cluster 0 e 13
  T_memory = c("CD27", "CD44"),
  
  T_CD4    = c("CD4", "STAT4","PRKCH","FYN","TOX","IFNG"), # cluster  4 e 15
  T_CD8    = c("CD8A", "GZMB","PRF1","NKG7","GNLY","GZMK"), # cluster 2 3 
  T_reg    = c("CD25"),
  
  NK_classic       = c("CD16","KLRD1","KLRF1", "CD127"), # cluster 3        
  NK_proliferating = c("MKI67","TOP2A","PCNA","STMN1")
)

TCD4_markers <-list( 
  TCD4_naive      = c("IL7R","CCR7","SELL"),  # cluster 0 e 13
  TCD4_memory     = c("CD27"), 
  TCD4_effector   = c("STAT4","PRKCH","FYN","TOX","IFNG") # cluster 4
)

TCD8_markers <-list(   
  TCD8_naive       = c("IL7R","CCR7","SELL"),  # cluster 0 e 13
  TCD8_memory      = c("CD27"),                   
  TCD8_effector    = c("GZMB","PRF1","NKG7","GNLY"), # cluster 3
  NK_classic       = c("CD16","KLRD1","KLRF1", "NKp46", "KIR2DL", "KIR3DL", "GZMK"), # cluster 2 
  NK_T             = c("CD3D","TRAC","CD8A","CD4","NCAM1"), # cluster 2
  NK_proliferating = c("MKI67","TOP2A","PCNA","STMN1")
)


DotPlot(object   = seurat_integrado,
        features = genes_markers_effector, 
        group.by = "SCT_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = TCD4_markers, 
        group.by = "SCT_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = TCD8_markers, 
        group.by = "SCT_snn_res.0.3")


genes_markers_B <- list(
  Bcell          = c("MS4A1","CD79A","CD79B"), # cluster 1 7 9 10 11 12 14
  Bcells_naive   = c("TCL1A","IGHD"), # cluster 1
  Bcells_memory  = c("CD27","BANK1","HLA-DRA"), # cluster 11 12 
  Plasmablasts   = c("MZB1","XBP1","JCHAIN","SDC1","TNFRSF17","IGHG1","IGKC") #cluster 6 7 11
)

genes_markers_myeloid <- list(
  Monocytes_classic     = c("LYZ","S100A8","S100A9","S100A12","CD14","VCAN","MNDA"), # cluster 5
  Monocytes_nonclassic  = c("FCGR3A","MS4A7","MS4A4A","CXCL10","APOBEC3A","CSF1R"), # cluster 8
  Monocytes_intermediate= c("HLA-DRA","CD74"), #cluster 16
  Dendritic_cDC2        = c("FCER1A","CLEC10A","CD1C"), 
  Dendritic_pDC         = c("LILRA4","IL3RA","IRF7"), # cluster 15
  Platelets             = c("PPBP","GP9","ITGA2B","TUBB1"),
  Erythroid             = c("HBB","HBA1","HBA2","ALAS2","SLC4A1")
)

DotPlot(object   = seurat_integrado,
        features = genes_markers_B, 
        group.by = "SCT_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = genes_markers_myeloid, 
        group.by = "SCT_snn_res.0.3")

DimPlot(
  seurat_integrado,
  reduction = "umap",
  label = TRUE,        # coloca o número dentro do cluster
  label.size = 6,      # tamanho da fonte
  repel = TRUE         # evita sobreposição
)

Idents(seurat_integrado) <- seurat_integrado$SCT_snn_res.0.3
seurat_integrado <- SCTransform(seurat_integrado, assay = "RNA", verbose = TRUE)
seurat_integrado <- PrepSCTFindMarkers(seurat_integrado, assay = "SCT", verbose = TRUE)
markers_cluster16 <- FindMarkers(seurat_integrado, ident.1 = 16, min.pct = 0.25)
head(markers_cluster16, n = 20)

head(markers_cluster6, n = 20)

markers_cluster7 <- FindMarkers(seurat_integrado, ident.1 = 7, min.pct = 0.25)
head(markers_cluster7, n = 20)

markers_cluster8 <- FindMarkers(seurat_integrado, ident.1 = 8, min.pct = 0.25)
head(markers_cluster8, n = 20)

markers_cluster10 <- FindMarkers(seurat_integrado, ident.1 = 10, min.pct = 0.25)
head(markers_cluster10, n = 20)

markers_cluster11 <- FindMarkers(seurat_integrado, ident.1 = 11, min.pct = 0.25)
head(markers_cluster11, n = 20)

# --- Mudando o nome dos clusters ---

cluster_names <- c(
  "0" = "T cell naive",
  "1" = "B cell naive",
  "2" = "NK classic",
  "3" = "T CD8+ effector",
  "4" = "T CD4+ effector",
  "5" = "Monocytes classic",
  "6" = "Plasmablasts",
  "7" = "Plasmablasts",
  "8" = "Monocytes non classic",
  "9" = "Plasmablasts",
  "10" = "Plasmablasts",
  "11" = "Plasmablasts",
  "12" = "Plasmablasts",
  "13" = "T cell naive",
  "15" = "Dendritic"
)

Idents(seurat_integrado) <- "SCT_snn_res.0.3"
names(cluster_names) <- levels(seurat_integrado)

seurat_integrado$celltype_global <- plyr::mapvalues(
  Idents(seurat_integrado),
  from = names(cluster_names),
  to   = cluster_names
)

# Cluster 4
DimPlot(sub, label = TRUE)
seurat_integrado$celltype_global <- as.character(seurat_integrado$celltype_global)
seurat_integrado$celltype_global[Cells(sub)[sub$refined_clusters == "0"]] <- "T cell naive"
seurat_integrado$celltype_global[Cells(sub)[sub$refined_clusters == "1"]] <- "T cell naive"
seurat_integrado$celltype_global[Cells(sub)[sub$refined_clusters == "2"]] <- "T cell naive"
seurat_integrado$celltype_global[Cells(sub)[sub$refined_clusters == "3"]] <- "Red Cell"

# cluster 16
seurat_integrado$celltype_global <- as.character(seurat_integrado$celltype_global)

seurat_integrado$celltype_global[Cells(sub1)[sub1$refined_clusters == "0"]] <- "T CD4+ effector"
seurat_integrado$celltype_global[Cells(sub1)[sub1$refined_clusters == "1"]] <- "T CD4+ effector"
seurat_integrado$celltype_global[Cells(sub1)[sub1$refined_clusters == "2"]] <- "T CD4+ effector"
seurat_integrado$celltype_global[Cells(sub1)[sub1$refined_clusters == "3"]] <- "Platelets"
seurat_integrado$celltype_global <- factor(seurat_integrado$celltype_global)


DimPlot(seurat_integrado, group.by = "celltype_global", label = TRUE)

DimPlot(seurat_integrado, group.by = "celltype_global", label = TRUE) + ggtitle("Anotação Manual dos Clusters")

# --- Umaps ---
DimPlot(seurat_integrado, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "celltype_global") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("PBMC cell types")

DimPlot(seurat_integrado, reduction = "umap", label = FALSE, pt.size = 0.6, group.by = "celltype_global",
        cols = c("T cell naive" = "#DE757BFF",
                 "NK classic" = "#80003FFF", 
                 "T CD8+ effector" = "#DC5FBDFF",
                 "T CD4+ effector" = "#40007FFF",
                 "Monocytes classic" =  "#60A8E0FF", 
                 "Monocytes non classic" = "#073F80FF",
                 "B cell naive" = "#0F8554FF",
                 "Plasmablasts" =  "#73AF48FF",
                 "Dendritic" = "#EDAD08FF",
                 "Platelets" = "#E17C05FF", 
                 "Red Cell" = "darkred")) +
  ggtitle("PBMC cell types")


DimPlot(seurat_integrado, reduction = "umap", label = TRUE, pt.size = 0.6, label.size = 7, group.by = "celltype_global",
        cols = c("T cell naive" = "#DE757BFF",
                 "NK classic" = "#80003FFF", 
                 "T CD8+ effector" = "#DC5FBDFF",
                 "T CD4+ effector" = "#40007FFF",
                 "Monocytes classic" =  "#60A8E0FF", 
                 "Monocytes non classic" = "#073F80FF",
                 "B cell naive" = "#0F8554FF",
                 "Plasmablasts" =  "#73AF48FF",
                 "Dendritic" = "#EDAD08FF",
                 "Platelets" = "#E17C05FF", 
                 "Red Cell" = "darkred")) +
  ggtitle("PBMC cell types")


saveRDS(seurat_integrado, file = "RDS/sct_harmony_merged.rds")

# --- Pseudobulk

genes_ig <- grep("^IG[HKL]", rownames(seurat_integrado), value = TRUE)

# calcular score IG por célula
seurat_integrado$IG_score <- colSums(
  GetAssayData(seurat_integrado, layer = "counts")[genes_ig, ]
)

# normalizar pelo total de counts
seurat_integrado$IG_score <- seurat_integrado$IG_score / seurat_integrado$nCount_RNA

seurat_integrado$macrogrupo <- case_when(
  seurat_integrado$cell_type %in% c("Monocytes classic", "Monocytes non classic", "Dendritic") ~ "Myeloid",
  seurat_integrado$cell_type %in% c("T cell naive", "T CD8+ effector", "T CD4+ effector", "NK classic") ~ "Lymphoid_T_NK",
  seurat_integrado$cell_type %in% c("B cell naive", "Plasmablasts") ~ "B_lineage",
  TRUE ~ "outros"
)


obj_myeloid <- subset(seurat_integrado, subset = macrogrupo == "Myeloid")
obj_tnk     <- subset(seurat_integrado, subset = macrogrupo == "Lymphoid_T_NK")
obj_b       <- subset(seurat_integrado, subset = macrogrupo == "B_lineage")

obj_myeloid <- subset(seurat_integrado,
                      subset = macrogrupo == "Myeloid" & IG_score < 0.05)

obj_tnk <- subset(seurat_integrado,
                  subset = macrogrupo == "Lymphoid_T_NK" & IG_score < 0.05)

obj_b <- subset(seurat_integrado,
                subset = macrogrupo == "B_lineage")


export_pseudobulk_by_celltype <- function(seurat_obj, macrogrupo_nome) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  counts <- GetAssayData(seurat_obj, slot = "counts")
  meta <- seurat_obj@meta.data
  
  # lista de cell types dentro do macrogrupo
  celltypes <- unique(meta$cell_type)
  
  for (ct in celltypes) {
    message("Processando: ", macrogrupo_nome, " - ", ct)
    
    meta_ct <- meta[meta$cell_type == ct, ]
    cells_by_sample <- split(meta_ct$cell_id, meta_ct$orig.ident)
    
    # pseudobulk counts
    pb_counts <- do.call(cbind, lapply(cells_by_sample, function(cells) {
      rowSums(counts[, cells, drop = FALSE])
    }))
    
    # metadata
    sample_meta <- meta_ct %>%
      dplyr::group_by(orig.ident, cell_type, macrogrupo) %>%
      dplyr::summarise(
        dengue_classification = ifelse(all(is.na(dengue_classification)), NA, na.omit(dengue_classification)[1]),
        .groups = "drop"
      ) %>%
      as.data.frame()
    
    rownames(sample_meta) <- sample_meta$orig.ident
    sample_meta <- sample_meta[colnames(pb_counts), , drop = FALSE]
    
    # exportar arquivos separados
    fname_counts <- paste0("pseudobulk_counts_", gsub(" ", "_", ct), ".csv")
    fname_meta   <- paste0("metadata_", gsub(" ", "_", ct), ".csv")
    
    write.csv(pb_counts, file = fname_counts)
    write.csv(sample_meta, file = fname_meta)
  }
}

export_pseudobulk_by_celltype(obj_myeloid, "Myeloid")
export_pseudobulk_by_celltype(obj_tnk, "Lymphoid_T_NK")
export_pseudobulk_by_celltype(obj_b, "B_lineage")

# --- exemplo para dendritic
pb_dendritic <- read.csv("pseudobulk_counts_Dendritic.csv", row.names = 1)
rowSums(pb_dendritic[grep("^IG", rownames(pb_dendritic)), ])








#--- anterior
meta <- seurat_integrado@meta.data
meta$cell_id <- rownames(meta)

# Lista de tipos celulares
tipos_celulares <- unique(meta$celltype_global)

# Lista para guardar os pseudobulks
pseudobulk_lista <- list()
metadata_lista <- list()

for (tipo in tipos_celulares) {
  message("Processando tipo celular: ", tipo)
  
  # Filtra células do tipo atual
  meta_tipo <- meta[meta$celltype_global == tipo, ]
  
  # Agrupa células por amostra
  cells_by_sample <- split(meta_tipo$cell_id, meta_tipo$orig.ident)
  
  # Gera pseudobulk
  pseudobulk_counts <- sapply(cells_by_sample, function(cells) {
    rowSums(counts[, cells, drop = FALSE])
  })
  
  # Cria metadados por amostra
  sample_meta <- meta_tipo %>%
    group_by(orig.ident) %>%
    summarise(
      dengue_classification = first(na.omit(dengue_classification)),
      age = first(na.omit(age)),
      virus = first(na.omit(virus)),
      sex = first(na.omit(sex)),
      group = first(na.omit(group)),
      infection = first(na.omit(infection)),
      dataset = first(na.omit(dataset)),
      timepoint = first(na.omit(timepoint)),
      disease = first(na.omit(disease))
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  
  
  
  rownames(sample_meta) <- sample_meta$orig.ident
  sample_meta <- sample_meta[colnames(pseudobulk_counts), , drop = FALSE]
  
  # Salva na lista
  pseudobulk_lista[[tipo]] <- pseudobulk_counts
  metadata_lista[[tipo]] <- sample_meta
}

for (tipo in names(pseudobulk_lista)) {
  write.csv(pseudobulk_lista[[tipo]], file = paste0("pseudobulk_", tipo, ".csv"))
  write.csv(metadata_lista[[tipo]], file = paste0("metadata_", tipo, ".csv"))
}
