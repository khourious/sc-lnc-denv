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
seurat_integrado <- readRDS("RDS/sct_harmony_annotation.rds")
integrated_sct_harmony <- seurat_integrado
table(seurat_integrado$RNA_snn_res.0.3)
seurat_integrado$RNA_snn_res.0.3

FeaturePlot(seurat_integrado, features = c("HBB","HBA1","HBA2")) # Red Cells

# Remover células que expressam HBB, HBA1 ou HBA2
seurat_integrado <- subset(
  seurat_integrado,
  subset = HBB < 1 & HBA1 < 1 & HBA2 < 1
)

FeaturePlot(seurat_integrado, features = c("PPBP","GP9","ITGA2B","TUBB1")) # plaletes

seurat_integrado <- subset(
  seurat_integrado,
  subset = PPBP < 1 & TUBB1 < 1 & ITGA2B < 1 & GP9 < 1
)

Idents(sub) <- sub$RNA_snn_res.0.3
sub <- SCTransform(sub, assay = "RNA", verbose = TRUE)
sub <- PrepSCTFindMarkers(sub, assay = "SCT", verbose = TRUE)
markers_cluster0<- FindMarkers(sub, ident.1 = 0, min.pct = 0.25)
head(markers_cluster0, n = 20)



# Verifique os clusters atuaisIdents(seurat_integrado) <- seurat_integrado$RNA_snn_res.0.3

# Exclui o cluster 
#seurat_integrado <- subset(seurat_integrado, idents = setdiff(levels(Idents(seurat_integrado)), "14"))

#sub$refined_clusters <- Idents(sub)
#seurat_integrado$refined_clusters <- NA
#seurat_integrado$refined_clusters[Cells(sub)] <- sub$refined_clusters

# Confirme que foi removido
table(Idents(seurat_integrado))

Idents(seurat_integrado) <- seurat_integrado$RNA_snn_res.0.3
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
        group.by = "RNA_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = TCD4_markers, 
        group.by = "RNA_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = TCD8_markers, 
        group.by = "RNA_snn_res.0.3")


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
        group.by = "RNA_snn_res.0.3")

DotPlot(object   = seurat_integrado,
        features = genes_markers_myeloid, 
        group.by = "RNA_snn_res.0.3")

DimPlot(
  seurat_integrado,
  reduction = "umap",
  label = TRUE,        # coloca o número dentro do cluster
  label.size = 6,      # tamanho da fonte
  repel = TRUE         # evita sobreposição
)

Idents(seurat_integrado) <- seurat_integrado$RNA_snn_res.0.3
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

Idents(seurat_integrado) <- "RNA_snn_res.0.3"
names(cluster_names) <- levels(seurat_integrado)

seurat_integrado$celltype_0.3  <- plyr::mapvalues(
  Idents(seurat_integrado),
  from = names(cluster_names),
  to   = cluster_names
)

# Cluster 0
DimPlot(sub, label = TRUE)
seurat_integrado$celltype_0.5_0.3  <- as.character(seurat_integrado$celltype_0.5_0.3 )

sub$refined_clusters <- Idents(sub)
seurat_integrado$celltype_0.3 [Cells(sub)[sub$refined_clusters == "0"]] <- "T cell naive"
seurat_integrado$celltype_0.3 [Cells(sub)[sub$refined_clusters == "1"]] <- "T cell naive"
seurat_integrado$celltype_0.3 [Cells(sub)[sub$refined_clusters == "2"]] <- "Red Cell"

# cluster 4
seurat_integrado$celltype_0.3  <- as.character(seurat_integrado$celltype_0.3 )

sub1$refined_clusters <- Idents(sub1)
seurat_integrado$celltype_0.3 [Cells(sub1)[sub1$refined_clusters == "0"]] <- "T CD4+ effector"
seurat_integrado$celltype_0.3 [Cells(sub1)[sub1$refined_clusters == "1"]] <- "T CD4+ effector"
seurat_integrado$celltype_0.3 [Cells(sub1)[sub1$refined_clusters == "2"]] <- "T CD4+ effector"
seurat_integrado$celltype_0.3 [Cells(sub1)[sub1$refined_clusters == "3"]] <- "T CD4+ effector"

seurat_integrado$celltype_0.3  <- factor(seurat_integrado$celltype_0.3 )


DimPlot(seurat_integrado, group.by = "celltype_0.3", label = TRUE)

DimPlot(seurat_integrado, group.by = "celltype_0.3", label = TRUE) + ggtitle("Anotação Manual dos Clusters")

saveRDS(seurat_integrado, file = "RDS/sct_harmony_annotation.rds")

# --- Umaps ---
DimPlot(seurat_integrado, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "celltype_0.3") +
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

DimPlot(seurat_integrado, reduction = "umap", label = FALSE, pt.size = 0.6, group.by = "celltype_0.3",
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


DimPlot(seurat_integrado, reduction = "umap", label = TRUE, pt.size = 0.6, label.size = 7, group.by = "celltype_0.3",
        cols = c("T cell Naive" = "#32994d",
                 "NK cell" = "#301c6f", 
                 "T CD4+ effector" = "#016c59",
                 "T CD8+ effector" = "#2f7298",
                 "Monocytes" =  "#a248ae", 
                 "B cell" = "#de6223",
                 "Plasmablasts" =  "#ed2a48",
                 "cDendritic" = "#ee894f",
                 "pDendritic" = "#c84659", 
                 "Pre-Plasmablasts" = "#f06772")) +
  ggtitle("PBMC cell types")

# --- Refinando anotação ----
Idents(seurat_integrado) <- seurat_integrado$RNA_snn_res.0.5
DimPlot(
  seurat_integrado,
  reduction = "umap",
  label = TRUE,        # coloca o número dentro do cluster
  label.size = 6,      # tamanho da fonte
  repel = TRUE         # evita sobreposição
)

# --- identificando os clusters ----

Idents(seurat_integrado) <- seurat_integrado$RNA_snn_res.0.5

DimPlot(seurat_integrado, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)

genes_markers_B <- list(
  Bcell          = c("MS4A1","CD79A","CD79B"),
  Bcells_naive   = c("TCL1A","IGHD"), 
  Bcells_memory  = c("CD27","BANK1","HLA-DRA"),
  Plasmablasts   = c("MZB1","XBP1","JCHAIN","SDC1","TNFRSF17","IGHG1","IGKC") 
)

genes_subset_B <- list(
  Bcell          = c("CD79A", "MS4A1", "CD19", "HLA-DRA"),
  Bcells_naive   = c("IGHD","IGHM"), 
  Bcells_memory  = c("CD27","IGHG1", "IGHA1","BANK1"), 
  Plasma   = c("JCHAIN","MZB1","CD38","SDC1","TNFRSF17","IGKC"),
  P_Plasma  = c("MKI67") 
)

genes_markers_myeloid <- list(
  Monocytes_classic     = c("LYZ","S100A8","S100A9","S100A12","CD14","VCAN","MNDA"), # cluster 5
  Monocytes_nonclassic  = c("FCGR3A","MS4A7","MS4A4A","CXCL10","APOBEC3A","CSF1R"), # cluster 8
  Monocytes_intermediate= c("HLA-DRA","CD74"), #cluster 16
  Dendritic_cDC2        = c("FCER1A","CLEC10A","CD1C"), 
  Dendritic_pDC         = c("LILRA4","IL3RA","IRF7"), # cluster 15
)

DotPlot(object   = seurat_integrado,
        features = genes_markers_B, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = genes_subset_B, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = genes_markers_myeloid, 
        group.by = "RNA_snn_res.0.5")

DimPlot(
  seurat_integrado,
  reduction = "umap",
  label = TRUE,        
  label.size = 6,      
  repel = TRUE         
)


# Monócitos
Monocyte_markers <- list(
  Classical     = c("LYZ"), 
  NonClassical  = c("CD14","FCGR3A"), 
  Intermediate  = c("C1QA")
)

# NK
NK_markers <- list(
  NK_cytotoxic = c("CD3D","NKG7","GNLY"), 
  NK_signaling = c("XCL1","XCL2")
)

# B cells
B_markers <- list(
  Naive     = c("CD79A","MS4A1","TCL1A"), 
  Memory    = c("IGHA1"), 
  Activated = c("IGHM","IGHG1")
)

# T CD4
TCD4_markers <- list(
  Naive   = c("CD3D","IL7R","TCF7"), 
  Memory  = c("CD27", "LAG3", "KLRB1"), 
  Effector= c("STAT4","PRKCH","FYN","TOX","IFNG")
)

# T CD8
TCD8_markers <- list(
  Naive     = c("CD3D","IL7R","TCF7","CD8A","CD8B"), 
  Effector  = c("GZMB","PRF1","NKG7","GNLY"), 
  Memory    = c("CD27"), 
  Exhausted = c("LAG3","H1FX")
)

# Tregs
Treg_markers <- list(
  Treg = c("FOXP3","IL2RA","CTLA4")
)

# MAIT
MAIT_markers <- list(
  MAIT = c("CD3D","IL7R","KLRB1")
)

DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)

# Exemplos de DotPlot
DotPlot(object   = seurat_integrado,
        features = Monocyte_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = NK_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = B_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = TCD4_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = TCD8_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = Treg_markers, 
        group.by = "RNA_snn_res.0.5")

DotPlot(object   = seurat_integrado,
        features = MAIT_markers, 
        group.by = "RNA_snn_res.0.5")


cluster_names <- c(
  "0" = "CD4+ Naive T",
  "1" = "Cytotoxit NK",
  "2" = "Monocytes",
  "3" = "Naive B",
  "4" = "CD8+ Naive T",
  "5" = "Plasmablast",
  "6" = "CD4+ Effector T",
  "7" = "CD4+ Memory T",
  "8" = "Memory B",
  "9" = "CD8+ Effector Memory T",
  "10" = "Plasmablast",
  "11" = "Pre-Plasmablast",
  "12" = "Activated B",
  "13" = "Signaling NK",
  "14" = "Cytotoxit NK",
  "15" = "Activated B",
  "16" = "pDC",
  "17" = "Plasmablast",
  "18" = "Monocytes",
  "19" = "CD4+ Naive T",
  "20" = "Pre-Plasmablast",
  "21" = "CD4+ Memory T",
  "22" = "cDC",
  "23" = "Monocytes"
)

Idents(seurat_integrado) <- "RNA_snn_res.0.5"
names(cluster_names) <- levels(seurat_integrado)

seurat_integrado$celltype_0.5 <- plyr::mapvalues(
  Idents(seurat_integrado),
  from = names(cluster_names),
  to   = cluster_names
)

DimPlot(seurat_integrado, group.by = "RNA_snn_res.0.5", label = TRUE)

DimPlot(seurat_integrado, group.by = "celltype_0.5", label = TRUE)

# --- ajustando clusters ---
red <- subset(seurat_integrado, idents = 4)
red <- FindVariableFeatures(red)
red <- NormalizeData(red)
red <- ScaleData(red)
red <- RunPCA(red)
red <- FindNeighbors(red, dims = 1:5)
red <- FindClusters(red, resolution = 0.1)
red <- RunUMAP(red, dims = 1:8)
DimPlot(red, label = TRUE)
Idents(red) <- red$RNA_snn_res.0.3

FeaturePlot(red, features = c("CD3D","IL7R","TCF7","CD8A","CD8B")) # CD8+ Naive T
FeaturePlot(red, features = c("LAG3","H1FX")) # CD8+ exhausted T



# cluster monocitos
plaq <- subset(seurat_integrado, idents = c(2,18,23))
plaq <- FindVariableFeatures(plaq)
plaq <- NormalizeData(plaq)
plaq <- ScaleData(plaq)
plaq <- RunPCA(plaq)
plaq <- FindNeighbors(plaq, dims = 1:3)
plaq <- FindClusters(plaq, resolution = 0.3)
plaq <- RunUMAP(plaq, dims = 1:3)
DimPlot(plaq, label = TRUE)

FeaturePlot(plaq, features = c("LYZ","S100A8","S100A9","S100A12","CD14","VCAN","MNDA")) # Monocytes_classic 2 3 4 5 6 7 9
FeaturePlot(plaq, features = c("FCGR3A","MS4A7","MS4A4A","CXCL10","APOBEC3A","CSF1R")) # Monocytes_nonclassic 0 1
FeaturePlot(plaq, features = c("HLA-DRA","CD74")) # Monocytes_intermediate 
FeaturePlot(plaq, features = c("FCER1A","CLEC10A","CD1C")) # Dendritic_cDC 8 
FeaturePlot(plaq, features = c("LILRA4","IL3RA","IRF7")) # Dendritic_pDC  

DimPlot(red, label = TRUE)

seurat_integrado$celltype_0.5 <- as.character(seurat_integrado$celltype_0.5)

red$refined_clusters <- Idents(red)
seurat_integrado$celltype_0.5[Cells(red)[red$refined_clusters == "0"]] <- "CD8+ exhausted T"
seurat_integrado$celltype_0.5[Cells(red)[red$refined_clusters == "1"]] <- "CD8+ exhausted T"
seurat_integrado$celltype_0.5[Cells(red)[red$refined_clusters == "2"]] <- "CD8+ Naive T"
seurat_integrado$celltype_0.5[Cells(red)[red$refined_clusters == "3"]] <- "CD8+ Naive T"
seurat_integrado$celltype_0.5[Cells(red)[red$refined_clusters == "4"]] <- "CD8+ Naive T"

plaq$refined_clusters <- Idents(plaq)
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "0"]] <- "Non Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "1"]] <- "Non Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "2"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "3"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "4"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "5"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "6"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "7"]] <- "Classical Monocytes"
seurat_integrado$celltype_0.5[Cells(plaq)[plaq$refined_clusters == "8"]] <- "cDC"

DimPlot(seurat_integrado, group.by = "celltype_0.5", label = TRUE)

celltype_colors <- c(
  # T cells (verdes e azuis diferenciados)
  "CD4+ Naive T"            = "#7fcdbb",
  "CD4+ Effector T"         = "#4eb3d3", 
  "CD4+ Memory T"           = "#1078c1",
  
  "CD8+ Naive T"            = "#c7e9c0",  
  "CD8+ Effector Memory T"  = "#60d256",
  "CD8+ exhausted T"        = "#006d2c",

  "Signaling NK"           = "#10adc1",   
  "Cytotoxit NK"           = "#045a8d", 
  
  # Monócitos (roxo)
  "Classical Monocytes"     = "#8073ac",  
  "Non Classical Monocytes"      = "#542788",   
  
  # Plasmablastos (rosa)
  "Plasmablast"             =  "#c51b7d",    
  "Pre-Plasmablast"         = "#de77ae",  
  
  # B cells (tons de marrom/laranja queimado bem distintos)
  "Naive B"                    = "#feb24c", 
  "Activated B"                = "#993404",   
  "Memory B"                   = "#d94801",  
  
   # Outros (tons neutros)
  "pDC"                     = "#d73027",
  "cDC"                     = "#d53e4f"
)

DimPlot(seurat_integrado, group.by = "celltype_0.5", cols = celltype_colors, 
        label = FALSE,  pt.size = 0.7)

saveRDS(seurat_integrado, file = "RDS/sct_harmony_annotation.rds")


# --- unificando anotação ---

seurat_integrado$celltype_0.5_fine <- dplyr::case_when(
  # T CD4+
  seurat_integrado$celltype_0.5 %in% c("CD4+ Naive T") ~ "CD4+ Naive T",
  seurat_integrado$celltype_0.5 %in% c("CD4+ Effector T", "CD4+ Memory T") ~ "CD4+ ActMem T",
  
  # T CD8+
  seurat_integrado$celltype_0.5 %in% c("CD8+ Naive T") ~ "CD8+ Naive T",
  seurat_integrado$celltype_0.5 %in% c("CD8+ Effector Memory T", "CD8+ exhausted T") ~ "CD8+ ActMem T",
  
  # B cells
  seurat_integrado$celltype_0.5 %in% c("Naive B") ~ "Naive B",
  seurat_integrado$celltype_0.5 %in% c("Memory B", "Activated B") ~ "ActMem B",
  
  # Plasmablastos
  seurat_integrado$celltype_0.5 %in% c("Plasmablast") ~ "Plasmablast",
  seurat_integrado$celltype_0.5 %in% c("Pre-Plasmablast") ~ "Pre-Plasmablast",
  
  # Monócitos
  seurat_integrado$celltype_0.5 %in% c("Classical Monocytes", "Non Classical Monocytes", "Monocytes") ~ "Monocytes",
  
  # NK
  seurat_integrado$celltype_0.5 %in% c("Signaling NK", "Cytotoxit NK") ~ "NK",
  
  # DCs
  seurat_integrado$celltype_0.5 %in% c("cDC") ~ "cDC",
  seurat_integrado$celltype_0.5 %in% c("pDC") ~ "pDC"
)

unique(seurat_integrado$celltype_0.5[is.na(seurat_integrado$celltype_0.5_fine)])


DimPlot(seurat_integrado, group.by = "celltype_0.5_fine", 
        label = TRUE,  pt.size = 0.7)

celltype_colors_fine <- c(
  # T cells (verdes e azuis diferenciados)
  "CD4+ Naive T"            = "#7fcdbb",
  "CD4+ ActMem T"           = "#1078c1",
  
  "CD8+ Naive T"            = "#c7e9c0",  
  "CD8+ ActMem T"           = "#006d2c",
  
  "NK"                      = "#045a8d", 
  
  # Monócitos (roxo)
  "Monocytes"              = "#a940d4",   
  
  # Plasmablastos (rosa)
  "Plasmablast"             =  "#c51b7d",    
  "Pre-Plasmablast"         = "#de77ae",  
  
  # B cells (tons de marrom/laranja queimado bem distintos)
  "Naive B"                    = "#feb24c", 
  "ActMem B"                   = "#d94801",  
  
  # Outros (tons neutros)
  "pDC"                     = "#fa120a",
  "cDC"                     = "#d53e4f"
)

DimPlot(seurat_integrado, group.by = "celltype_0.5_fine", 
        cols = celltype_colors_fine, label = TRUE,  pt.size = 0.7, label.size = 4)

# --- Pseudobulk

integrated_sct_harmony <- JoinLayers(seurat_integrado, assay = "RNA")

pseudobulk <- AggregateExpression(
  integrated_sct_harmony,
  group.by = c("celltype_0.5_fine", "orig.ident"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

counts_mat <- pseudobulk$RNA
dim(counts_mat) 
colnames(counts_mat)    
rownames(counts_mat)[1:10] 

write.csv(counts_mat, "pseudobulk_counts.csv")
