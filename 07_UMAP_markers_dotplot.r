library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)

seurat_obj <- readRDS("seurat_obj_reclustered.rds")

# --- testar UMAP ---

# UMAP com paleta refinada e estilo clean
library(Seurat)
library(ggplot2)

# Defina uma paleta de cores para cada população
cluster_colors <- c(
  "TCD4"        = "#1f78b4",
  "TCD8"        = "#e31a1c",
  "NK"          = "#6a3d9a",
  "Treg"        = "#33a02c",
  "Tfh"         = "#fb9a99",
  "B"           = "#ff7f00",
  "Plasmablast" = "#b15928",
  "Plasma"      = "#a6cee3",
  "Dendritic"   = "#b2df8a",
  "Mono"        = "#cab2d6"
)

DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type",
        pt.size = 0.4, label = TRUE, repel = TRUE) +
  scale_color_manual(values = cluster_colors) +
  theme_void(base_size = 14) +
  theme(legend.position = "right", legend.title = element_blank())



# --- marcadores

T CD4	CD3D, CD4, IL7R, CCR7
T CD8	CD3D, CD8A, CD8B, GZMA, GZMB, PRF1
NK	NKG7, KLRD1, GNLY, CD56 (NCAM1), PRF1
Tregs	FOXP3, IL2RA (CD25), CTLA4
T helper (folicular)	BCL6, CXCR5
Células B	MS4A1 (CD20), CD79A, CD79B
Plasmablasts	MZB1, XBP1, TNFRSF17 (BCMA)
Plasma cells	SDC1 (CD138), MZB1, XBP1
Dendríticas	FCER1A, CLEC9A, CD1C, LYZ
Monócitos	CD14, FCGR3A, LYZ


FeaturePlot(seurat_obj, features = c("CD3D","CD4","CD8A","CD56","PRF1","GZMA","MS4A1","MZB1","XBP1","BCL6"),
            cols = c("lightgrey","red"), pt.size = 0.4, ncol = 3) +
  theme_void(base_size = 14)

library(Seurat)
library(ggplot2)

População	Marcadores principais
T CD4	CD3D, CD4, IL7R, CCR7
T CD8	CD3D, CD8A, CD8B, GZMA, GZMB, PRF1
NK	NCAM1 (CD56), NKG7, GNLY, PRF1, GZMA, GZMB, KLRD1
Treg	FOXP3, IL2RA (CD25), CTLA4
Tfh (T helper folicular)	BCL6, CXCR5
Células B	MS4A1 (CD20), CD79A, CD79B
Plasmablasts	MZB1, XBP1, TNFRSF17 (BCMA)
Plasma cells	SDC1 (CD138), MZB1, XBP1
Dendríticas	FCER1A, CLEC9A, CD1C, LYZ
Monócitos	CD14, FCGR3A, LYZ

# Lista de marcadores
markers <- c(
  "CD3D","CD4","CD8A","CD8B",
  "NCAM1","NKG7","GNLY","PRF1","GZMA","GZMB","KLRD1",
  "FOXP3","IL2RA","CTLA4",
  "BCL6","CXCR5",
  "MS4A1","CD79A","CD79B",
  "MZB1","XBP1","TNFRSF17","SDC1",
  "FCER1A","CLEC9A","CD1C","LYZ",
  "CD14","FCGR3A"
)

seurat_obj$cell_type <- factor(seurat_obj$cell_type,
  levels = c("TCD4","TCD8","Treg","Tfh","NK","B","Plasmablast","Plasma","Dendritic","Mono"))

# DotPlot agrupado por tipo celular
DotPlot(seurat_obj, features = markers, group.by = "cell_type") +
  RotatedAxis() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  ggtitle("Markers by Cell Type")

ggsave("DotPlot_markers.pdf", width = 12, height = 6)
ggsave("DotPlot_markers.png", width = 12, height = 6, dpi = 300)