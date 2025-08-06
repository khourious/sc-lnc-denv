BiocManager::install("scDblFinder", force = T)
install.packages("ggplot2")
install.packages("dplyr")
install.packages("cluster")
install.packages("igraph")
install.packages("devtools")
devtools::install_github("immunogenomics/LISI")

# Carregar pacotes
library(igraph)
library(lisi)
library(scDblFinder)
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(cluster)

setwd("~/denv")

# --- lista DENV ---
arquivos_h5 <- list(
  dt1_control = "dataset1/denv_control/filtered_feature_bc_matrix.h5",
  dt1_DF = "dataset1/denv_DF/filtered_feature_bc_matrix.h5",
  dt1_DHF = "dataset1/denv_DHF/filtered_feature_bc_matrix.h5",
  dt2_DF_1  = "dataset2/denv_nat01/filtered_feature_bc_matrix.h5",
  dt2_DF_2 = "dataset2/denv_nat02/filtered_feature_bc_matrix.h5",
  dt3_primary = "dataset3/denv_primary/filtered_feature_bc_matrix.h5",
  dt3_secundary = "dataset3/denv_secundary/filtered_feature_bc_matrix.h5",
  dt4_control = "dataset4/denv_control/filtered_feature_bc_matrix.h5",
  dt4_DF = "dataset4/denv_DF/filtered_feature_bc_matrix.h5",
  dt4_DWS = "dataset4/denv_DWS/filtered_feature_bc_matrix.h5",
  dt4_SD = "dataset4/denv_SD/filtered_feature_bc_matrix.h5"
)

# --- barcode DENV ---
amostras_por_arquivo <- list(
  dt1_control = "Healthy_Control_run1",
  dt1_DF = c("DF_Day_minus_1_run1","DF_Day_minus_1_run2","DF_Day_minus_2_run1","DF_Def_run1",
             "DF_Def_run2","DF_Wk2_run1"),
  dt1_DHF = c("DHF_Day_minus_1_run1","DHF_Day_minus_1_run2","DHF_Day_minus_2_run1",
              "DHF_Def_run1","DHF_Def_run2","DHF_Wk2_run1"),
  dt2_DF_1  = c("SRR12215051","SRR12215052","SRR12215053"),
  dt2_DF_2 = c("SRR12215054","SRR12215055","SRR12215056"),
  dt3_primary = c("SRR11088622_Primary1_D1","SRR11088623_Primary1_D3","SRR11088624_Primary2_D1",
                  "SRR11088625_Primary2_D5","SRR11088626_Primary3_D1","SRR11088627_Primary3_D2"),
  dt3_secundary = c("SRR11088628_Secondary1_D0","SRR11088629_Secondary1_D1","SRR11088630_Secondary1_D7",
                    "SRR11088631_Secondary2_D1","SRR11088632_Secondary2_D1","SRR11088633_Secondary2_D5",
                    "SRR11088634_Secondary3_D0","SRR11088635_Secondary3_D1","SRR11088636_Secondary3_D5"),
  dt4_control = c("SRR22739533","SRR22739534","SRR22739543","SRR22739544","SRR22739551","SRR22739552"),
  dt4_DF = c("SRR22739530","SRR22739535","SRR22739536","SRR22739537","SRR22739539","SRR22739541",
             "SRR22739545","SRR22739547","SRR22739549","SRR22739550","SRR22739555"),
  dt4_DWS = c("SRR22739525","SRR22739527","SRR22739529","SRR22739538"),
  dt4_SD = c("SRR22739526","SRR22739528","SRR22739531","SRR22739532","SRR22739540","SRR22739542",
             "SRR22739546","SRR22739548","SRR22739553","SRR22739554")
)

# --- limiar mt
limiares_mt <- list(
  "Healthy_Control_run1"= 10,
  "DF_Day_minus_1_run1"= 10,
  "DF_Day_minus_1_run2"= 10,
  "DF_Day_minus_2_run1"= 10,
  "DF_Def_run1"= 10,
  "DF_Def_run2"= 10,
  "DF_Wk2_run1"= 10,
  "DHF_Day_minus_1_run1"= 10,
  "DHF_Day_minus_1_run2"= 10,
  "DHF_Day_minus_2_run1"= 10,
  "DHF_Def_run1"= 10,
  "DHF_Def_run2"= 10,
  "DHF_Wk2_run1"= 10,
  "SRR12215051"= 10,
  "SRR12215052"= 10,
  "SRR12215053"= 10,
  "SRR12215054"= 10,
  "SRR12215055"= 10,
  "SRR12215056"= 10,
  "SRR11088622_Primary1_D1"= 10,
  "SRR11088623_Primary1_D3"= 10,
  "SRR11088624_Primary2_D1"= 10,
  "SRR11088625_Primary2_D5"= 10,
  "SRR11088626_Primary3_D1"= 10,
  "SRR11088627_Primary3_D2"= 10,
  "SRR11088628_Secondary1_D0"= 10,
  "SRR11088629_Secondary1_D1"= 10,
  "SRR11088630_Secondary1_D7"= 10,
  "SRR11088631_Secondary2_D1"= 10,
  "SRR11088632_Secondary2_D1"= 10,
  "SRR11088633_Secondary2_D5"= 10,
  "SRR11088634_Secondary3_D0"= 10,
  "SRR11088635_Secondary3_D1"= 10,
  "SRR11088636_Secondary3_D5"= 10,
  "SRR22739533"= 10,
  "SRR22739534"= 10,
  "SRR22739543"= 10,
  "SRR22739544"= 10,
  "SRR22739551"= 10,
  "SRR22739552"= 10,
  "SRR22739530"= 10,
  "SRR22739535"= 10,
  "SRR22739536"= 10,
  "SRR22739537"= 10,
  "SRR22739539"= 10,
  "SRR22739541"= 10,
  "SRR22739545"= 10,
  "SRR22739547"= 10,
  "SRR22739549"= 10,
  "SRR22739550"= 10,
  "SRR22739555"= 10,
  "SRR22739525"= 10,
  "SRR22739527"= 10,
  "SRR22739529"= 10,
  "SRR22739538"= 10,
  "SRR22739526"= 10,
  "SRR22739528"= 10,
  "SRR22739531"= 10,
  "SRR22739532"= 10,
  "SRR22739540"= 10,
  "SRR22739542"= 10,
  "SRR22739546"= 10,
  "SRR22739548"= 10,
  "SRR22739553"= 10,
  "SRR22739554"= 10
  )
## --- Função para processar amostras ---
processar_amostras <- function(arquivos_h5, amostras_por_arquivo, limiares_mt) {
  objs <- list()
  for (nome in names(arquivos_h5)) {
    counts <- Read10X_h5(arquivos_h5[[nome]])
    total <- ncol(counts)
    amostras <- amostras_por_arquivo[[nome]]
    n <- length(amostras)
    por_amostra <- floor(total / n)
    
    for (i in seq_along(amostras)) {
      ini <- (i - 1) * por_amostra + 1
      fim <- if (i == n) total else i * por_amostra
      barcodes <- colnames(counts)[ini:fim]
      sub <- counts[, barcodes]
      colnames(sub) <- paste0(amostras[i], "_", barcodes)
      
      seu <- CreateSeuratObject(sub, project = amostras[i])
      seu$orig.ident <- amostras[i]
      
      # Calcular percent.mt
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
      
      # Filtrar por limiar
      if (amostras[i] %in% names(limiares_mt)) {
        limiar <- limiares_mt[[amostras[i]]]
        seu <- subset(seu, subset = percent.mt < limiar)
        message("🧹 ", amostras[i], ": filtrado com percent.mt <", limiar)
      }
      
      # Identificar doublets
      sce <- as.SingleCellExperiment(seu)
      sce <- scDblFinder(sce)
      seu <- seu[, sce$scDblFinder.class == "singlet"]
      
      objs[[amostras[i]]] <- seu
    }
  }
  return(objs)
}

# --- Leitura e pré-processamento ---
seurat_list <- processar_amostras(arquivos_h5, amostras_por_arquivo, limiares_mt)


####### - PCA + Harmony
# --- Normalização padrão ---
harmony_merged_pca <- lapply(harmony_merged_pca, NormalizeData)
harmony_merged_pca <- lapply(harmony_merged_pca, FindVariableFeatures)

# --- Merge dos objetos para integração ---
harmony_merged_pca <- merge(harmony_merged_pca[[1]], y = harmony_merged_pca[-1], add.cell.ids = names(harmony_merged_pca))
harmony_merged_pca <- ScaleData(harmony_merged_pca, vars.to.regress = "percent.mt")
harmony_merged_pca <- RunUMAP(harmony_merged_pca, reduction = "harmony", dims = 1:30, n.neighbors = 50)


# --- Harmony ----
harmony_merged_pca <- RunHarmony(harmony_merged_pca, group.by.vars = "orig.ident")

# --- clustering com Harmony embeddings
harmony_merged_pca <- FindNeighbors(harmony_merged_pca, reduction = "harmony", dims = 1:30)
harmony_merged_pca <- FindClusters(harmony_merged_pca, resolution = 0.5)

harmony_merged_pca@meta.data
# --- Visualizar UMAP com clusters
DimPlot(harmony_merged_pca, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(harmony_merged_pca, reduction = "umap", group.by = "orig.ident", label = FALSE)

# --- Identificar marcadores de cluster
harmony_merged_pca <- JoinLayers(harmony_merged_pca)
markers <- FindAllMarkers(harmony_merged_pca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# --- Salvar objeto final
saveRDS(harmony_merged_pca, file = "harmony_merged_pca.rds")
write.csv(markers, file = "harmony_pca_markers.csv", row.names = FALSE)


############
# SCT + seurat 

# --- Normalização com SCTransform por amostra ---
seurat_list_sct <- lapply(seurat_list, function(obj) {
  SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
})

# --- Integração com anchors ---
features <- SelectIntegrationFeatures(seurat_list_sct, nfeatures = 3000)
seurat_list_sct <- PrepSCTIntegration(seurat_list_sct, anchor.features = features)

# --- Anchors ---
anchors <- FindIntegrationAnchors(seurat_list_sct, normalization.method = "SCT", anchor.features = features)

# --- integração ---
merged_sct_seurat <- IntegrateData(anchors, normalization.method = "SCT")

# --- PCA + UMAP ---
merged_sct_seurat<- RunPCA(merged_sct_seurat)
merged_sct_seurat <- RunUMAP(merged_sct_seurat, dims = 1:30)
merged_sct_seurat <- FindNeighbors(merged_sct_seurat, dims = 1:30)

# --- Visualizar UMAP com clusters
DimPlot(merged_sct_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(merged_sct_seurat, reduction = "umap", group.by = "orig.ident", label = FALSE)

# --- Identificar marcadores de cluster
markers_sct_seurat <- FindAllMarkers(merged_sct_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# --- Salvar objeto final
saveRDS(merged_sct_seurat, file = "sct_seurat_merged.rds")
write.csv(markers_sct_seurat, file = "seurat_sct_markers.csv", row.names = FALSE)

############
######### sct + harmony

# --- Merge após SCT ---
merged_sct_harmony <- merge(seurat_list_sct[[1]], y = seurat_list_sct[-1], add.cell.ids = names(seurat_list_sct))

# --- PCA e Harmony ---
merged_sct_harmony <- RunPCA(merged_sct_harmony)
merged_sct_harmony <- RunHarmony(merged_sct_harmony, group.by.vars = "orig.ident")
merged_sct_harmony <- RunUMAP(merged_sct_harmony, reduction = "harmony", dims = 1:30)

# -- vizinhos e clusters ---
merged_sct_harmony <- FindNeighbors(merged_sct_harmony, reduction = "harmony", dims = 1:30)
merged_sct_harmony <- FindClusters(merged_sct_harmony)

# --- Visualização ---
DimPlot(harmony_merged_pca group.by = "orig.ident") + ggtitle("PCA + Harmony")
DimPlot(merged_sct_seurat, group.by = "orig.ident") + ggtitle("SCT + Seurat")
DimPlot(merged_sct_harmony, group.by = "orig.ident") + ggtitle("SCT + Harmony")


# Silhouette Score

sil_pca <- silhouette(as.numeric(Idents(merged_pca)), dist(Embeddings(merged_pca, "harmony")))
sil_sct_seurat <- silhouette(as.numeric(Idents(merged_sct_seurat)), dist(Embeddings(merged_sct_seurat, "pca")))
sil_sct_harmony <- silhouette(as.numeric(Idents(merged_sct_harmony)), dist(Embeddings(merged_sct_harmony, "harmony")))

merged_sct_seurat <- FindClusters(merged_sct_seurat)

