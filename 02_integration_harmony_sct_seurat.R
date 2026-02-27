BiocManager::install("scDblFinder", force = T)
install.packages("ggplot2")
install.packages("dplyr")
install.packages("cluster")
install.packages("harmony")
install.packages("devtools")
install.packages("hdf5r")
install.packages("SeuratDisk")
BiocManager::install("glmGamPoi")

devtools::install_github("immunogenomics/LISI")

# Carregar pacotes
library(igraph)
library(lisi)
library(scDblFinder)
library(Seurat)
library(hdf5r)
library(harmony)
library(ggplot2)
library(dplyr)
library(cluster)
library(patchwork)
library(glmGamPoi)


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
  #dt1_control = "Healthy_Control_run1",
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
  #"Healthy_Control_run1"= 10,
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
processar_amostras <- function(arquivos_h5, amostras_por_arquivo, limiares_mt, metadata_samples, metadata_denv) {
  objs <- list()
  
  for (nome_arquivo in names(arquivos_h5)) {
    counts <- Read10X_h5(arquivos_h5[[nome_arquivo]])
    total <- ncol(counts)
    amostras <- amostras_por_arquivo[[nome_arquivo]]
    n <- length(amostras)
    por_amostra <- floor(total / n)
    
    # 🔗 Metadados do arquivo H5
    meta_denv <- metadata_denv[metadata_denv$sample_id == nome_arquivo, ]
    
    for (i in seq_along(amostras)) {
      ini <- (i - 1) * por_amostra + 1
      fim <- if (i == n) total else i * por_amostra
      barcodes <- colnames(counts)[ini:fim]
      sub <- counts[, barcodes]
      colnames(sub) <- paste0(amostras[i], "_", barcodes)
      
      seu <- CreateSeuratObject(sub, project = amostras[i])
      seu$orig.ident <- amostras[i]
      
      # 🔗 Metadados da amostra individual
      meta_sample <- metadata_samples[metadata_samples$sample_id == amostras[i], ]
      
      if (nrow(meta_sample) == 1) {
        for (col in colnames(meta_sample)) {
          if (col != "sample_id") {
            seu[[col]] <- meta_sample[[col]]
          }
        }
      } else {
        warning("⚠️ Metadados individuais não encontrados para ", amostras[i])
      }
      
      # 🔗 Metadados do arquivo H5 (aplicado a todas as células da amostra)
      if (nrow(meta_denv) == 1) {
        for (col in colnames(meta_denv)) {
          if (col != "sample_id") {
            seu[[col]] <- meta_denv[[col]]
          }
        }
      } else {
        warning("⚠️ Metadados DENV não encontrados para arquivo ", nome_arquivo)
      }
      
      message("📦 ", amostras[i], ": ", ncol(seu), " células iniciais")
      
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
      
      if (amostras[i] %in% names(limiares_mt)) {
        limiar <- limiares_mt[[amostras[i]]]
        seu <- subset(seu, subset = percent.mt < limiar)
        message("🧹 ", amostras[i], ": filtrado com percent.mt <", limiar, " → ", ncol(seu), " células restantes")
      }
      
      if (ncol(seu) == 0) {
        message("⚠️ ", amostras[i], ": nenhuma célula após filtro mitocondrial. Ignorado.")
        next
      }
      
      objs[[amostras[i]]] <- seu
    }
  }
  
  return(objs)
}


names(arquivos_h5)
names(amostras_por_arquivo)

# --- Leitura e pré-processamento ---
seurat_list <- processar_amostras(
  arquivos_h5 = arquivos_h5,
  amostras_por_arquivo = amostras_por_arquivo,
  limiares_mt = limiares_mt,
  metadata_samples = metadata_samples,
  metadata_denv = metadata_denv
)



####### - PCA + Harmony
# --- Normalização padrão ---
seurat_list <- lapply(seurat_list, NormalizeData)
seurat_list <- lapply(seurat_list, FindVariableFeatures)
seurat_list <- lapply(seurat_list, ScaleData)
seurat_list <- lapply(seurat_list, RunPCA)

# --- Merge dos objetos para integração ---
harmony_merged_pca <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
harmony_merged_pca <- ScaleData(harmony_merged_pca, vars.to.regress = "percent.mt")
harmony_merged_pca <- RunPCA(harmony_merged_pca)

# --- Harmony ----

harmony_merged_pca <- RunHarmony(harmony_merged_pca,
                                 "orig.ident",
                                 "pca",
                                 dims.use = 1:30)



# --- UMAP e clustering usando Harmony ----
harmony_merged_pca <- RunUMAP(harmony_merged_pca, reduction = "harmony", dims = 1:30, n.neighbors = 50)
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
saveRDS(harmony_merged_pca, file = "harmony_merged_pca_novo.rds")
write.csv(markers, file = "harmony_pca_markers_novo.csv", row.names = FALSE)

metadata_denv
############
# SCT + seurat 

# --- Normalização com SCTransform por amostra ---
# --- SCTransform por amostra ---
seurat_list_sct <- lapply(seurat_list, function(obj) {
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = TRUE)
  
  # Split do assay RNA por orig.ident
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
  
  return(obj)
})

# --- Merge dos objetos SCT ---
merged_sct <- merge(
  x = seurat_list_sct[[1]],
  y = seurat_list_sct[-1],
  add.cell.ids = names(seurat_list_sct)
)

# --- Definir assay ativo como SCT ---
DefaultAssay(merged_sct) <- "SCT"
merged_sct  <- NormalizeData(merged_sct)
merged_sct  <- FindVariableFeatures(merged_sct)
merged_sct  <- ScaleData(merged_sct)
merged_sct  <- RunPCA(merged_sct)

# --- Vizinhos e clusters ---
merged_sct  <- FindNeighbors(merged_sct, dims = 1:30, reduction = "pca")
merged_sct  <- FindClusters(merged_sct, resolution = 2, cluster.name = "unintegrated_clusters")

# --- integração ---
merged_sct_seurat <- IntegrateData(anchors, normalization.method = "SCT")

# --- PCA + UMAP ---
merged_sct_seurat <- RunPCA(merged_sct_seurat)
merged_sct_seurat <- RunUMAP(merged_sct_seurat, dims = 1:30)
merged_sct_seurat <- FindNeighbors(merged_sct_seurat, dims = 1:30)


# --- Visualizar UMAP com clusters
DimPlot(merged_sct_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Identificar marcadores de cluster
markers_sct_seurat <- FindAllMarkers(merged_sct_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# --- Salvar objeto final
saveRDS(merged_sct_seurat, file = "merged_sct_seurat.rds")
write.csv(markers_sct_seurat, file = "seurat_sct_markers.csv", row.names = FALSE)


############
######### sct + harmony

# --- Harmony ---

merged_sct <- IntegrateLayers(
  object = merged_sct, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT",
  verbose = FALSE
)

merged_sct <- FindNeighbors(merged_sct, reduction = "harmony", dims = 1:30)
merged_sct <- FindClusters(merged_sct, resolution = 2, cluster.name = "harmony_clusters2")

merged_sct <- RunUMAP(merged_sct, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(
  merged_sct,
  reduction = "umap.harmony",
  group.by = c("orig.ident")
)

merged_sct_harmony <- merged_sct
head(merged_sct_harmony@meta.data)
# --- visulização não integrada ---
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "orig.ident")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "group")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "disease")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "timepoint")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "virus")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "age")
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "sex") 
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "infection") 
DimPlot(merged_sct_harmony, reduction = "umap.harmony", group.by = "dataset") 

# --- integração e UMAP ---
integrated_sct_harmony <- RunHarmony(merged_sct_harmony, group.by.vars = "dataset")
integrated_sct_harmony <- RunUMAP(integrated_sct_harmony, reduction = "harmony", dims = 1:30)

# -- vizinhos e clusters ---
integrated_sct_harmony <- FindNeighbors(integrated_sct_harmony, reduction = "harmony", dims = 1:30)
integrated_sct_harmony <- FindClusters(integrated_sct_harmony, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1))

head(integrated_sct_harmony@meta.data)
# --- Visualização Integrada ---
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "group")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "disease")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "timepoint")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "virus")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "age")
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "sex") 
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "infection") 
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "dataset") 
DimPlot(integrated_sct_harmony, reduction = "umap", group.by = "SCT_snn_res.0.2", label = TRUE) + ggtitle("Seurat Clusters - resolution = 0.2")

# --- Identificar marcadores de cluster
integrated_sct_harmony <- JoinLayers(integrated_sct_harmony)
integrated_sct_harmony <- PrepSCTFindMarkers(integrated_sct_harmony)
markers <- FindAllMarkers(integrated_sct_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# --- Visualização ---
DimPlot(harmony_merged_pca, group.by = "orig.ident") + ggtitle("PCA + Harmony")
DimPlot(merged_sct_seurat, group.by = "orig.ident") + ggtitle("SCT + Seurat")
DimPlot(integrated_sct_harmony, group.by = "orig.ident") + ggtitle("SCT + Harmony")
write.csv(markers , file = "sct_harmony_markers.csv", row.names = FALSE)
saveRDS(integrated_sct_harmony, file = "sct_harmony_merged_novo.rds")



library(ggplot2)
library(dplyr)

# Agrupar e contar células por timepoint
df <- merged_sct_harmony@meta.data %>%
  dplyr::count(timepoint)

# Converter timepoint para fator ordenado (se necessário)
df$timepoint <- factor(df$timepoint, levels = sort(unique(df$timepoint)))

# Criar gráfico de linha
ggplot(df, aes(x = timepoint, y = n, group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 2) +
  theme_minimal() +
  labs(
    title = "Quantidade de células ao longo do tempo",
    x = "Timepoint clínico",
    y = "Número de células"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Amostras com tempo relativo à defervescência
amostras_defervescente <- unlist(amostras_por_arquivo[c("dt1_control", "dt1_DF", "dt1_DHF", "dt2_DF_1", "dt2_DF_2", "dt3_primary", "dt3_secundary")])

# Amostras com tempo absoluto de febre
amostras_febre_absoluta <- unlist(amostras_por_arquivo[c("dt4_control", "dt4_DF", "dt4_DWS", "dt4_SD")])

library(ggplot2)
library(dplyr)

# Agrupar e contar células por timepoint
df <- merged_sct_harmony@meta.data %>%
  dplyr::count(timepoint)

# Converter timepoint para fator ordenado (se necessário)
df$timepoint <- factor(df$timepoint, levels = sort(unique(df$timepoint)))

# Criar gráfico de linha
ggplot(df, aes(x = timepoint, y = n, group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 2) +
  theme_minimal() +
  labs(
    title = "Quantidade de células ao longo do tempo",
    x = "Timepoint clínico",
    y = "Número de células"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Amostras com tempo relativo à defervescência
amostras_defervescente <- unlist(amostras_por_arquivo[c("dt1_control", "dt1_DF", "dt1_DHF", "dt2_DF_1", "dt2_DF_2", "dt3_primary", "dt3_secundary")])

# Amostras com tempo absoluto de febre
amostras_febre_absoluta <- unlist(amostras_por_arquivo[c("dt4_control", "dt4_DF", "dt4_DWS", "dt4_SD")])

# Extrair metadados
meta <-  merged_sct_harmony@meta.data

# Defervescente
df_def <- meta %>%
  filter(orig.ident %in% amostras_defervescente) %>%
  count(timepoint) %>%
  mutate(timepoint = factor(timepoint, levels = c("control","-5","-4", "-3", "-2", "-1", "0", 
                                                  "1", "2", "3", "5", "7", "14", "180"))) %>%
  arrange(timepoint)


# Febre absoluta
df_febre <- meta %>%
  filter(orig.ident %in% amostras_febre_absoluta) %>%
  count(timepoint) %>%
  mutate(timepoint = factor(timepoint, levels = c("control", "1", "2", "3", "4", "5", "6", "7", "T"))) %>%
  arrange(timepoint)


# Gráfico 1: tempo relativo à defervescência
p1 <- ggplot(df_def, aes(x = timepoint, y = n ,group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 2) +
  theme_minimal() +
  labs(title = "Células ao longo do tempo (defervescência)", x = "Dia relativo", y = "Número de células")

# Gráfico 2: tempo absoluto de febre
p2 <- ggplot(df_febre, aes(x = timepoint, y = n, group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 2) +
  theme_minimal() +
  labs(title = "Células ao longo do tempo (dias de febre)", x = "Dia de febre", y = "Número de células")

p1 + p2
