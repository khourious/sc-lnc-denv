# Instalar pacotes se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")

# Carregar pacotes
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

# --- lista HFpEF ---
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

# --- Lista final ---
seurat_list <- list()

# -------------------------------
# 2. Processar cada arquivo .h5
# -------------------------------

for (nome_arquivo in names(arquivos_h5)) {
  caminho <- arquivos_h5[[nome_arquivo]]
  amostras <- amostras_por_arquivo[[nome_arquivo]]
  
  # Lê o arquivo
  counts <- Read10X_h5(caminho)
  
  total_cells <- ncol(counts)
  n_amostras <- length(amostras)
  cells_por_amostra <- floor(total_cells / n_amostras)
  
  for (i in seq_along(amostras)) {
    start <- (i - 1) * cells_por_amostra + 1
    end <- if (i == n_amostras) total_cells else i * cells_por_amostra
    
    # Seleciona subconjunto de células
    barcodes <- colnames(counts)[start:end]
    sub_counts <- counts[, barcodes]
    
    # Renomear as células com prefixo da amostra
    colnames(sub_counts) <- paste0(amostras[i], "_", barcodes)
    
    # Criar objeto Seurat
    seu <- CreateSeuratObject(counts = sub_counts, project = amostras[i])
    seu$amostra <- amostras[i]
    
    seurat_list[[amostras[i]]] <- seu
  }
}

if (!dir.exists("qc_doublets")) {
  dir.create("qc_doublets")
}

for (amostra in names(seurat_list)) {
  cat("Processando:", amostra, "\n")
  
  seu <- seurat_list[[amostra]]
  
  # Pré-processamento
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:30)
  
  # Converter para SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)
  
  # Rodar scDblFinder
  sce <- scDblFinder(sce)
  
  # Adicionar resultados ao Seurat
  seu$doublet_class <- sce$scDblFinder.class
  seu$doublet_score <- sce$scDblFinder.score
  
  # Scatterplot UMAP: singlet vs doublet
  p <- DimPlot(seu, group.by = "doublet_class", reduction = "umap", pt.size = 0.5) +
    ggtitle(paste("Doublets -", amostra)) +
    theme_minimal()
  
  # Salvar gráfico
  ggsave(
    filename = paste0("qc_doublets/", amostra, "_doublets_umap.png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  # Filtrar singlets e doublets
  singlets <- subset(seu, subset = doublet_class == "singlet")
  doublets <- subset(seu, subset = doublet_class == "doublet")
  
  # Salvar objetos, mesmo que estejam vazios
  saveRDS(singlets, file = paste0("qc_doublets/", amostra, "_singlets.rds"))
  saveRDS(doublets, file = paste0("qc_doublets/", amostra, "_doublets.rds"))
}
