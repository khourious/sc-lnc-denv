# Script R para Anotação de PBMC com CellTypist e Visualização com CellxGene
# Autor: Assistente Manus
# Data: 2025

# ============================================================================
# INSTALAÇÃO DE PACOTES
# ============================================================================

# Função para instalar pacotes se não estiverem disponíveis
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Instalar pacotes necessários
cat("Instalando pacotes R necessários...\n")
install_if_missing("reticulate")
install_if_missing("Seurat")
install_if_missing("ggplot2")

# Configurar Python e instalar CellTypist
cat("Configurando ambiente Python...\n")
library(reticulate)

# Instalar CellTypist no Python (descomente se necessário)
# py_install("celltypist")

# ============================================================================
# CARREGAMENTO DE DADOS
# ============================================================================

cat("Carregando dados de PBMC...\n")

# OPÇÃO 1: Carregar seus próprios dados
# Descomente e modifique as linhas abaixo para usar seus dados reais:
# 
# # Para dados em formato 10X:
# library(Seurat)
# pbmc.data <- Read10X(data.dir = "path/to/your/10x/data/")
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC_Annotation")
#
# # Para dados em formato RDS:
# pbmc <- readRDS("path/to/your/pbmc_data.rds")
#
# # Para dados em formato CSV:
# expr_matrix <- read.csv("path/to/your/expression_matrix.csv", row.names = 1)
# pbmc <- CreateSeuratObject(counts = expr_matrix, project = "PBMC_Annotation")

# OPÇÃO 2: Dados de exemplo (para demonstração)
# Remova esta seção quando usar seus dados reais
set.seed(123)
n_genes <- 2000
n_cells <- 500

# Simular genes comuns em PBMC
gene_names <- c(
  paste0("CD", sample(1:300, 50)),  # Marcadores CD
  paste0("IL", sample(1:50, 30)),   # Interleucinas
  paste0("TNF", sample(1:20, 10)),  # TNF
  paste0("IFNG", sample(1:10, 5)),  # Interferon gamma
  paste0("Gene_", 1:(n_genes-95))   # Genes genéricos
)

expr_matrix <- matrix(
  rpois(n_genes * n_cells, lambda = 3), 
  nrow = n_genes, 
  ncol = n_cells
)
rownames(expr_matrix) <- gene_names[1:n_genes]
colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)

# Criar objeto Seurat
library(Seurat)
pbmc <- CreateSeuratObject(
  counts = expr_matrix, 
  project = "PBMC_Annotation_Demo",
  min.cells = 3,
  min.features = 200
)

cat("Dados carregados com sucesso!\n")
cat("Número de células:", ncol(pbmc), "\n")
cat("Número de genes:", nrow(pbmc), "\n")

# ============================================================================
# PRÉ-PROCESSAMENTO DOS DATAS
# ============================================================================

cat("Iniciando pré-processamento dos dados...\n")

# Calcular métricas de qualidade
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filtrar células de baixa qualidade (ajuste conforme necessário)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Normalização
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identificar genes altamente variáveis
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Escalonamento dos dados
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Análise de componentes principais (PCA)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP para visualização
pbmc <- RunUMAP(pbmc, dims = 1:10)

cat("Pré-processamento concluído!\n")

# ============================================================================
# ANOTAÇÃO COM CELLTYPIST
# ============================================================================

cat("Iniciando anotação com CellTypist...\n")

# Importar CellTypist
celltypist <- import("celltypist")

# Preparar dados para CellTypist
# CellTypist espera uma matriz (células x genes) com dados normalizados
expression_data <- as.matrix(t(GetAssayData(pbmc, slot = "data")))

cat("Carregando modelo CellTypist...\n")

# Carregar modelo pré-treinado
# Modelos disponíveis: https://www.celltypist.org/models
# Para PBMC, recomenda-se: "Immune_All_Low.pkl" ou "Immune_All_High.pkl"
tryCatch({
  model <- celltypist$models$Model$load(model_name = "Immune_All_Low.pkl")
  cat("Modelo \'Immune_All_Low.pkl\' carregado com sucesso!\n")
}, error = function(e) {
  cat("Erro ao carregar modelo. Tentando baixar...\n")
  # Baixar modelo se não estiver disponível
  celltypist$models$download_models(model = "Immune_All_Low.pkl")
  model <- celltypist$models$Model$load(model_name = "Immune_All_Low.pkl")
})

# Realizar predição
cat("Realizando predições...\n")
predictions <- celltypist$annotate(
  X = expression_data, 
  model = model, 
  majority_voting = TRUE
)

# Adicionar predições ao objeto Seurat
pbmc$celltypist_predicted_labels <- predictions$predicted_labels$majority_voting
pbmc$celltypist_conf_score <- predictions$predicted_labels$conf_score

cat("Anotação com CellTypist concluída!\n")

# ============================================================================
# VISUALIZAÇÃO DOS RESULTADOS
# ============================================================================

cat("Gerando visualizações...\n")

# Visualizar clusters originais
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters Seurat")

# Visualizar anotações CellTypist
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "celltypist_predicted_labels", label = TRUE) +
  ggtitle("Anotações CellTypist") +
  theme(legend.position = "bottom")

# Salvar plots
ggsave("seurat_clusters.png", p1, width = 10, height = 8)
ggsave("celltypist_annotations.png", p2, width = 12, height = 10)

# Mostrar plots
print(p1)
print(p2)

# ============================================================================
# EXPORTAR DADOS PARA CELLXGENE
# ============================================================================

cat("Preparando dados para CellxGene...\n")

# Salvar objeto Seurat
saveRDS(pbmc, file = "pbmc_annotated.rds")

# Para usar com CellxGene, você pode exportar para formato h5ad
# Instale SeuratDisk se necessário:
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")

tryCatch({
  library(SeuratDisk)
  
  # Salvar como h5seurat
  SaveH5Seurat(pbmc, filename = "pbmc_annotated.h5seurat", overwrite = TRUE)
  
  # Converter para h5ad (formato do scanpy/cellxgene)
  Convert("pbmc_annotated.h5seurat", dest = "h5ad", overwrite = TRUE)
  
  cat("Dados exportados para \'pbmc_annotated.h5ad\'\n")
  cat("Para visualizar no CellxGene, execute:\n")
  cat("cellxgene launch pbmc_annotated.h5ad\n")
  
}, error = function(e) {
  cat("SeuratDisk não está disponível. Dados salvos apenas como RDS.\n")
  cat("Para instalar SeuratDisk:\n")
  cat("install.packages(\'remotes\')\n")
  cat("remotes::install_github(\'mojaveazure/seurat-disk\')\n")
})

# ============================================================================
# RELATÓRIO FINAL
# ============================================================================

cat("\n" %+% "="*60 %+% "\n")
cat("RELATÓRIO DE ANOTAÇÃO CONCLUÍDO\n")
cat("="*60 %+% "\n")
cat("Número total de células:", ncol(pbmc), "\n")
cat("Número de genes:", nrow(pbmc), "\n")
cat("Tipos de células identificados:\n")

# Contar tipos de células
cell_types <- table(pbmc$celltypist_predicted_labels)
for(i in 1:length(cell_types)) {
  cat("  -", names(cell_types)[i], ":", cell_types[i], "células\n")
}

cat("\nArquivos gerados:\n")
cat("  - pbmc_annotated.rds (objeto Seurat)\n")
cat("  - seurat_clusters.png (visualização clusters)\n")
cat("  - celltypist_annotations.png (visualização anotações)\n")

if(file.exists("pbmc_annotated.h5ad")) {
  cat("  - pbmc_annotated.h5ad (para CellxGene)\n")
}

cat("\nPara visualizar no CellxGene:\n")
cat("1. Instale CellxGene: pip install cellxgene\n")
cat("2. Execute: cellxgene launch pbmc_annotated.h5ad\n")
cat("3. Abra o navegador no endereço mostrado\n")

cat("\n" %+% "="*60 %+% "\n")

