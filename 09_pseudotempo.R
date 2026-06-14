if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('Seurat')
install.packages("remotes")
install.packages("R.utils")
devtools::install_github("satijalab/seurat-wrappers")
library(Seurat)
library(monocle3)
library(SeuratWrappers)

setwd("~/denv")

# --- importar RDS ---
seurat_obj <- readRDS("RDS/sct_harmony_annotation.rds")

# --- Selecionar subgrupo ---
subset_obj <- subset(seurat_obj, subset = group %in% 
                       c("control", "acute", "defervescence"))

# --- Subset por tipos celulares da linhagem B ---
b_cells <- subset(subset_obj, subset = celltype_0.5 %in% 
                    c("Naive B", "ActMem B", "Pre-Plasmablast", "Plasmablast"))

DefaultAssay(b_cells) <- "RNA"

# --- Converter para Monocle3 ---
cds <- as.cell_data_set(b_cells)

# --- Pré-processamento e redução de dimensionalidade ---
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# --- Se quiser usar os clusters do Seurat, pode fazer ---
# cds$celltype <- b_cells$celltype_0.5

# --- Clusterização ---
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# --- Definir raiz ---
root_cells <- WhichCells(b_cells, expression = celltype_0.5 == "Naive B")
cds <- order_cells(cds, root_cells = root_cells)

# --- Selecionar células naive como raiz ---
root_cells <- WhichCells(b_cells, expression = celltype_0.5 == "Naive B")
if(length(root_cells) == 0) stop("Nenhuma célula naive encontrada")
cds <- order_cells(cds, root_cells = root_cells)

# --- Organizar rownome em outra coluna ---
rowData(cds)$gene_short_name <- rownames(cds)

# --- Visualizar trajetória colorida por pseudotempo ---
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)

plot_cells(cds, color_cells_by = "dengue_classification", label_cell_groups = FALSE) +
  facet_wrap(~condition)  # ou usar uma única figura com cores

# --- Identificar mRNAs que variam significativamente ao longo da trajetória ---
pr_test <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# --- Identificar TODOS os genes significativos ---
pr_test <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
sig_genes <- subset(pr_test, morans_I > 0.1 & q_value < 0.05 & status == "OK")

# 2. Classificar os genes como mRNA ou lncRNA (usando anotação do seu objeto Seurat)
#    Supondo que você tenha um vetor `lncRNA_symbols` com os símbolos dos lncRNAs
#    (extraído da sua anotação original). Se não tiver, crie a partir do GTF.
#    Por enquanto, vou assumir que você tem uma lista `is_lncRNA` no rowData do cds.
#    Se não tiver, você pode adicionar:
#    rowData(cds)$is_lncRNA <- rownames(cds) %in% seus_lncRNAs_conhecidos

# --- Matriz de expressão dos genes significativos, ordenada por pseudotempo ---
sig_ids <- rownames(sig_genes)
expr_sig <- as.matrix(exprs(cds)[sig_ids, ])

# --- Ordenar células por pseudotempo ---
cell_order <- order(pseudotime(cds))
expr_ordered <- expr_sig[, cell_order]

# --- Clustering não supervisionado dos perfis --- 
set.seed(123)
kmeans_clust <- kmeans(expr_ordered, centers = 3)
cluster_assign <- kmeans_clust$cluster

# --- Ordenar genes por cluster para o heatmap --- 
gene_order <- order(cluster_assign)
expr_ordered <- expr_ordered[gene_order, ]

# --- Anotação para o heatmap --- 
annotation_col <- data.frame(
  condition = colData(cds)$dengue_classification[cell_order],
  pseudotime = pseudotime(cds)[cell_order]
)
rownames(annotation_col) <- colnames(expr_ordered)

# --- Adicionar anotação de tipo de gene (mRNA vs lncRNA) nas linhas --- 
gene_type <- ifelse(rownames(expr_ordered) %in% lncRNA_symbols, "lncRNA", "mRNA")
annotation_row <- data.frame(type = gene_type)
rownames(annotation_row) <- rownames(expr_ordered)

# --- Heatmap final --- 
pheatmap(expr_ordered, 
         scale = "row", 
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         gaps_row = cumsum(table(cluster_assign)[order(unique(cluster_assign))]),
         main = "Temporal clusters of genes (mRNAs and lncRNAs)")

# --- distribuição de lncRNAs por cluster e condição --- 
cluster_assignment_df <- data.frame(
  gene = rownames(expr_ordered),
  cluster = cluster_assign[gene_order],
  type = gene_type
)

# --- quantos lncRNAs em cada cluster? E mRNAs? --- 
table(cluster_assignment_df$cluster, cluster_assignment_df$type)

# --- Para entender a dinâmica da expressão em DF vs DHF, podemos calcular o perfil médio
#    de expressão de cada cluster separadamente por condição. Primeiro, precisamos dos
#    valores normalizados (log-counts) e agrupar por condição e percentil de pseudotempo.

# --- Dividir o pseudotempo em 10 bins --- 
pseudotime_bins <- cut(pseudotime(cds), breaks = 10, labels = 1:10)

# Para cada condição (DF, DHF, control), calcular a média de expressão de cada cluster
# dentro de cada bin. Isso pode ser feito com um loop ou com dplyr.

library(dplyr)

# Criar um dataframe longo: célula, pseudotime_bin, condition, e a expressão de cada gene
# Isso é mais complexo; vou simplificar: calcular a média por cluster e condição ao longo do pseudotime.

# Primeiro, calcular a expressão média de cada cluster (média dos genes do cluster)
# para cada célula.
expr_cluster <- matrix(NA, nrow = length(unique(cluster_assign)), ncol = ncol(expr_ordered))
for (cl in unique(cluster_assign)) {
  genes_cl <- rownames(expr_ordered)[cluster_assign[gene_order] == cl]
  if (length(genes_cl) > 0) {
    expr_cluster[cl, ] <- colMeans(expr_ordered[genes_cl, , drop = FALSE])
  }
}
rownames(expr_cluster) <- paste0("Cluster", unique(cluster_assign))

# Agora, associar cada célula à sua condição e pseudotime bin
cell_metadata <- data.frame(
  cell = colnames(expr_ordered),
  condition = colData(cds)$dengue_classification[cell_order],
  pseudotime_bin = pseudotime_bins[cell_order]
)

# Para cada cluster, calcular a média por condição e bin
# (Aqui usaremos um reshape simples)
library(reshape2)
expr_cluster_df <- melt(expr_cluster, varnames = c("cluster", "cell"), value.name = "expression")
expr_cluster_df <- merge(expr_cluster_df, cell_metadata, by = "cell")

# Agora, agregar a média por cluster, condition e pseudotime_bin
cluster_summary <- expr_cluster_df %>%
  group_by(cluster, condition, pseudotime_bin) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

# Plotar as curvas (exemplo usando ggplot2)
ggplot(cluster_summary, aes(x = as.numeric(pseudotime_bin), y = mean_expr, color = condition)) +
  geom_line() + geom_point() +
  facet_wrap(~cluster, scales = "free_y") +
  labs(title = "Temporal dynamics of gene clusters by condition",
       x = "Pseudotime bin", y = "Mean expression")
