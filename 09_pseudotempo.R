
library(Seurat)
library(monocle3)
library(SeuratWrappers)

subset_obj <- subset(seurat_obj, subset = group %in% 
                    c("control", "acute", "defervescence"))

# Subset por tipos celulares da linhagem B
b_cells <- subset(subset_obj, subset = celltype_0.5 %in% 
                    c("Naive B", "ActMem B", "Pre-Plasmablast", "Plasmablast"))

DefaultAssay(b_cells) <- "RNA"

# Converter para Monocle3
cds <- as.cell_data_set(b_cells)

# Pré-processamento e redução de dimensionalidade
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Se você quiser usar os clusters do Seurat, pode fazer:
# cds@clusters$UMAP$clusters <- subset_obj$celltype_0.5

# Clusterização (opcional, útil para learn_graph)
cds <- cluster_cells(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# Aprender o gráfico da trajetória
cds <- learn_graph(cds)

# Selecionar células naive como raiz
root_cells <- WhichCells(b_cells, expression = celltype_0.5 == "Naive B")
if(length(root_cells) == 0) stop("Nenhuma célula naive encontrada")
cds <- order_cells(cds, root_cells = root_cells)

# Visualizar trajetória colorida por pseudotempo
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)

# Visualizar expressão de um lncRNA específico (se existir)
if("LNC_DHF_01" %in% rownames(cds)) {
  plot_genes_in_pseudotime(cds["LNC_DHF_01", ], color_cells_by = "condition")
}

# Identificar lncRNAs que variam significativamente ao longo da trajetória
pr_test <- graph_test(cds, neighbor_graph = "principal_graph")
seus_lncRNAs <- # lista de lncRNA
sig_lncRNAs <- pr_test[pr_test$gene_short_name %in% seus_lncRNAs & pr_test$q_value < 0.05, ]