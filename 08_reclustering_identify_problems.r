library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)

T_cells <- subset(seurat_obj, subset = cell_type %in% c("T_CD4", "T_CD8"))
T_cells <- subset(seurat_obj, idents = "T_cells")  # ou pelo marcador CD3+


T_cells <- NormalizeData(T_cells)
T_cells <- FindVariableFeatures(T_cells)
T_cells <- ScaleData(T_cells)
T_cells <- RunPCA(T_cells)
T_cells <- RunUMAP(T_cells, dims = 1:20)
T_cells <- FindNeighbors(T_cells, dims = 1:20)
T_cells <- FindClusters(T_cells, resolution = 0.5)  # ajuste resolução


# CD4, CD8A/CD8B → separar T helper vs T citotóxicos

# FOXP3 → Tregs

# Evitar clusters com IGH/IGK/IGL → provavelmente contaminantes B/plasmablastos

# suposições: `T_cells` é o Seurat object com as T já reclusterizadas
# e meta contém coluna "orig.ident" ou "sample" com o indivíduo.
# ajuste nomes conforme seu objeto.

# 1) imunoglobulinas
ig_genes <- c("IGHG1","IGHG2","IGHM","IGHA1","IGHA2","IGKC","IGLC2","IGLC3","IGLC1","IGK","IGL") 
ig_genes <- intersect(ig_genes, rownames(T_cells))

# 2) Calcular expressão agregada de IG por célula
expr_norm <- GetAssayData(T_cells, slot = "data")  # log-normalized
if(length(ig_genes) == 0) stop("Nenhum gene IG encontrado nos nomes das linhas")

T_cells$IG_score <- Matrix::colSums(expr_norm[ig_genes, , drop = FALSE])

# 3) threshold para IG+
quant_thresh <- quantile(T_cells$IG_score, 0.99)   # top 1%
fixed_thresh <- 0.5                                # ajuste se fizer sentido
# escolha: usar quantil
T_cells$IG_pos <- T_cells$IG_score >= quant_thresh
table(T_cells$IG_pos)

# 4) UMAP colorido por indivíduo (orig.ident) e por IG_pos
p1 <- DimPlot(T_cells, group.by = "orig.ident", pt.size = 0.6) + ggtitle("UMAP por indivíduo")
p2 <- FeaturePlot(T_cells, features = "IG_score", cols = c("lightgrey","red"), pt.size = 0.6) + ggtitle("IG score")
p3 <- DimPlot(T_cells, group.by = "IG_pos", cols = c("grey","blue"), pt.size = 0.6) + ggtitle("Células IG+ (azul)")

plot_grid(p1, p2, p3, ncol = 3)

# 5) Tabela por indivíduo
meta <- T_cells@meta.data %>% 
  mutate(orig.ident = as.character(orig.ident),
         cluster = as.character(Idents(T_cells)))

summary_df <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    total_cells = n(),
    IG_pos_cells = sum(IG_pos),
    frac_IG = IG_pos_cells / total_cells
  ) %>%
  arrange(desc(frac_IG))

# 6) Tabela por indivíduo e cluster
by_ind_cluster <- meta %>%
  group_by(orig.ident, cluster) %>%
  summarise(
    n_cells = n(),
    n_IG = sum(IG_pos),
    frac_IG = n_IG / n_cells
  ) %>%
  arrange(orig.ident, desc(frac_IG))

# visualizar
print(summary_df)
head(by_ind_cluster)

# 7) Barplot: fração IG+ por indivíduo
ggplot(summary_df, aes(x = reorder(orig.ident, -frac_IG), y = frac_IG)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = scales::percent(frac_IG, accuracy = 0.1)), vjust = -0.5, size = 3) +
  labs(x = "Indivíduo", y = "Fração de células IG+", title = "Fração de células IG+ por indivíduo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8) Heatmap / tile das frações por indivíduo x cluster (opcional)
library(viridis)
heat_df <- by_ind_cluster %>% mutate(frac_IG = ifelse(is.na(frac_IG), 0, frac_IG))
ggplot(heat_df, aes(x = cluster, y = orig.ident, fill = frac_IG)) +
  geom_tile() +
  scale_fill_viridis(name = "frac IG+") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Indivíduo", title = "Fração IG+ por indivíduo e cluster")

# 9) Violin do IG_score por indivíduo (para ver distribuição)
ggplot(meta, aes(x = orig.ident, y = IG_score)) +
  geom_violin(fill = "lightgrey") +
  geom_jitter(height = 0, width = 0.2, size = 0.5, alpha = 0.4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Indivíduo", y = "IG score", title = "Distribuição de IG score por indivíduo")

# 10) Teste simples: Fisher / proporção para cada indivíduo vs resto
test_results <- summary_df %>%
  rowwise() %>%
  mutate(
    other_IG = sum(summary_df$IG_pos_cells) - IG_pos_cells,
    other_total = sum(summary_df$total_cells) - total_cells,
    fisher_p = fisher.test(matrix(c(IG_pos_cells, total_cells - IG_pos_cells,
                                    other_IG, other_total - other_IG), nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(fisher_p, method = "BH")) %>%
  arrange(p_adj)

test_results

# 11) Verificar se IG+ estão concentradas em poucos barcodes (doublets)
ig_cells <- WhichCells(T_cells, expression = IG_pos == TRUE)
length(ig_cells)  # número de células IG+

# 12) Verificar percentuais por cluster e por indivíduo
table(meta[ig_cells, "orig.ident"])
table(meta[ig_cells, "cluster"])

# 13) Se suspeitar de doublets, rodar DoubletFinder (exemplo)
# library(DoubletFinder)
# siga o workflow do pacote: pK, pN, nExp estimativas e remover doublets

# --- Unir seurats ---

all(colnames(T_cells) %in% colnames(seurat_obj))

# Pegar os Idents do objeto reclusterizado
new_clusters <- Idents(T_cells)

# Criar coluna nova no objeto original
seurat_obj$T_subcluster <- as.character(Idents(seurat_obj))  # backup dos clusters originais
seurat_obj$T_subcluster[names(new_clusters)] <- as.character(new_clusters)

Idents(seurat_obj) <- seurat_obj$T_subcluster

DimPlot(seurat_obj, group.by = "T_subcluster")   # mostra tudo junto

# salvar objeto Seurat com os novos clusters
saveRDS(seurat_obj, file = "seurat_obj_reclustered.rds")

# carregar depois
seurat_obj <- readRDS("seurat_obj_reclustered.rds")




