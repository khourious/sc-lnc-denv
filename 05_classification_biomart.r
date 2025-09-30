library(dplyr)
library(Matrix)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(biomaRt)

# Monocitos
# --- Heatmap dos 30 genes mais expressos ---
top_genes <- res_df %>% arrange(padj) %>% head(30) %>% pull(gene) # extrai os 30 genes mais significativos
heatmap_data <- pseudobulk_counts[top_genes, ] # filtra os dados de contagem para esses genes

pheatmap(log2(heatmap_data + 1), cluster_rows = TRUE, cluster_cols = TRUE, # log2 transformation
         annotation_col = sample_meta["dengue_classification"], # adiciona anotação
         main = "Top 30 Genes - Monocytes") # título do heatmap

top_genes <- resdf_df %>% arrange(padj) %>% head(30) %>% pull(gene) # extrai os 30 genes mais significativos
heatmap_data <- pseudobulk_counts[top_genes, ] # filtra os dados de contagem para esses genes

pheatmap(log2(heatmap_data + 1), cluster_rows = TRUE, cluster_cols = TRUE, # log2 transformation
         annotation_col = sample_meta["dengue_classification"], # adiciona anotação
         main = "Top 30 Genes - Monocytes") # título do heatmap

top_genes <- resdhf_df %>% arrange(padj) %>% head(30) %>% pull(gene) # extrai os 30 genes mais significativos
heatmap_data <- pseudobulk_counts[top_genes, ] # filtra os dados de contagem para esses genes

pheatmap(log2(heatmap_data + 1), cluster_rows = TRUE, cluster_cols = TRUE, # log2 transformation
         annotation_col = sample_meta["dengue_classification"], # adiciona anotação
         main = "Top 30 Genes - Monocytes") # título do heatmap

# --- Identificação de lncRNAs e genes codificadores de proteínas ---
# Usando biomaRt para anotação
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # humano

# Anotação dos genes diferencialmente expressos
annot <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = res_df$gene, # DF vs DHF
  mart = mart
)

annot_df <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = resdf_df$gene, # DF vs control
  mart = mart
)

annot_dhf <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = resdf_df$gene, # DHF vs control
  mart = mart
)

# Mesclando a anotação com os resultados
res_annot <- merge(res_df, annot, by.x = "gene", by.y = "external_gene_name")
res_annot_df <- merge(resdf_df, annot, by.x = "gene", by.y = "external_gene_name")
res_annot_dhf <- merge(resdhf_df, annot, by.x = "gene", by.y = "external_gene_name")

# --- Salvar resultados anotados ---
write.csv(res_annot, "DEG_results/monocytes/DF_DHF_biomart.csv", row.names = FALSE)
write.csv(res_annot_df, "DEG_results/monocytes/DF_control_biomart.csv", row.names = FALSE)
write.csv(res_annot_dhf, "DEG_results/monocytes/DHF_control_biomart.csv", row.names = FALSE)
