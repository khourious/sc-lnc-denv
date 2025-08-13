remotes::install_github("satijalab/azimuth")
remotes::install_github("Jacobog02/azimuth_stopgap")
update.packages(oldPkgs = c("withr", "rlang"))

library(Seurat)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(Azimuth)

# Azimuth HubMap. Ele é um dataset de 10 mil células PBMC gerado pela 10X Genomics

# Instalar Azimuth se necessário
remotes::install_github("satijalab/azimuth")

# Carregar Azimuth
library(Azimuth)

# Carregar referência PBMC oficial
reference <- LoadReference("pbmc")

# Rodar anotação automática
seurat_annotado <- RunAzimuth(seurat_integrado, reference = reference)

# Visualizar os tipos celulares anotados
DimPlot(seurat_annotado, group.by = "predicted.celltype.l1", label = TRUE) +
  ggtitle("Tipos celulares gerais (Azimuth)")


DimPlot(seurat_annotado, group.by = "predicted.celltype.l2", label = TRUE) +
  ggtitle("Subtipos celulares refinados (Azimuth)")


###############################
#### --- com seurat

install.packages("SeuratData")
library(SeuratData)

# Baixar e carregar referência PBMC
InstallData("pbmc3k")
reference <- LoadData("pbmc3k")

# https://seurat.nygenome.org/src/contrib/pbmc3k.SeuratData_3.1.4.tar.gz
install.packages("pbmc3k.SeuratData_3.1.4.tar.gz", repos = NULL, type = "source")

library(pbmc3k.SeuratData)
pbmc3k <- LoadData("pbmc3k")

DefaultAssay(seurat_integrado) <- "SCT"

anchors <- FindTransferAnchors(
  reference = reference,
  query = seurat_integrado,
  dims = 1:30,
  reference.assay = "RNA",
  query.assay = "SCT",  
  reduction = "pcaproject"
)

seurat_integrado <- TransferData(
  anchorset = anchors,
  refdata = reference$celltype,
  dims = 1:30
)

DimPlot(seurat_integrado, group.by = "predicted.id", label = TRUE) +
  ggtitle("Anotação automática com referência Seurat")



