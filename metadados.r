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

amostras_por_arquivo <- list(
  dt1_control = "Healthy_Control_run1",
  dt1_DF = c("DF_Day_minus_1_run1","DF_Day_minus_1_run2","DF_Day_minus_2_run1","DF_Def_run1",
             "DF_Def_run2","DF_Wk2_run1"),
  dt1_DHF = c("DHF_Day_minus_1_run1","DHF_Day_minus_1_run2","DHF_Day_minus_2_run1",
              "DHF_Def_run1","DHF_Def_run2","DHF_Wk2_run1"),
  dt2_DF_1  = c("SRR12215051","SRR12215052","SRR12215053"),
  dt2_DF_2 = c("SRR12215054","SRR12215055","SRR12215056"),
  dt3_DF = c("SRR11088622_Primary1_D1","SRR11088623_Primary1_D3","SRR11088624_Primary2_D1",
              "SRR11088625_Primary2_D5","SRR11088626_Primary3_D1","SRR11088627_Primary3_D2",
              "SRR11088634_Secondary3_D0","SRR11088635_Secondary3_D1","SRR11088636_Secondary3_D5"),
  dt3_DHF = c("SRR11088628_Secondary1_D0","SRR11088629_Secondary1_D1","SRR11088630_Secondary1_D7",
              "SRR11088631_Secondary2_D1","SRR11088632_Secondary2_D1","SRR11088633_Secondary2_D5"),
  dt4_control = c("SRR22739533","SRR22739534","SRR22739543","SRR22739544","SRR22739551","SRR22739552"),
  dt4_DF = c("SRR22739530","SRR22739535","SRR22739536","SRR22739537","SRR22739539","SRR22739541",
             "SRR22739545","SRR22739547","SRR22739549","SRR22739550","SRR22739555"),
  dt4_DWS = c("SRR22739525","SRR22739527","SRR22739529","SRR22739538"),
  dt4_SD = c("SRR22739526","SRR22739528","SRR22739531","SRR22739532","SRR22739540","SRR22739542",
             "SRR22739546","SRR22739548","SRR22739553","SRR22739554")
)

# --- metadata DENV ---

metadata_denv <- data.frame(
  sample_id = c("dt1_control", "dt1_DF", "dt1_DHF", "dt2_DF_1", "dt2_DF_2", "dt3_DF", "dt3_DHF", "dt4_control", "dt4_DF", "dt4_DWS", "dt4_SD"),
  age = c("adult", "adult", "adult", "child", "child", "child", "child", "child", "child", "child", "child"),
  virus = c("control", "DENV-3", "DENV-3", "DENV-1", "DENV-1", "DENV-1", "DENV-3/1", "control", "NONE", "NONE", "NONE")
)

metadata_samples <- data.frame(
  sample_id = c("Healthy_Control_run1","DF_Day_minus_1_run1","DF_Day_minus_1_run2","DF_Day_minus_2_run1","DF_Def_run1",
             "DF_Def_run2","DF_Wk2_run1",
             "DHF_Day_minus_1_run1","DHF_Day_minus_1_run2","DHF_Day_minus_2_run1",
              "DHF_Def_run1","DHF_Def_run2","DHF_Wk2_run1",
              "SRR12215051","SRR12215052","SRR12215053",
              "SRR12215054","SRR12215055","SRR12215056",
              "SRR11088622_Primary1_D1","SRR11088623_Primary1_D3",
              "SRR11088624_Primary2_D1","SRR11088625_Primary2_D5",
              "SRR11088626_Primary3_D1","SRR11088627_Primary3_D2",
              "SRR11088634_Secondary3_D0","SRR11088635_Secondary3_D1","SRR11088636_Secondary3_D5",
              "SRR11088628_Secondary1_D0","SRR11088629_Secondary1_D1","SRR11088630_Secondary1_D7",
              "SRR11088631_Secondary2_D1","SRR11088632_Secondary2_D1","SRR11088633_Secondary2_D5",
              "SRR22739533","SRR22739534","SRR22739543","SRR22739544","SRR22739551","SRR22739552", 
              "SRR22739530","SRR22739535","SRR22739536","SRR22739537","SRR22739539","SRR22739541",
             "SRR22739545","SRR22739547","SRR22739549","SRR22739550","SRR22739555",
             "SRR22739525","SRR22739527","SRR22739529","SRR22739538",
             "SRR22739526","SRR22739528","SRR22739531","SRR22739532","SRR22739540","SRR22739542",
             "SRR22739546","SRR22739548","SRR22739553","SRR22739554"),

  sex = c("female", "male", "male", "male", "male", "male", "male",
              "male", "male", "male", "male", "male", "male",
              "male", "male", "male", 
              "male", "male", "male",
              "male", "male",
              "male", "male",
              "male", "male",
              "female", "female",
              "male", "male", "male",
              "male", "male", "male",
              "male", "male", "male",
              "female", "female","female",
              "female", "female","female",
              "male", "NONE" ,"male", "female","female", "male",
              "female", "male", "female", "female", "male", "female", 
              "female", "male", "male", "male", "male", 
              "male", "male", "male", "female",
              "female", "female", "female", "male", "male", "female",
              "none", "female", "female", "male"),

  group = c("control", "acute", "acute", "acute", "defervescence", "defervescence", "convalescent",
            "acute", "acute", "acute", "defervescence", "defervescence", "convalescent",
            "acute", "acute", "acute",
            "acute", "acute", "acute",
            "acute", "post-defervescence",
            "acute", "post-defervescence",
            "acute", "post-defervescence",
            "defervescence", "post-defervescence", "post-defervescence",
            "defervescence", "post-defervescence", "post-defervescence",
            "post-defervescence", "post-defervescence", "post-defervescence",
            "control", "control", "control", "control", "control", "control",
            "acute", "acute", "acute", "acute", "acute", "acute",
            "acute", "acute", "acute", "acute",
            "acute", "acute", "acute", "acute",
            "acute", "acute", "acute", "acute", "acute", "acute",
            "acute", "acute", "acute", "acute"),
  infection = c("control", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY",
               "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY",
               "PRIMARY", "PRIMARY", "PRIMARY", 
               "PRIMARY", "PRIMARY", "PRIMARY",
               "PRIMARY", "PRIMARY",
              "PRIMARY", "PRIMARY",
              "PRIMARY", "PRIMARY",
              "SECUNDARY", "SECONDARY", "SECONDARY",
              "SECUNDARY", "SECUNDARY", "SECUNDARY",
              "control", "control", "control", "control", "control", "control",
              "SECUNDARY","SECUNDARY", "SECUNDARY", "SECUNDARY", "PRIMARY", "indeterminate",
              "SECUNDARY", "PRIMARY", "SECUNDARY", "SECUNDARY", "SECUNDARY", 
              "PRIMARY", "equiv", "PRIMARY", "PRIMARY",
              "PRIMARY", "SECUNDARY", "PRIMARY", "SECUNDARY", "SECUNDARY", "SECUNDARY",
              "none", "SECONDARY", "SECONDARY","SECONDARY")


# Inicializar lista para armazenar os metadados expandidos
metadados_expandidos <- list()

# Loop por cada entrada em amostras_por_arquivo
for (nome_arquivo in names(amostras_por_arquivo)) {
  amostras <- amostras_por_arquivo[[nome_arquivo]]
  
  # Buscar os metadados correspondentes
  metadado <- metadata_denv[metadata_denv$sample_id == nome_arquivo, ]
  
  # Criar data.frame para cada amostra individual
  df <- data.frame(
    sample_id = amostras,
    arquivo = nome_arquivo,
    age = metadado$age,
    virus = metadado$virus,
    stringsAsFactors = FALSE
  )
  
  # Adicionar à lista
  metadados_expandidos[[nome_arquivo]] <- df
}

# Unir tudo em um único data.frame
metadata_final <- do.call(rbind, metadados_expandidos)


library(dplyr)

seurat_integrado <- readRDS("caminho/para/seu_arquivo.rds")


seurat_integrado@meta.data <- seurat_integrado@meta.data %>%
  left_join(metadata_final, by = "sample_id")


DimPlot(seurat_integrado, reduction = "umap", group.by = "virus", label = TRUE) +
  ggtitle("UMAP por Sorotipo Viral")

DimPlot(seurat_integrado, reduction = "umap", group.by = "age", label = TRUE) +
  ggtitle("UMAP por idade")

seurat_integrado@meta.data <- seurat_integrado@meta.data %>%
  left_join(metadata_samples, by = "sample_id")

DimPlot(seurat_integrado, reduction = "umap", group.by = "group", label = TRUE) +
  ggtitle("UMAP por Grupo Clínico")

DimPlot(seurat_integrado, reduction = "umap", group.by = "infection", label = TRUE) +
  ggtitle("UMAP por Infecção")

  