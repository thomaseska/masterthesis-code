# Cleaned up script with data exploration steps in order

library(Seurat)
library(data.table)
library(ggplot2)

library(dplyr)

#FUNCTIONS
# create a Seurat object from a dgcMatrix and metadata data table
# compute variable features, scale data, compute PCA and UMAP
setUpSeurat <- function(count_dgc_matrix, metadata, nvar_features=3000,
                        pca = T, umap = T, dims = 15, res = 0.5) {
  # create object
  seur_obj <- CreateSeuratObject(count_dgc_matrix)
  
  # add metadata
  seur_obj <- AddMetaData(seur_obj, metadata)
  
  # find variable features
  seur_obj <- FindVariableFeatures(seur_obj,
                                   selection.method = "vst",
                                   nfeatures = nvar_features)
  
  # scale data
  seur_obj <- ScaleData(seur_obj)
  
  # PCA
  if (pca) {
    seur_obj <- RunPCA(seur_obj)
    if (umap) {
      seur_obj <- FindNeighbors(seur_obj, dims = 1:dims)
      seur_obj <- FindClusters(seur_obj, resolution = res)
      seur_obj <- RunUMAP(seur_obj, dims = 1:dims, n.neighbors = 26)
    }
  }
  
  return(seur_obj)
  
}

# Add tcell subtype annotations to overall metadata and prepare for Seurat::AddMetaData()
add_tcell_anno <- function(metadata_dt, tcell_metadata_dt){
  
  tcell_metadata_lookup <- tcell_metadata_dt[,c("Cell", "cellSubType")]
  tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)
  
  inds <- match(metadata_dt$Cell, tcell_metadata_lookup$Cell)
  metadata_dt[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]
  
  # Seurat works with data.frames
  # if needed convert back with as.data.table(df, keep.rownames="Cell")
  metadata_df <- as.data.frame(metadata_dt)
  rownames(metadata_df) <- metadata_df$Cell
  metadata_df$Cell <- NULL
  
  return(metadata_df)
}

# load data
count_mat_cohort1 <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1863_counts_cells_cohort1.rds")
metadata_cohort1 <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")

coh1_metadata <- add_tcell_anno(metadata_cohort1, tcell_metadata)
coh1_seur_obj <- setUpSeurat(count_mat_cohort1, coh1_metadata)

coh1_seur_obj <- SetIdent(coh1_seur_obj, value = coh1_seur_obj@meta.data$cellType)
DimPlot(coh1_seur_obj, reduction = "umap", raster = FALSE)

pretty_names <- data.table(data_name=c("B_cell", "Cancer_cell", "Endothelial_cell", "Fibroblast",
                                       "Mast_cell", "Myeloid_cell", "pDC", "TCELL_CD4_EM",
                                       "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating", "TCELL_CD4_N",
                                       "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                                       "TCELL_CD8_EMRA", "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating",
                                       "TCELL_CD8_N", "TCELL_CD8_RM", "TCELL_gdT", "TCELL_NK_CYTO", 
                                       "TCELL_NK_REST", "TCELL_Vg9Vd2_gdT"),
                           plot_name=as.factor(c("B cell", "Cancer cell", "Endothelial cell", "Fibroblast",
                                                 "Mast cell", "Myeloid cell", "pDC", "(T) CD4+ EM",
                                                 "(T) CD4+ EX", "(T) CD4+ EX Proliferating", "(T) CD4+ Naive",
                                                 "(T) CD4+ Reg", "(T) CD4+ Reg Proliferating", "(T) CD8+ EM",
                                                 "(T) CD8+ EMRA", "(T) CD8+ EX", "(T) CD8+ EX Proliferating",
                                                 "(T) CD8+ Naive", "(T) CD8+ RM", "(T) gdT", "(T) NK CYTO", 
                                                 "(T) NK REST", "(T) Vg9Vd2 gdT")),
                           type=c("other", "other", "other", "other",
                                  "other", "other", "other", "T cell",
                                  "T cell", "T cell", "T cell",
                                  "T cell", "T cell", "T cell",
                                  "T cell", "T cell", "T cell",
                                  "T cell", "T cell", "T cell", "T cell", 
                                  "T cell", "T cell"),
                           relevant=c(FALSE, FALSE, FALSE, FALSE,
                                      FALSE, FALSE, FALSE, TRUE,
                                      TRUE, TRUE, FALSE,
                                      TRUE, TRUE, TRUE,
                                      FALSE, TRUE, TRUE,
                                      FALSE, TRUE, FALSE, FALSE, 
                                      FALSE, FALSE))                      


coh1_metadata <- as.data.table(coh1_metadata, keep.rownames = "Cell")
reformed_meta <- merge(coh1_metadata, pretty_names, by.x = "cellType", by.y = "data_name", all.x = TRUE)

ggplot(reformed_meta, aes(factor(plot_name, levels= pretty_names$plot_name), color = BC_type))+geom_bar(position = "dodge")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("cell type")+
  ylab("number of cells")+
  scale_color_discrete(name = "T cell")

ggplot(coh1_metadata[!is.na(plot_name)], aes(reorder(plot_name, plot_name, function(x) -length(x)), fill=type))+geom_bar(position = "dodge")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("Number of Cells")+
  scale_fill_discrete("Type")

tcell_object <- subset(coh1_seur_obj, idents = c("TCELL_CD4_EM", "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating", "TCELL_CD4_N",
                                                 "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                                                 "TCELL_CD8_EMRA", "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating",
                                                 "TCELL_CD8_N", "TCELL_CD8_RM", "TCELL_gdT", "TCELL_NK_CYTO", 
                                                 "TCELL_NK_REST", "TCELL_Vg9Vd2_gdT"))
# rerun reductions
tcell_object <- ScaleData(tcell_object)
tcell_object <- RunPCA(tcell_object)
tcell_object <- FindNeighbors(tcell_object, dims = 1:15)
tcell_object <- FindClusters(tcell_object, resolution = 0.5)
tcell_object <- RunUMAP(tcell_object, dims = 1:15, min.dist = 0.5, n.neighbors = 5)

coh1_metadata_df <- as.data.frame(coh1_metadata)
rownames(coh1_metadata_df) <- coh1_metadata_df$Cell
coh1_metadata_df$Cell <- NULL
tcell_object <- AddMetaData(tcell_object, coh1_metadata_df)

tcell_object <- SetIdent(tcell_object, value = tcell_object@meta.data$plot_name)
p2 <- DimPlot(tcell_object, reduction = "umap", label = TRUE, repel = TRUE)


metadata_cohort1 <- merge(metadata_cohort1, pretty_names_short, by.x = "cellType", by.y="data_name")
metadata_cohort1 <- as.data.frame(metadata_cohort1)
rownames(metadata_cohort1) <- metadata_cohort1$Cell
metadata_cohort1$Cell <- NULL

coh1_seur_obj <- AddMetaData(coh1_seur_obj, metadata_cohort1)
coh1_seur_obj <- SetIdent(coh1_seur_obj, value=coh1_seur_obj@meta.data$plot_name)

coh1_seur_obj <- FindNeighbors(coh1_seur_obj, dims = 1:15)
coh1_seur_obj <- FindClusters(coh1_seur_obj, resolution = 0.5)
coh1_seur_obj <- RunUMAP(coh1_seur_obj, dims = 1:15)

p1 <- DimPlot(coh1_seur_obj, reduction="umap", label=TRUE, raster=FALSE)


ggarrange(p1, p2, labels = "AUTO")
