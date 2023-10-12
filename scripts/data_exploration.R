# Basic Dataset statistics

library(Seurat)
library(data.table)
library(ggplot2)

#FUNCTIONS
# create a Seurat object from a dgcMatrix and metadata data table
# compute variable features, scale data, compute PCA and UMAP
setUpSeurat <- function(count_dgc_matrix, metadata_dt, nvar_features=3000,
                        pca = T, umap = T, dims = 15, res = 0.5) {
  # create object
  seur_obj <- CreateSeuratObject(count_dgc_matrix)
  rm(count_dgc_matrix)
  # add metadata
  metadata_df <- as.data.frame(metadata_dt)
  rownames(metadata_df) <- metadata_df$Cell
  metadata_df$Cell <- NULL
  seur_obj <- AddMetaData(seur_obj, metadata_df)
  
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

# COHORT 1
count_mat_cohort1 <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1863_counts_cells_cohort1.rds")
metadata_cohort1 <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
cohort1_seur_obj <- setUpSeurat(count_mat_cohort1, metadata_cohort1)
rm(count_mat_cohort1)


ggplot(metadata_cohort1, aes(cellType))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  labs(title = "Number of Cells per Type in Cohort 2 (Bassez et al)")+
  xlab("Cell Type")

LabelPoints(VariableFeaturePlot(cohort1_seur_obj),
            points = head(VariableFeatures(cohort1_seur_obj), 10),
            repel = TRUE)

VlnPlot(cohort1_seur_obj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident")

DimPlot(cohort1_seur_obj, reduction = "pca", label = TRUE)
ElbowPlot(cohort1_seur_obj)
DimPlot(cohort1_seur_obj, reduction = "umap", group.by = "cellType", label = TRUE, repel = TRUE)
DimPlot(cohort1_seur_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)



# assume counts as dgcMatrix
count_mat <- readRDS("data/bassez/1864-counts_tcell_cohort1.rds")
metadata <- fread("data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")

# number of cells per cell type
metadata[, .N, by=cellSubType]
# Barplot per celltype
ggplot(metadata, aes(cellSubType))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  labs(title = "Number of Cells per Subtype in T-Cell Dataset (Bassez et al)")+
  xlab("Cell Subtype")


# create a Seurat object
tcell_seur_obj <- CreateSeuratObject(count_mat)

# add metadata 
metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$Cell
#metadata_df$Cell <- NULL
tcell_seur_obj <- AddMetaData(tcell_seur_obj, metadata_df)

# set ident to patient id
tcell_seur_obj <- SetIdent(tcell_seur_obj, value = "patient_id")
# set ident to subcelltype
tcell_seur_obj <- SetIdent(tcell_seur_obj, value = "cellSubType")

VlnPlot(tcell_seur_obj, features = c("nFeature_RNA", "nCount_RNA"))

tcell_seur_obj <- FindVariableFeatures(tcell_seur_obj, selection.method = "vst", nfeatures = 3000)
var_feat_plt <- VariableFeaturePlot(tcell_seur_obj)
LabelPoints(var_feat_plt,
                    points = head(VariableFeatures(tcell_seur_obj), 10),
                    repel = TRUE)
# outlier HBB: hemoglobin subunit beta

# scale data
tcell_seur_obj <- ScaleData(tcell_seur_obj)

# PCA
tcell_seur_obj <- RunPCA(tcell_seur_obj)
DimPlot(tcell_seur_obj, reduction = "pca", label = TRUE)

# cluster
ElbowPlot(tcell_seur_obj) # --> to find good number of clustering dimensions
tcell_seur_obj <- FindNeighbors(tcell_seur_obj, dims = 1:15)
tcell_seur_obj <- FindClusters(tcell_seur_obj, resolution = 0.5)

# UMAP
tcell_seur_obj <- RunUMAP(tcell_seur_obj, dims = 1:15)
um1 <- DimPlot(tcell_seur_obj, reduction = "umap", group.by = "cellSubType", label = TRUE, repel = TRUE)
um2 <- DimPlot(tcell_seur_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
um1 | um2

DimPlot(tcell_seur_obj, reduction = "umap", group.by = "patient_id", label = TRUE, repel = TRUE)


# put tcell descriptors into cohort1

tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
coh1_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")

tcell_metadata_lookup <- tcell_metadata[,c("Cell", "cellSubType")]
tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)

inds <- match(coh1_metadata$Cell, tcell_metadata_lookup$Cell)
coh1_metadata[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]


# randomly sample cells: keep all t cells, take 1000 of all other types
coh1_metadata[, .N, by=cellType]
celltypes <- coh1_metadata$cellType %>% unique()

sampled_cells <- sapply(celltypes, function(type){
  if (startsWith(type, "TCELL_") || type == "T_cell") {
    size <- coh1_metadata[cellType==type]$Cell %>% length()
  } else {
    size <-min(1000, coh1_metadata[cellType==type]$Cell %>% length())
  }
  sampled <- sample(coh1_metadata[cellType==type]$Cell, size = size, replace = FALSE)
  return(sampled)
})
sampled_cells <- unlist(sampled_cells)

sampled_metadata <- coh1_metadata[Cell %in% sampled_cells]
fwrite(sampled_metadata, "data/bassez/sampled_cells_metadata_cohort1.csv")

# sample the count matrix

coh1_count_mat <- readRDS("data/bassez/1863-counts_cells_cohort1.rds")
coh1_count_mat <- coh1_count_mat[,colnames(coh1_count_mat) %in% sampled_cells]
saveRDS(coh1_count_mat, "data/bassez/sampled_counts_cohort1.rds")


# scDEED

check_umap <- scDEED(coh1_count_mat, num_pc = 17, use_method = "umap", visualization = TRUE)
saveRDS(check_umap, "/nfs/data/tcell_deconvolution/data/bassez/scDEED/scDEED_obj.rds")

tcell_count_mat <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1864-counts_tcell_cohort1.rds")
check_tcell_umap <- scDEED(tcell_count_mat, num_pc = 16, use_method = "umap", visualization = TRUE)

tcell_seur_obj <- CreateSeuratObject(tcell_count_mat)
tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
metadata_df <- as.data.frame(tcell_metadata)
rownames(metadata_df) <- metadata_df$Cell
tcell_seur_obj <- AddMetaData(tcell_seur_obj, metadata_df)
tcell_seur_obj <- FindVariableFeatures(tcell_seur_obj, selection.method = "vst", nfeatures = 3000)

tcell_seur_obj <- ScaleData(tcell_seur_obj)

# PCA
tcell_seur_obj <- RunPCA(tcell_seur_obj)
DimPlot(tcell_seur_obj, reduction = "pca", label = TRUE)

tcell_seur_obj <- FindNeighbors(tcell_seur_obj, dims = 1:16)
tcell_seur_obj <- FindClusters(tcell_seur_obj, resolution = 0.5)

# UMAP
tcell_seur_obj <- RunUMAP(tcell_seur_obj, dims = 1:16, n.neighbors = 5, min.dist = 0.9)

tcell_seur_obj <- SetIdent(tcell_seur_obj, value = tcell_seur_obj@meta.data$cellSubType)

DimPlot(tcell_seur_obj, reduction = "umap", label = TRUE)


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

ggplot(reformed_meta, aes(factor(plot_name, levels= pretty_names$plot_name), fill = type))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("cell type")+
  ylab("number of cells")+
  scale_fill_discrete(name = "T cell")

