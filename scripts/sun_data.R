library(Seurat)
library(data.table)

sun_mtx <- ReadMtx(mtx = "data/sun/matrix.mtx.gz",
                   cells = "data/sun/barcodes.tsv.gz",
                   features = "data/sun/features.tsv.gz")
sun_meta <- fread("data/sun/OMIX001073-20-04.tsv")

# on server
sun_mtx <- ReadMtx(mtx = "/nfs/data/tcell_deconvolution/data/sun/matrix.mtx.gz",
                   cells = "/nfs/data/tcell_deconvolution/data/sun/barcodes.tsv.gz",
                   features = "/nfs/data/tcell_deconvolution/data/sun/features.tsv.gz")
sun_meta <- fread("/nfs/data/tcell_deconvolution/data/sun/OMIX001073-20-04.tsv")
# 

sun_obj <- CreateSeuratObject(counts = sun_mtx)
rm(sun_mtx)

sun_meta <- as.data.frame(sun_meta)
rownames(sun_meta) <- sun_meta$V1
sun_meta$V1 <- NULL

sun_obj<- AddMetaData(sun_obj, sun_meta)
sun_obj <- SetIdent(sun_obj, value = sun_obj@meta.data$cluster)

sun_obj <- FindVariableFeatures(sun_obj, selection.method = "vst", nfeatures = 2000)
sun_obj <- ScaleData(sun_obj)

sun_obj <- RunPCA(sun_obj)
sun_obj <- FindNeighbors(sun_obj, dims = 1:30)
sun_obj <- FindClusters(sun_obj, resolution = 0.5)
sun_obj <- RunUMAP(sun_obj, dims = 1:30)

sun_obj <- SetIdent(sun_obj, value = sun_obj@meta.data$cluster)
DimPlot(sun_obj, reduction = "umap", label = T)


# load bassez data on Server
tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
coh1_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
tcell_metadata_lookup <- tcell_metadata[,c("Cell", "cellSubType")]
tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)
inds <- match(coh1_metadata$Cell, tcell_metadata_lookup$Cell)
coh1_metadata[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]

coh1_count_mat <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1863_counts_cells_cohort1.rds")

bass_obj <- CreateSeuratObject(coh1_count_mat)
rm(coh1_count_mat)

coh1_metadata <- as.data.frame(coh1_metadata)
rownames(coh1_metadata) <- coh1_metadata$Cell
coh1_metadata$Cell <- NULL

bass_obj <- AddMetaData(bass_obj, coh1_metadata)

# integration

## normalize the data
bass_obj <- NormalizeData(bass_obj)
sun_obj <- NormalizeData(sun_obj)


## find features
bass_obj <- FindVariableFeatures(bass_obj, selection.method = "vst", nfeatures = 2000)
sun_obj <- FindVariableFeatures(sun_obj, selection.method = "vst", nfeatures = 2000)

bass_obj <- SCTransform(bass_obj)
bass_obj <- RunPCA(bass_obj)

sun_obj <- SCTransform(sun_obj, conserve.memory = TRUE)
saveRDS(sun_obj, "/nfs/data/tcell_deconvolution/data/sun/seurat_obj_sct.rds")

obj_list <- list(bass_obj, sun_obj)

## select features
features <- SelectIntegrationFeatures(obj_list)

obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, verbose = F, features = features)


anchor_list <- FindIntegrationAnchors(obj_list, anchor.features = features, reduction = "rpca", normalization.method = "SCT")
saveRDS(anchor_list, "/nfs/data/tcell_deconvolution/data/integration_anchor_list.rds")


integrated_obj <- IntegrateData(anchorset = anchor_list, normalization.method = "SCT")
saveRDS(integrated_obj, "/nfs/data/tcell_deconvolution/data/integrated_seurat_obj.rds")

# data origin metadata
data_origin <- data.frame("Cell"=sun_meta$V1, "Origin"="Sun")
data_origin <- rbind(data_origin, data.frame("Cell"=coh1_metadata$Cell, "Origin"="Bassez"))
rownames(data_origin) <- c(sun_meta$V1, coh1_metadata$Cell)
integrated_obj <- AddMetaData(integrated_obj, data_origin)
integrated_obj <- SetIdent(integrated_obj, value = integrated_obj@meta.data$Origin)

integrated_obj <- RunPCA(integrated_obj, assay = "integrated")
PCAPlot(object = integrated_obj, split.by = "Origin", raster=TRUE)

integrated_obj <- RunUMAP(integrated_obj, dims = 1:50)
DimPlot(integrated_obj, split.by = "Origin", raster = FALSE)


# label transfer only

transfer_anchors <- FindTransferAnchors(reference = bass_obj, query = sun_obj, reference.assay = "SCT", query.assay = "SCT")
predictions <- TransferData(anchorset = transfer_anchors, refdata = bass_obj$cellType, dims = 1:30)

sun_obj <- AddMetaData(sun_obj, metadata = predictions)
sun_obj <- SetIdent(sun_obj, value = sun_obj@meta.data$predicted.id)
DimPlot(sun_obj, reduction = "umap", label = TRUE, raster = FALSE)
saveRDS(sun_obj, "/nfs/data/tcell_deconvolution/data/sun/seurat_obj_sct_with_bass_labels.rds")


tcell_list <- c("TCELL_CD8_RM", "TCELL_NK_REST", "TCELL_CD4_EM", "TCELL_CD8_EM", "TCELL_gdT", "TCELL_CD4_N",
                "TCELL_Vg9Vd2_gdT", "TCELL_CD4_REG", "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating", "TCELL_CD8_EX", 
                "TCELL_CD8_EX_Proliferating")

sampled_sun_obj <- subset(sun_obj, idents = tcell_list)
sampled_sun_obj <- RunPCA(sampled_sun_obj)
sampled_sun_obj <- FindNeighbors(sampled_sun_obj, dims=1:30)
sampled_sun_obj <- FindClusters(sampled_sun_obj)
sampled_sun_obj <- RunUMAP(sampled_sun_obj, dims = 1:30)

sampled_sun_obj <- SetIdent(sampled_sun_obj, value = sampled_sun_obj@meta.data$predicted.id)
DimPlot(sampled_sun_obj, reduction = "umap", label=TRUE, raster=FALSE)
saveRDS(sampled_sun_obj, "/nfs/data/tcell_deconvolution/data/sun/seurat_obj_sct_with_bass_labels_only_tcells.rds")



# label transfer for tcells
# subset bass_obj

bass_obj <- SetIdent(bass_obj, value = bass_obj@meta.data$cellType)
bass_obj <- subset(bass_obj, idents = c("TCELL_gdT", "TCELL_NK_REST", "TCELL_CD4_N", "TCELL_CD8_EM",
                                        "TCELL_CD8_N", "TCELL_CD8_RM", "TCELL_CD4_EM", "TCELL_CD4_REG",
                                        "TCELL_CD4_EX", "TCELL_CD8_EX_Proliferating",
                                        "TCELL_CD8_EX", "TCELL_CD8_EMRA", "TCELL_Vg9Vd2_gdT",
                                        "TCELL_NK_CYTO", "TCELL_CD4_EX_Proliferating", 
                                        "TCELL_CD4_REG_Proliferating"))

# subset sun_obj
sun_obj <- subset(sun_obj, idents = "T cells & NK cells")


# transfer
transfer_anchors <- FindTransferAnchors(reference = bass_obj, query = sun_obj, reference.assay = "SCT", query.assay = "SCT", normalization.method = "SCT")
predictions <- TransferData(anchorset = transfer_anchors, refdata = bass_obj$cellType, dims = 1:30)


sun_obj <- AddMetaData(sun_obj, metadata = predictions)
sun_obj <- SetIdent(sun_obj, value = sun_obj@meta.data$predicted.id)
DimPlot(sun_obj, reduction = "umap", label = TRUE, raster = FALSE)


