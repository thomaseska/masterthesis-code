# hierachichal clustering
# --> cell type relations, which can we distinguish?

library(HGC)
library(Seurat)
library(data.table)        

sampled_meta <- fread("data/bassez/sampled_cells_metadata_cohort1.csv")
sampled_counts <- readRDS("data/bassez/sampled_counts_cohort1.rds")


sampled_seurat <- setUpSeurat(coh1_count_mat, coh1_metadata)
saveRDS(sampled_seurat, "data/bassez/sampled_seurat.rds")
sampled_seurat <- FindClusteringTree(sampled_seurat)
HGC.PlotDendrogram(sampled_seurat@graphs$ClusteringTree,
                   labels = sampled_seurat@meta.data$cellType,
                   k=20)
sampled_seurat <- SetIdent(sampled_seurat, value = "cellType")
sampled_seurat <- BuildClusterTree(sampled_seurat, assay = "RNA", slot = "counts")

PlotClusterTree(sampled_seurat)

seur_dist_mat <- sampled_seurat@tools$BuildClusterTree %>% cophenetic()
seur_dist_mat <- melt(as.matrix(seur_dist_mat), varnames = c("row", "col"))
seur_dist_mat <- as.data.table(seur_dist_mat)

dist_mat <- seur_dist_mat

signature <- readRDS("data/bassez/deconvolution_results/DWLS/signature.rds")
dist_mat <- hclust(dist(t(signature))) %>% cophenetic()
dist_mat <- melt(as.matrix(dist_mat), varnames = c("row", "col"))
dist_mat <- as.data.table(dist_mat)

obv_types <- c("B_cell", "Myeloid_cell", "Fibroblast", "Cancer_cell", "Endothelial_cell")
ggplot(dist_mat[!row %in% obv_types & !col %in% obv_types], aes(row, col, fill=value))+geom_tile()

clust1 <- c("TCELL_CD4_EM", "TCELL_CD8_RM", "TCELL_CD8_EM", "TCELL_gdT", "TCELL_Vg9Vd2_gdT")
ggplot(dist_mat[row %in% clust1 & col %in% clust1], aes(row, col, fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value, 3)))
clust2 <- c("TCELL_CD4_EX", "TCELL_CD8_EX")
ggplot(dist_mat[row %in% clust2 & col %in% clust2], aes(row, col, fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value, 3)))
clust3 <- c("TCELL_CD8_N", "TCELL_CD4_N")
ggplot(dist_mat[row %in% clust3 & col %in% clust3], aes(row, col, fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value, 3)))
clust4 <- c("TCELL_CD8_EMRA", "TCELL_NK_CYTO", "TCELL_NK_REST")
ggplot(dist_mat[row %in% clust4 & col %in% clust4], aes(row, col, fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value, 3)))
clust5 <- c("TCELL_CD4_EX_Proliferating", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EX_Proliferating")
ggplot(dist_mat[row %in% clust5 & col %in% clust5], aes(row, col, fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value, 3)))
