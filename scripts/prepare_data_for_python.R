library(anndata)
library(MatrixExtra)


# load data
tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
coh1_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
tcell_metadata_lookup <- tcell_metadata[,c("Cell", "cellSubType")]
tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)
inds <- match(coh1_metadata$Cell, tcell_metadata_lookup$Cell)
coh1_metadata[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]

coh1_count_mat <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1863_counts_cells_cohort1.rds")

# sample the data
# pick 1000 cells per cell type, if there are less than 1000 take all
sampled_cells <- sapply(unique(coh1_metadata$cellType), function(type){
  if (startsWith(type, "TCELL_")) {
    size <- min(coh1_metadata[cellType==type]$Cell %>% length(), 1000)
  } else if (type == "T_cell"){
    size <- 0
  } else {
    size <-min(1000, coh1_metadata[cellType==type]$Cell %>% length())
  }
  sampled <- sample(coh1_metadata[cellType==type]$Cell, size = size, replace = FALSE)
  return(sampled)
})
sampled_cells <- unlist(sampled_cells)
sampled_meta <- coh1_metadata[Cell %in% sampled_cells]
sampled_count_mat <- coh1_count_mat[,colnames(coh1_count_mat) %in% sampled_cells]

t_count_mat <- t_shallow(sampled_count_mat)

testset <- t_count_mat[1:10, 1:10]
testcells <- coh1_metadata[Cell %in% rownames(testset)]
testcells <- testcells[order(match(Cell, rownames(testset)))]
ad <- AnnData(X = testset, obs =data.frame(group = testcells$cellType, row.names = rownames(testset)), var = data.frame(row.names = colnames(testset)))
ad


# real
sampled_meta <- sampled_meta[Cell %in% rownames(t_count_mat)]
sampled_meta <- sampled_meta[order(match(Cell, rownames(t_count_mat)))]
ad <- AnnData(X = t_count_mat,
              obs = data.frame(groups = sampled_meta$cellType, 
                               row.names = rownames(t_count_mat)),
              var = data.frame(row.names = colnames(t_count_mat)))

write_h5ad(ad, "/nfs/data/tcell_deconvolution/data/bassez/nsforest_prepared_obj_sampled.h5ad")
