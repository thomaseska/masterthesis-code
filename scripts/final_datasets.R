# create single cell (bayesprism & dwls) and bulk (even & mirror) datasets
# for sc draw from dataset, then remove drawn cells from overall data
# for bulk create datasets with SimBu (count+cpm)

# load data
tcell_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
coh1_metadata <- fread("/nfs/data/tcell_deconvolution/data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
tcell_metadata_lookup <- tcell_metadata[,c("Cell", "cellSubType")]
tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)
inds <- match(coh1_metadata$Cell, tcell_metadata_lookup$Cell)
coh1_metadata[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]
# remove ambigous T cells
coh1_metadata <- coh1_metadata[cellType!="T_cell"]

coh1_count_mat <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/1863_counts_cells_cohort1.rds")

celltypes <- coh1_metadata$cellType %>% unique()

# sc sample function
cell_sampler <- function(type, metadata, method) {
  if (method == "bayesprism") {
    maxsize <- 1000
    size <- min(maxsize, 0.5 * (metadata[cellType==type]$Cell %>% length()))
  } else if (method == "dwls") {
    maxsize <- 100
    size <- min(maxsize,(metadata[cellType==type]$Cell %>% length()))
  } else {
    return(NA)
  }
  

  sampled <- sample(metadata[cellType==type]$Cell, size = size, replace = FALSE)
  return(sampled)
}

# sample sc
sampled_cells <- sapply(celltypes, cell_sampler, metadata = coh1_metadata, method = "bayesprism")
sampled_cells <- unlist(sampled_cells)
sc1_bayes_sampled_metadata <- coh1_metadata[Cell %in% sampled_cells]

sampled_cells <- sapply(celltypes, cell_sampler, metadata = sc1_bayes_sampled_metadata, method = "dwls")
sampled_cells <- unlist(sampled_cells)
sc1_dwls_sampled_metadata <- coh1_metadata[Cell %in% sampled_cells]

updated_metadata <- coh1_metadata[!(Cell %in% sc1_bayes_sampled_metadata$Cell)]

write.csv2(sc1_bayes_sampled_metadata,
           "/nfs/data/tcell_deconvolution/data/bassez/final_data/sc_bayes_sampled_metadata.csv")
write.csv2(sc1_dwls_sampled_metadata,
           "/nfs/data/tcell_deconvolution/data/bassez/final_data/sc_dwls_sampled_metadata.csv")

# create sc pipeline data
create_sc_data <- function(count_mat, metadata, folder_path) {
  sampled_sc_count_mat <- count_mat[,colnames(count_mat) %in% metadata$Cell]
  saveRDS(sampled_sc_count_mat, file.path(folder_path, "matrix_counts.rds"))
  # print(dim(sc_count_mat))
  
  batch <- metadata[, c("Cell", "patient_id")]
  batch <- batch[order(match(Cell, colnames(sampled_sc_count_mat)))]
  batch <- batch$patient_id
  saveRDS(batch, file.path(folder_path, "batch.rds"))
  
  celltype_annotation <- metadata[, c("Cell", "cellType")]
  celltype_annotation <- celltype_annotation[order(match(Cell, 
                                                         colnames(sampled_sc_count_mat))
  )]
  celltype_annotation <- celltype_annotation$cellType
  saveRDS(celltype_annotation, file.path(folder_path, "celltype_annotations.rds"))
  
}

# ...for bayesprism
folder_path <- "/nfs/data/tcell_deconvolution/data/bassez/final_data/singleCell/finalBayes"
create_sc_data(coh1_count_mat, sc1_bayes_sampled_metadata, folder_path)

# ...for dwls
folder_path <- "/nfs/data/tcell_deconvolution/data/bassez/final_data/singleCell/finalDWLS"
create_sc_data(coh1_count_mat, sc1_dwls_sampled_metadata, folder_path)

# create bulk pipeline data
create_bulk_data <- function(count_mat, metadata, folder_path, dataset_name, scenario){
  
  count_mat <- count_mat[,colnames(count_mat) %in% metadata$Cell]
  
  annotation <- data.frame("ID"=metadata$Cell,
                           "cell_type"=metadata$cellType,
                           row.names = metadata$Cell)
  sim_ds <- SimBu::dataset(annotation = annotation,
                           count_matrix = count_mat,
                           tpm_matrix = counts_to_cpm(count_mat),
                           name = "bassez",
                           scale_tpm = FALSE)
  simulation <- SimBu::simulate_bulk(data = sim_ds,
                                     scenario = scenario,
                                     scaling_factor = "NONE",
                                     ncells = 10000,
                                     nsamples = 10,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                     run_parallel = FALSE)
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]],
          file = file.path(folder_path, paste0(dataset_name, "_tpm.rds")))
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
          file = file.path(folder_path, paste0(dataset_name, "_counts.rds")))
  
  facs <- simulation$cell_fractions %>% t()
  saveRDS(facs, file.path(folder_path, paste0(dataset_name, "_facs.rds")))
  
  
}

# ...with even fractions
folder_path <- "/nfs/data/tcell_deconvolution/data/bassez/final_data/PBMC/finalEven/"
create_bulk_data(coh1_count_mat, updated_metadata, folder_path, "finalEven", "even")

# ...with mirror fractions
folder_path <- "/nfs/data/tcell_deconvolution/data/bassez/final_data/PBMC/finalMirror/"
create_bulk_data(coh1_count_mat, updated_metadata, folder_path, "finalMirror", "mirror_db")



# simulations
# underrep

lapply(celltypes, function(type){
  anti_type_meta <- updated_metadata[cellType!=type]
  
  dataset_name <- paste0("final_no_", type, "even")
  folder_path <- paste0("/nfs/data/tcell_deconvolution/data/bassez/final_data/PBMC/final_no_", type, "even")
  dir.create(folder_path)
  create_bulk_data(coh1_count_mat, anti_type_meta, folder_path, dataset_name, "even")
  
  dataset_name <- paste0("final_no_", type, "mirror")
  folder_path <- paste0("/nfs/data/tcell_deconvolution/data/bassez/final_data/PBMC/final_no_", type, "mirror")
  dir.create(folder_path)
  create_bulk_data(coh1_count_mat, anti_type_meta, folder_path, dataset_name, "mirror_db")
  
  
  
})

# overrep

create_bulk_data_overrep <- function(count_mat, metadata, folder_path, dataset_name, weighted_type){
  
  count_mat <- count_mat[,colnames(count_mat) %in% metadata$Cell]
  
  annotation <- data.frame("ID"=metadata$Cell,
                           "cell_type"=metadata$cellType,
                           row.names = metadata$Cell)
  sim_ds <- SimBu::dataset(annotation = annotation,
                           count_matrix = count_mat,
                           tpm_matrix = counts_to_cpm(count_mat),
                           name = "bassez",
                           scale_tpm = FALSE)
  simulation <- SimBu::simulate_bulk(data = sim_ds,
                                     scenario = "weighted",
                                     weighted_cell_type = weighted_type,
                                     weighted_amount = 0.95,
                                     scaling_factor = "NONE",
                                     ncells = 10000,
                                     nsamples = 10,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                     run_parallel = FALSE)
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]],
          file = file.path(folder_path, paste0(dataset_name, "_tpm.rds")))
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
          file = file.path(folder_path, paste0(dataset_name, "_counts.rds")))
  
  facs <- simulation$cell_fractions %>% t()
  saveRDS(facs, file.path(folder_path, paste0(dataset_name, "_facs.rds")))
  
  
}

lapply(celltypes, function(type){
  dataset_name <- paste0("final_only_", type)
  folder_path <- paste0("/nfs/data/tcell_deconvolution/data/bassez/final_data/PBMC/final_only_", type)
  dir.create(folder_path)
  
  create_bulk_data_overrep(coh1_count_mat, updated_metadata, folder_path, dataset_name, type)
})
