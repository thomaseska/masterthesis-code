#!/usr/bin/Rscript

library(docopt)
library(data.table)
library(dplyr)
library(SimBu)
library(Matrix)
source("~/masterthesis/scripts/helpers.R")

"Usage:
  create_datasets.R <sc_count_mat> <metadata> <folder_path> <dataset_name> <method> <spillover> <overrep>
Options:
<sc_count_mat> path to sc matrix
<metadata> path to metadata of sc matrix
<folder_path> path to folder to save output to
<dataset_name> name of the dataset
<method> create datasets for dwls or bayesprism
<spillover> logical, whether to create spillover datasets
<overrep> logical, whether to create overrepresentation datasets" -> doc



args <- docopt::docopt(doc)



prepare_benchmark_pipeline_data <- function(metadata, sc_count_mat,
                                            folder_path, dataset_name,
                                            simbu_ncells=10000, simbu_nsamples=10,
                                            sc=TRUE){
  # set.seed(123)
  # randomly sample cells: keep all t cells, take 1000 of all other types
  # t cell subtypes are preceded by "TCELL_". Subtype is not always available

  metadata <- metadata[cellType != "T_cell"]
  sc_count_mat <- sc_count_mat[, colnames(sc_count_mat) %in% metadata$Cell]
  celltypes <- metadata$cellType %>% unique()
  
  if (sc) {

    sampled_cells <- sapply(celltypes, function(type){
        if (startsWith(type, "TCELL_")) {
        size <- metadata[cellType==type]$Cell %>% length()
        } else if (type == "T_cell"){
        size <- 0
        } else {
        size <-min(1000, metadata[cellType==type]$Cell %>% length())
        }
        sampled <- sample(metadata[cellType==type]$Cell, size = size, replace = FALSE)
        return(sampled)
    })
    sampled_cells <- unlist(sampled_cells)
    
    sampled_metadata <- metadata[Cell %in% sampled_cells]
    dir.create(folder_path)
    fwrite(sampled_metadata, file.path(folder_path, "sampled_cells_metadata.csv"))
    
    print(sampled_metadata)
    
    
    # prepare pipeline data
    # - sampled_metadata and sampled_counts are created in data_exploration.R
    # - bulk data is created in pseudobulk.R
    
    # Single Cell
    # need matrix_counts.rds, batch.rds and celltype_annotation.rds
    
    dir.create(file.path(folder_path, "singleCell"))
    dir.create(file.path(folder_path, "singleCell", dataset_name))
    
    sampled_sc_count_mat <- sc_count_mat[,colnames(sc_count_mat) %in% sampled_cells]
    saveRDS(sampled_sc_count_mat, file.path(folder_path, "singleCell", dataset_name, 
                                    "matrix_counts.rds"))
    # print(dim(sc_count_mat))
    
    batch <- sampled_metadata[, c("Cell", "patient_id")]
    batch <- batch[order(match(Cell, colnames(sampled_sc_count_mat)))]
    batch <- batch$patient_id
    saveRDS(batch, file.path(folder_path, "singleCell", dataset_name, "batch.rds"))
    
    celltype_annotation <- sampled_metadata[, c("Cell", "cellType")]
    celltype_annotation <- celltype_annotation[order(match(Cell, 
                                                            colnames(sampled_sc_count_mat))
                                                    )]
    celltype_annotation <- celltype_annotation$cellType
    saveRDS(celltype_annotation, file.path(folder_path, "singleCell", 
                                            dataset_name, "celltype_annotations.rds"))
  }
  # Bulk
  # need simulated bulk data and data_facs.rds (celltype fractions)
  
  dir.create(file.path(folder_path, "PBMC"))
  dir.create(file.path(folder_path, "PBMC", dataset_name))
  
  annotation <- data.frame("ID"=metadata$Cell,
                           "cell_type"=metadata$cellType,
                           row.names = metadata$Cell)
  sim_ds <- SimBu::dataset(annotation = annotation,
                           count_matrix = sc_count_mat,
                           name = "bassez")
  simulation <- SimBu::simulate_bulk(data = sim_ds,
                                     scenario = "even",
                                     scaling_factor = "NONE",
                                     ncells = simbu_ncells,
                                     nsamples = simbu_nsamples,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                     run_parallel = FALSE)
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
          file = file.path(folder_path, "PBMC", dataset_name,
                           paste0(dataset_name, "_counts.rds")))
  
  
  facs <- simulation$cell_fractions %>% t()
  saveRDS(facs, file.path(folder_path, "PBMC", dataset_name,
                          paste0(dataset_name, "_facs.rds")))
  
  
}

prepare_benchmark_pipeline_data_DWLS <- function(metadata, sc_count_mat,
                                            folder_path, dataset_name,
                                            simbu_ncells=10000, simbu_nsamples=10,
                                            sc=TRUE){
  # set.seed(123)
  # randomly sample cells: keep all t cells, take 1000 of all other types
  # t cell subtypes are preceded by "TCELL_". Subtype is not always available
  metadata <- metadata[cellType != "T_cell"]
  celltypes <- metadata$cellType %>% unique()
  sc_count_mat <- sc_count_mat[, colnames(sc_count_mat) %in% metadata$Cell]
  if (sc) {

    sampled_cells <- sapply(celltypes, function(type){
      if (startsWith(type, "TCELL_")) {
        size <- min(metadata[cellType==type]$Cell %>% length(), 100)
      } else if (type == "T_cell"){
        size <- 0
      } else {
        size <-min(100, metadata[cellType==type]$Cell %>% length())
      }
      sampled <- sample(metadata[cellType==type]$Cell, size = size, replace = FALSE)
      return(sampled)
    })
    sampled_cells <- unlist(sampled_cells)
    
    sampled_metadata <- metadata[Cell %in% sampled_cells]
    dir.create(folder_path)
    fwrite(sampled_metadata, file.path(folder_path, "sampled_cells_metadata.csv"))
    
    print(sampled_metadata)
    
    
    # prepare pipeline data
    # - sampled_metadata and sampled_counts are created in data_exploration.R
    # - bulk data is created in pseudobulk.R
    
    # Single Cell
    # need matrix_counts.rds, batch.rds and celltype_annotation.rds
    
    dir.create(file.path(folder_path, "singleCell"))
    dir.create(file.path(folder_path, "singleCell", dataset_name))
    
    sampled_sc_count_mat <- sc_count_mat[,colnames(sc_count_mat) %in% sampled_cells]
    saveRDS(sampled_sc_count_mat, file.path(folder_path, "singleCell", dataset_name, 
                                    "matrix_counts.rds"))
    # print(dim(sc_count_mat))
    
    batch <- sampled_metadata[, c("Cell", "patient_id")]
    batch <- batch[order(match(Cell, colnames(sampled_sc_count_mat)))]
    batch <- batch$patient_id
    saveRDS(batch, file.path(folder_path, "singleCell", dataset_name, "batch.rds"))
    
    celltype_annotation <- sampled_metadata[, c("Cell", "cellType")]
    celltype_annotation <- celltype_annotation[order(match(Cell, 
                                                          colnames(sampled_sc_count_mat))
    )]
    celltype_annotation <- celltype_annotation$cellType
    saveRDS(celltype_annotation, file.path(folder_path, "singleCell", 
                                          dataset_name, "celltype_annotations.rds"))
    
  }

  # Bulk
  # need simulated bulk data and data_facs.rds (celltype fractions)
  
  dir.create(file.path(folder_path, "PBMC"))
  dir.create(file.path(folder_path, "PBMC", dataset_name))
  

  
  annotation <- data.frame("ID"=metadata$Cell,
                           "cell_type"=metadata$cellType,
                           row.names = metadata$Cell)
  sim_ds <- SimBu::dataset(annotation = annotation,
                           count_matrix = sc_count_mat,
                           tpm_matrix = counts_to_cpm(sc_count_mat),
                           name = "bassez")
  simulation <- SimBu::simulate_bulk(data = sim_ds,
                                     scenario = "even",
                                     scaling_factor = "NONE",
                                     ncells = simbu_ncells,
                                     nsamples = simbu_nsamples,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                     run_parallel = FALSE)
  
  saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]],
          file = file.path(folder_path, "PBMC", dataset_name,
                           paste0(dataset_name, "_tpm.rds")))
  
  
  facs <- simulation$cell_fractions %>% t()
  saveRDS(facs, file.path(folder_path, "PBMC", dataset_name,
                          paste0(dataset_name, "_facs.rds")))
  
  
}

create_spillover_datasets <- function(metadata, sc_count_mat,
                                      folder_path, dataset_name, method,
                                      simbu_ncells=10000, simbu_nsamples=10){
  
  tcelltypes <- metadata[startsWith(cellType, "TCELL_")]$cellType %>% unique()
  
  sapply(tcelltypes, function(tcell){
    print(paste("Dataset for celltype", tcell))
    new_met <- metadata[cellType != tcell]
    new_cells <- new_met$Cell
    new_count <- sc_count_mat[,colnames(sc_count_mat) %in% new_cells]
    
    new_dataset_name <- paste(dataset_name, "minus", tcell, sep = "_")
    
    if (method == "dwls") {
      prepare_benchmark_pipeline_data_DWLS(metadata = new_met,
                                    sc_count_mat = new_count,
                                    folder_path = folder_path,
                                    dataset_name = new_dataset_name,
                                    simbu_ncells = simbu_ncells,
                                    simbu_nsamples = simbu_nsamples,
                                    sc=FALSE)
    } else if (method == "bayesprism") {
       prepare_benchmark_pipeline_data(metadata = new_met,
                                    sc_count_mat = new_count,
                                    folder_path = folder_path,
                                    dataset_name = new_dataset_name,
                                    simbu_ncells = simbu_ncells,
                                    simbu_nsamples = simbu_nsamples,
                                    sc=FALSE)
    }
    
  })
}

create_overrep_datasets <- function(metadata, sc_count_mat,
                                      folder_path, dataset_name, method,
                                      simbu_ncells=10000, simbu_nsamples=10){
  
  tcelltypes <- metadata[startsWith(cellType, "TCELL_")]$cellType %>% unique()
  
  sapply(tcelltypes, function(tcell){
    print(paste("Dataset for celltype", tcell))
    new_dataset_name <- paste(dataset_name, "overrep", tcell, sep = "_")
    
    prepare_overrep_data(metadata = metadata,
                         sc_count_mat = sc_count_mat,
                         folder_path = folder_path,
                         dataset_name = new_dataset_name,
                         method = method,
                         overrep_celltype = tcell,
                         simbu_ncells = simbu_ncells,
                         simbu_nsamples = simbu_nsamples)
    
  })
}


prepare_overrep_data <- function(metadata,
                                 sc_count_mat,
                                 folder_path,
                                 dataset_name,
                                 method,
                                 overrep_celltype,
                                 simbu_ncells = 10000,
                                 simbu_nsamples = 10) {
  
  dir.create(file.path(folder_path, "PBMC"))
  dir.create(file.path(folder_path, "PBMC", dataset_name))
  
  annotation <- data.frame("ID"=metadata$Cell,
                           "cell_type"=metadata$cellType,
                           row.names = metadata$Cell)
  if (method == "bayesprism") {
    sim_ds <- SimBu::dataset(annotation = annotation,
                             count_matrix = sc_count_mat,
                             name = "bassez")
    simulation <- SimBu::simulate_bulk(data = sim_ds,
                                     scenario = "weighted",
                                     weighted_cell_type = overrep_celltype,
                                     weighted_amount = 0.99,
                                     scaling_factor = "NONE",
                                     ncells = simbu_ncells,
                                     nsamples = simbu_nsamples,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                     run_parallel = FALSE)
  
    saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
          file = file.path(folder_path, "PBMC", dataset_name,
                           paste0(dataset_name, "_counts.rds")))
  } else if (method == "dwls") {
    sim_ds <- SimBu::dataset(annotation = annotation,
                             count_matrix = sc_count_mat,
                             tpm_matrix = counts_to_cpm(sc_count_mat),
                             name = "bassez")
    simulation <- SimBu::simulate_bulk(data = sim_ds,
                                       scenario = "weighted",
                                       weighted_cell_type = overrep_celltype,
                                       weighted_amount = 0.8,
                                       scaling_factor = "NONE",
                                       ncells = simbu_ncells,
                                       nsamples = simbu_nsamples,
                                       BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                       run_parallel = FALSE)
    
    saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]],
            file = file.path(folder_path, "PBMC", dataset_name,
                             paste0(dataset_name, "_tpm.rds")))
  }
  
  facs <- simulation$cell_fractions %>% t()
  saveRDS(facs, file.path(folder_path, "PBMC", dataset_name,
                          paste0(dataset_name, "_facs.rds")))
  
  
}



sc_count_mat <- readRDS(args$sc_count_mat)
metadata <- fread(args$metadata)

if(args$method == "bayesprism") {
    prepare_benchmark_pipeline_data(metadata = metadata,
                                    sc_count_mat = sc_count_mat,
                                    folder_path = args$folder_path,
                                    dataset_name = args$dataset_name)
} else if (args$method == "dwls") {
    prepare_benchmark_pipeline_data_DWLS(metadata = metadata,
                                        sc_count_mat = sc_count_mat,
                                        folder_path = args$folder_path,
                                        dataset_name = args$dataset_name)
}

if (args$spillover) {
  create_spillover_datasets(metadata = metadata,
                            sc_count_mat = sc_count_mat,
                            folder_path = args$folder_path,
                            dataset_name = args$dataset_name,
                            method = args$method)
}

if (args$overrep) {
  create_overrep_datasets(metadata = metadata,
                          sc_count_mat = sc_count_mat,
                          folder_path = args$folder_path,
                          dataset_name = args$dataset_name,
                          method = args$method)
}
