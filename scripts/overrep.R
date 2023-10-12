

create_100percent_datasets <- function(metadata, sc_count_mat,
                                       folder_path, dataset_name, method,
                                       simbu_ncells=10000, simbu_nsamples=10){
  tcelltypes <- metadata[startsWith(cellType, "TCELL_")]$cellType %>% unique()
  
  sapply(tcelltypes, function(tcell){
    print(paste("Dataset for celltype", tcell))
    new_dataset_name <- paste(dataset_name, "only", tcell, sep = "_")

    new_meta <- metadata[cellType == tcell]
    new_sc_count <- sc_count_mat[,colnames(sc_count_mat) %in% new_meta$Cell]
    

    # Bulk
    # need simulated bulk data and data_facs.rds (celltype fractions)
    
    dir.create(file.path(folder_path, "PBMC"))
    dir.create(file.path(folder_path, "PBMC", new_dataset_name))
    print("XX")
    annotation <- data.frame("ID"=new_meta$Cell,
                             "cell_type"=new_meta$cellType,
                             row.names = new_meta$Cell)
    print(dim(new_sc_count))
    print(dim(annotation))
    sim_ds <- SimBu::dataset(annotation = annotation,
                             count_matrix = new_sc_count,
                             tpm_matrix = counts_to_cpm(new_sc_count),
                             name = "bassez")
    print("ZZ")
    simulation <- SimBu::simulate_bulk(data = sim_ds,
                                       scenario = "even",
                                       scaling_factor = "NONE",
                                       ncells = simbu_ncells,
                                       nsamples = simbu_nsamples,
                                       BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                       run_parallel = FALSE)
    
    saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
            file = file.path(folder_path, "PBMC", new_dataset_name,
                             paste0(new_dataset_name, "_counts.rds")))
    saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]],
            file = file.path(folder_path, "PBMC", new_dataset_name,
                             paste0(new_dataset_name, "_tpm.rds")))
    
    
    facs <- simulation$cell_fractions %>% t()
    saveRDS(facs, file.path(folder_path, "PBMC", new_dataset_name,
                            paste0(new_dataset_name, "_facs.rds")))
    
    
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
                                     weighted_amount = 0.8,
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


folders <- list.dirs(path = "data/bassez/deconvolution_results/bayesprism_bassezEven_overrep/",
                     recursive = FALSE, full.names = TRUE)

fractions <- lapply(folders, function(folder){
  deconv <- readRDS(file.path(folder, "deconvolution.rds"))
  metrics <- readRDS(file.path(folder, "results_metric.rds"))
  # runtime_deconv <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_deconvolution.rds")
  # runtime_signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_signature.rds")
  # signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/signature.rds")
  if (grepl("_overrep_", folder, fixed = TRUE)) {
    namebuild <- strsplit(folder, "/")[[1]]
    namebuild <- namebuild[length(namebuild)]
    namebuild <- strsplit(namebuild, "_", fixed = TRUE)[[1]]
    
#    celltype <- namebuild[(which(namebuild=="overrep")+1):(length(namebuild)-3)]
#    celltype <- paste(celltype, collapse = "_")
#    to_add <- data.table("cell_type"=celltype,
#                         "sample"=paste0("weighted_sample", 1:10),
#                         "fraction"=0)
    to_add <- setNames(data.table(matrix(nrow = 0, ncol = 3)), c("cell_type", "sample", "fraction"))
    
    namebuild <- namebuild[(which(namebuild=="overrep")+1):(length(namebuild)-3)]
    name <- paste("dataset_with_80%", paste(namebuild, collapse = "_"), sep = "_")
  } else {
    name <- "dataset_full"
    to_add <- setNames(data.table(matrix(nrow = 0, ncol = 3)), c("cell_type", "sample", "fraction"))
  }
  
  
  # compare predicted and true fractions
  ground_truth <- metrics$facs_groud_truth %>% as.data.table(keep.rownames = "cell_type")
  ground_truth <- data.table::melt(ground_truth, id.vars = "cell_type",
                                   variable.name = "sample",
                                   value.name = "fraction")
  
  ground_truth <- rbind(ground_truth, to_add)
  
  deconv <- as.data.table(deconv, keep.rownames = "sample" )
  predicted <- data.table::melt(deconv, id.vars = "sample",
                                variable.name = "cell_type",
                                value.name = "fraction")
  
  fractions <- merge(ground_truth, predicted,
                     by = c("cell_type", "sample"), all = TRUE,
                     suffixes = c(".true", ".predicted"))
  fractions$diff <- fractions$fraction.true - fractions$fraction.predicted
  
  fractions$dataset <- name
  
  print(fractions)
})
fractions <- rbindlist(fractions)
fractions <- fractions[cell_type!="T_cell"]
mean_fractions <- fractions[, lapply(.SD, mean),
                            by=c("cell_type", "dataset"),
                            .SDcols = c("fraction.true", "fraction.predicted", "diff")]
sd_fractions <- fractions[, lapply(.SD, sd),
                          by=c("cell_type", "dataset"),
                          .SDcols = c("fraction.true", "fraction.predicted", "diff")]


ggplot(mean_fractions, aes(dataset, cell_type, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "True Fractions - \nPredicted Fractions")+
  labs(title = "Difference between True and Deconvoluted Fractions in Overrepresentation Datasets
                Bassez, Bayesprism, Even Fractions, 10000 cell/sample")

ggplot(mean_fractions[dataset=="dataset_with_80% _TCELL_CD4_EX"], aes(cell_type, fraction.true))+
  geom_col()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  labs(title = "Cell input Fractions")

