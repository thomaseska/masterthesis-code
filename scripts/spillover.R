# Spillover Analysis:
# use all Cell types to compute signatures
# use all cell types except X as data to deconvolute for all X
# check how much of cell type X was detected (=spillover from other cell types)

create_spillover_datasets <- function(metadata, sc_count_mat,
                                      folder_path, dataset_name,
                                      simbu_ncells=10000, simbu_nsamples=10){
  
  tcelltypes <- metadata[startsWith(cellType, "TCELL_")]$cellType %>% unique()
  
  sapply(tcelltypes, function(tcell){
    print(paste("Dataset for celltype", tcell))
    new_met <- metadata[cellType != tcell]
    new_cells <- new_met$Cell
    new_count <- sc_count_mat[,colnames(sc_count_mat) %in% new_cells]
    
    new_dataset_name <- paste(dataset_name, "minus", tcell, sep = "_")
    
    prepare_benchmark_pipeline_data_DWLS(metadata = new_met,
                                    sc_count_mat = new_count,
                                    folder_path = folder_path,
                                    dataset_name = new_dataset_name,
                                    simbu_ncells = simbu_ncells,
                                    simbu_nsamples = simbu_nsamples)
    
  })
}

create_spillover_datasets(metadata = coh1_metadata,
                          sc_count_mat = coh1_count_mat,
                          folder_path = "data/bassez/benchmark_data",
                          dataset_name = "bassezEvenDWLS")


relevant_types <- c("TCELL_CD4_EM", "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating",
                    "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                    "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating", "TCELL_CD8_RM")

folders <- list.dirs(path = "data/bassez/deconvolution_results/dwls_bassez_even_spillover/",
                     recursive = FALSE, full.names = TRUE)

fractions <- lapply(folders, function(folder){
  deconv <- readRDS(file.path(folder, "deconvolution.rds"))
  metrics <- readRDS(file.path(folder, "results_metric.rds"))
  # runtime_deconv <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_deconvolution.rds")
  # runtime_signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_signature.rds")
  # signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/signature.rds")
  if (grepl("_minus_", folder, fixed = TRUE)) {
    namebuild <- strsplit(folder, "/")[[1]]
    namebuild <- namebuild[length(namebuild)]
    namebuild <- strsplit(namebuild, "_", fixed = TRUE)[[1]]
    
    celltype <- namebuild[(which(namebuild=="minus")+1):(length(namebuild)-3)]
    celltype <- paste(celltype, collapse = "_")
    to_add <- data.table("cell_type"=celltype,
                         "sample"=paste0("even_sample", 1:10),
                         "fraction"=0)
    
    namebuild <- namebuild[which(namebuild=="minus"):(length(namebuild)-3)]
    name <- paste("dataset", paste(namebuild, collapse = "_"), sep = "_")
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
  labs(title = "Difference between True and Deconvoluted Fractions in spillover Datasets
                Bassez, DWLS, Even Fractions, 10000 cell/sample")
