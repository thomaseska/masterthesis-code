library(omnideconv)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(SimBu)
library(gridExtra)
library(goseq)

prepare_benchmark_pipeline_data <- function(metadata, sc_count_mat,
                                            folder_path, dataset_name,
                                            simbu_ncells=10000, simbu_nsamples=10){
  # set.seed(123)
  # randomly sample cells: keep all t cells, take 1000 of all other types
  # t cell subtypes are preceded by "TCELL_". Subtype is not always available

  metadata <- metadata[cellType != "T_cell"]
  celltypes <- metadata$cellType %>% unique()
  
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
                                            simbu_ncells=10000, simbu_nsamples=10){
  set.seed(123)
  # randomly sample cells: keep all t cells, take 1000 of all other types
  # t cell subtypes are preceded by "TCELL_". Subtype is not always available
  
  celltypes <- metadata$cellType %>% unique()
  
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
# put tcell descriptors into cohort1
tcell_metadata <- fread("data/bassez/1870-BIOKEY_metaData_tcells_cohort1_web.xls")
coh1_metadata <- fread("data/bassez/1872-BIOKEY_metaData_cohort1_web.xls")
tcell_metadata_lookup <- tcell_metadata[,c("Cell", "cellSubType")]
tcell_metadata_lookup$cellSubType <- paste0("TCELL_", tcell_metadata_lookup$cellSubType)
inds <- match(coh1_metadata$Cell, tcell_metadata_lookup$Cell)
coh1_metadata[!is.na(inds)]$cellType <- tcell_metadata_lookup$cellSubType[na.omit(inds)]

coh1_count_mat <- readRDS("data/bassez/1863-counts_cells_cohort1.rds")

prepare_benchmark_pipeline_data(metadata = coh1_metadata,
                                sc_count_mat = coh1_count_mat,
                                folder_path = "data/bassez/benchmark_data",
                                dataset_name = "bassezEven")

prepare_benchmark_pipeline_data_DWLS(metadata = coh1_metadata,
                                sc_count_mat = coh1_count_mat,
                                folder_path = "data/bassez/benchmark_data",
                                dataset_name = "bassezEvenDWLS")

# analysis of omnideconv results

deconv <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalEven_counts_ct0_rep0/deconvolution.rds")
metrics <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalEven_counts_ct0_rep0/results_metric.rds")
# runtime_deconv <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_deconvolution.rds")
# runtime_signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_signature.rds")
# signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/signature.rds")


# compare predicted and true fractions
ground_truth <- metrics$facs_groud_truth %>% as.data.table(keep.rownames = "cell_type")
ground_truth <- data.table::melt(ground_truth, id.vars = "cell_type",
                  variable.name = "sample",
                  value.name = "fraction")

# for missing TCELL_CD4_EX
to_add <- data.table("cell_type"="TCELL_CD4_EX",
                     "sample"=paste0("mirror_db_sample", 1:10),
                     "fraction"=0)
ground_truth <- rbind(ground_truth, to_add)

deconv <- as.data.table(deconv, keep.rownames = "sample" )
predicted <- data.table::melt(deconv, id.vars = "sample",
                  variable.name = "cell_type",
                  value.name = "fraction")

fractions <- merge(ground_truth, predicted,
                   by = c("cell_type", "sample"), all = TRUE,
                   suffixes = c(".true", ".predicted"))

ggplot(fractions, aes(fraction.true, fraction.predicted))+
  geom_point(aes(color=cell_type))+
  stat_cor()+
  geom_abline(slope = 1)+
  #coord_equal(ratio = 1)+
  labs(title = "Predicted Fractions VS. True Fraction per Cell Type
       Bassez, DWLS, all cell types, 10000 cell/sample")

ggplot(fractions, aes(cell_type, fraction.true))+
  geom_col(aes(color = "True"))+
  geom_col(aes(y=fraction.predicted, color = "Predicted"), alpha = 0.5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_color_manual(breaks = c("True", "Predicted"), values = c("red", "blue")) #+
  #facet_wrap(vars(sample))

p1 <- ggplot(fractions, aes(sample, fraction.true, fill = cell_type))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p2 <- ggplot(fractions, aes(sample, fraction.predicted, fill = cell_type))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

grid.arrange(p1, p2)

rmse_by_type <- metrics$rmse_cell_type %>% as.data.table()
ggplot(rmse_by_type, aes(cell_type, RMSE))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

fractions$diff <- fractions$fraction.true - fractions$fraction.predicted

mean_diff_by_type <- fractions[, mean(diff), by = cell_type]
setnames(mean_diff_by_type, "V1", "mean_diff")
ggplot(mean_diff_by_type, aes(cell_type, mean_diff))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(fractions, aes(sample, cell_type, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "True Fractions - \nPredicted Fractions")+
  labs(title = "Difference between True and Deconvoluted Fractions per Sample and Cell Type
                Bassez, DWLS, All cell types, Even Fractions, 10000 cell/sample")
  


# no spillover

deconv_ALL <- readRDS("data/bassez/deconvolution_results/DWLS/deconvolution.rds")
metrics_ALL <- readRDS("data/bassez/deconvolution_results/DWLS/results_metric.rds")
# runtime_deconv <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_deconvolution.rds")
# runtime_signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/runtime_signature.rds")
# signature <- readRDS("data/bassez/deconvolution_results/bayesprism2/signature.rds")


# compare predicted and true fractions
ground_truth_ALL <- metrics_ALL$facs_groud_truth %>% as.data.table(keep.rownames = "cell_type")
ground_truth_ALL <- data.table::melt(ground_truth_ALL, id.vars = "cell_type",
                     variable.name = "sample",
                     value.name = "fraction")
deconv_ALL <- as.data.table(deconv_ALL, keep.rownames = "sample" )
predicted_ALL <- data.table::melt(deconv_ALL, id.vars = "sample",
                  variable.name = "cell_type",
                  value.name = "fraction")

fractions_ALL <- merge(ground_truth_ALL, predicted_ALL,
                   by = c("cell_type", "sample"), all = TRUE,
                   suffixes = c(".true", ".predicted"))

fractions_ALL$diff <- (fractions_ALL$fraction.true - fractions_ALL$fraction.predicted)# * fractions_ALL$fraction.true

ggplot(fractions_ALL, aes(sample, cell_type, fill=diff))+geom_tile()+
  #geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "True Fractions - \nPredicted Fractions")+
  labs(title = "Difference between True and Deconvoluted Fractions per Sample and Cell Type
                Bassez, DWLS, all cell types, 10000 cell/sample")

ggplot(fractions_ALL, aes(fraction.true, fraction.predicted))+
  geom_point(aes(color=cell_type))+
  stat_cor()+
  # geom_smooth(method = "lm")+
  geom_abline(slope = 1)+
  #coord_equal(ratio = 1)+
  labs(title = "Predicted Fractions VS. True Fraction per Cell Type
       Bassez, DWLS, all cell types, 10000 cell/sample")

ggplot(fractions_ALL[, mean(fraction.true), by=cell_type], aes(cell_type, V1))+
  geom_col()+
  #geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
