library(SimBu)
library(data.table)
library(dplyr)

# create Seurat object in data_exploration.R

sim_ds <- SimBu::dataset_seurat(seurat_obj = tcell_seur_obj,
                                count_assay = "RNA",
                                cell_id_col = "Cell",
                                cell_type_col = "cellSubType")

# via counts
annotation <- data.frame("ID"=sampled_metadata$Cell, "cell_type"=sampled_metadata$cellType,
                         row.names = sampled_metadata$Cell)
sim_ds <- SimBu::dataset(annotation = annotation,
                         count_matrix = coh1_count_mat,
                         name = "bassez")
simulation <- SimBu::simulate_bulk(data = sim_ds,
                                   scenario = "mirror_db",
                                   scaling_factor = "NONE",
                                   ncells = 100,
                                   nsamples = 10,
                                   BPPARAM = BiocParallel::MulticoreParam(workers = 3),
                                   run_parallel = TRUE)
SimBu::plot_simulation(simulation)

saveRDS(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]],
        file = "data/bassez/sampled_simbu_bulk_counts.rds")

