# SimBu Tutorial
# http://omnideconv.org/SimBu/articles/simulator_documentation.html

library(SimBu)

counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

colnames(counts) <- paste0("cell_",rep(1:300))
colnames(tpm) <- paste0("cell_",rep(1:300))
rownames(counts) <- paste0("gene_",rep(1:1000))
rownames(tpm) <- paste0("gene_",rep(1:1000))

annotation <- data.frame("ID"=paste0("cell_",rep(1:300)), 
                         "cell_type"=c(rep("T cells CD4",50), 
                                       rep("T cells CD8",50),
                                       rep("Macrophages",100),
                                       rep("NK cells",10),
                                       rep("B cells",70),
                                       rep("Monocytes",20)))

ds <- SimBu::dataset(annotation = annotation,
                     count_matrix = counts,
                     tpm_matrix = tpm,
                     name = "test_dataset")

simulation <- SimBu::simulate_bulk(data = ds,
                                   scenario = "random", 
                                   scaling_factor = "NONE", 
                                   ncells=100, 
                                   nsamples = 10, 
                                   BPPARAM = BiocParallel::MulticoreParam(workers = 4),  #this will use 4 threads to run the simulation
                                   run_parallel = TRUE)                                 #multi-threading to TRUE
dim(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]])