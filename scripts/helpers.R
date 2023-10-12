# miscellaneous functions

library(biomaRt)
library(Matrix)


# generate tpms from counts
counts_to_tpm <- function(count_mat) {
  # get gene lengths
  genes <- rownames(count_mat)
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  positions <- getBM(attributes = c("hgnc_symbol", "start_position", "end_position"),
                     filters = c("hgnc_symbol"), values = genes, mart = mart)
  positions$gene_length <- positions$end_position - positions$start_position
  positions <- positions[!duplicated(positions$hgnc_symbol),]
  print(dim(count_mat))
  print(length(positions$hgnc_symbol))
  # rpk: divide counts by gene length in kb
  count_mat <- count_mat / (positions$gene_length / 1000)
  # scaling factor: sum rpk per cell, divide by 1,000,000
  scal_fac <- colSums(count_mat) / 1e6
  # tpm: divide rpk by scaling factor
  # dgcMatrix@x stores none-zero matrix values in 1d
  # dgcMatrix@p stores cumulative number of non-zero element per column
  # --> only scale non-zero values, keep sparse format
  count_mat@x <- count_mat@x / rep.int(scal_fac, diff(count_mat@p))
  print(dim(count_mat))
  return(count_mat)
}

counts_to_cpm <- function(count_mat) {
  # scaling factor: sum counts per cell, divide by 1,000,000
  scal_fac <- colSums(count_mat) / 1e6
  # cpm: divide counts by scaling factor
  # dgcMatrix@x stores none-zero matrix values in 1d
  # dgcMatrix@p stores cumulative number of non-zero element per column
  # --> only scale non-zero values, keep sparse format
  count_mat@x <- count_mat@x / rep.int(scal_fac, diff(count_mat@p))
  return(count_mat)
}

pretty_names <- data.table(data_name=c("B_cell", "Cancer_cell", "Endothelial_cell", "Fibroblast",
                                       "Mast_cell", "Myeloid_cell", "pDC", "TCELL_CD4_EM",
                                       "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating", "TCELL_CD4_N",
                                       "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                                       "TCELL_CD8_EMRA", "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating",
                                       "TCELL_CD8_N", "TCELL_CD8_RM", "TCELL_gdT", "TCELL_NK_CYTO", 
                                       "TCELL_NK_REST", "TCELL_Vg9Vd2_gdT"),
                           plot_name=as.factor(c("B cell", "Cancer cell", "Endothelial cell", "Fibroblast",
                                                 "Mast cell", "Myeloid cell", "pDC", "(T) CD4+ EM",
                                                 "(T) CD4+ EX", "(T) CD4+ EX Proliferating", "(T) CD4+ Naive",
                                                 "(T) CD4+ Reg", "(T) CD4+ Reg Proliferating", "(T) CD8+ EM",
                                                 "(T) CD8+ EMRA", "(T) CD8+ EX", "(T) CD8+ EX Proliferating",
                                                 "(T) CD8+ Naive", "(T) CD8+ RM", "(T) gdT", "(T) NK CYTO", 
                                                 "(T) NK REST", "(T) Vg9Vd2 gdT")),
                           type=c("Non T Cell", "Non T Cell", "Non T Cell", "Non T Cell",
                                  "Non T Cell", "Non T Cell", "Non T Cell", "T Cell CD4",
                                  "T Cell CD4", "T Cell CD4", "T Cell CD4",
                                  "T Cell CD4", "T Cell CD4", "T Cell CD8",
                                  "T Cell CD8", "T Cell CD8", "T Cell CD8",
                                  "T Cell CD8", "T Cell CD8", "T Cell Other", "T Cell Other", 
                                  "T Cell Other", "T Cell Other"),
                           relevant=c(FALSE, FALSE, FALSE, FALSE,
                                      FALSE, FALSE, FALSE, TRUE,
                                      TRUE, TRUE, FALSE,
                                      TRUE, TRUE, TRUE,
                                      FALSE, TRUE, TRUE,
                                      FALSE, TRUE, FALSE, FALSE, 
                                      FALSE, FALSE))                      

pretty_names_short <- data.table(data_name=c("B_cell", "Cancer_cell", "Endothelial_cell", "Fibroblast",
                                       "Mast_cell", "Myeloid_cell", "pDC", "T_cell"),
                           plot_name=as.factor(c("B cell", "Cancer cell", "Endothelial cell", "Fibroblast",
                                                 "Mast cell", "Myeloid cell", "pDC", "T & NK cell")))                      


sample_lookup <- data.table(sample=c(paste0("even_sample", 1:10), paste0("mirror_db_sample", 1:10)),
                            sample_new=factor(c(rep(paste0("sample", 1:10), 2)), levels = paste0("sample", 1:10)))
